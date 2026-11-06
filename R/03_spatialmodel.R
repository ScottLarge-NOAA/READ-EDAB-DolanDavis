# ==============================================================================
# R/03_spatialmodel.R
# Purpose: Analytical pipeline to identify the optimal spatial mesh cutoff for 
#          each Dolan-Davis storm intensity class (1-5). Fits full models, 
#          runs spatial CV, retains the best model per class (highest ELPD), 
#          cleans up suboptimal files, and extracts area-integrated indices.
# ==============================================================================

# --- 1. SETUP & LIBRARIES -----------------------------------------------------
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(sf)
library(sdmTMB)
library(sdmTMBextra)
library(future.apply)
library(R.utils)
library(RhpcBLASctl)

# Establish directories
out_dir <- here::here("outputs", "models", "dolan_davis")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Load analytical data
analytical_model_df <- readRDS(here::here("data", "processed", "analytical_model_df.rds"))

# Generate spatially independent folds on the FULL dataset using k-means
set.seed(42)
k_folds <- 5
analytical_model_df$fold_id <- kmeans(
  analytical_model_df[, c("X_km", "Y_km")], 
  centers = k_folds
)$cluster

# --- 2. BARRIER & PRIOR CONFIGURATION -----------------------------------------
# Prepare EPU boundaries in UTM Zone 19N
epu_water <- NEFSCspatial::epu_sf |> 
  st_make_valid() |> 
  st_transform(32619)

# Prepare North American coastline for spatial barriers
land_sf <- rnaturalearth::ne_countries(
  scale = 10, 
  continent = "North America", 
  returnclass = "sf"
) |>
  st_make_valid() |> 
  st_transform(32619)

# Crop landmass to EPU bounding box to optimize mesh calculations
epu_bbox <- st_as_sfc(st_bbox(epu_water))
coastline_sf <- st_intersection(land_sf, epu_bbox) |> 
  st_make_valid()

# Calculate domain range for PC Priors
x_range <- diff(range(analytical_model_df$X_km, na.rm = TRUE))
y_range <- diff(range(analytical_model_df$Y_km, na.rm = TRUE))
max_dist_km <- sqrt(x_range^2 + y_range^2)

# Formulate Penalized Complexity (PC) Priors
spatial_priors <- sdmTMBpriors(
  matern_s  = pc_matern(range_gt = max_dist_km * 0.1, sigma_lt = 5, range_prob = 0.05, sigma_prob = 0.05),
  matern_st = pc_matern(range_gt = max_dist_km * 0.1, sigma_lt = 5, range_prob = 0.05, sigma_prob = 0.05)
)

# --- 3. PARALLEL SETUP & GRID -------------------------------------------------
plan(multisession, workers = 14L)
mesh_cutoffs  <- c(25, 50, 75, 100, 125)
storm_classes <- c(1:5)
st_types      <- c("iid", "ar1", "rw")

registry <- expand_grid(
  mesh_cutoffs = mesh_cutoffs, 
  storm_classes = storm_classes, 
  st_types = st_types
) |>
  mutate(model_id = row_number())

message(sprintf("Starting tournament for %d mesh/class combinations...", nrow(registry)))

# --- 4. EVALUATION PIPELINE ---------------------------------------------------
results_list <- future_lapply(seq_len(nrow(registry)), function(i) {
  
  # Prevent nested parallelization RAM bloat
  TMB::openmp(n = 1) 
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
  
  dd_c       <- registry$storm_classes[i]
  cutoff_val <- registry$mesh_cutoffs[i]
  m_id       <- registry$model_id[i]
  st_t       <- registry$st_types[i]
  
  # Subset data to align fold_ids properly during CV
  sub_dat <- analytical_model_df |> filter(dd_class == dd_c)
  
  # Construct candidate mesh and barrier
  mesh_i <- make_mesh(data = sub_dat, xy_cols = c("X_km", "Y_km"), cutoff = cutoff_val)
  bmesh_i <- add_barrier_mesh(
    spde_obj = mesh_i, 
    barrier_sf = coastline_sf, 
    range_fraction = 0.2, 
    proj_scaling = 1000
  )
  
  # 3-Tier Convergence Check
  check_convergence <- function(m) {
    if (is.null(m)) return(FALSE)
    if (!isTRUE(m$pos_def_hessian) || !isTRUE(suppressWarnings(max(abs(m$gradients), na.rm = TRUE)) < 0.001)) return(FALSE)
    sane <- tryCatch({ all(unlist(suppressMessages(suppressWarnings(sdmTMB::sanity(m))))) }, error = function(e) FALSE)
    if (!isTRUE(sane)) return(FALSE)
    ran_pars <- tryCatch(broom.mixed::tidy(m, effects = "ran_pars"), error = function(e) NULL)
    if (!is.null(ran_pars) && isTRUE(any(ran_pars$estimate[ran_pars$term %in% c("sigma_E", "sigma_O")] > 10, na.rm = TRUE))) return(FALSE)
    return(TRUE)
  }
  
  temp_file <- file.path(out_dir, paste0("temp_class_", dd_c, "_cutoff_", cutoff_val, "_", st_t, ".rds"))
  
  # 4.1 Fit Full Model
  fit <- tryCatch({
    R.utils::withTimeout({
      sdmTMB(
        formula = n_storms ~ 1,
        data = sub_dat,
        mesh = bmesh_i,
        time = "year",
        spatial = "on",
        spatiotemporal = st_t,
        family = nbinom2(link = "log"),
        offset = "log_effort",
        priors = spatial_priors,
        control = sdmTMBcontrol(multiphase = TRUE, newton_loops = 1, eval.max = 2000, iter.max = 2000, get_joint_precision = TRUE)
      )
    }, timeout = 3600, onTimeout = "error")
  }, error = function(e) NULL)
  
  if (!is.null(fit) && !check_convergence(fit)) {
    fit <- tryCatch(run_extra_optimization(fit, newton_loops = 1, nlminb_loops = 1), error = function(e) fit)
  }
  
  # 4.2 Cross-Validation
  if (check_convergence(fit)) {
    saveRDS(fit, temp_file)
    
    cv_res <- tryCatch({
      R.utils::withTimeout({
        sdmTMB_cv(
          formula = n_storms ~ 1,
          data = sub_dat,
          mesh = bmesh_i,
          time = "year",
          spatial = "on",
          spatiotemporal = st_t,
          family = nbinom2(link = "log"),
          offset = "log_effort",
          fold_ids = sub_dat$fold_id, 
          k_folds = k_folds,
          priors = spatial_priors,
          control = sdmTMBcontrol(newton_loops = 1)
        )
      }, timeout = 7200, onTimeout = "error")
    }, error = function(e) NULL)
    
    cv_elpd <- if(!is.null(cv_res)) cv_res$sum_loglik else NA_real_
    mod_aic <- AIC(fit)
    
  } else {
    cv_elpd <- NA_real_
    mod_aic <- NA_real_
  }
  
  # 4.3 Return metrics
  tibble(
    model_id = m_id,
    dd_class = dd_c,
    st_type = st_t,
    cutoff = cutoff_val,
    nodes = mesh_i$mesh$n,
    converged = check_convergence(fit),
    aic = mod_aic,
    elpd = cv_elpd,
    file_path = temp_file
  )
}, future.seed = TRUE, future.packages = c("sdmTMB", "dplyr", "sf", "TMB", "RhpcBLASctl", "broom.mixed"))

plan(sequential)

# --- 5. MODEL SELECTION & DISK CLEANUP ----------------------------------------
message("Tournament complete. Identifying best models and cleaning up...")

all_results <- bind_rows(results_list)

# Retain highest ELPD model per DD class
best_models <- all_results |> 
  filter(converged == TRUE, !is.na(elpd)) |> 
  group_by(dd_class) |> 
  slice_max(elpd, n = 1, with_ties = FALSE) |> 
  ungroup() |>
  mutate(final_file = file.path(out_dir, paste0("best_fit_class_", dd_class, ".rds")))

# Promote winning files and purge losers
for (i in seq_len(nrow(all_results))) {
  current_file <- all_results$file_path[i]
  
  if (file.exists(current_file)) {
    if (all_results$model_id[i] %in% best_models$model_id) {
      target_name <- best_models$final_file[best_models$model_id == all_results$model_id[i]]
      file.rename(current_file, target_name)
    } else {
      file.remove(current_file)
    }
  }
}

write_csv(all_results, file.path(out_dir, "tournament_results_all.csv"))
write_csv(best_models, file.path(out_dir, "tournament_results_best.csv"))

# --- 6. STANDARDIZED INDEX EXTRACTION -----------------------------------------
message("Step 6: Calculating Area-Integrated EPU Footprints from Winning Models...")

# Construct the 10x10 km marine prediction grid
spatial_grid <- st_make_grid(epu_water, cellsize = c(10000, 10000), what = "centers") |> st_as_sf()
marine_grid_sf <- st_join(spatial_grid, epu_water, join = st_intersects, left = FALSE)
coords <- st_coordinates(marine_grid_sf)

grid_base <- marine_grid_sf |>
  st_drop_geometry() |>
  mutate(
    X_km = coords[, "X"] / 1000,
    Y_km = coords[, "Y"] / 1000,
    cell_area_km2 = 100
  ) |>
  select(EPU, X_km, Y_km, cell_area_km2)

pred_grid <- expand_grid(
  grid_base,
  year = seq(min(analytical_model_df$year), max(analytical_model_df$year))
) |>
  mutate(log_effort = 0) # Standardize predictions to 1 unit of effort

epu_names <- unique(pred_grid$EPU)

# Iterate over winning models and extract basin indices using Area Masking
basin_indices_list <- lapply(seq_len(nrow(best_models)), function(k) {
  
  target_class <- best_models$dd_class[k]
  fit_file <- best_models$final_file[k]
  
  if (!file.exists(fit_file)) return(NULL)
  
  fit_obj <- readRDS(fit_file)
  
  pred_obj <- predict(
    fit_obj, 
    newdata = pred_grid, 
    return_tmb_object = TRUE,
    offset = pred_grid$log_effort
  )
  
  class_epu_indices <- map_dfr(epu_names, function(target_epu) {
    target_areas <- ifelse(pred_grid$EPU == target_epu, pred_grid$cell_area_km2, 0)
    idx <- get_index(obj = pred_obj, area = target_areas, bias_correct = TRUE)
    idx |> mutate(EPU = target_epu, dd_class = target_class)
  })
  
  rm(fit_obj, pred_obj); gc()
  return(class_epu_indices)
})

basin_class_indices <- bind_rows(basin_indices_list)
write_rds(basin_class_indices, file.path(out_dir, "dolan_davis_raw_basin_indices.rds"))

# Save dependencies required for the downstream visualization script
save(
  epu_water, coastline_sf, pred_grid, best_models, 
  file = file.path(out_dir, "viz_dependencies.RData")
)
message("Pipeline complete: Models and indices saved.")