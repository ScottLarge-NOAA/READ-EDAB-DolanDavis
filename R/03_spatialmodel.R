# ==============================================================================
# R/02_sdm_dolan_davis_tournament.R
# Purpose: Analytical pipeline to identify the optimal spatial mesh cutoff for 
#          each Dolan-Davis storm intensity class (1-5).
#          Fits full models, runs spatial CV, retains the best model per class
#          (highest ELPD), and cleans up suboptimal model files from disk.
# ==============================================================================

library(dplyr)
library(tidyr)
library(sf)
library(sdmTMB)
library(sdmTMBextra)
library(future.apply)
library(R.utils)
library(RhpcBLASctl)

# --- 1. SETUP & DIRECTORIES ---
out_dir <- here::here("outputs", "models", "dolan_davis")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

analytical_model_df <- readRDS(here::here("data", "processed", "analytical_model_df.rds"))

# Generate spatially independent folds on the FULL dataset
set.seed(42)
k_folds <- 5
analytical_model_df$fold_id <- kmeans(
  analytical_model_df[, c("X_km", "Y_km")], 
  centers = k_folds
)$cluster

# --- 2. BARRIER & PRIOR CONFIGURATION ---
epu_water <- NEFSCspatial::epu_sf |> st_make_valid() |> st_transform(32619)

land_sf <- rnaturalearth::ne_countries(scale = 10, continent = "North America", returnclass = "sf") |>
  st_make_valid() |> st_transform(32619)

epu_bbox <- st_as_sfc(st_bbox(epu_water))
coastline_sf <- st_intersection(land_sf, epu_bbox) |> st_make_valid()

x_range <- diff(range(analytical_model_df$X_km, na.rm = TRUE))
y_range <- diff(range(analytical_model_df$Y_km, na.rm = TRUE))
max_dist_km <- sqrt(x_range^2 + y_range^2)

spatial_priors <- sdmTMBpriors(
  matern_s = pc_matern(range_gt = max_dist_km * 0.1, sigma_lt = 5, range_prob = 0.05, sigma_prob = 0.05),
  matern_st = pc_matern(range_gt = max_dist_km * 0.1, sigma_lt = 5, range_prob = 0.05, sigma_prob = 0.05)
)

# --- 3. PARALLEL SETUP & GRID ---
plan(multisession, workers = 14L)
mesh_cutoffs <- c(25, 50, 75, 100, 125)
storm_classes <- c(1:5)
st_types <- c("iid", "ar1", "rw")

registry <- expand_grid(mesh_cutoffs, storm_classes, st_types) |>
  mutate(model_id = row_number())

message(sprintf("Starting tournament for %d mesh/class combinations...", nrow(registry)))

# --- 4. EVALUATION PIPELINE ---
results_list <- future_lapply(seq_len(nrow(registry)), function(i) {
  
  TMB::openmp(n = 1) 
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
  
  dd_c <- registry$storm_classes[i]
  cutoff_val <- registry$mesh_cutoffs[i]
  m_id <- registry$model_id[i]
  st_t <- registry$st_types[i]
  
  # Subset data first, so fold_ids match the row count during CV
  sub_dat <- analytical_model_df |> filter(dd_class == dd_c)
  
  mesh_i <- make_mesh(data = sub_dat, xy_cols = c("X_km", "Y_km"), cutoff = cutoff_val)
  bmesh_i <- add_barrier_mesh(spde_obj = mesh_i, barrier_sf = coastline_sf, range_fraction = 0.2, proj_scaling = 1000)
  
  check_convergence <- function(m) {
    if (is.null(m)) return(FALSE)
    if (!isTRUE(m$pos_def_hessian) || !isTRUE(suppressWarnings(max(abs(m$gradients), na.rm = TRUE)) < 0.001)) return(FALSE)
    sane <- tryCatch({ all(unlist(suppressMessages(suppressWarnings(sdmTMB::sanity(m))))) }, error = function(e) FALSE)
    if (!isTRUE(sane)) return(FALSE)
    ran_pars <- tryCatch(tidy(m, effects = "ran_pars"), error = function(e) NULL)
    if (!is.null(ran_pars) && isTRUE(any(ran_pars$estimate[ran_pars$term %in% c("sigma_E", "sigma_O")] > 10, na.rm = TRUE))) return(FALSE)
    return(TRUE)
  }
  
  # Temporary file path for this specific run
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
  
  # 4.2 Cross-Validation (Only if full model converges)
  if (check_convergence(fit)) {
    saveRDS(fit, temp_file) # Save the converged full model to disk
    
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
    
    # # --- RMSE CALCULATION ADDITION ---
    # if (!is.null(cv_res) && "cv_predicted" %in% names(cv_res$data)) {
    #   # Exponentiate the log-link predictions back to response space (counts)
    #   out_of_sample_preds <- exp(cv_res$data$cv_predicted)
    #   # Calculate Root Mean Square Error against the observed raw counts
    #   cv_rmse <- sqrt(mean((cv_res$data$n_storms - out_of_sample_preds)^2, na.rm = TRUE))
    # } else {
    #   cv_rmse <- NA_real_
    # }
    
    mod_aic <- AIC(fit)
    
  } else {
    cv_elpd <- NA_real_
    # cv_rmse <- NA_real_
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
    # rmse = cv_rmse,
    file_path = temp_file
  )
}, future.seed = TRUE, future.packages = c("sdmTMB", "dplyr", "sf", "TMB", "RhpcBLASctl"))

plan(sequential)

# --- 5. MODEL SELECTION & DISK CLEANUP ---
message("Tournament complete. Identifying best models and cleaning up...")

all_results <- bind_rows(results_list)

# Find the best cutoff for each DD class based on highest ELPD
best_models <- all_results |> 
  filter(converged == TRUE, !is.na(elpd)) |> 
  group_by(dd_class) |> 
  slice_max(elpd, n = 1, with_ties = FALSE) |> 
  ungroup() |>
  mutate(final_file = file.path(out_dir, paste0("best_fit_class_", dd_class, ".rds")))

# Rename the winning files and delete the losers
for (i in seq_len(nrow(all_results))) {
  current_file <- all_results$file_path[i]
  
  if (file.exists(current_file)) {
    # If this row is a winner, rename the file to the final production name
    if (all_results$model_id[i] %in% best_models$model_id) {
      target_name <- best_models$final_file[best_models$model_id == all_results$model_id[i]]
      file.rename(current_file, target_name)
    } else {
      # If it's a loser, delete it to save disk space
      file.remove(current_file)
    }
  }
}

write_csv(all_results, file.path(out_dir, "tournament_results_all.csv"))
write_csv(best_models, file.path(out_dir, "tournament_results_best.csv"))
# ==============================================================================
# --- 6. STANDARDIZED INDEX EXTRACTION ---
# ==============================================================================
message("Step 6: Calculating Area-Integrated EPU Footprints from Winning Models...")

# 6.1 Construct the 10x10 km marine prediction grid
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
  mutate(log_effort = 0) # Standardize to 1 unit of effort

epu_names <- unique(pred_grid$EPU)

# 6.2 Iterate over the winning models and extract indices using Area Masking
basin_indices_list <- lapply(seq_len(nrow(best_models)), function(k) {
  
  target_class <- best_models$dd_class[k]
  fit_file <- best_models$final_file[k]

  
  if (!file.exists(fit_file)) return(NULL)
  
  fit_obj <- readRDS(fit_file)
  
  # Predict onto the standard grid
  pred_obj <- predict(
    fit_obj, 
    newdata = pred_grid, 
    return_tmb_object = TRUE,
    offset = pred_grid$log_effort
  )
  
  # Extract area-integrated index for each basin
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


# ==============================================================================
# --- 7. ABSOLUTE FOOTPRINT PLOTTING ---
# ==============================================================================

message("Step 7: Generating Visualizations...")

fig_dir <- here::here("outputs", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

plot_data <- basin_class_indices |>
  mutate(dd_class = as.factor(dd_class))

# Construct the time series plot
storm_plot <- ggplot(plot_data |> filter(EPU != "SS"), aes(x = year, y = est, color = dd_class, fill = dd_class)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  facet_wrap(EPU ~ ., scales = "free_y") +
  scale_color_viridis_d(name = "DD Class", option = "plasma") +
  scale_fill_viridis_d(name = "DD Class", option = "plasma") +
  # OPTIONAL: Uncomment the line below if Class 1 visually compresses the higher classes
  scale_y_log10() +
  labs(
    x = "Year",
    y = "Area-Integrated Storm Footprint",
    title = "Standardized Storm Intensity Indices by Basin"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey95", color = "grey50"),
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank()
)

print(storm_plot)
ggsave(
  filename = file.path(fig_dir, "dolan_davis_epu_indices.png"),
  plot = storm_plot, width = 10, height = 7, dpi = 300
)
# ==============================================================================
# --- 8. CLIMATE ANOMALY (Z-SCORE) PLOTTING ---
# ==============================================================================
# 8.1 Calculate the NOAA 30-year climate normal baseline (1991-2020)
ref_stats <- basin_class_indices |>
  filter(year >= 1991 & year <= 2020, EPU != "SS") |>
  group_by(EPU, dd_class) |>
  summarize(
    mu_ref = mean(est, na.rm = TRUE),
    sd_ref = sd(est, na.rm = TRUE),
    .groups = "drop"
)

# 8.2 Apply the baseline to the entire time series
anomaly_df <- basin_class_indices |>
  filter(EPU != "SS") |> 
  left_join(ref_stats, by = c("EPU", "dd_class")) |>
  mutate(
    z_score = (est - mu_ref) / sd_ref,
    anomaly_sign = ifelse(z_score > 0, "Positive", "Negative")
)

write_csv(anomaly_df, file.path(out_dir, "dolan_davis_zscore_anomalies_1991_2020.csv"))

anomaly_plot <- ggplot(anomaly_df |> filter(EPU != "SS"), aes(x = year, y = z_score, fill = anomaly_sign)) +
  geom_col() +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  annotate("rect", xmin = 1991, xmax = 2020, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "black") +
  facet_grid(dd_class ~ EPU) +
  scale_fill_manual(values = c("Positive" = "#b2182b", "Negative" = "#2166ac")) +
  labs(
    x = "Year",
    y = "Storm Footprint Anomaly (Z-Score)",
    title = "Spatiotemporal Storm Anomalies (Dolan-Davis Index)",
    subtitle = "Standardized against the 1991-2020 NOAA Climate Normal (shaded). \nRows: Intensity Class | Columns: EPU Basin"
  ) +
  theme_bw() +
  theme(
    legend.position = "none", 
    strip.background = element_rect(fill = "grey95", color = "grey50"),
    strip.text = element_text(face = "bold", size = 9),
    panel.grid.minor = element_blank()
  )
print(anomaly_plot)

ggsave(
  filename = file.path(fig_dir, "dolan_davis_anomalies.png"),
  plot = anomaly_plot, width = 12, height = 9, dpi = 300
)


# 1. Prepare a single-year spatial prediction grid
# The spatial field (omega_s) does not change by year, so one year is sufficient
spatial_grid_single <- pred_grid |>
  filter(year == min(year)) |>
  select(EPU, X_km, Y_km, year, log_effort)

# 2. Extract omega_s from the winning model for each class
omega_list <- lapply(seq_len(nrow(best_models)), function(k) {
  
  target_class <- best_models$dd_class[k]
  fit_file <- best_models$final_file[k]
  
  if (!file.exists(fit_file)) return(NULL)
  
  fit_obj <- readRDS(fit_file)
  
  # Predict onto the single-year grid
  pred_obj <- predict(fit_obj, newdata = spatial_grid_single)
  
  # Extract coordinates and the spatial random effect column (omega_s)
  pred_obj |>
    select(X_km, Y_km, omega_s) |>
    mutate(dd_class = target_class)
})

omega_df <- bind_rows(omega_list) |>
  mutate(dd_class_label = paste("Class", dd_class))

# 3. Convert prediction coordinates back to an sf object for mapping
omega_sf <- omega_df |>
  # Convert km back to meters to match the EPSG:32619 CRS of your polygons
  mutate(X_m = X_km * 1000, Y_m = Y_km * 1000) |>
  st_as_sf(coords = c("X_m", "Y_m"), crs = 32619)

# 4. Generate the spatial map
omega_plot <- ggplot() +
  # Add the gridded spatial field predictions
  geom_sf(data = omega_sf, aes(color = omega_s), size = 1.2, shape = 15) +
  # Overlay the EPU boundaries for geographic reference
  geom_sf(data = epu_water, fill = NA, color = "black", linewidth = 0.5) +
  # Overlay the terrestrial barrier mesh to visualize the boundaries
  geom_sf(data = coastline_sf, fill = "grey30", color = NA) +
  # Facet by storm class horizontally
  facet_wrap(~ dd_class_label, nrow = 1) +
  # Use a divergent color palette strictly centered on zero
  scale_color_gradient2(
    low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
    name = "Spatial\nEffect (ω)"
  ) +
  labs(
    title = "Time-Invariant Spatial Fields (ω) by Storm Intensity Class",
    subtitle = "Red: Chronic baseline storm hotspots | Blue: Chronic baseline coldspots",
    x = NULL, 
    y = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = "grey95", color = "grey50"),
    strip.text = element_text(face = "bold", size = 11),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

print(omega_plot)

# Save the high-resolution map
ggsave(
  filename = file.path(fig_dir, "dolan_davis_omega_maps.png"),
  plot = omega_plot, 
  width = 16, 
  height = 5, 
  dpi = 300
)