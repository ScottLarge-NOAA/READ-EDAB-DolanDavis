# ==============================================================================
# R/04_visualizations.R
# Purpose: Visualizes outputs from the Dolan-Davis spatial modeling pipeline.
#          Generates absolute footprint timelines, standardized climate anomaly 
#          (Z-score) charts, and maps of time-invariant spatial fields (omega).
# ==============================================================================

# --- 1. SETUP & DATA INGESTION ------------------------------------------------
library(dplyr)
library(readr)
library(ggplot2)
library(sf)
library(sdmTMB)

out_dir <- here::here("outputs", "models", "dolan_davis")
fig_dir <- here::here("outputs", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Load indices and spatial components prepared in 03_spatialmodel.R
basin_class_indices <- readRDS(file.path(out_dir, "dolan_davis_raw_basin_indices.rds"))
load(file.path(out_dir, "viz_dependencies.RData"))

# --- 2. ABSOLUTE FOOTPRINT PLOTTING -------------------------------------------
message("Generating Absolute Footprint Plot...")

plot_data <- basin_class_indices |>
  mutate(dd_class = as.factor(dd_class))

storm_plot <- ggplot(plot_data |> filter(EPU != "SS"), aes(x = year, y = est, color = dd_class, fill = dd_class)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  facet_wrap(EPU ~ ., scales = "free_y") +
  scale_color_viridis_d(name = "DD Class", option = "plasma") +
  scale_fill_viridis_d(name = "DD Class", option = "plasma") +
  scale_y_log10() + # Uncomment if Class 1 compresses higher classes
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

ggsave(
  filename = file.path(fig_dir, "dolan_davis_epu_indices.png"),
  plot = storm_plot, width = 10, height = 7, dpi = 300
)

# --- 3. CLIMATE ANOMALY (Z-SCORE) PLOTTING ------------------------------------
message("Generating Climate Anomaly Plot...")

# 3.1 Calculate the NOAA 30-year climate normal baseline (1991-2020)
ref_stats <- basin_class_indices |>
  filter(year >= 1991 & year <= 2020, EPU != "SS") |>
  group_by(EPU, dd_class) |>
  summarize(
    mu_ref = mean(est, na.rm = TRUE),
    sd_ref = sd(est, na.rm = TRUE),
    .groups = "drop"
  )

# 3.2 Standardize against the baseline
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

ggsave(
  filename = file.path(fig_dir, "dolan_davis_anomalies.png"),
  plot = anomaly_plot, width = 12, height = 9, dpi = 300
)

# --- 4. SPATIAL HOTSPOT MAPPING (TIME-AVERAGED) -------------------------------
message("Generating Time-Averaged Spatial Hotspot Maps...")

# Predict across ALL years to ensure ar1, rw, and iid fields are fully captured
omega_list <- lapply(seq_len(nrow(best_models)), function(k) {
  target_class <- best_models$dd_class[k]
  fit_file <- best_models$final_file[k]
  if (!file.exists(fit_file)) return(NULL)
  fit_obj <- readRDS(fit_file)
  
  # Predict onto the full multi-year grid
  pred_full <- predict(fit_obj, newdata = pred_grid)
  
  # Calculate the true spatial average over time for each cell
  pred_full |>
    group_by(X_km, Y_km) |>
    summarize(
      # Average the log-link estimates across all years
      mean_log_est = mean(est, na.rm = TRUE), 
      .groups = "drop"
    ) |>
    mutate(
      dd_class = target_class,
      # Mean-center the estimates within the class to create a divergent +/- scale 
      # (Mimics the omega_s visualization exactly, but accurately captures AR1 fields)
      relative_hotspot = mean_log_est - mean(mean_log_est, na.rm = TRUE)
    )
})

omega_df <- bind_rows(omega_list) |>
  mutate(dd_class_label = paste("Class", dd_class))

# Format coordinates for sf plotting
omega_sf <- omega_df |>
  mutate(X_m = X_km * 1000, Y_m = Y_km * 1000) |>
  st_as_sf(coords = c("X_m", "Y_m"), crs = 32619)

omega_plot <- ggplot() +
  geom_sf(data = omega_sf, aes(color = relative_hotspot), size = 1.2, shape = 15) +
  geom_sf(data = epu_water, fill = NA, color = "black", linewidth = 0.5) +
  geom_sf(data = coastline_sf, fill = "grey30", color = NA) +
  facet_wrap(~ dd_class_label, nrow = 1) +
  scale_color_gradient2(
    low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
    name = "Relative\nHotspot"
  ) +
  labs(
    title = "Time-Averaged Spatial Hotspots by Storm Intensity Class",
    subtitle = "Red: Above-average chronic footprint | Blue: Below-average chronic footprint",
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

ggsave(
  filename = file.path(fig_dir, "dolan_davis_hotspot_maps.png"),
  plot = omega_plot, 
  width = 16, 
  height = 5, 
  dpi = 300
)

message("Visualization pipeline complete.")