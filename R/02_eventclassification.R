library(tidyverse)
library(arrow)
library(sf)
library(lubridate)
library(here)

# A. Spatial Metadata: Intersect Buoys with EPUs
epu <- read_sf(file.path("data/data-raw/EPU_NOESTUARIES.shp")) %>% 
  st_transform(4326) %>% st_make_valid()

# Project buoy metadata to match the EPU CRS strictly before spatial joining
buoy_locations <- buoydata::buoy_data |> 
  dplyr::as_tibble() |> 
  sf::st_as_sf(coords = c("LON", "LAT"), crs = 4326, remove = FALSE) |> 
  sf::st_transform(crs = sf::st_crs(epu)) |> 
  sf::st_make_valid() |> 
  rename(buoy_id = ID)

# Spatial inner join: strictly retain buoys inside EPU polygons and extract metadata
epu_wave_buoys <- buoy_locations |> 
  sf::st_join(epu, join = sf::st_intersects, left = FALSE) |> 
  sf::st_drop_geometry() |> 
  dplyr::select(station = buoy_id, EPU, LON, LAT) # Ensure column names match time series

# B. Load Time Series and Merge Spatial Metadata
buoy_data <- arrow::open_dataset(here::here("data", "buoy_cache")) |> 
  dplyr::collect() |> 
  # Ensure strict chronological order per buoy for run-length encoding
  dplyr::arrange(station, time) |> 
  # Filter to EPU buoys and attach spatial metadata
  dplyr::inner_join(epu_wave_buoys, by = "station") |> 
  mutate(year = lubridate::year(time))

# Calculate annual observation effort per buoy
# Rationale: Standardizing by true year length prevents leap-year aliasing in the log-offset
buoy_effort <- buoy_data |> 
  dplyr::mutate(
    # Truncate time to the top of the hour to handle high-frequency data
    time_hr = lubridate::floor_date(time, "hour")
  ) |> 
  dplyr::group_by(station, year, EPU, LON, LAT) |> 
  dplyr::summarize(
    # Count unique hours observed rather than total rows
    n_unique_hours = dplyr::n_distinct(time_hr),
    yr_length_hrs = dplyr::if_else(lubridate::leap_year(year[1]), 366 * 24, 365 * 24),
    obs_effort = n_unique_hours / yr_length_hrs,
    .groups = "drop"
  )

# C. Calculate Observation Effort
# Standardizing by true year length prevents leap-year aliasing in the log-offset
pot_threshold_m <- 2.0
# Define maximum allowable gap in hours before splitting an event
max_gap_hrs <- 12 

storm_events <- buoy_data |> 
  dplyr::group_by(station, year) |> 
  dplyr::arrange(station, time) |> 
  dplyr::mutate(
    is_storm = WVHT >= pot_threshold_m,
    # Calculate time difference between consecutive observations
    time_gap_hrs = as.numeric(difftime(time, dplyr::lag(time, default = min(time)), units = "hours")),
    # Trigger a new event if the wave height drops below threshold OR if there is a massive data gap
    event_break = !is_storm | (is_storm & time_gap_hrs > max_gap_hrs),
    event_id = cumsum(event_break)
  ) |> 
  dplyr::filter(is_storm) |> 
  dplyr::group_by(station, year, EPU, event_id) |> 
  dplyr::summarize(
    # Explicit time-based duration: max(time) - min(time) + 1 hour for the bounding hour
    duration_hrs = as.numeric(difftime(max(time), min(time), units = "hours")) + 1,
    max_wvht_m = max(WVHT, na.rm = TRUE),
    n_observations = dplyr::n(), # Track actual data points vs temporal duration
    .groups = "drop"
  ) |> 
  dplyr::mutate(
    # Dolan-Davis power index
    storm_power = (max_wvht_m^2) * duration_hrs,
    
    # DD92 Classification
    dd_class = dplyr::case_when(
      storm_power >= 4700 ~ 5, 
      storm_power >= 1500 ~ 4, 
      storm_power >= 370  ~ 3, 
      storm_power >= 110  ~ 2, 
      storm_power >  0    ~ 1, 
      TRUE ~ NA_real_
    )
  )

head(storm_events)

# Aggregate distinct storms by buoy, year, and intensity class
storm_counts <- storm_events |> 
  dplyr::group_by(station, year, dd_class) |> 
  dplyr::summarize(n_storms = dplyr::n(), .groups = "drop")

# Create a complete factorial grid to retain true zeros for unobserved event classes
grid_base <- tidyr::expand_grid(
  dd_class = 1:5,
  buoy_effort |> dplyr::select(station, year, EPU, LON, LAT, obs_effort)
)

final_analytical_data <- grid_base |> 
  dplyr::left_join(storm_counts, by = c("station", "year", "dd_class")) |> 
  dplyr::mutate(
    n_storms = tidyr::replace_na(n_storms, 0),
    log_effort = log(obs_effort) # Log-offset representing exposure duration
  )

# Enforce geometry integrity and project coordinates to Cartesian km
analytical_sf <- final_analytical_data |> 
  sf::st_as_sf(coords = c("LON", "LAT"), crs = 4326) |> 
  sf::st_transform(crs = 32619) |>  
  sf::st_make_valid()

# Extract projected X/Y arrays directly into the dataframe for sdmTMB formula mapping
analytical_model_df <- analytical_sf |> 
  dplyr::mutate(
    X_km = sf::st_coordinates(analytical_sf)[, 1] / 1000,
    Y_km = sf::st_coordinates(analytical_sf)[, 2] / 1000
  ) |> 
  sf::st_drop_geometry()

# 1. Ensure categorical predictors are correctly typed for sdmTMB
analytical_model_df <- analytical_model_df |> 
  dplyr::mutate(
    # sdmTMB requires factors for categorical fixed effects
    dd_class = factor(dd_class, levels = 1:5, ordered = TRUE),
    EPU = as.factor(EPU),
    year_f = as.factor(year) # Useful if modeling year as an IID random intercept
  )

# 2. Create output directory if it doesn't exist
fs::dir_create(here::here("data", "processed"))

# 3. Save artifacts
saveRDS(
  analytical_model_df, 
  here::here("data", "processed", "analytical_model_df.rds")
)

saveRDS(
  storm_events, 
  here::here("data", "processed", "storm_events_pot.rds")
)
