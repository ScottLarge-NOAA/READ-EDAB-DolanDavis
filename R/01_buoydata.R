# ==============================================================================
# DOLAN-DAVIS STORM INDEX: SPATIAL SUBSETTING & DATA EXTRACTION
# ==============================================================================

library(tidyverse)
library(sf)
library(arrow)
library(lubridate)

# 1. Define Paths & Endpoints
cache_dir <- here::here("data", "buoy_cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

coastwatch_base <- "https://coastwatch.pfeg.noaa.gov/erddap/tabledap/cwwcNDBCMet.csv"
current_year <- lubridate::year(Sys.Date())

# ==============================================================================
# 2. SPATIAL SUBSETTING
# ==============================================================================
message("Mapping spatial boundaries and identifying target buoys...")

epu <- read_sf(file.path("data/data-raw/EPU_NOESTUARIES.shp")) %>% 
  st_transform(4326) %>% st_make_valid()

buoy_locations <- buoydata::buoy_data %>% 
  as_tibble() %>%
  st_as_sf(coords = c("LON", "LAT"), crs = 4326, remove = FALSE) %>%
  st_transform(crs = st_crs(epu)) |> 
  rename(buoy_id = ID)

# Inner join to retain only buoys within EPU polygons
epu_wave_buoys <- buoy_locations %>%
  st_join(epu, join = st_intersects, left = FALSE) %>%
  filter(grepl("^[0-9]{5}$", buoy_id))

target_stations <- unique(epu_wave_buoys$buoy_id)
message(sprintf("Identified %d target wave buoys within EPU spatial footprints.", length(target_stations)))

# ==============================================================================
# 3. PRIMARY EXTRACTION: COASTWATCH REST API
# ==============================================================================
message("Executing primary extraction via CoastWatch REST API...")

purrr::walk(target_stations, function(current_buoy) {
  output_file <- file.path(cache_dir, paste0("buoy_", current_buoy, ".parquet"))
  if (file.exists(output_file)) return()
  
  query_url <- paste0(coastwatch_base, "?station,time,WVHT&station=%22", current_buoy, "%22&WVHT%3E=0")
  
  max_retries <- 3
  attempt <- 1
  backoff <- c(10, 30, 60)
  
  while (attempt <= max_retries) {
    res <- tryCatch({
      suppressMessages(
        readr::read_csv(
          query_url,
          skip = 2, 
          col_names = c("station", "time", "WVHT"),
          col_types = "cTc" 
        )
      )
    }, error = function(e) {
      if (grepl("404|Not Found", conditionMessage(e))) return("SKIP")
      return(conditionMessage(e))
    })
    
    if (is.data.frame(res)) {
      if (nrow(res) > 0) {
        clean_data <- res %>%
          mutate(
            station = as.character(station),
            WVHT = as.numeric(WVHT)
          ) %>%
          filter(!is.na(WVHT)) 
        
        if (nrow(clean_data) > 0) {
          arrow::write_parquet(clean_data, output_file)
          message(sprintf("  [API Success] %s: Cached %d observations.", current_buoy, nrow(clean_data)))
        }
      }
      break
    }
    
    if (identical(res, "SKIP")) break
    
    Sys.sleep(backoff[attempt])
    attempt <- attempt + 1
  }
})

# ==============================================================================
# 4. SECONDARY EXTRACTION: NDBC FLAT-FILE RESCUE
# ==============================================================================
# Identify stations that failed the ERDDAP pull (due to 503s or being dropped from the mirror)
cached_files <- list.files(cache_dir, pattern = "\\.parquet$")
cached_stations <- gsub("buoy_|\\.parquet", "", cached_files)
missing_stations <- setdiff(target_stations, cached_stations)

if (length(missing_stations) > 0) {
  message(sprintf("Initiating NDBC historical flat-file rescue for %d missing stations...", length(missing_stations)))
  
  get_ndbc_year <- function(station, yr) {
    url <- sprintf("https://www.ndbc.noaa.gov/data/historical/stdmet/%sh%d.txt.gz", station, yr)
    lines <- tryCatch(readLines(url, warn = FALSE), error = function(e) NULL)
    if (is.null(lines) || length(lines) < 3) return(NULL)
    
    df <- tryCatch(read.table(text = lines, header = TRUE, fill = TRUE, comment.char = ""), error = function(e) NULL)
    if (is.null(df)) return(NULL)
    
    names(df) <- toupper(gsub("^#", "", names(df)))
    if (grepl("[A-Za-z]", df[1, 1])) df <- df[-1, ]
    if (!"WVHT" %in% names(df)) return(NULL)
    
    time_cols <- df %>% select(any_of(c("YY", "YYYY", "MM", "DD", "HH", "MIN")))
    yr_val <- if ("YYYY" %in% names(time_cols)) as.numeric(time_cols$YYYY) else {
      y <- as.numeric(time_cols$YY)
      ifelse(y < 50, y + 2000, y + 1900)
    }
    
    tibble(
      station = as.character(station),
      time = ISOdatetime(yr_val, as.numeric(time_cols$MM), as.numeric(time_cols$DD), 
                         as.numeric(time_cols$HH), if("MIN" %in% names(time_cols)) as.numeric(time_cols$MIN) else 0, 
                         0, tz = "UTC"),
      WVHT = as.numeric(df$WVHT)
    ) %>% filter(WVHT < 99, !is.na(WVHT), !is.na(time))
  }
  
  purrr::walk(missing_stations, function(current_buoy) {
    output_file <- file.path(cache_dir, paste0("buoy_", current_buoy, ".parquet"))
    buoy_data_list <- list()
    
    # Scrape from inception of NDBC program to current year
    for (yr in 1970:current_year) {
      yr_data <- get_ndbc_year(current_buoy, yr)
      if (!is.null(yr_data) && nrow(yr_data) > 0) buoy_data_list[[as.character(yr)]] <- yr_data
    }
    
    if (length(buoy_data_list) > 0) {
      clean_data <- bind_rows(buoy_data_list)
      arrow::write_parquet(clean_data, output_file)
      message(sprintf("  [Rescue Success] %s: Cached %d observations.", current_buoy, nrow(clean_data)))
    } else {
      message(sprintf("  [Confirmed Missing] %s: No historical wave data exists on primary server.", current_buoy))
    }
  })
} else {
  message("All target stations successfully extracted via primary API.")
}

message("Data extraction pipeline complete.")