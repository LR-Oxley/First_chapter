library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(zoo)

terraOptions(memfrac = 0.4, progress = 1)

# ============================================================
# Temperature analysis for one SSP
# Outputs:
# 1. Arctic mean soil temperature 1850–2099
# 2. Temperature profiles with depth
# 3. Thawed months per year
# 4. Monthly seasonal cycle for selected periods
# 5. Combined all-SSP plot from saved CSVs
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Please provide SSP as argument, e.g. 370")
}

ssp <- as.character(args[1])

base_dir <- "monthly_mineralised"
out_dir <- file.path(base_dir, "temperature_analysis_60N")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

temp_file <- file.path(
  paste0("60N/mean_verttemp_ssp", ssp, "_60deg.nc")
)

thawed_file <- file.path(
  "60N",
  paste0("arctic_total_thawed_", ssp, "_no_lim_mean.nc")
)

if (!file.exists(temp_file)) {
  stop("Missing temperature file: ", temp_file)
}

if (!file.exists(thawed_file)) {
  stop("Missing thawed file: ", thawed_file)
}

ssp_labels <- c(
  "126" = "SSP1-2.6",
  "245" = "SSP2-4.5",
  "370" = "SSP3-7.0",
  "585" = "SSP5-8.5"
)

ssp_colors_raw <- c(
  "126" = "blue",
  "245" = "orange",
  "370" = "#D73027",
  "585" = "#7B3294"
)

ssp_colors_label <- c(
  "SSP1-2.6" = "blue",
  "SSP2-4.5" = "orange",
  "SSP3-7.0" = "#D73027",
  "SSP5-8.5" = "#7B3294"
)


ssp_label <- ssp_labels[[ssp]]
ssp_color <- ssp_colors_raw[[ssp]]

years <- 1850:2099

depth_vals <- c(
  0.01, 0.035, 0.075, 0.15,
  0.25, 0.4, 0.6, 0.85,
  1.25, 2.0, 3.0, 4.25
)

n_depths_expected <- length(depth_vals)

# ============================================================
# Depth thickness weights
# ============================================================

get_depth_bounds <- function(depth_vals) {
  
  z_top <- numeric(length(depth_vals))
  z_bottom <- numeric(length(depth_vals))
  
  z_top[1] <- 0
  z_bottom[length(depth_vals)] <- 5
  
  mids <- (depth_vals[-1] + depth_vals[-length(depth_vals)]) / 2
  
  z_bottom[-length(depth_vals)] <- mids
  z_top[-1] <- mids
  
  data.frame(
    depth = depth_vals,
    z_top = z_top,
    z_bottom = z_bottom,
    thickness = z_bottom - z_top
  )
}

depth_bounds <- get_depth_bounds(depth_vals)
depth_weights <- depth_bounds$thickness / sum(depth_bounds$thickness)

# ============================================================
# Read input rasters
# ============================================================
# Arctic extent
arctic_ext <- ext(-180, 180, 60, 90)


cat("Reading temperature file:", temp_file, "\n")
temp_full <- rast(temp_file)

cat("Reading thawed file:", thawed_file, "\n")
thawed_total <- rast(thawed_file)
thawed_total <- crop(thawed_total, arctic_ext)
crs(temp_full) <- crs(thawed_total)
cat("Initial geometry check:\n")
print(compareGeom(temp_full[[1]], thawed_total[[1]], stopOnError = FALSE))


depth_metadata <- as.numeric(depth(temp_full))
time_metadata <- as.Date(time(temp_full))

cat("First 24 layer depths:\n")
print(depth_metadata[1:min(24, length(depth_metadata))])

cat("First 24 layer dates:\n")
print(time_metadata[1:min(24, length(time_metadata))])


# ============================================================
# Create analysis mask on temperature grid
# ============================================================

cat("Creating 60N analysis mask on temperature grid...\n")

# Temperature-grid template north of 60°N
template_60N <- crop(temp_full[[1]], arctic_ext)

# Thawed-N data north of 60°N
thawed_total_crop <- crop(thawed_total, arctic_ext)

# Create numeric mask:
# land = 1
# outside analysis domain = NA
analysis_mask_thawgrid <- ifel(
  !is.na(thawed_total_crop[[nlyr(thawed_total_crop)]]),
  1,
  NA
)

# Transfer mask to temperature grid
analysis_mask_tempgrid <- resample(
  analysis_mask_thawgrid,
  template_60N,
  method = "near"
)

# Force any non-land values to NA
analysis_mask_tempgrid <- ifel(
  analysis_mask_tempgrid == 1,
  1,
  NA
)

names(analysis_mask_tempgrid) <- "analysis_mask"

# Cell area on temperature grid
area_rast <- cellSize(template_60N, unit = "m")

# Keep only analysis-mask cells
area_rast <- mask(
  area_rast,
  analysis_mask_tempgrid
)

area_total <- global(
  area_rast,
  "sum",
  na.rm = TRUE
)[1, 1]

cat("Total analysis area [m2]:", area_total, "\n")
cat("Total analysis area [million km2]:",
    area_total / 1e12, "\n")

compareGeom(
  area_rast,
  analysis_mask_tempgrid,
  stopOnError = TRUE
)

# -----------------------------
# Plot mask
# -----------------------------

df_landmask <- as.data.frame(
  analysis_mask_tempgrid,
  xy = TRUE,
  na.rm = TRUE
)

names(df_landmask) <- c(
  "x",
  "y",
  "analysis_mask"
)

# p_mask <- ggplot(
#   df_landmask,
#   aes(x = x, y = y)
# ) +
#   geom_tile(fill = "darkgreen") +
#   coord_equal() +
#   theme_bw() +
#   labs(
#     title = paste0(
#       "Analysis mask north of 60°N – ",
#       ssp_label
#     ),
#     x = "Longitude",
#     y = "Latitude"
#   )
# 
# print(p_mask)
# 
# ggsave(
#   file.path(out_dir, paste0("ssp", ssp, "_analysis_mask.png")),
#   p_mask,
#   width = 8,
#   height = 5,
#   dpi = 300
# )

# ============================================================
# Remove duplicate monthly blocks
# ============================================================

n_depths_full <- length(unique(as.numeric(depth(temp_full))))

if (n_depths_full != n_depths_expected) {
  stop(
    "Depth mismatch: expected ", n_depths_expected,
    " depths but file has ", n_depths_full
  )
}

dates_layer <- as.Date(time(temp_full))

dates_month_all <- dates_layer[seq(1, length(dates_layer), by = n_depths_full)]

cat("Monthly blocks before cleaning:", length(dates_month_all), "\n")
cat("Unique monthly dates:", length(unique(dates_month_all)), "\n")
cat("Duplicate monthly blocks:", sum(duplicated(dates_month_all)), "\n")

keep_month_blocks <- which(!duplicated(dates_month_all))

keep_layers <- unlist(lapply(keep_month_blocks, function(j) {
  ((j - 1) * n_depths_full + 1):(j * n_depths_full)
}))

temp_full <- temp_full[[keep_layers]]

dates_layer <- as.Date(time(temp_full))
dates_month_all <- dates_layer[seq(1, length(dates_layer), by = n_depths_full)]

cat("Layers after duplicate removal:", nlyr(temp_full), "\n")
cat("Duplicate months after cleaning:", sum(duplicated(dates_month_all)), "\n")

# ============================================================
# Keep only 1850–2099
# ============================================================

years_month_all <- as.integer(format(dates_month_all, "%Y"))

keep_month_blocks <- which(years_month_all >= 1850 & years_month_all <= 2099)

keep_layers <- unlist(lapply(keep_month_blocks, function(j) {
  ((j - 1) * n_depths_full + 1):(j * n_depths_full)
}))

temp_full <- temp_full[[keep_layers]]

dates_layer <- as.Date(time(temp_full))
dates_month_all <- dates_layer[seq(1, length(dates_layer), by = n_depths_full)]

cat("Layers after keeping 1850–2099:", nlyr(temp_full), "\n")
cat("Final monthly blocks:", length(dates_month_all), "\n")
cat("Expected monthly blocks:", length(years) * 12, "\n")
cat("Final start:", min(dates_month_all), "\n")
cat("Final end:", max(dates_month_all), "\n")

n_layers <- nlyr(temp_full)
n_months <- n_layers / n_depths_full
n_depths <- n_depths_full

if (n_months != floor(n_months)) {
  stop("Number of layers is not divisible by n_depths")
}

if (n_months != length(years) * 12) {
  warning(
    "Expected ", length(years) * 12,
    " months but found ", n_months
  )
}


# # ============================================================
# # Mean surface soil temperature map, SSP370, 2000–2014
# # Masked to total_thawed region
# # ============================================================
# 
# library(tidyterra)
# library(sf)
# library(rnaturalearth)
# library(rnaturalearthdata)
# 
# ssp <- "370"
# 
# temp_file <- file.path(
#   paste0("60N/mean_verttemp_ssp", ssp, "_60deg.nc")
# )
# 
# thawed_file <- file.path(
#   "total_thawed_extended",
#   paste0("arctic_total_thawed_", ssp, "_no_lim_mean.nc")
# )
# 
# temp <- rast(temp_file)
# thawed <- rast(thawed_file)
# 
# n_depths <- length(unique(as.numeric(depth(temp))))
# 
# n_depths <- length(unique(as.numeric(depth(temp))))
# 
# # surface layer = first depth layer
# surface_layers <- seq(3, nlyr(temp), by = n_depths)
# # surface layer = first depth layer
# surface_dates <- as.Date(time(temp))[surface_layers]
# surface_years <- as.integer(format(surface_dates, "%Y"))
# surface_months <- as.integer(format(surface_dates, "%m"))
# 
# summer_months <- c(6, 7, 8)
# 
# keep_surface_layers <- surface_layers[
#   surface_years >= 2000 &
#     surface_years <= 2014 
#     # surface_months %in% summer_months
# ]
# 
# surface_temp_2000_2014 <- mean(
#   temp[[keep_surface_layers]] - 273.15,
#   na.rm = TRUE
# )
# 
# names(surface_temp_2000_2014) <- "surface_temp_C"
# 
# # mask from total_thawed region
# thawed_mask <- !is.na(thawed[[nlyr(thawed)]])
# 
# thawed_mask_tempgrid <- resample(
#   thawed_mask,
#   surface_temp_2000_2014,
#   method = "near"
# )
# 
# surface_temp_2000_2014 <- mask(
#   surface_temp_2000_2014,
#   thawed_mask_tempgrid
# )
# 
# # Arctic polar projection
# crs_arctic <- "+proj=laea +lat_0=90 +lon_0=0 +datum=WGS84 +units=m +no_defs"
# 
# surface_temp_laea <- project(
#   surface_temp_2000_2014,
#   crs_arctic,
#   method = "bilinear"
# )
# 
# # crop plotting extent tightly to your masked area
# 
# surface_ext_laea <- ext(surface_temp_laea)
# 
# # optional: add small buffer around your region
# 
# buffer_m <- 300000
# 
# surface_ext_laea_buffered <- ext(
#   
#   xmin(surface_ext_laea) - buffer_m,
#   
#   xmax(surface_ext_laea) + buffer_m,
#   
#   ymin(surface_ext_laea) - buffer_m,
#   
#   ymax(surface_ext_laea) + buffer_m
#   
# )
# 
# surface_temp_laea_crop <- crop(surface_temp_laea, surface_ext_laea_buffered)
# 
# coast_laea_crop <- st_crop(
#   
#   coast_laea,
#   
#   xmin = xmin(surface_ext_laea_buffered),
#   
#   xmax = xmax(surface_ext_laea_buffered),
#   
#   ymin = ymin(surface_ext_laea_buffered),
#   
#   ymax = ymax(surface_ext_laea_buffered)
#   
# )
# coast <- ne_coastline(scale = "medium", returnclass = "sf")
# coast_laea <- st_transform(coast, crs_arctic)
# 
# p_surface_temp <- ggplot() +
#   geom_spatraster(data = surface_temp_laea_crop) +
#   geom_sf(
#     data = coast_laea_crop,
#     colour = "black",
#     linewidth = 0.25,
#     fill = NA
#   ) +
#   coord_sf(
#     xlim = c(xmin(surface_ext_laea_buffered), xmax(surface_ext_laea_buffered)),
#     ylim = c(ymin(surface_ext_laea_buffered), ymax(surface_ext_laea_buffered)),
#     expand = FALSE
#   ) +
#   scale_fill_gradientn(
#     colours = c(
#       "#08306B",
#       "#2171B5",
#       "#6BAED6",
#       "#C6DBEF",
#       "#FEE090",
#       "#FDAE61",
#       "#F46D43",
#       "#D73027"
#     ),
#     
#     na.value = "white",
#     
#     name = expression("Soil temp. ("*degree*C*")")
#     
#   )+
#   #theme_void(base_size = 15) +
#   labs(
#     title = "Mean JJA surface soil temperature, SSP370",
#     subtitle = "2000–2014 mean"
#   )
# 
# ggsave(
#   file.path(out_dir, "ssp370_surface_temperature_2000_2014_total_thawed_mask.png"),
#   p_surface_temp,
#   width = 10,
#   height = 10,
#   dpi = 300
# )

# ============================================================
# Monthly depth-weighted temperature
# ============================================================

# cat("Calculating monthly depth-weighted temperature...\n")
# 
# monthly_depth_weighted_C <- rep(NA_real_, n_months)
# temp_depth_monthly_list <- vector("list", n_months * n_depths)
# 
# counter <- 1L

# for (i in seq_len(n_months)) {
#   
#   depth_layers <- ((i - 1L) * n_depths + 1L):(i * n_depths)
#   
#   temp_month_C <- crop(
#     temp_full[[depth_layers]],
#     arctic_ext
#   ) - 273.15
#   
#   temp_month_C <- mask(
#     temp_month_C,
#     analysis_mask_tempgrid
#   )
#   
#   compareGeom(
#     temp_month_C[[1]],
#     area_rast,
#     stopOnError = TRUE
#   )
#   
#   # Require all depths to be present for depth-weighted temperature
#   complete_depth_mask <- app(
#     !is.na(temp_month_C),
#     sum
#   ) == n_depths
#   
#   # Initialize weighted sum
#   temp_depth_weighted <- temp_month_C[[1]] * 0
#   
#   for (d in seq_len(n_depths)) {
#     
#     temp_depth <- temp_month_C[[d]]
#     
#     temp_depth_weighted <-
#       temp_depth_weighted +
#       temp_depth * depth_weights[d]
#     
#     numerator_d <- global(
#       temp_depth * area_rast,
#       "sum",
#       na.rm = TRUE
#     )[1, 1]
#     
#     valid_area_d <- global(
#       ifel(!is.na(temp_depth), area_rast, NA),
#       "sum",
#       na.rm = TRUE
#     )[1, 1]
#     
#     temp_depth_monthly_list[[counter]] <- data.frame(
#       SSP_raw = ssp,
#       SSP = ssp_label,
#       Year = as.integer(format(dates_month_all[i], "%Y")),
#       Month = as.integer(format(dates_month_all[i], "%m")),
#       Depth = depth_vals[d],
#       Temperature_C = ifelse(
#         valid_area_d > 0,
#         numerator_d / valid_area_d,
#         NA_real_
#       )
#     )
#     
#     counter <- counter + 1L
#   }
#   
#   temp_depth_weighted <- ifel(
#     complete_depth_mask,
#     temp_depth_weighted,
#     NA
#   )
#   
#   numerator <- global(
#     temp_depth_weighted * area_rast,
#     "sum",
#     na.rm = TRUE
#   )[1, 1]
#   
#   valid_area <- global(
#     ifel(!is.na(temp_depth_weighted), area_rast, NA),
#     "sum",
#     na.rm = TRUE
#   )[1, 1]
#   
#   monthly_depth_weighted_C[i] <- ifelse(
#     valid_area > 0,
#     numerator / valid_area,
#     NA_real_
#   )
#   
#   if (i %% 100 == 0 || i == n_months) {
#     cat(
#       "Processed month", i, "of", n_months,
#       "| Date:", as.character(dates_month_all[i]), "\n"
#     )
#   }
#   
#   rm(
#     temp_month_C,
#     temp_depth,
#     temp_depth_weighted,
#     complete_depth_mask
#   )
#   
#   if (i %% 50 == 0) {
#     gc()
#   }
# }
# 
# temp_depth_monthly_df <- bind_rows(temp_depth_monthly_list)
# 
# write_csv(
#   temp_depth_monthly_df,
#   file.path(out_dir, paste0("ssp", ssp, "_depth_monthly_temperature.csv"))
# )

# ============================================================
# 1. Arctic mean soil temperature, depth-weighted
# ============================================================

# temp_monthly_mean_df <- data.frame(
#   SSP_raw = ssp,
#   SSP = ssp_label,
#   Year = as.integer(format(dates_month_all, "%Y")),
#   Month = as.integer(format(dates_month_all, "%m")),
#   Temperature_C = monthly_depth_weighted_C
# )
# 
# temp_annual_mean_df <- temp_monthly_mean_df %>%
#   group_by(SSP_raw, SSP, Year) %>%
#   summarise(
#     mean_temperature_C = mean(Temperature_C, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   arrange(Year) %>%
#   mutate(
#     mean_temperature_20yr_C = rollmean(
#       mean_temperature_C,
#       20,
#       fill = NA,
#       align = "center"
#     ),
#     ref_1880_1900_C = mean(
#       mean_temperature_C[Year >= 1880 & Year <= 1900],
#       na.rm = TRUE
#     ),
#     anomaly_1880_1900_C = mean_temperature_C - ref_1880_1900_C,
#     anomaly_1880_1900_20yr_C = rollmean(
#       anomaly_1880_1900_C,
#       20,
#       fill = NA,
#       align = "center"
#     )
#   )
# 
# write_csv(
#   temp_monthly_mean_df,
#   file.path(out_dir, paste0("ssp", ssp, "_depth_weighted_monthly_mean_temperature.csv"))
# )
# 
# write_csv(
#   temp_annual_mean_df,
#   file.path(out_dir, paste0("ssp", ssp, "_depth_weighted_annual_mean_temperature.csv"))
# )
# 
# p1 <- ggplot(
#   temp_annual_mean_df,
#   aes(x = Year, y = mean_temperature_C)
# ) +
#   geom_line(color = ssp_color, linewidth = 0.4, alpha = 0.5) +
#   geom_line(
#     aes(y = mean_temperature_20yr_C),
#     color = "black",
#     linewidth = 1.1,
#     na.rm = TRUE
#   ) +
#   geom_vline(xintercept = 2015, linetype = "dashed") +
#   theme_bw(base_size = 15) +
#   labs(
#     title = paste0("Depth-weighted Arctic soil temperature - ", ssp_label),
#     subtitle = "Colored line = annual mean, black line = 20-year running mean",
#     x = "Year",
#     y = expression("Soil temperature ("*degree*C*")")
#   )
# 
# ggsave(
#   file.path(out_dir, paste0("ssp", ssp, "_01_arctic_mean_soil_temperature.png")),
#   p1,
#   width = 8,
#   height = 5,
#   dpi = 300
# )

# ============================================================
# 2. mean pan-Arctic temperature at 0.85 m depth
# ============================================================

target_depth_m <- 1.0

depth_index_target <- which.min(abs(depth_vals - target_depth_m))

cat(
  "Requested depth:", target_depth_m, "m — using nearest available layer:",
  depth_vals[depth_index_target], "m (index", depth_index_target, ")\n"
)


monthly_target_depth_C <- rep(NA_real_, n_months)
temp_depth_monthly_list <- vector("list", n_months * n_depths)

counter <- 1L

for (i in seq_len(n_months)) {
  
  depth_layers <- ((i - 1L) * n_depths + 1L):(i * n_depths)
  
  temp_month_C <- crop(
    temp_full[[depth_layers]],
    arctic_ext
  ) - 273.15
  
  temp_month_C <- mask(
    temp_month_C,
    analysis_mask_tempgrid
  )
  
  compareGeom(
    temp_month_C[[1]],
    area_rast,
    stopOnError = TRUE
  )
  
  # Per-depth diagnostic CSV (unchanged) -- still records every depth
  for (d in seq_len(n_depths)) {
    
    temp_depth <- temp_month_C[[d]]
    
    numerator_d <- global(
      temp_depth * area_rast,
      "sum",
      na.rm = TRUE
    )[1, 1]
    
    valid_area_d <- global(
      ifel(!is.na(temp_depth), area_rast, NA),
      "sum",
      na.rm = TRUE
    )[1, 1]
    
    temp_depth_monthly_list[[counter]] <- data.frame(
      SSP_raw = ssp,
      SSP = ssp_label,
      Year = as.integer(format(dates_month_all[i], "%Y")),
      Month = as.integer(format(dates_month_all[i], "%m")),
      Depth = depth_vals[d],
      Temperature_C = ifelse(
        valid_area_d > 0,
        numerator_d / valid_area_d,
        NA_real_
      )
    )
    
    counter <- counter + 1L
  }
  
  # Pan-Arctic mean at the single target depth layer only
  temp_target_depth <- temp_month_C[[depth_index_target]]
  
  numerator <- global(
    temp_target_depth * area_rast,
    "sum",
    na.rm = TRUE
  )[1, 1]
  
  valid_area <- global(
    ifel(!is.na(temp_target_depth), area_rast, NA),
    "sum",
    na.rm = TRUE
  )[1, 1]
  
  monthly_target_depth_C[i] <- ifelse(
    valid_area > 0,
    numerator / valid_area,
    NA_real_
  )
  
  if (i %% 100 == 0 || i == n_months) {
    cat(
      "Processed month", i, "of", n_months,
      "| Date:", as.character(dates_month_all[i]), "\n"
    )
  }
  
  rm(
    temp_month_C,
    temp_depth,
    temp_target_depth
  )
  
  if (i %% 50 == 0) {
    gc()
  }
}

temp_monthly_mean_df <- data.frame(
  SSP_raw = ssp,
  SSP = ssp_label,
  Year = as.integer(format(dates_month_all, "%Y")),
  Month = as.integer(format(dates_month_all, "%m")),
  Depth_m = depth_vals[depth_index_target],
  Temperature_C = monthly_target_depth_C
)



write_csv(
  temp_monthly_mean_df,
  file.path(out_dir, paste0("ssp", ssp, "_1m_depth_monthly_mean_temperature.csv"))
)

write_csv(
  temp_annual_mean_df,
  file.path(out_dir, paste0("ssp", ssp, "_1m_depth_annual_mean_temperature.csv"))
)

p1 <- ggplot(
    temp_annual_mean_df,
    aes(x = Year, y = mean_temperature_C)
  ) +
    geom_line(color = ssp_color, linewidth = 0.4, alpha = 0.5) +
    geom_line(
      aes(y = mean_temperature_20yr_C),
      color = "black",
      linewidth = 1.1,
      na.rm = TRUE
    ) +
    geom_vline(xintercept = 2015, linetype = "dashed") +
    theme_bw(base_size = 15) +
    labs(
      title = paste0("Arctic soil temperature at ~1 m depth - ", ssp_label),
     subtitle = "Colored line = annual mean, black line = 20-year running mean",
     x = "Year",
     y = expression("Soil temperature ("*degree*C*")")
    )

ggsave(
  file.path(out_dir, paste0("ssp", ssp, "_01_arctic_mean_soil_temperature.png")),
  p1,
  width = 8,
  height = 5,
  dpi = 300
  )



monthly_cycle_periods <- temp_monthly_mean_df %>%
  mutate(
    Period = case_when(
      Year >= 1880 & Year <= 1900 ~ "1880-1900",
      Year >= 2000 & Year <= 2020 ~ "2000-2020",
      Year >= 2080 & Year <= 2099 ~ "2080-2099",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Period)) %>%
  group_by(SSP_raw, SSP, Period, Month) %>%
  summarise(
    mean_temperature_C = mean(Temperature_C, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Period = factor(
      Period,
      levels = c("1880-1900", "2000-2020", "2080-2099")
    )
  )

write_csv(
  monthly_cycle_periods,
  file.path(out_dir, paste0("ssp", ssp, "_monthly_temperature_cycle_periods.csv"))
)

p4 <- ggplot(
  monthly_cycle_periods,
  aes(x = Month, y = mean_temperature_C, color = Period)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  scale_x_continuous(
    breaks = 1:12,
    labels = month.abb
  ) +
  theme_bw(base_size = 15) +
  labs(
    title = paste0("Monthly soil temperature at ~1 m depth - ", ssp_label),
    x = "Month",
    y = expression("Soil temperature ("*degree*C*")"),
    color = "Period"
  )

ggsave(
  file.path(out_dir, paste0("ssp", ssp, "_04_monthly_temperature_cycle_periods.png")),
  p4,
  width = 8,
  height = 5,
  dpi = 300
)



# 
# # ============================================================
# # 2. Temperature profiles with depth
# # ============================================================
# 
# profile_periods <- temp_depth_monthly_df %>%
#   mutate(
#     Period = case_when(
#       Year >= 1880 & Year <= 1900 ~ "1880-1900",
#       Year >= 2000 & Year <= 2020 ~ "2000-2020",
#       Year >= 2080 & Year <= 2099 ~ "2080-2099",
#       TRUE ~ NA_character_
#     )
#   ) %>%
#   filter(!is.na(Period)) %>%
#   group_by(SSP_raw, SSP, Period, Depth) %>%
#   summarise(
#     mean_temperature_C = mean(Temperature_C, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   mutate(
#     Period = factor(
#       Period,
#       levels = c("1880-1900", "2000-2020", "2080-2099")
#     )
#   )
# 
# write_csv(
#   profile_periods,
#   file.path(out_dir, paste0("ssp", ssp, "_temperature_profiles_periods.csv"))
# )
# 
# p2 <- ggplot(
#   profile_periods,
#   aes(x = mean_temperature_C, y = Depth, color = Period)
# ) +
#   geom_line(linewidth = 1.1) +
#   geom_point(size = 2) +
#   scale_y_reverse() +
#   theme_bw(base_size = 15) +
#   labs(
#     title = paste0("Arctic mean soil temperature profiles - ", ssp_label),
#     x = expression("Soil temperature ("*degree*C*")"),
#     y = "Depth [m]",
#     color = "Period"
#   )
# 
# ggsave(
#   file.path(out_dir, paste0("ssp", ssp, "_02_temperature_profiles_depth.png")),
#   p2,
#   width = 8,
#   height = 6,
#   dpi = 300
# )
# 
# # ============================================================
# # 3. Thawed months per year
# # ============================================================
# 
# thawed_months_df <- temp_monthly_mean_df %>%
#   mutate(
#     thawed_month = Temperature_C > 0
#   ) %>%
#   group_by(SSP_raw, SSP, Year) %>%
#   summarise(
#     thawed_months_per_year = sum(thawed_month, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   arrange(Year) %>%
#   mutate(
#     thawed_months_20yr = rollmean(
#       thawed_months_per_year,
#       20,
#       fill = NA,
#       align = "center"
#     )
#   )
# 
# write_csv(
#   thawed_months_df,
#   file.path(out_dir, paste0("ssp", ssp, "_thawed_months_per_year.csv"))
# )
# 
# p3 <- ggplot(
#   thawed_months_df,
#   aes(x = Year, y = thawed_months_per_year)
# ) +
#   geom_line(color = ssp_color, linewidth = 0.4, alpha = 0.5) +
#   geom_line(
#     aes(y = thawed_months_20yr),
#     color = "black",
#     linewidth = 1.1,
#     na.rm = TRUE
#   ) +
#   geom_vline(xintercept = 2015, linetype = "dashed") +
#   theme_bw(base_size = 15) +
#   labs(
#     title = paste0("Thawed months per year - ", ssp_label),
#     subtitle = "Calculated from depth-weighted Arctic mean soil temperature > 0 °C",
#     x = "Year",
#     y = "Thawed months per year"
#   )
# 
# ggsave(
#   file.path(out_dir, paste0("ssp", ssp, "_03_thawed_months_per_year.png")),
#   p3,
#   width = 8,
#   height = 5,
#   dpi = 300
# )
# 
# # ============================================================
# # 4. Monthly temperature cycle
# # ============================================================
# 
# monthly_cycle_periods <- temp_monthly_mean_df %>%
#   mutate(
#     Period = case_when(
#       Year >= 1880 & Year <= 1900 ~ "1880-1900",
#       Year >= 2000 & Year <= 2020 ~ "2000-2020",
#       Year >= 2080 & Year <= 2099 ~ "2080-2099",
#       TRUE ~ NA_character_
#     )
#   ) %>%
#   filter(!is.na(Period)) %>%
#   group_by(SSP_raw, SSP, Period, Month) %>%
#   summarise(
#     mean_temperature_C = mean(Temperature_C, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   mutate(
#     Period = factor(
#       Period,
#       levels = c("1880-1900", "2000-2020", "2080-2099")
#     )
#   )
# 
# write_csv(
#   monthly_cycle_periods,
#   file.path(out_dir, paste0("ssp", ssp, "_monthly_temperature_cycle_periods.csv"))
# )
# 
# p4 <- ggplot(
#   monthly_cycle_periods,
#   aes(x = Month, y = mean_temperature_C, color = Period)
# ) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
#   geom_line(linewidth = 1.1) +
#   geom_point(size = 2) +
#   scale_x_continuous(
#     breaks = 1:12,
#     labels = month.abb
#   ) +
#   theme_bw(base_size = 15) +
#   labs(
#     title = paste0("Monthly depth-weighted soil temperature - ", ssp_label),
#     x = "Month",
#     y = expression("Soil temperature ("*degree*C*")"),
#     color = "Period"
#   )
# 
# ggsave(
#   file.path(out_dir, paste0("ssp", ssp, "_04_monthly_temperature_cycle_periods.png")),
#   p4,
#   width = 8,
#   height = 5,
#   dpi = 300
# )
# 
# # ============================================================
# # 5. Combined all-SSP annual plot from saved CSVs
# # ============================================================
# 
# cat("Checking whether all SSP annual CSVs exist for combined plot...\n")
# 
# ssps <- c("126", "245", "370", "585")
# 
# annual_csvs <- file.path(
#   out_dir,
#   paste0("ssp", ssps, "_depth_weighted_annual_mean_temperature.csv")
# )
# 
# names(annual_csvs) <- ssps
# 
# if (all(file.exists(annual_csvs))) {
#   
#   annual_all_df <- bind_rows(lapply(annual_csvs, read_csv, show_col_types = FALSE))
#   
#   p_all <- ggplot(
#     annual_all_df,
#     aes(x = Year, y = mean_temperature_C, color = SSP)
#   ) +
#     geom_line(linewidth = 0.5, alpha = 0.5) +
#     geom_line(
#       aes(y = mean_temperature_20yr_C),
#       linewidth = 1.1,
#       na.rm = TRUE
#     ) +
#     scale_color_manual(values = ssp_colors_label) +
#     theme_bw(base_size = 15) +
#     theme(legend.position = "bottom") +
#     labs(
#       title = "Depth-weighted Arctic soil temperature north of 60N",
#       subtitle = "Thin lines = annual means; thick lines = 20-year running means",
#       x = "Year",
#       y = expression("Soil temperature ("*degree*C*")"),
#       color = NULL
#     )
#   
#   ggsave(
#     file.path(out_dir, "all_ssps_depth_weighted_annual_soil_temperature_60N.png"),
#     p_all,
#     width = 8,
#     height = 5,
#     dpi = 300
#   )
#   
#   write_csv(
#     annual_all_df,
#     file.path(out_dir, "all_ssps_depth_weighted_annual_soil_temperature_60N.csv")
#   )
#   
# } else {
#   
#   cat("Not all annual CSVs exist yet. Skipping combined all-SSP plot.\n")
#   print(annual_csvs[!file.exists(annual_csvs)])
# }
# 
# cat("Finished SSP", ssp, "\n")
# 
# 
# library(dplyr)
# library(ggplot2)
# library(readr)
# library(zoo)
# 
# out_dir <- "monthly_mineralised/temperature_analysis_60N"
# 
# ssps <- c("126", "245", "370", "585")
# 
# ssp_colors <- c(
#   "SSP1-2.6" = "blue",
#   "SSP2-4.5" = "orange",
#   "SSP3-7.0" = "#D73027",
#   "SSP5-8.5" = "#7B3294"
# )
# 
# surface_df <- bind_rows(
#   lapply(ssps, function(ssp) {
#     
#     df <- read_csv(
#       file.path(
#         out_dir,
#         paste0("ssp", ssp, "_depth_monthly_temperature.csv")
#       ),
#       show_col_types = FALSE
#     )
#     
#     df %>%
#       filter(Depth == min(Depth)) %>%      # surface layer (0.01 m)
#       group_by(SSP, Year) %>%
#       summarise(
#         surface_temp_C = mean(Temperature_C),
#         .groups = "drop"
#       ) %>%
#       arrange(Year) %>%
#       mutate(
#         surface_temp_20yr = zoo::rollmean(
#           surface_temp_C,
#           20,
#           fill = NA,
#           align = "center"
#         )
#       )
#   })
# )
# 
# ggplot(surface_df, aes(x = Year, colour = SSP)) +
#   geom_line(
#     aes(y = surface_temp_C),
#     linewidth = 0.35,
#     alpha = 0.35
#   ) +
#   geom_line(
#     aes(y = surface_temp_20yr),
#     linewidth = 1.1,
#     na.rm = TRUE
#   ) +
#   scale_colour_manual(values = ssp_colors) +
#   theme_bw(base_size = 15) +
#   theme(
#     legend.position = "bottom"
#   ) +
#   labs(
#     x = "Year",
#     y = expression("Surface soil temperature (" * degree * "C)"),
#     colour = NULL
#   )






# ============================================================
# Analyse all SSPs: 1 m depth temperature
# - Annual pan-Arctic mean, 1850-2099
# - Monthly mean cycle for 1880-1900, 2000-2020, 2080-2099
# Reads the per-SSP CSVs written by the monthly-temperature script
# ============================================================

library(dplyr)
library(readr)
library(ggplot2)
library(zoo)

ssps <- c("126", "245", "370", "585")

ssp_labels <- c(
  "126" = "SSP1-2.6",
  "245" = "SSP2-4.5",
  "370" = "SSP3-7.0",
  "585" = "SSP5-8.5"
)

ssp_colors <- c(
  "SSP1-2.6" = "blue",
  "SSP2-4.5" = "orange",
  "SSP3-7.0" = "#D73027",
  "SSP5-8.5" = "#7B3294"
)

base_dir <- "monthly_mineralised"
out_dir <- file.path(base_dir, "temperature_analysis_60N")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Read and combine per-SSP monthly CSVs
# ------------------------------------------------------------

temp_monthly_list <- lapply(ssps, function(ssp) {
  
  in_file <- file.path(
    out_dir,
    paste0("ssp", ssp, "_1m_depth_monthly_mean_temperature.csv")
  )
  
  if (!file.exists(in_file)) {
    stop("Missing file: ", in_file)
  }
  
  df <- read_csv(in_file, show_col_types = FALSE)
  
  df$SSP_raw <- ssp
  df$SSP <- unname(ssp_labels[ssp])
  
  df
})

temp_monthly_all <- bind_rows(temp_monthly_list)

temp_monthly_all$SSP <- factor(
  temp_monthly_all$SSP,
  levels = unname(ssp_labels)
)

cat("Combined monthly rows:", nrow(temp_monthly_all), "\n")
print(head(temp_monthly_all))

# ------------------------------------------------------------
# 2. Annual pan-Arctic mean temperature, per SSP
# ------------------------------------------------------------

temp_annual_mean_all <- temp_monthly_all %>%
  group_by(SSP_raw, SSP, Year) %>%
  summarise(
    mean_temperature_C = mean(Temperature_C, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(SSP, Year) %>%
  group_by(SSP) %>%
  mutate(
    mean_temperature_20yr_C = rollmean(
      mean_temperature_C,
      20,
      fill = NA,
      align = "center"
    ),
    ref_2000_2020_C = mean(
      mean_temperature_C[Year >= 2000 & Year <= 2020],
      na.rm = TRUE
    ),
    anomaly_2000_2020_C = mean_temperature_C - ref_2000_2020_C,
    anomaly_2000_2020_20yr_C = rollmean(
      anomaly_2000_2020_C,
      20,
      fill = NA,
      align = "center"
    )
  ) %>%
  ungroup()

write_csv(
  temp_annual_mean_all,
  file.path(out_dir, "all_ssps_1m_depth_annual_mean_temperature.csv")
)

# Combined annual mean plot, all SSPs
p_annual_all <- ggplot(
  temp_annual_mean_all,
  aes(x = Year, y = anomaly_2000_2020_C, color = SSP)
) +
  geom_line(linewidth = 1.1, na.rm = TRUE) +
  geom_vline(xintercept = 2015, linetype = "dashed") +
  scale_color_manual(values = ssp_colors) +
  theme_bw(base_size = 15) +
  labs(
    title = "Annual soil temperature at 1 m depth",
    x = "Year",
    y = expression("Soil temperature ("*degree*C*")"),
    color = "Scenario"
  )

print(p_annual_all)

ggsave(
  file.path(out_dir, "all_ssps_1m_annual_mean_temperature.png"),
  p_annual_all,
  width = 9,
  height = 5.5,
  dpi = 300
)

# Optional: separate plot per SSP (matches earlier single-SSP style)
for (ssp in ssps) {
  
  df_ssp <- temp_annual_mean_all %>% filter(SSP_raw == ssp)
  ssp_label <- unname(ssp_labels[ssp])
  ssp_color <- unname(ssp_colors[ssp_label])
  
  p_ssp <- ggplot(
    df_ssp,
    aes(x = Year, y = mean_temperature_C)
  ) +
    geom_line(color = ssp_color, linewidth = 0.4, alpha = 0.5) +
    geom_line(
      aes(y = mean_temperature_20yr_C),
      color = "black",
      linewidth = 1.1,
      na.rm = TRUE
    ) +
    geom_vline(xintercept = 2015, linetype = "dashed") +
    theme_bw(base_size = 15) +
    labs(
      title = paste0("Annual pan-Arctic soil temperature at ~1 m depth - ", ssp_label),
      subtitle = "Colored line = annual mean, black line = 20-year running mean",
      x = "Year",
      y = expression("Soil temperature ("*degree*C*")")
    )
  
  ggsave(
    file.path(out_dir, paste0("ssp", ssp, "_1m_annual_mean_temperature.png")),
    p_ssp,
    width = 8,
    height = 5,
    dpi = 300
  )
}

# ------------------------------------------------------------
# 3. Monthly mean temperature cycle for selected periods, per SSP
# ------------------------------------------------------------

monthly_cycle_periods_all <- temp_monthly_all %>%
  mutate(
    Period = case_when(
      Year >= 1880 & Year <= 1900 ~ "1880-1900",
      Year >= 2000 & Year <= 2020 ~ "2000-2020",
      Year >= 2080 & Year <= 2099 ~ "2080-2099",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Period)) %>%
  group_by(SSP_raw, SSP, Period, Month) %>%
  summarise(
    mean_temperature_C = mean(Temperature_C, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Period = factor(
      Period,
      levels = c("1880-1900", "2000-2020", "2080-2099")
    )
  )

write_csv(
  monthly_cycle_periods_all,
  file.path(out_dir, "all_ssps_1m_monthly_temperature_cycle_periods.csv")
)

# Faceted plot: one panel per SSP, periods as colored lines
p_monthly_facet <- ggplot(
  monthly_cycle_periods_all,
  aes(x = Month, y = mean_temperature_C, color = Period)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  facet_wrap(~SSP, ncol = 2) +
  scale_x_continuous(
    breaks = 1:12,
    labels = month.abb
  ) +
  theme_bw(base_size = 13) +
  labs(
    title = "Monthly soil temperature at 1 m depth",
    x = "Month",
    y = expression("Soil temperature ("*degree*C*")"),
    color = "Period"
  )

print(p_monthly_facet)

ggsave(
  file.path(out_dir, "all_ssps_1m_monthly_temperature_cycle_periods_facets.png"),
  p_monthly_facet,
  width = 10,
  height = 8,
  dpi = 300
)

# Optional: separate plot per SSP (matches earlier single-SSP style)
for (ssp in ssps) {
  
  df_ssp <- monthly_cycle_periods_all %>% filter(SSP_raw == ssp)
  ssp_label <- unname(ssp_labels[ssp])
  
  p_ssp <- ggplot(
    df_ssp,
    aes(x = Month, y = mean_temperature_C, color = Period)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_line(linewidth = 1.1) +
    geom_point(size = 2) +
    scale_x_continuous(
      breaks = 1:12,
      labels = month.abb
    ) +
    theme_bw(base_size = 15) +
    labs(
      title = paste0("Monthly soil temperature at ~1 m depth - ", ssp_label),
      x = "Month",
      y = expression("Soil temperature ("*degree*C*")"),
      color = "Period"
    )
  
  ggsave(
    file.path(out_dir, paste0("ssp", ssp, "_1m_monthly_temperature_cycle_periods.png")),
    p_ssp,
    width = 8,
    height = 5,
    dpi = 300
  )
}

cat("Done. Outputs saved in:\n", out_dir, "\n")