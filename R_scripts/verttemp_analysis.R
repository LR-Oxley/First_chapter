# vertical soil temperature analysis: 
library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(zoo)
library(readr)

# --------------------------------------------------
# Settings
# --------------------------------------------------

base_dir <- "monthly_mineralised"
out_dir <- file.path(base_dir, "temperature_analysis_fastest")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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

years <- 1850:2099
n_depths <- 12

depth_vals <- c(
  0.01, 0.035, 0.075, 0.15,
  0.25, 0.4, 0.6, 0.85,
  1.25, 2.0, 3.0, 4.25
)

# --------------------------------------------------
# Depth weights
# --------------------------------------------------

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

# --------------------------------------------------
# Fastest SSP function
# --------------------------------------------------

process_temperature_ssp_fastest <- function(ssp) {
  
  temp_file <- file.path(
    base_dir,
    paste0("mean_verttemp_ssp", ssp, "_final.nc")
  )
  
  if (!file.exists(temp_file)) {
    warning("Missing file: ", temp_file)
    return(NULL)
  }
  
  cat("Reading:", temp_file, "\n")
  
  temp_full <- rast(temp_file)
  
  n_layers <- nlyr(temp_full)
  n_months <- n_layers / n_depths
  
  if (n_months != floor(n_months)) {
    stop("Layer number is not divisible by n_depths for SSP ", ssp)
  }
  
  # area raster
  area_rast <- cellSize(temp_full[[1]], unit = "m")
  
  cat("Calculating area-weighted mean for all layers...\n")
  
  # numerator: sum(temp * area)
  numerator <- global(
    temp_full * area_rast,
    "sum",
    na.rm = TRUE
  )[, 1]
  
  # denominator: valid area per layer
  valid_area <- global(
    (!is.na(temp_full)) * area_rast,
    "sum",
    na.rm = TRUE
  )[, 1]
  
  layer_mean_raw <- numerator / valid_area
  
  # Put layer means into matrix:
  # rows = depth, columns = month
  temp_matrix <- matrix(
    layer_mean_raw,
    nrow = n_depths,
    ncol = n_months
  )
  
  # depth-weighted monthly mean
  monthly_depth_weighted_raw <- as.numeric(
    colSums(temp_matrix * depth_weights)
  )
  
  data.frame(
    SSP_raw = ssp,
    SSP = ssp_labels[[ssp]],
    Year = rep(years, each = 12)[seq_len(n_months)],
    Month = rep(1:12, length.out = n_months),
    Temperature_raw = monthly_depth_weighted_raw
  )
}

# --------------------------------------------------
# Process all SSPs
# --------------------------------------------------

temperature_monthly_df <- bind_rows(
  lapply(ssps, process_temperature_ssp_fastest)
) %>%
  mutate(
    SSP = factor(
      SSP,
      levels = c("SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5")
    ),
    
    # Use this if raw temperature is Kelvin
    Temperature_C = Temperature_raw - 273.15
    
    # Use this instead if raw temperature is Fahrenheit
    # Temperature_C = (Temperature_raw - 32) * 5 / 9
  )

write_csv(
  temperature_monthly_df,
  file.path(out_dir, "all_ssp_depth_weighted_monthly_temperature.csv")
)

# --------------------------------------------------
# Annual mean
# --------------------------------------------------

temperature_annual_df <- temperature_monthly_df %>%
  group_by(SSP, Year) %>%
  summarise(
    mean_temperature_C = mean(Temperature_C, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(SSP) %>%
  arrange(Year) %>%
  mutate(
    mean_temperature_20yr_C = rollmean(
      mean_temperature_C,
      20,
      fill = NA,
      align = "center"
    ),
    ref_1880_1900_C = mean(
      mean_temperature_C[Year >= 1880 & Year <= 1900],
      na.rm = TRUE
    ),
    anomaly_1880_1900_C = mean_temperature_C - ref_1880_1900_C,
    anomaly_1880_1900_20yr_C = rollmean(
      anomaly_1880_1900_C,
      20,
      fill = NA,
      align = "center"
    )
  ) %>%
  ungroup()

write_csv(
  temperature_annual_df,
  file.path(out_dir, "all_ssp_depth_weighted_annual_temperature.csv")
)

# --------------------------------------------------
# Plot absolute temperature
# --------------------------------------------------

p_temp_abs <- ggplot(
  temperature_annual_df,
  aes(x = Year, y = mean_temperature_C, colour = SSP)
) +
  geom_line(linewidth = 0.4, alpha = 0.35) +
  geom_line(
    aes(y = mean_temperature_20yr_C),
    linewidth = 1.1,
    na.rm = TRUE
  ) +
  geom_vline(xintercept = 2015, linetype = "dashed") +
  scale_colour_manual(values = ssp_colors, drop = FALSE) +
  theme_bw(base_size = 15) +
  theme(legend.position = "bottom") +
  labs(
    x = "Year",
    y = expression("Depth-weighted soil temperature ("*degree*C*")"),
    colour = "SSP scenario"
  )

ggsave(
  file.path(out_dir, "all_ssp_depth_weighted_annual_temperature.png"),
  p_temp_abs,
  width = 8,
  height = 5,
  dpi = 300
)

print(p_temp_abs)

# --------------------------------------------------
# Plot anomaly
# --------------------------------------------------

p_temp_anom <- ggplot(
  temperature_annual_df,
  aes(x = Year, y = anomaly_1880_1900_C, colour = SSP)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_line(linewidth = 0.4, alpha = 0.35) +
  geom_line(
    aes(y = anomaly_1880_1900_20yr_C),
    linewidth = 1.1,
    na.rm = TRUE
  ) +
  geom_vline(xintercept = 2015, linetype = "dashed") +
  scale_colour_manual(values = ssp_colors, drop = FALSE) +
  theme_bw(base_size = 15) +
  theme(legend.position = "bottom") +
  labs(
    x = "Year",
    y = expression(Delta*" depth-weighted soil temperature ("*degree*C*")"),
    colour = "SSP scenario"
  )

ggsave(
  file.path(out_dir, "all_ssp_depth_weighted_temperature_anomaly_1880_1900.png"),
  p_temp_anom,
  width = 8,
  height = 5,
  dpi = 300
)

print(p_temp_anom)















temp_file <- "60N/mean_verttemp_ssp370_60deg.nc"

temp_full <- rast(temp_file)
plot(temp_full[[2]])
library(terra)
library(dplyr)
library(ggplot2)
library(tidyr)

# Example: temp_full is your monthly depth-resolved soil temperature raster
# units assumed: Kelvin
# layers ordered as: month1_depth1, month1_depth2, ..., month1_depth12, month2_depth1, ...

n_depths_full <- length(unique(as.numeric(depth(temp_full))))
depth_vals <- as.numeric(depth(temp_full))[1:n_depths_full]

temp_dates_raw <- as.Date(time(temp_full))

temp_dates_all <- temp_dates_raw[
  seq(1, length(temp_dates_raw), by = n_depths_full)
]

# Optional: restrict years
years_keep <- 1950:2099
temp_years <- as.integer(format(temp_dates_all, "%Y"))
month_keep <- which(temp_years %in% years_keep)

# Build table of monthly mean temperature by depth
temp_cycle_list <- list()

for (d in seq_len(n_depths_full)) {
  
  layer_idx <- ((month_keep - 1) * n_depths_full) + d
  
  temp_d <- temp_full[[layer_idx]]
  
  temp_mean_C <- global(temp_d, "mean", na.rm = TRUE)[, 1] - 273.15
  
  temp_cycle_list[[d]] <- data.frame(
    Date = temp_dates_all[month_keep],
    Year = as.integer(format(temp_dates_all[month_keep], "%Y")),
    Month = as.integer(format(temp_dates_all[month_keep], "%m")),
    Depth_m = depth_vals[d],
    Temp_C = temp_mean_C
  )
}

temp_cycle_df <- bind_rows(temp_cycle_list)

# Mean monthly cycle across selected years
temp_climatology <- temp_cycle_df %>%
  group_by(Month, Depth_m) %>%
  summarise(
    Temp_C = mean(Temp_C, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(temp_climatology,
       aes(x = Month, y = Temp_C, colour = factor(Depth_m))) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  scale_x_continuous(breaks = 1:12) +
  theme_bw() +
  labs(
    title = "Mean monthly soil temperature cycle",
    subtitle = paste0(min(years_keep), "–", max(years_keep)),
    x = "Month",
    y = "Soil temperature (°C)",
    colour = "Depth (m)"
  )



library(terra)
library(dplyr)
library(ggplot2)

# -----------------------------
# Settings
# -----------------------------

years_keep <- 1950:2099

n_depths_full <- length(unique(as.numeric(depth(temp_full))))
depth_vals <- as.numeric(depth(temp_full))[1:n_depths_full]

temp_dates_raw <- as.Date(time(temp_full))

temp_dates_all <- temp_dates_raw[
  seq(1, length(temp_dates_raw), by = n_depths_full)
]

temp_years <- as.integer(format(temp_dates_all, "%Y"))
month_keep <- which(temp_years %in% years_keep)

# -----------------------------
# Extract monthly Arctic mean temperature per depth
# -----------------------------

temp_monthly_list <- list()

for (d in seq_len(n_depths_full)) {
  
  layer_idx <- ((month_keep - 1) * n_depths_full) + d
  
  temp_d <- temp_full[[layer_idx]]
  
  temp_mean_C <- global(temp_d, "mean", na.rm = TRUE)[, 1] - 273.15
  
  temp_monthly_list[[d]] <- data.frame(
    Date = temp_dates_all[month_keep],
    Year = as.integer(format(temp_dates_all[month_keep], "%Y")),
    Month = as.integer(format(temp_dates_all[month_keep], "%m")),
    Depth_m = depth_vals[d],
    Temp_C = temp_mean_C
  )
}

temp_monthly_df <- bind_rows(temp_monthly_list)

# -----------------------------
# Plot full monthly time series
# -----------------------------

ggplot(temp_monthly_df,
       aes(x = Date, y = Temp_C, colour = factor(Depth_m))) +
  geom_line(linewidth = 0.35) +
  theme_bw() +
  labs(
    title = "Monthly Arctic mean soil temperature by depth",
    subtitle = "1950–2099",
    x = "Year",
    y = "Soil temperature (°C)",
    colour = "Depth (m)"
  )

selected_depths <- c(0.01, 0.15, 0.6, 1.15, 2.75, 4.25)

temp_monthly_df_sel <- temp_monthly_df %>%
  filter(Depth_m %in% selected_depths)

ggplot(temp_monthly_df_sel,
       aes(x = Date, y = Temp_C, colour = factor(Depth_m))) +
  geom_line(linewidth = 0.45) +
  theme_bw() +
  labs(
    title = "Monthly Arctic mean soil temperature by selected depths",
    subtitle = "1950–2099",
    x = "Year",
    y = "Soil temperature (°C)",
    colour = "Depth (m)"
  )





####################

# --------------------------------------------------
# Mean temperature time series by depth
# --------------------------------------------------

n_depths_full <- length(unique(as.numeric(depth(temp_full))))

depth_ids <- seq_len(n_depths_full)

temp_ts_list <- list()

for (d in depth_ids) {
  
  idx <- seq(d, nlyr(temp_full), by = n_depths_full)
  
  temp_depth <- temp_full[[idx]]
  
  temp_mean <- global(
    temp_depth,
    "mean",
    na.rm = TRUE
  )[,1] - 273.15
  
  temp_ts_list[[d]] <- data.frame(
    Date = as.Date(time(temp_depth)),
    Temp_C = temp_mean,
    Depth = round(as.numeric(depth(temp_depth)[1]), 3)
  )
}

temp_ts_df <- bind_rows(temp_ts_list)

plot_depths <- c(1, 4, 8, 12)

ggplot(
  temp_ts_df %>%
    filter(Depth %in% unique(temp_ts_df$Depth)[plot_depths]),
  aes(Date, Temp_C, colour = factor(Depth))
) +
  geom_line(linewidth = 0.4) +
  theme_bw() +
  geom_vline(
    xintercept = as.Date("2015-01-01"),
    linetype = "dashed"
  ) +
  labs(
    colour = "Depth (m)",
    y = "Temperature (°C)"
  )

dates_all <- as.Date(time(temp_full))

length(dates_all)
length(unique(dates_all))

surface_idx <- seq(1, nlyr(temp_full), by = n_depths_full)

surface_temp <- temp_full[[surface_idx]]

df <- data.frame(
  Date = as.Date(time(surface_temp)),
  Temp = global(surface_temp, "mean", na.rm = TRUE)[,1] - 273.15
)

ggplot(df, aes(Date, Temp)) +
  geom_line() +
  theme_bw() +
  geom_vline(
    xintercept = as.Date("2015-01-01"),
    colour = "red"
  )


head(time(temp_full), 30)
unique(depth(temp_full))[1:12]
data.frame(
  
  layer = 1:30,
  
  date = as.Date(time(temp_full))[1:30],
  
  depth = as.numeric(depth(temp_full))[1:30]
  
)



n_depths_full <- length(unique(as.numeric(depth(temp_full))))

dates_layer <- as.Date(time(temp_full))

dates_month <- dates_layer[seq(1, length(dates_layer), by = n_depths_full)]

length(dates_month)

length(unique(dates_month))

sum(duplicated(dates_month))

range(dates_month)


surface_idx <- seq(1, nlyr(temp_full), by = n_depths_full)

surface_temp <- temp_full[[surface_idx]]

df_surface <- data.frame(
  
  Date = dates_month,
  
  Temp_C = global(surface_temp, "mean", na.rm = TRUE)[, 1] - 273.15
  
)

ggplot(df_surface, aes(Date, Temp_C)) +
  
  geom_line(linewidth = 0.3) +
  
  geom_vline(xintercept = as.Date("2015-01-01"), colour = "red") +
  
  theme_bw() +
  
  labs(
    
    title = "Mean surface soil temperature",
    
    subtitle = "Shallowest depth layer",
    
    x = "Year",
    
    y = "Temperature (°C)"
    
  )





####### check if k_t is correct
# ============================================================
# Diagnostic: check kT depth values for one month
# ============================================================

dup_dates <- temp_month_dates[duplicated(temp_month_dates)]

print(length(dup_dates))

print(head(dup_dates, 20))

print(tail(dup_dates, 20))

temp_dates_raw <- as.Date(time(temp_full))

temp_dates_all <- temp_dates_raw[seq(1, length(temp_dates_raw), by = n_depths_full)]

sum(duplicated(temp_dates_all))

range(temp_dates_all)

check_kT_depths <- function(check_year,
                            check_month,
                            temp_full,
                            depth_vals,
                            n_depths_full,
                            depth_keep,
                            area_rast,
                            land_mask,
                            Ea,
                            Ed,
                            t_opt_C,
                            k_T_ref_scalar,
                            out_dir,
                            ssp) {
  
  dates_layer <- as.Date(time(temp_full))
  
  dates_month_all <- dates_layer[seq(1, length(dates_layer), by = n_depths_full)]
  
  month_block_all <- seq_along(dates_month_all)
  
  month_table <- data.frame(
    block = month_block_all,
    date = dates_month_all,
    year = as.integer(format(dates_month_all, "%Y")),
    month = as.integer(format(dates_month_all, "%m"))
  )
  
  cat("Monthly blocks:", nrow(month_table), "\n")
  cat("Unique monthly dates:", length(unique(month_table$date)), "\n")
  cat("Duplicated monthly dates:", sum(duplicated(month_table$date)), "\n")
  
  # Keep first occurrence of each date
  month_table_unique <- month_table[!duplicated(month_table$date), ]
  
  j_row <- which(
    month_table_unique$year == check_year &
      month_table_unique$month == check_month
  )
  
  if (length(j_row) != 1) {
    print(month_table_unique[month_table_unique$year == check_year, ])
    stop("Could not find exactly one unique monthly date.")
  }
  
  j <- month_table_unique$block[j_row]
  
  cat("Using monthly block:", j, "\n")
  cat("Using date:", as.character(month_table_unique$date[j_row]), "\n")
  
  depth_layers <- ((j - 1) * n_depths_full) + depth_keep
  
  cat("Checking kT for:", check_year, "-", sprintf("%02d", check_month), "\n")
  cat("Monthly index:", j, "\n")
  cat("Date:", as.character(dates_month[j]), "\n")
  cat("Selected raster layers:\n")
  print(depth_layers)
  
  temp_K <- temp_full[[depth_layers]]
  temp_K <- terra::mask(temp_K, land_mask)
  temp_C <- temp_K - 273.15
  
  kT_raw <- peaked_arrhenius(
    temp_C,
    Ea = Ea,
    Ed = Ed,
    t_opt_C = t_opt_C
  )
  
  kT_raw <- terra::ifel(temp_C < 0, 0, kT_raw)
  kT_raw <- terra::ifel(is.na(temp_C), NA, kT_raw)
  
  kT_factor <- kT_raw / k_T_ref_scalar
  
  # Temperature profile
  temp_mean <- terra::global(
    temp_C,
    "mean",
    weights = area_rast,
    na.rm = TRUE
  )[, 1]
  
  temp_min <- terra::global(
    temp_C,
    "min",
    na.rm = TRUE
  )[, 1]
  
  temp_max <- terra::global(
    temp_C,
    "max",
    na.rm = TRUE
  )[, 1]
  
  # Raw kT profile
  kT_raw_mean <- terra::global(
    kT_raw,
    "mean",
    weights = area_rast,
    na.rm = TRUE
  )[, 1]
  
  kT_raw_min <- terra::global(
    kT_raw,
    "min",
    na.rm = TRUE
  )[, 1]
  
  kT_raw_max <- terra::global(
    kT_raw,
    "max",
    na.rm = TRUE
  )[, 1]
  
  # Normalised kT profile
  kT_factor_mean <- terra::global(
    kT_factor,
    "mean",
    weights = area_rast,
    na.rm = TRUE
  )[, 1]
  
  kT_factor_min <- terra::global(
    kT_factor,
    "min",
    na.rm = TRUE
  )[, 1]
  
  kT_factor_max <- terra::global(
    kT_factor,
    "max",
    na.rm = TRUE
  )[, 1]
  
  check_df <- data.frame(
    depth_m = depth_vals,
    temp_min_C = temp_min,
    temp_mean_C = temp_mean,
    temp_max_C = temp_max,
    kT_raw_min = kT_raw_min,
    kT_raw_mean = kT_raw_mean,
    kT_raw_max = kT_raw_max,
    kT_factor_min = kT_factor_min,
    kT_factor_mean = kT_factor_mean,
    kT_factor_max = kT_factor_max
  )
  
  print(check_df)
  
  write.csv(
    check_df,
    file.path(
      out_dir,
      paste0(
        "check_kT_depths_",
        ssp, "_",
        check_year, "_",
        sprintf("%02d", check_month),
        ".csv"
      )
    ),
    row.names = FALSE
  )
  
  p_temp <- ggplot(check_df, aes(x = temp_mean_C, y = depth_m)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_y_reverse() +
    theme_bw() +
    labs(
      title = paste0("Mean soil temperature profile: ", check_year, "-", sprintf("%02d", check_month)),
      x = "Temperature (°C)",
      y = "Depth (m)"
    )
  
  ggsave(
    file.path(
      out_dir,
      paste0("check_temperature_profile_", ssp, "_", check_year, "_", sprintf("%02d", check_month), ".png")
    ),
    p_temp,
    width = 6,
    height = 6,
    dpi = 300
  )
  
  p_kT <- ggplot(check_df, aes(x = kT_factor_mean, y = depth_m)) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_y_reverse() +
    theme_bw() +
    labs(
      title = paste0("Mean normalised kT profile: ", check_year, "-", sprintf("%02d", check_month)),
      x = "kT / kT_ref",
      y = "Depth (m)"
    )
  
  ggsave(
    file.path(
      out_dir,
      paste0("check_kT_profile_", ssp, "_", check_year, "_", sprintf("%02d", check_month), ".png")
    ),
    p_kT,
    width = 6,
    height = 6,
    dpi = 300
  )
  
  invisible(check_df)
}

check_kT_depths(
  check_year = 2010,
  check_month = 7,
  temp_full = temp_full,
  depth_vals = depth_vals,
  n_depths_full = n_depths_full,
  depth_keep = depth_keep,
  area_rast = area_rast,
  land_mask = land_mask,
  Ea = Ea_nitrification,
  Ed = Ed_nitrification,
  t_opt_C = t_opt_C,
  k_T_ref_scalar = k_T_ref_scalar,
  out_dir = out_dir,
  ssp = ssp
)


check_kT_depths(
  
  check_year = 2010,
  
  check_month = 1,
  
  temp_full = temp_full,
  
  depth_vals = depth_vals,
  
  n_depths_full = n_depths_full,
  
  depth_keep = depth_keep,
  
  area_rast = area_rast,
  
  land_mask = land_mask,
  
  Ea = Ea_nitrification,
  
  Ed = Ed_nitrification,
  
  t_opt_C = t_opt_C,
  
  k_T_ref_scalar = k_T_ref_scalar,
  
  out_dir = out_dir,
  
  ssp = ssp
  
)

temp_mean_C

kT_raw_mean

kT_factor_mean




# test for soil temperature
# # choose depth layer, e.g. shallowest depth
depth_id <- 2
n_depths <- 12

# extract only this depth for every month
idx <- seq(depth_id, nlyr(temp_full), by = n_depths)
temp_depth <- temp_full[[idx]]

# spatial mean per month
temp_mean <- global(temp_depth, "mean", na.rm = TRUE)[,1]

df_temp <- data.frame(
  Date = as.Date(time(temp_depth)),
  Temp_K = temp_mean,
  Temp_C = temp_mean - 273.15
)

# plot full time series
ggplot(df_temp, aes(Date, Temp_C)) +
  geom_line(linewidth = 0.4) +
  theme_bw() +
  labs(
    title = "Monthly mean soil temperature",
    subtitle = paste("Depth:", depth(temp_full)[depth_id], "m"),
    y = "Temperature (°C)",
    x = "Year"
  )

df_zoom <- df_temp %>%

  filter(Date >= as.Date("2005-01-01"),

         Date <= as.Date("2025-12-31"))

ggplot(df_zoom, aes(Date, Temp_C)) +

  geom_line() +

  theme_bw()




years <- as.numeric(format(time(temp_full)[seq(1, nlyr(temp_full), by = 12)], "%Y"))

temp_ts <- extract(

  temp_full[[seq(1, nlyr(temp_full), by = 12)]],

  cell

)[1,-1]

plot(years, temp_ts - 273.15, type="l")

abline(v = 2015, col="red")


length(time(temp_full))

nlyr(temp_full)

length(unique(depth(temp_full)))

depth1 <- temp_full[[seq(12, nlyr(temp_full), by = 12)]]

df <- data.frame(
  Date = as.Date(time(depth1)),
  Temp = global(depth1, "mean", na.rm = TRUE)[,1] - 273.15
)

df %>%
  filter(Date >= as.Date("2010-01-01"),
         Date <= as.Date("2020-12-31")) %>%
  print(n = 30)


library(dplyr)

df_year <- df %>%

  mutate(Year = as.integer(format(Date, "%Y"))) %>%

  group_by(Year) %>%

  summarise(T = mean(Temp))

plot(df_year$Year, df_year$T, type="l", )


abline(v=2015, col="red")

