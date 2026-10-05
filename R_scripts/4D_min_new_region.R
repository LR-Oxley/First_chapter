# ============================================================
# Monthly mineralised and bioavailable N
# Years: 1850–2099
# Dynamic depth-pool version, with soil moisture
#
# Main logic:
# Thawed permafrost N is calculated relative to the
# 1850–1900 mean baseline at each grid cell.
#
# baseline-relative thawed N =
# total thawed N in current year - mean total thawed N from 1850–1900
##
# Newly thawed N is calculated as the positive year-to-year
# increase in this baseline-relative thawed N pool.
#
# This avoids adding the initial active-layer N pool in 1850
# and prevents an artificial bioavailable N peak at the start.
#
# Newly thawed organic N is distributed vertically according to
# the current ALD increment/depth profile and then added to
# dynamic organic depth pools.
#
# Organic N is mineralised monthly using:
# - temperature modifier k_T
# - soil moisture modifier k_moisture
# - base mineralisation rate
#
# Bioavailable N =
# monthly mineralised organic N
# + rapid inorganic fraction of newly thawed N
#
# ============================================================
terra::gdal(drivers = TRUE)

library(terra)
library(dplyr)
library(ggplot2)


# -----------------------------
# 0. Setup
# -----------------------------

args <- commandArgs(trailingOnly = TRUE)
ssp <- args[1]

if (is.na(ssp) || length(ssp) == 0) {
  ssp <- "370"
}

out_dir <- "monthly_mineralised/whole_region_mean_12deg"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_dir_yearly <- file.path(out_dir, paste0("_w_temp_w_sm_yearly_nc_", ssp))
dir.create(out_dir_yearly, recursive = TRUE, showWarnings = FALSE)


log_dir <- "logs"
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(
  log_dir,
  paste0(
    "region_mineralisation_w_temp_w_sm_mean_",
    ssp,
    "_",
    format(Sys.time(), "%Y%m%d_%H%M%S"),
    ".log"
  )
)

log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, append = TRUE, type = "message")

start_time <- Sys.time()

cat("====================================================\n")
cat("STARTED SCRIPT\n")
cat("Time:", as.character(start_time), "\n")
cat("SSP:", ssp, "\n")
cat("====================================================\n")

on.exit({
  end_time <- Sys.time()
  cat("====================================================\n")
  cat("SCRIPT FINISHED\n")
  cat("End:", as.character(end_time), "\n")
  cat("Runtime:", as.character(end_time - start_time), "\n")
  cat("Log file:", log_file, "\n")
  cat("====================================================\n")
  sink(type = "message")
  sink()
  close(log_con)
}, add = TRUE)

# -----------------------------
# 1. Parameters
# -----------------------------

k_base <- 0.002

Ea_nitrification <- 80000
Ed_nitrification <- 200000
t_opt_C <- 28

f_inorg_rapid <- 0.1136
f_org <- 1 - f_inorg_rapid

average_arctic_thaw_length <- 0.3

#test_years <- 1850:2014
test_years <- 1850:2099
start_year <- min(test_years)

baseline_start <- 1850
baseline_end <- 1900

ref_start <- 2000
ref_end <- 2014

profile_depth_max <- 5

common_extent <- ext(-179.95, 179.95, 45, 90)

# -----------------------------
# 2. File paths
# -----------------------------

ald_file <- paste0("mean_ssp", ssp, "_fillnans.nc")

thawed_file <- paste0(
  "total_thawed_extended/arctic_total_thawed_",
  ssp,
  "_no_lim_mean.nc"
)

temp_file <- paste0(
  "monthly_mineralised/mean_verttemp_ssp",
  ssp,
  "_final.nc"
)

# # soil moisture data
theta_file <- paste0("soil_moisture/sm_", ssp, "_present.nc")



# -----------------------------
# 3. Load data
# -----------------------------

cat("Loading rasters...\n")

ALD <- rast(ald_file)
ALD <- crop(ALD, common_extent)

thawed_total <- rast(thawed_file)
thawed_total <- crop(thawed_total, common_extent)
#plot(thawed_total[[2]])
land_mask <- !is.na(thawed_total[[nlyr(thawed_total)]])
plot(land_mask)
temp_full <- rast(temp_file)
crs(temp_full) <- crs(thawed_total)
#plot(temp_full[[2]])
#temp_full<-mask(temp_full, land_mask)

theta <- rast(theta_file)
theta <- crop(theta, common_extent)
theta <- resample(theta, thawed_total, method = "bilinear")
crs(theta) <- crs(thawed_total)
#plot(theta[[2]])

LC <- rast("LC_remapnn_corr.nc")
LC <- crop(LC, common_extent)
LC <- resample(LC, thawed_total[[1]], method = "near")

cat("Initial raster info:\n")
print(ALD)
print(thawed_total)
print(temp_full)
print(LC)

cat("Geometry check temp vs thawed:\n")
print(compareGeom(LC[[1]], thawed_total[[1]], stopOnError = FALSE))
cat("Geometry check temp_full vs thawed:\n")
print(compareGeom(temp_full[[1]], thawed_total[[1]], stopOnError = FALSE))
cat("Geometry check ALD vs thawed:\n")
print(compareGeom(ALD[[1]], thawed_total[[1]], stopOnError = FALSE))
cat("Geometry check theta vs thawed:\n")
print(compareGeom(theta[[1]], thawed_total[[1]], stopOnError = FALSE))

####################################################################################
# temperature: make sure there are no duplicate years
####################################################################################
n_depths_full <- length(unique(as.numeric(depth(temp_full))))
dates_layer <- as.Date(time(temp_full))
dates_month_all <- dates_layer[seq(1, length(dates_layer), by = n_depths_full)]

cat("Monthly blocks before cleaning:", length(dates_month_all), "\n")
cat("Unique monthly dates:", length(unique(dates_month_all)), "\n")
cat("Duplicates:", sum(duplicated(dates_month_all)), "\n")

keep_month_blocks <- which(!duplicated(dates_month_all))

keep_layers <- unlist(lapply(keep_month_blocks, function(j) {
  ((j - 1) * n_depths_full + 1):(j * n_depths_full)
}))

temp_full <- temp_full[[keep_layers]]

cat("Layers after removing duplicate months:", nlyr(temp_full), "\n")

dates_layer <- as.Date(time(temp_full))

dates_month_all <- dates_layer[seq(1, length(dates_layer), by = n_depths_full)]

sum(duplicated(dates_month_all))

####################################################################################

# ============================================================
# Optional test mode
# ============================================================

test_mode <- TRUE

if (test_mode) {
  
  cat("RUNNING SMALL REGION TEST MODE\n")
  
  test_ext <- ext(90, 90.5, 60, 60.5)
  
  ALD <- crop(ALD, test_ext)
  thawed_total <- crop(thawed_total, test_ext)
  temp_full <- crop(temp_full, test_ext)
  LC <- crop(LC, test_ext)
  
  cat("Small-region raster checks:\n")
  print(ALD)
  print(thawed_total)
  print(temp_full)
  print(LC)
}

# -----------------------------
# 4. Subset thawed_total and ALD to 1850–2014
# -----------------------------

thawed_years <- as.numeric(format(time(thawed_total), "%Y"))

years_to_keep <- which(
  thawed_years >= min(test_years) &
    thawed_years <= max(test_years)
)

if (length(years_to_keep) == 0) {
  stop("No thawed_total layers found for 1850–2014.")
}

thawed_total <- thawed_total[[years_to_keep]]

ALD_years <- as.numeric(format(time(ALD), "%Y"))

ALD_idx <- which(
  ALD_years >= min(test_years) &
    ALD_years <= max(test_years)
)

if (length(ALD_idx) == 0) {
  stop("No ALD layers found for 1850–2014.")
}

ALD <- ALD[[ALD_idx]]

if (nlyr(ALD) != nlyr(thawed_total)) {
  stop("ALD and thawed_total do not have the same number of yearly layers.")
}

cat("Subset thawed_total to", nlyr(thawed_total), "layers.\n")
cat("Subset ALD to", nlyr(ALD), "layers.\n")

# ============================================================
# 4b. Baseline correction of thawed permafrost N
# ============================================================

cat("Correcting thawed N relative to 1850–1900 mean baseline...\n")

thawed_years <- as.numeric(format(time(thawed_total), "%Y"))

baseline_idx <- which(
  thawed_years >= baseline_start &
    thawed_years <= baseline_end
)

if (length(baseline_idx) == 0) {
  stop("No baseline years found in thawed_total.")
}

land_mask <- !is.na(thawed_total[[nlyr(thawed_total)]])

baseline_thawed_N <- app(
  terra::mask(thawed_total[[baseline_idx]], land_mask),
  mean,
  na.rm = TRUE
)


baseline_thawed_N <- terra::mask(baseline_thawed_N, land_mask)

thawed_permafrost_N <- thawed_total - baseline_thawed_N

thawed_permafrost_N <- terra::mask(
  thawed_permafrost_N,
  land_mask
)

#plot(baseline_thawed_N)

cat("Baseline thawed N range:\n")
print(global(baseline_thawed_N, range, na.rm = TRUE))

cat("Corrected thawed permafrost N range:\n")
print(global(thawed_permafrost_N, range, na.rm = TRUE))

area_rast <- cellSize(thawed_permafrost_N[[1]], unit="m")

thawed_permafrost_pg <- terra::global(
  thawed_permafrost_N * area_rast,
  "sum",
  na.rm = TRUE
)[, 1] / 1e12

thawed_permafrost_df <- data.frame(
  Year = as.integer(format(time(thawed_total), "%Y")),
  thawed_permafrost_Pg = thawed_permafrost_pg
)

p_thawed_permafrost <- ggplot(
  thawed_permafrost_df,
  aes(x = Year, y = thawed_permafrost_Pg)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.2) +
  theme_bw() +
  labs(
    title = "Baseline-corrected thawed permafrost N",
    subtitle = "Relative to 1850–1900 mean thawed N",
    x = "Year",
    y = "Thawed permafrost N (Pg N)"
  )

print(p_thawed_permafrost)

library(dplyr)
library(ggplot2)

thawed_permafrost_df <- thawed_permafrost_df %>%
  arrange(Year) %>%
  mutate(
    annual_change_Pg = thawed_permafrost_Pg - lag(thawed_permafrost_Pg),
    annual_positive_change_Pg = pmax(annual_change_Pg, 0)
  )

p_annual_change <- ggplot(
  thawed_permafrost_df,
  aes(x = Year, y = annual_change_Pg)
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_col() +
  theme_bw() +
  labs(
    title = "Annual change in baseline-corrected thawed permafrost N",
    subtitle = "Year-to-year difference in thawed permafrost N pool",
    x = "Year",
    y = expression("Annual change (Pg N yr"^{-1}*")")
  )

print(p_annual_change)

p_positive_change <- ggplot(
  thawed_permafrost_df,
  aes(x = Year, y = annual_positive_change_Pg)
) +
  geom_col() +
  theme_bw() +
  labs(
    title = "Annual newly thawed permafrost N increment",
    subtitle = "Only positive year-to-year increases",
    x = "Year",
    y = expression("Positive annual increment (Pg N yr"^{-1}*")")
  )

print(p_positive_change)
# -----------------------------
# 5. Temperature temporal/depth subset
# -----------------------------

n_depths_full <- length(unique(as.numeric(depth(temp_full))))
all_depth_vals <- as.numeric(depth(temp_full))[1:n_depths_full]

depth_keep <- which(all_depth_vals <= profile_depth_max)

depth_vals <- all_depth_vals[depth_keep]
n_depths <- length(depth_vals)

temp_dates_raw <- as.Date(time(temp_full))
plot(temp_full[[2]])
temp_dates_all <- temp_dates_raw[
  seq(1, length(temp_dates_raw), by = n_depths_full)
]

temp_years_all <- as.integer(format(temp_dates_all, "%Y"))

temp_idx_months <- which(temp_years_all %in% test_years)

temp_month_dates <- temp_dates_all[temp_idx_months]
n_months <- length(temp_month_dates)

cat("Using", n_depths, "depth layers and", n_months, "months.\n")
print(depth_vals)

#plot(temp_full[[2]])

############## check



# check_temp_month <- function(year_check, month_check) {
#   
#   j <- which(
#     as.integer(format(temp_month_dates, "%Y")) == year_check &
#       as.integer(format(temp_month_dates, "%m")) == month_check
#   )
#   
#   depth_layers <- get_temp_layers_from_full(
#     month_index_selected = j,
#     temp_idx_months = temp_idx_months,
#     n_depths_full = n_depths_full,
#     depth_keep = depth_keep
#   )
#   
#   temp_C <- temp_full[[depth_layers]] - 273.15
#   names(temp_C) <- paste0(depth_vals, "m")
#   
#   print(global(temp_C, c("mean", "min", "max"), na.rm = TRUE))
#   
#   invisible(temp_C)
# }
# 
# jan2000 <- check_temp_month(2000, 1)
# feb2000 <- check_temp_month(2000, 2)
# dec2000 <- check_temp_month(2000, 12)
# 
# check_unfrozen_fraction <- function(temp_C) {
#   
#   unfrozen <- temp_C > 0
#   
#   global(unfrozen, "mean", na.rm = TRUE)
# }
# 
# check_unfrozen_fraction(jan2000)
# check_unfrozen_fraction(feb2000)
# check_unfrozen_fraction(dec2000)
# 
# plot(global(jan2000, "mean", na.rm = TRUE)[,1],
#      type = "b",
#      xaxt = "n",
#      xlab = "Depth layer",
#      ylab = "Mean soil temperature (°C)",
#      main = "January 2000 mean temperature by depth")
# 
# axis(1, at = seq_along(depth_vals), labels = depth_vals)
# abline(h = 0, lty = 2)
# 
# 
# unfrozen_fraction <- global(temp_full > 273.15, "mean", na.rm = TRUE)
# data.frame(
#   depth_m = depth_vals,
#   unfrozen_fraction = unfrozen_fraction[,1]
# )
# 
# plot(temp_C[[12]], main = "January 2000 temperature at 4.25 m")
# plot(temp_C[[12]] > 0, main = "Unfrozen cells at 4.25 m, January 2000")
#   
# xy <- c(100, 67.5)
# 
# terra::extract(
#   
#   temp_full[[layers]] - 273.15,
#   
#   matrix(xy, ncol = 2)
#   
# )
# 
# 
# 
# plot(temp_full[[layers[1]]] - 273.15,
#      
#      main="Surface January 2000")
# 
# plot(temp_full[[layers[1] + 6]] - 273.15,
#      
#      main="Surface July 2000")
# 
# 
# 
# year_check <- 2000
# month_check <- 1
# 
# j <- which(
#   as.integer(format(temp_month_dates, "%Y")) == year_check &
#     as.integer(format(temp_month_dates, "%m")) == month_check
# )
# 
# depth_layers <- get_temp_layers_from_full(
#   month_index_selected = j,
#   temp_idx_months = temp_idx_months,
#   n_depths_full = n_depths_full,
#   depth_keep = depth_keep
# )
# 
# temp_jan2000 <- temp_full[[depth_layers]] - 273.15
# 
# temp_jan2000_60N <- crop(
#   temp_jan2000,
#   ext(xmin(temp_jan2000), xmax(temp_jan2000), 60, ymax(temp_jan2000))
# )
# 
# sapply(1:nlyr(temp_jan2000_60N), function(i) {
#   global(temp_jan2000_60N[[i]], "mean", na.rm = TRUE)[1, 1]
# })
# 
# global(
#   temp_jan2000_60N,
#   c("mean", "min", "max"),
#   na.rm = TRUE
# )
# 
# 
# plot(temp_jan2000_60N[[2]])
# 
# global(temp_jan2000_60N > 0, "mean", na.rm = TRUE)
# 


# ============================================================
# 6. Helper functions
# ============================================================

peaked_arrhenius <- function(temp_C, Ea, Ed, t_opt_C = 28) {
  R_gas <- 8.314462618
  
  temp_K <- temp_C + 273.15
  t_opt_K <- t_opt_C + 273.15
  
  hlp1 <- temp_K - t_opt_K
  hlp2 <- temp_K * t_opt_K * R_gas
  
  numerator <- Ed * exp(Ea * hlp1 / hlp2)
  denominator <- Ed - Ea * (1 - exp(Ed * hlp1 / hlp2))
  
  numerator / denominator
}

get_temp_layers_from_full <- function(month_index_selected,
                                      temp_idx_months,
                                      n_depths_full,
                                      depth_keep) {
  
  month_index_full <- temp_idx_months[month_index_selected]
  
  ((month_index_full - 1) * n_depths_full) + depth_keep
}

get_depth_bounds <- function(depth_mid, max_depth) {
  depth_bounds <- numeric(length(depth_mid) + 1)
  depth_bounds[1] <- 0
  
  for (i in 2:length(depth_mid)) {
    depth_bounds[i] <- (depth_mid[i - 1] + depth_mid[i]) / 2
  }
  
  depth_bounds[length(depth_bounds)] <-
    depth_mid[length(depth_mid)] +
    (depth_mid[length(depth_mid)] - depth_bounds[length(depth_bounds) - 1])
  
  depth_bounds[depth_bounds > max_depth] <- max_depth
  depth_bounds
}

profile_integral_raster <- function(z1, z2, a, b, k) {
  a * (z2 - z1) + (b / k) * (exp(-k * z1) - exp(-k * z2))
}

get_month_depth_layers <- function(month_index, n_depths) {
  ((month_index - 1) * n_depths + 1):(month_index * n_depths)
}

prepare_yearly_raster <- function(r, start_year = 1850) {
  n <- nlyr(r)
  yrs <- start_year:(start_year + n - 1)
  names(r) <- as.character(yrs)
  time(r) <- as.Date(paste0(yrs, "-07-01"))
  r
}

# ============================================================
# 7. Land-cover parameters and dynamic depth weights
# ============================================================

taiga_classes <- c(1, 2, 3, 4, 5, 8, 9)
tundra_classes <- c(6, 7, 10, 12, 14)
wetlands_classes <- 11
barren_classes <- c(13, 15, 16)

params <- list(
  taiga  = c(a = 0.007, b = 0.097,  k = 2.7),
  tundra = c(a = 0.010, b = 0.017,  k = 1.9),
  barren = c(a = 0.000, b = 0.0161, k = 1.6)
)

depth_bounds <- get_depth_bounds(depth_vals, max_depth = profile_depth_max)

cat("Depth bounds used:\n")
print(depth_bounds)

depth_thickness <- diff(depth_bounds)

depth_weights <- depth_thickness / sum(depth_thickness)

weighted_depth_mean <- function(r_stack, weights) {
  
  if (nlyr(r_stack) != length(weights)) {
    stop("Number of raster layers and depth weights do not match.")
  }
  
  weighted_sum <- r_stack[[1]] * 0
  weight_sum <- r_stack[[1]] * 0
  
  for (d in seq_len(nlyr(r_stack))) {
    
    valid <- !is.na(r_stack[[d]])
    
    weighted_sum <- weighted_sum + terra::ifel(
      valid,
      r_stack[[d]] * weights[d],
      0
    )
    
    weight_sum <- weight_sum + terra::ifel(
      valid,
      weights[d],
      0
    )
  }
  
  out <- weighted_sum / weight_sum
  out <- terra::ifel(weight_sum > 0, out, NA)
  
  out
}

cat("Depth thicknesses used for temperature weighting:\n")

print(depth_thickness)

cat("Depth weights used for temperature weighting:\n")

print(depth_weights)

make_depth_weights_ALD_increment <- function(previous_max_ALD,
                                             current_ALD,
                                             LC,
                                             depth_bounds,
                                             n_depths,
                                             params,
                                             profile_depth_max,
                                             land_mask) {
  
  previous_max_ALD <- terra::ifel(
    previous_max_ALD > profile_depth_max,
    profile_depth_max,
    previous_max_ALD
  )
  
  current_ALD <- terra::ifel(
    current_ALD > profile_depth_max,
    profile_depth_max,
    current_ALD
  )
  
  taiga_mask <- LC %in% taiga_classes
  tundra_mask <- LC %in% tundra_classes
  wetlands_mask <- LC %in% wetlands_classes
  barren_mask <- LC %in% barren_classes
  
  layer_mass_list <- vector("list", n_depths)
  
  for (d in seq_len(n_depths)) {
    
    z1_layer <- depth_bounds[d]
    z2_layer <- depth_bounds[d + 1]
    
    z1 <- terra::ifel(previous_max_ALD > z1_layer, previous_max_ALD, z1_layer)
    z2 <- terra::ifel(current_ALD < z2_layer, current_ALD, z2_layer)
    
    valid_overlap <- z2 > z1
    
    taiga_mass <- profile_integral_raster(
      z1, z2,
      params$taiga["a"],
      params$taiga["b"],
      params$taiga["k"]
    )
    
    tundra_mass <- profile_integral_raster(
      z1, z2,
      params$tundra["a"],
      params$tundra["b"],
      params$tundra["k"]
    )
    
    barren_mass <- profile_integral_raster(
      z1, z2,
      params$barren["a"],
      params$barren["b"],
      params$barren["k"]
    )
    
    wetland_mass <- z2 - z1
    
    layer_mass <- terra::ifel(
      taiga_mask, taiga_mass,
      terra::ifel(
        tundra_mask, tundra_mass,
        terra::ifel(
          wetlands_mask, wetland_mass,
          terra::ifel(
            barren_mask, barren_mass,
            NA
          )
        )
      )
    )
    
    layer_mass <- terra::ifel(valid_overlap, layer_mass, 0)
    layer_mass <- terra::mask(layer_mass, land_mask)
    
    layer_mass_list[[d]] <- layer_mass
  }
  
  layer_mass_stack <- rast(layer_mass_list)
  names(layer_mass_stack) <- paste0("depth_", seq_len(n_depths))
  
  total_mass <- app(layer_mass_stack, sum, na.rm = TRUE)
  
  weights <- layer_mass_stack / total_mass
  
  weights <- terra::ifel(
    total_mass > 0,
    weights,
    0
  )
  
  weights <- terra::mask(weights, land_mask)
  names(weights) <- paste0("depth_", seq_len(n_depths))
  
  weights
}


# ============================================================
# Standalone depth distribution diagnostic
# ============================================================
# 
# test_year <- 2050
# 
# cat("\n====================================================\n")
# cat("DEPTH DISTRIBUTION DIAGNOSTIC\n")
# cat("Year:", test_year, "\n")
# cat("====================================================\n")
# 
# years_ALD <- as.numeric(format(time(ALD), "%Y"))
# test_i <- which(years_ALD == test_year)
# 
# if (length(test_i) != 1) {
#   stop("test_year not found exactly once in ALD time axis.")
# }
# 
# if (test_i <= 1) {
#   stop("test_year must not be the first year.")
# }
# 
# # -----------------------------
# # Current ALD and previous maximum ALD
# # -----------------------------
# 
# current_ALD_test <- ALD[[test_i]]
# 
# previous_max_ALD_test <- app(
#   ALD[[1:(test_i - 1)]],
#   max,
#   na.rm = TRUE
# )
# 
# previous_max_ALD_test <- terra::mask(previous_max_ALD_test, land_mask)
# current_ALD_test <- terra::mask(current_ALD_test, land_mask)
# 
# # Cap ALD to profile depth max
# previous_max_ALD_cap <- terra::ifel(
#   previous_max_ALD_test > profile_depth_max,
#   profile_depth_max,
#   previous_max_ALD_test
# )
# 
# current_ALD_cap <- terra::ifel(
#   current_ALD_test > profile_depth_max,
#   profile_depth_max,
#   current_ALD_test
# )
# 
# # Cells with new ALD depth interval inside 0–5 m
# active_ALD_increment <- current_ALD_cap > previous_max_ALD_cap
# active_ALD_increment <- terra::mask(active_ALD_increment, land_mask)
# 
# # -----------------------------
# # Newly thawed N increment
# # -----------------------------
# 
# current_thawed_test <- thawed_permafrost_N[[test_i]]
# 
# previous_thawed_max_test <- app(
#   thawed_permafrost_N[[1:(test_i - 1)]],
#   max,
#   na.rm = TRUE
# )
# 
# previous_thawed_max_test <- terra::mask(previous_thawed_max_test, land_mask)
# 
# thawed_increment_test <- current_thawed_test - previous_thawed_max_test
# thawed_increment_test <- terra::ifel(thawed_increment_test < 0, 0, thawed_increment_test)
# thawed_increment_test <- terra::mask(thawed_increment_test, land_mask)
# 
# active_N_increment <- thawed_increment_test > 0
# active_N_increment <- terra::mask(active_N_increment, land_mask)
# 
# # Cells where both conditions are true
# active_test <- active_N_increment & active_ALD_increment
# active_test <- terra::mask(active_test, land_mask)
# 
# active_test_mask <- terra::ifel(active_test, 1, NA)
# 
# NContentLayer_active <- terra::mask(
#   NContentLayer_test,
#   active_test_mask
# )
# 
# weight_sum <- app(NContentLayer_active, sum, na.rm = TRUE)
# weight_sum <- terra::mask(weight_sum, active_test_mask)
# cat("Active cells summary:\n")
# cat("Cells with positive thawed N increment:\n")
# print(freq(active_N_increment))
# 
# cat("Cells with positive ALD increment within profile depth:\n")
# print(freq(active_ALD_increment))
# 
# cat("Cells with both N increment and ALD increment:\n")
# print(freq(active_test))
# 
# # -----------------------------
# # Compute depth weights
# # -----------------------------
# 
# NContentLayer_test <- make_depth_weights_ALD_increment(
#   previous_max_ALD = previous_max_ALD_test,
#   current_ALD = current_ALD_test,
#   LC = LC,
#   depth_bounds = depth_bounds,
#   n_depths = n_depths,
#   params = params,
#   profile_depth_max = profile_depth_max,
#   land_mask = land_mask
# )
# 
# NContentLayer_active <- terra::mask(
#   NContentLayer_test,
#   active_test,
#   maskvalues = FALSE
# )
# 
# # -----------------------------
# # Check weight sums
# # -----------------------------
# 
# weight_sum <- app(NContentLayer_active, sum, na.rm = TRUE)
# 
# weight_sum <- terra::mask(weight_sum, active_test_mask)
# cat("Weight sum range:\n")
# 
# print(global(weight_sum, range, na.rm = TRUE))
# 
# cat("Weight sum mean:\n")
# 
# print(global(weight_sum, mean, na.rm = TRUE))
# 
# weight_ok <- abs(weight_sum - 1) < 1e-6
# 
# print(freq(weight_ok))
# print(freq(weight_ok))
# 
# cat("\nNumber of cells with weight_sum not close to 1:\n")
# weight_bad <- active_test & abs(weight_sum - 1) >= 1e-6
# weight_bad <- terra::mask(weight_bad, active_test, maskvalues = FALSE)
# print(freq(weight_bad))
# 
# bad_cells <- weight_sum == 0
# bad_cells_mask <- terra::ifel(bad_cells, 1, NA)
# 
# cat("LC classes where weight_sum = 0:\n")
# print(freq(mask(LC, bad_cells_mask)))
# 
# 
# # -----------------------------
# # Mean depth-weight profile
# # -----------------------------
# 
# mean_weights <- global(NContentLayer_active, mean, na.rm = TRUE)
# 
# depth_mid <- depth_vals[seq_len(n_depths)]
# 
# mean_weight_df <- data.frame(
#   Depth_m = depth_mid,
#   Mean_weight = mean_weights[, 1]
# )
# 
# mean_weight_df$Mean_weight_norm <- mean_weight_df$Mean_weight /
#   sum(mean_weight_df$Mean_weight, na.rm = TRUE)
# 
# cat("\nMean depth weights:\n")
# print(mean_weight_df)
# 
# cat("\nSum of mean weights:\n")
# print(sum(mean_weight_df$Mean_weight, na.rm = TRUE))
# 
# cat("\nSum of normalized mean weights:\n")
# print(sum(mean_weight_df$Mean_weight_norm, na.rm = TRUE))
# 
# # -----------------------------
# # Conservation check
# # -----------------------------
# 
# newly_thawed_org_test <- thawed_increment_test * f_org * NContentLayer_active
# 
# newly_sum <- app(newly_thawed_org_test, sum, na.rm = TRUE)
# newly_sum <- terra::mask(newly_sum, active_test, maskvalues = FALSE)
# 
# expected_org <- thawed_increment_test * f_org
# expected_org <- terra::mask(expected_org, active_test, maskvalues = FALSE)
# 
# conservation_error <- newly_sum - expected_org
# 
# cat("\nConservation check: sum(depth organic N) - thawed_increment * f_org\n")
# print(global(conservation_error, range, na.rm = TRUE))
# 
# cat("\nMean conservation error:\n")
# print(global(conservation_error, mean, na.rm = TRUE))
# 
# ratio <- newly_sum / expected_org
# 
# global(ratio, range, na.rm = TRUE)
# 
# global(ratio, mean, na.rm = TRUE)
# 
# # -----------------------------
# # Plots
# # -----------------------------
# df <- data.frame(
#   Expected = values(expected_org),
#   Distributed = values(newly_sum)
# )
# 
# df <- na.omit(df)
# 
# ggplot(df, aes(total_thawed_N_201, sum)) +
#   geom_point(alpha = 0.2, size = 0.9) +
#   geom_abline(
#     slope = 1,
#     intercept = 0,
#     colour = "red",
#     linewidth = 0.01
#   ) +
#   coord_equal() +
#   theme_bw() +
#   labs(
#     x = "Expected organic N",
#     y = "Distributed organic N"
#   )
# 
# plot(
#   weight_sum,
#   main = paste("Sum of depth weights", test_year)
# )
# 
# plot(
#   conservation_error,
#   main = paste("Conservation error", test_year)
# )
# 
# plot(
#   NContentLayer_active,
#   main = paste("Depth weights", test_year)
# )
# 
# p_depth_profile <- ggplot(
#   mean_weight_df,
#   aes(x = Depth_m, y = Mean_weight_norm)
# ) +
#   geom_col() +
#   coord_flip() +
#   scale_x_reverse() +
#   theme_bw() +
#   labs(
#     title = paste("Mean depth distribution of newly thawed organic N", test_year),
#     subtitle = "Weights averaged only over cells with positive N and ALD increment",
#     x = "Depth [m]",
#     y = "Normalized mean depth weight"
#   )
# 
# print(p_depth_profile)
# 
# ggsave(
#   filename = file.path(
#     out_dir,
#     paste0("depth_distribution_diagnostic_", ssp, "_", test_year, ".png")
#   ),
#   plot = p_depth_profile,
#   width = 6,
#   height = 5,
#   dpi = 300
# )
# 
# write.csv(
#   mean_weight_df,
#   file.path(
#     out_dir,
#     paste0("depth_distribution_diagnostic_", ssp, "_", test_year, ".csv")
#   ),
#   row.names = FALSE
# )
# 
# # ============================================================
# # Compare model depth weights with expected Palmtag distribution
# # ============================================================
# 
# library(tidyr)
# 
# # Make sure active mask is 1/NA
# active_test_mask <- terra::ifel(active_test, 1, NA)
# 
# # Use only valid active cells
# NContentLayer_active <- terra::mask(NContentLayer_test, active_test_mask)
# 
# # -----------------------------
# # Helper: independent Palmtag expected weights
# # -----------------------------
# 
# make_palmtag_expected_weights <- function(previous_max_ALD,
#                                           current_ALD,
#                                           biome_mask,
#                                           biome_params,
#                                           depth_bounds,
#                                           n_depths,
#                                           profile_depth_max,
#                                           active_mask,
#                                           is_wetland = FALSE) {
#   
#   previous_max_ALD <- terra::ifel(
#     previous_max_ALD > profile_depth_max,
#     profile_depth_max,
#     previous_max_ALD
#   )
#   
#   current_ALD <- terra::ifel(
#     current_ALD > profile_depth_max,
#     profile_depth_max,
#     current_ALD
#   )
#   
#   mass_list <- vector("list", n_depths)
#   
#   for (d in seq_len(n_depths)) {
#     
#     z1_layer <- depth_bounds[d]
#     z2_layer <- depth_bounds[d + 1]
#     
#     z1 <- terra::ifel(previous_max_ALD > z1_layer, previous_max_ALD, z1_layer)
#     z2 <- terra::ifel(current_ALD < z2_layer, current_ALD, z2_layer)
#     
#     valid_overlap <- z2 > z1
#     
#     if (is_wetland) {
#       mass <- z2 - z1
#     } else {
#       mass <- profile_integral_raster(
#         z1, z2,
#         biome_params["a"],
#         biome_params["b"],
#         biome_params["k"]
#       )
#     }
#     
#     mass <- terra::ifel(valid_overlap, mass, 0)
#     mass <- terra::mask(mass, biome_mask)
#     mass <- terra::mask(mass, active_mask)
#     
#     mass_list[[d]] <- mass
#   }
#   
#   mass_stack <- rast(mass_list)
#   names(mass_stack) <- paste0("depth_", seq_len(n_depths))
#   
#   total_mass <- app(mass_stack, sum, na.rm = TRUE)
#   
#   expected_weights <- mass_stack / total_mass
#   expected_weights <- terra::ifel(total_mass > 0, expected_weights, NA)
#   expected_weights <- terra::mask(expected_weights, active_mask)
#   
#   names(expected_weights) <- paste0("depth_", seq_len(n_depths))
#   
#   expected_weights
# }
# 
# # -----------------------------
# # Biome masks
# # -----------------------------
# 
# biome_masks <- list(
#   Taiga = LC %in% taiga_classes,
#   Tundra = LC %in% tundra_classes,
#   Barren = LC %in% barren_classes,
#   Wetlands = LC %in% wetlands_classes
# )
# 
# biome_param_list <- list(
#   Taiga = params$taiga,
#   Tundra = params$tundra,
#   Barren = params$barren,
#   Wetlands = NA
# )
# 
# # -----------------------------
# # Build comparison table
# # -----------------------------
# 
# comparison_list <- list()
# 
# for (biome in names(biome_masks)) {
#   
#   cat("Processing biome:", biome, "\n")
#   
#   biome_mask <- terra::mask(biome_masks[[biome]], active_test_mask)
#   biome_mask <- terra::ifel(biome_mask, 1, NA)
#   
#   model_biome <- terra::mask(NContentLayer_active, biome_mask)
#   
#   palmtag_biome <- make_palmtag_expected_weights(
#     previous_max_ALD = previous_max_ALD_test,
#     current_ALD = current_ALD_test,
#     biome_mask = biome_mask,
#     biome_params = biome_param_list[[biome]],
#     depth_bounds = depth_bounds,
#     n_depths = n_depths,
#     profile_depth_max = profile_depth_max,
#     active_mask = active_test_mask,
#     is_wetland = biome == "Wetlands"
#   )
#   
#   model_mean <- global(model_biome, mean, na.rm = TRUE)[, 1]
#   palmtag_mean <- global(palmtag_biome, mean, na.rm = TRUE)[, 1]
#   
#   # Normalize means so each profile sums to 1
#   model_mean <- model_mean / sum(model_mean, na.rm = TRUE)
#   palmtag_mean <- palmtag_mean / sum(palmtag_mean, na.rm = TRUE)
#   
#   comparison_list[[biome]] <- data.frame(
#     Biome = biome,
#     Depth_m = depth_vals,
#     Model = model_mean,
#     Palmtag_expected = palmtag_mean
#   )
# }
# 
# comparison_df <- bind_rows(comparison_list)
# 
# comparison_long <- comparison_df %>%
#   pivot_longer(
#     cols = c(Model, Palmtag_expected),
#     names_to = "Source",
#     values_to = "Weight"
#   )
# 
# # -----------------------------
# # Plot comparison
# # -----------------------------
# 
# p_compare <- ggplot(
#   comparison_long,
#   aes(x = factor(Depth_m), y = Weight, fill = Source)
# ) +
#   geom_col(
#     position = position_dodge(width = 0.9),
#     width = 0.9
#   ) +
#   coord_flip() +
#   scale_x_discrete(limits = rev(levels(factor(comparison_long$Depth_m)))) +
#   facet_wrap(~ Biome, ncol = 2) +
#   theme_bw(base_size = 13) +
#   theme(
#     panel.grid.major.y = element_blank(),
#     legend.position = "bottom"
#   ) +
#   labs(
#     title = paste0("Model depth weights vs Palmtag distribution (", test_year, ")"),
#     x = "Depth [m]",
#     y = "Fraction of newly thawed organic N",
#     fill = NULL
#   )
# 
# print(p_compare)
# 
# ggsave(
#   filename = file.path(
#     out_dir,
#     paste0("palmtag_vs_model_depth_distribution_", ssp, "_", test_year, ".png")
#   ),
#   plot = p_compare,
#   width = 9,
#   height = 6,
#   dpi = 300
# )
# 
# write.csv(
#   comparison_df,
#   file.path(
#     out_dir,
#     paste0("palmtag_vs_model_depth_distribution_", ssp, "_", test_year, ".csv")
#   ),
#   row.names = FALSE
# )
# 
# cat("\nDiagnostic files saved to:\n")
# cat(out_dir, "\n")
# ============================================================
#  Moisture reference
# ============================================================

theta <- crop(theta, ext(thawed_permafrost_N))
theta <- resample(theta, thawed_permafrost_N[[1]], method = "bilinear")
theta <- terra::mask(theta, land_mask)

cat("Final geometry check theta vs thawed_permafrost_N:\n")
print(compareGeom(theta[[1]], thawed_permafrost_N[[1]], stopOnError = FALSE))

f_theta_present <- clamp(1 - theta, 0, 1)
f_theta_present <- mask(f_theta_present, land_mask)

area_theta <- cellSize(theta[[1]], unit = "m")

f_theta_ref_scalar <- global(
  f_theta_present,
  "mean",
  weights = area_theta,
  na.rm = TRUE
)[1, 1]

cat("f_theta_ref_scalar:", f_theta_ref_scalar, "\n")

if (!is.finite(f_theta_ref_scalar) || f_theta_ref_scalar <= 0) {
  stop("f_theta_ref_scalar is invalid.")
}


get_temp_layers_from_full <- function(month_index_selected,
                                      temp_idx_months,
                                      n_depths_full,
                                      depth_keep) {
  
  month_index_full <- temp_idx_months[month_index_selected]
  
  ((month_index_full - 1) * n_depths_full) + depth_keep
}

# ============================================================
# 8. k_T reference scalar
# ============================================================
compute_k_T_ref_scalar_surface_summer <- function(temp_full,
                                                  temp_month_dates,
                                                  temp_idx_months,
                                                  n_depths_full,
                                                  depth_keep,
                                                  Ea,
                                                  Ed,
                                                  t_opt_C = 28,
                                                  ref_start = 2000,
                                                  ref_end = 2014,
                                                  ref_months_use = 6:8,
                                                  land_mask,
                                                  use_fixed_T_ref = TRUE,
                                                  fixed_T_ref_C = 8) {
  
  if (use_fixed_T_ref) {
    
    T_ref_C <- fixed_T_ref_C
    
  } else {
    
    temp_years  <- as.integer(format(temp_month_dates, "%Y"))
    temp_months <- as.integer(format(temp_month_dates, "%m"))
    
    ref_months <- which(
      temp_years >= ref_start &
        temp_years <= ref_end &
        temp_months %in% ref_months_use
    )
    
    if (length(ref_months) == 0) {
      stop("No reference temperature months found.")
    }
    
    area_rast <- terra::cellSize(temp_full[[1]], unit = "m")
    area_rast <- terra::mask(area_rast, land_mask)
    
    # get surface layer
    surface_depth_keep <- depth_keep[1]
    
    ref_temps <- numeric(length(ref_months))
    
    for (ii in seq_along(ref_months)) {
      
      j <- ref_months[ii]
      
      cat(
        "Surface summer reference month",
        ii, "of", length(ref_months),
        "date:", as.character(temp_month_dates[j]), "\n"
      )
      
      surface_layer <- get_temp_layers_from_full(
        month_index_selected = j,
        temp_idx_months = temp_idx_months,
        n_depths_full = n_depths_full,
        depth_keep = surface_depth_keep
      )
      
      temp_surface_C <- temp_full[[surface_layer]] - 273.15
      temp_surface_C <- terra::mask(temp_surface_C, land_mask)
      
      ref_temps[ii] <- terra::global(
        temp_surface_C,
        "mean",
        weights = area_rast,
        na.rm = TRUE
      )[1, 1]
      
      rm(temp_surface_C)
      gc()
    }
    
    T_ref_C <- mean(ref_temps, na.rm = TRUE)
  }
  
  k_T_ref_scalar <- peaked_arrhenius(
    temp_C = T_ref_C,
    Ea = Ea,
    Ed = Ed,
    t_opt_C = t_opt_C
  )
  
  cat("Surface summer T_ref_C:", T_ref_C, "\n")
  cat("Surface summer k_T_ref_scalar:", k_T_ref_scalar, "\n")
  
  if (!is.finite(k_T_ref_scalar) || k_T_ref_scalar <= 0) {
    stop("Surface-summer k_T_ref_scalar is invalid.")
  }
  
  k_T_ref_scalar
}


# old version: calculates the ref k t 0-5 m depth over all months;
# compute_k_T_ref_scalar_depth_pools <- function(temp_full,
#                                                temp_month_dates,
#                                                temp_idx_months,
#                                                n_depths_full,
#                                                depth_keep,
#                                                Ea,
#                                                Ed,
#                                                t_opt_C = 28,
#                                                ref_start = 2000,
#                                                ref_end = 2014,
#                                                land_mask) {
#
#   temp_years  <- as.integer(format(temp_month_dates, "%Y"))
#   temp_months <- as.integer(format(temp_month_dates, "%m"))
#
#   ref_months <- which(
#     temp_years >= ref_start &
#       temp_years <= ref_end &
#       temp_months %in% 1:12
#   )
#
#   if (length(ref_months) == 0) {
#     stop("No temperature months found for reference period.")
#   }
#
#   area_rast <- cellSize(temp_full[[1]], unit = "m")
#   area_rast <- terra::mask(area_rast, land_mask)
#
#   ref_means <- numeric(length(ref_months))
#
#   for (ii in seq_along(ref_months)) {
#
#     j <- ref_months[ii]
#
#     cat(
#       "Reference month",
#       ii, "of", length(ref_months),
#       "date:", as.character(temp_month_dates[j]), "\n"
#     )
#
#     depth_layers <- get_temp_layers_from_full(
#       month_index_selected = j,
#       temp_idx_months = temp_idx_months,
#       n_depths_full = n_depths_full,
#       depth_keep = depth_keep
#     )
#
#     temp_month_depths <- temp_full[[depth_layers]] - 273.15
#
#     k_T_depths <- peaked_arrhenius(
#       temp_month_depths,
#       Ea = Ea,
#       Ed = Ed,
#       t_opt_C = t_opt_C
#     )
#
#     k_T_depths <- terra::ifel(temp_month_depths < 0, 0, k_T_depths)
#     k_T_depths <- terra::ifel(is.na(temp_month_depths), NA, k_T_depths)
#     k_T_depths <- terra::mask(k_T_depths, land_mask)
#
#     k_T_depth_mean <- weighted_depth_mean(k_T_depths, depth_weights)
#     k_T_depth_mean <- terra::mask(k_T_depth_mean, land_mask)
#
#     ref_means[ii] <- terra::global(
#       k_T_depth_mean,
#       "mean",
#       weights = area_rast,
#       na.rm = TRUE
#     )[1, 1]
#
#     rm(temp_month_depths, k_T_depths, k_T_depth_mean)
#     gc()
#   }
#
#   k_T_ref_scalar <- mean(ref_means, na.rm = TRUE)
#
#   cat("Final k_T_ref_scalar:", k_T_ref_scalar, "\n")
#
#   if (!is.finite(k_T_ref_scalar) || k_T_ref_scalar <= 0) {
#     stop("k_T_ref_scalar is invalid.")
#   }
#
#   k_T_ref_scalar
# }



###############################################################################
###  helper functions:

use_temperature_scaling <- TRUE

# with temperature:
compute_k_factor_depths <- function(temp_full,
                                    month_index_selected,
                                    temp_idx_months,
                                    n_depths_full,
                                    depth_keep,
                                    Ea,
                                    Ed,
                                    t_opt_C,
                                    k_T_ref_scalar,
                                    land_mask,
                                    use_temperature_scaling = TRUE) {
  
  depth_layers <- get_temp_layers_from_full(
    month_index_selected = month_index_selected,
    temp_idx_months = temp_idx_months,
    n_depths_full = n_depths_full,
    depth_keep = depth_keep
  )
  
  if (!use_temperature_scaling) {
    out <- rast(
      replicate(length(depth_keep), land_mask * 1, simplify = FALSE)
    )
    out <- terra::mask(out, land_mask)
    names(out) <- paste0("depth_", seq_along(depth_keep))
    return(out)
  }
  
  temp_month_depths <- temp_full[[depth_layers]] - 273.15
  
  k_T_depths <- peaked_arrhenius(
    temp_month_depths,
    Ea = Ea,
    Ed = Ed,
    t_opt_C = t_opt_C
  )
  
  k_T_depths <- terra::ifel(temp_month_depths < 0, 0, k_T_depths)
  k_T_depths <- terra::ifel(is.na(temp_month_depths), NA, k_T_depths)
  
  k_factor_depths <- k_T_depths / k_T_ref_scalar
  k_factor_depths <- terra::mask(k_factor_depths, land_mask)
  
  names(k_factor_depths) <- paste0("depth_", seq_along(depth_keep))
  
  k_factor_depths
}


mineralise_depth_pools <- function(organic_pool_depths,
                                   k_factor_depths,
                                   base_mineralisation_rate_monthly) {
  
  k_depths <- base_mineralisation_rate_monthly * k_factor_depths
  
  mineralised_depths <- organic_pool_depths * k_depths
  
  mineralised_depths <- terra::ifel(
    mineralised_depths > organic_pool_depths,
    organic_pool_depths,
    mineralised_depths
  )
  
  mineralised_depths <- terra::ifel(mineralised_depths < 0, 0, mineralised_depths)
  
  organic_pool_depths <- organic_pool_depths - mineralised_depths
  organic_pool_depths <- terra::ifel(organic_pool_depths < 0, 0, organic_pool_depths)
  
  mineralised_this_month <- app(mineralised_depths, sum, na.rm = TRUE)
  
  list(
    mineralised_this_month = mineralised_this_month,
    organic_pool_depths = organic_pool_depths
  )
}

cat("Before k_T_ref_scalar\n")


k_T_ref_scalar <- compute_k_T_ref_scalar_surface_summer(
  temp_full = temp_full,
  temp_month_dates = temp_month_dates,
  temp_idx_months = temp_idx_months,
  n_depths_full = n_depths_full,
  depth_keep = depth_keep,
  Ea = Ea_nitrification,
  Ed = Ed_nitrification,
  t_opt_C = t_opt_C,
  ref_start = ref_start,
  ref_end = ref_end,
  ref_months_use = 6:8,
  land_mask = land_mask,
  use_fixed_T_ref = TRUE,
  fixed_T_ref_C = 8
)

# old:
# k_T_ref_scalar <- compute_k_T_ref_scalar_depth_pools(
#   temp_full = temp_full,
#   temp_month_dates = temp_month_dates,
#   temp_idx_months = temp_idx_months,
#   n_depths_full = n_depths_full,
#   depth_keep = depth_keep,
#   Ea = Ea_nitrification,
#   Ed = Ed_nitrification,
#   t_opt_C = t_opt_C,
#   ref_start = ref_start,
#   ref_end = ref_end,
#   land_mask = land_mask
# )

#k_T_ref_scalar <- 1




# ============================================================
# 9. Main calculation

# version that includes all months for mineralisation
compute_flux_monthly_depth_pools <- function(
    permafrost_total_thawed,
    ALD,
    temp_full,
    temp_month_dates,
    k_base_mineralisation,
    Ea,
    Ed,
    theta,
    k_T_ref_scalar,
    f_theta_ref_scalar,
    start_year = 1850,
    t_opt_C = 28,
    f_inorg_rapid = 0.1136,
    n_depths,
    average_arctic_thaw_length = 0.3,
    depth_bounds,
    params,
    profile_depth_max,
    LC,
    write_yearly_rasters = FALSE,
    use_temperature_scaling = TRUE
) {
  
  f_org <- 1 - f_inorg_rapid
  
  n_years <- nlyr(permafrost_total_thawed)
  years <- start_year:(start_year + n_years - 1)
  temp_years <- as.integer(format(temp_month_dates, "%Y"))
  
  land_mask <- !is.na(permafrost_total_thawed[[n_years]])
  
  if (nlyr(ALD) != n_years) {
    stop("ALD and permafrost_total_thawed do not have the same number of layers.")
  }
  
  permafrost_total_thawed <- terra::mask(permafrost_total_thawed, land_mask)
  ALD <- terra::mask(ALD, land_mask)
  
  baseline_idx <- which(years >= 1850 & years <= 1900)
  
  baseline_thawed_mean <- app(
    permafrost_total_thawed[[baseline_idx]],
    mean,
    na.rm = TRUE
  )
  
  baseline_thawed_mean <- terra::mask(baseline_thawed_mean, land_mask)
  
  diagnostic_rows <- list()
  
  organic_pool_yearly_list <- vector("list", n_years)
  thawed_increment_yearly_list <- vector("list", n_years)
  new_thawed_organic_N_yearly_list <- vector("list", n_years)
  new_thawed_total_N_yearly_list <- vector("list", n_years)
  max_ALD_yearly_list <- vector("list", n_years)
  net_organic_pool_change_yearly_list <- vector("list", n_years)
  rapid_bioavailable_pool_yearly_list <- vector("list", n_years)
  
  base_mineralisation_rate_monthly <-
    k_base_mineralisation * average_arctic_thaw_length
  
  cat("Base mineralisation rate per month:",
      base_mineralisation_rate_monthly, "\n")
  
  zero_layer <- terra::mask(permafrost_total_thawed[[1]] * 0, land_mask)
  
  previous_baseline_relative_thawed <- zero_layer
  
  previous_max_ALD <- terra::mask(
    ALD[[1]],
    land_mask
  )
  
  organic_pool_depths <- rast(
    replicate(n_depths, zero_layer, simplify = FALSE)
  )
  
  names(organic_pool_depths) <- paste0("depth_", seq_len(n_depths))
  
  organic_pool_depths <- terra::mask(
    organic_pool_depths,
    land_mask
  )
  
  for (i in seq_len(n_years)) {
    
    yr <- years[i]
    cat("Processing year", yr, "(", i, "of", n_years, ")\n")
    
    mineralised_year_list <- list()
    bioavailable_year_list <- list()
    k_T_year_list <- list()
    k_combined_year_list <- list()
    
    current_thawed <- terra::mask(
      permafrost_total_thawed[[i]],
      land_mask
    )
    
    current_baseline_relative_thawed <- current_thawed - baseline_thawed_mean
    
    current_baseline_relative_thawed <- terra::ifel(
      current_baseline_relative_thawed < 0,
      0,
      current_baseline_relative_thawed
    )
    
    current_baseline_relative_thawed <- terra::mask(
      current_baseline_relative_thawed,
      land_mask
    )
    
    current_ALD <- terra::mask(
      ALD[[i]],
      land_mask
    )
    
    new_thawed_total_N <- current_baseline_relative_thawed -
      previous_baseline_relative_thawed
    
    # should negative new_thawed_total_N be set to 0?
    # new_thawed_total_N <- terra::ifel(
    #   new_thawed_total_N < 0,
    #   0,
    #   new_thawed_total_N
    # )
    
    new_thawed_total_N <- terra::mask(
      new_thawed_total_N,
      land_mask
    )
    
    new_thawed_organic_N <- new_thawed_total_N * f_org
    new_thawed_organic_N <- terra::mask(new_thawed_organic_N, land_mask)
    
    thawed_increment <- new_thawed_total_N
    
    if (i == 1) {
      
      NContentLayer_year <- rast(
        replicate(n_depths, zero_layer, simplify = FALSE)
      )
      names(NContentLayer_year) <- paste0("depth_", seq_len(n_depths))
      NContentLayer_year <- terra::mask(NContentLayer_year, land_mask)
      
    } else {
      
      NContentLayer_year <- make_depth_weights_ALD_increment(
        previous_max_ALD = previous_max_ALD,
        current_ALD = current_ALD,
        LC = LC,
        depth_bounds = depth_bounds,
        n_depths = n_depths,
        params = params,
        profile_depth_max = profile_depth_max,
        land_mask = land_mask
      )
    }
    
    new_thawed_total_N_yearly_list[[i]] <- new_thawed_total_N
    new_thawed_organic_N_yearly_list[[i]] <- new_thawed_organic_N
    
    newly_thawed_org_depths <- new_thawed_organic_N * NContentLayer_year
    newly_thawed_org_depths <- terra::mask(newly_thawed_org_depths, land_mask)
    
    
    
    organic_pool_depths <- organic_pool_depths + newly_thawed_org_depths
    organic_pool_depths <- terra::mask(organic_pool_depths, land_mask)
    
    inorg_rapid_yearly <- new_thawed_total_N * f_inorg_rapid
    
    inorg_rapid_yearly <- terra::mask(inorg_rapid_yearly, land_mask)
    rapid_bioavailable_pool_yearly_list[[i]] <- inorg_rapid_yearly
    
    
    theta_year <- theta[[min(i, nlyr(theta))]]
    theta_year <- terra::mask(theta_year, land_mask)
    
    f_theta <- terra::clamp(1 - theta_year, lower = 0, upper = 1)
    k_moisture <- f_theta / f_theta_ref_scalar
    
    k_moisture <- terra::ifel(
      land_mask,
      terra::ifel(is.na(k_moisture), 1, k_moisture),
      NA
    )
    k_moisture <- terra::mask(k_moisture, land_mask)
    
    month_idx_year <- which(temp_years == yr)
    
    if (length(month_idx_year) == 0) {
      stop(paste("No temperature months found for year", yr))
    }
    
    for (j in month_idx_year) {
      
      month_number <- as.integer(format(temp_month_dates[j], "%m"))
      
      
      k_T_depths <- compute_k_factor_depths(
        temp_full = temp_full,
        month_index_selected = j,
        temp_idx_months = temp_idx_months,
        n_depths_full = n_depths_full,
        depth_keep = depth_keep,
        Ea = Ea,
        Ed = Ed,
        t_opt_C = t_opt_C,
        k_T_ref_scalar = k_T_ref_scalar,
        land_mask = land_mask,
        use_temperature_scaling = use_temperature_scaling
      )
      
      k_T_depths <- terra::mask(k_T_depths, land_mask)
      
      k_T_depth_mean <- weighted_depth_mean(k_T_depths, depth_weights)
      k_T_depth_mean <- terra::mask(k_T_depth_mean, land_mask)
      
      organic_pool_sum <- app(organic_pool_depths, sum, na.rm = TRUE)
      
      k_T_pool_weighted <- app(
        k_T_depths * organic_pool_depths,
        sum,
        na.rm = TRUE
      ) / organic_pool_sum
      
      k_T_pool_weighted <- terra::ifel(
        organic_pool_sum > 0,
        k_T_pool_weighted,
        NA
      )
      
      k_T_pool_weighted <- terra::mask(k_T_pool_weighted, land_mask)
      
      k_combined_depths <- k_T_depths * k_moisture
      k_combined_depths <- terra::mask(k_combined_depths, land_mask)
      
      k_combined_pool_weighted <- app(
        k_combined_depths * organic_pool_depths,
        sum,
        na.rm = TRUE
      ) / organic_pool_sum
      
      k_combined_pool_weighted <- terra::ifel(
        organic_pool_sum > 0,
        k_combined_pool_weighted,
        NA
      )
      
      k_combined_pool_weighted <- terra::mask(
        k_combined_pool_weighted,
        land_mask
      )
      
      k_combined_depth_mean <- weighted_depth_mean(
        k_combined_depths,
        depth_weights
      )
      k_combined_depth_mean <- terra::mask(k_combined_depth_mean, land_mask)
      
      if (yr == 1850) {
        
        cat("Organic pool before mineralisation:\n")
        
        print(
          global(
            app(organic_pool_depths, sum, na.rm = TRUE),
            c("sum", "mean", "max"),
            na.rm = TRUE
          )
        )
        
      }
      
      mineralisation_result <- mineralise_depth_pools(
        organic_pool_depths = organic_pool_depths,
        k_factor_depths = k_combined_depths,
        base_mineralisation_rate_monthly = base_mineralisation_rate_monthly
      )
      
      mineralised_this_month <- terra::mask(
        mineralisation_result$mineralised_this_month,
        land_mask
      )
      
      if (yr == 1850) {
        
        cat("Mineralised this month:\n")
        
        print(
          global(
            mineralised_this_month,
            c("sum", "mean", "max"),
            na.rm = TRUE
          )
        )
        
      }
      
      organic_pool_depths <- terra::mask(
        mineralisation_result$organic_pool_depths,
        land_mask
      )
      
      inorg_rapid_this_month <- zero_layer
      
      if (month_number %in% 6:8) {
        inorg_rapid_this_month <- inorg_rapid_yearly / 3
      }
      
      bioavailable_this_month <- terra::mask(
        mineralised_this_month + inorg_rapid_this_month,
        land_mask
      )
      
      bioavailable_this_month <- terra::ifel(
        bioavailable_this_month < 0,
        0,
        bioavailable_this_month
      )
      
      area_rast <- terra::cellSize(k_T_depth_mean, unit = "m")
      
      diagnostic_rows[[length(diagnostic_rows) + 1]] <- data.frame(
        Date = temp_month_dates[j],
        Year = yr,
        Month = month_number,
        mineralised_pg_monthly = terra::global(
          mineralised_this_month * area_rast,
          "sum",
          na.rm = TRUE
        )[1, 1] / 1e12,
        bioavailable_pg_monthly = terra::global(
          bioavailable_this_month * area_rast,
          "sum",
          na.rm = TRUE
        )[1, 1] / 1e12,
        k_T_mean_monthly = terra::global(
          k_T_depth_mean,
          "mean",
          weights = area_rast,
          na.rm = TRUE
        )[1, 1],
        k_moisture_mean_monthly = terra::global(
          k_moisture,
          "mean",
          weights = area_rast,
          na.rm = TRUE
        )[1, 1],
        k_combined_pool_weighted_mean = terra::global(
          k_combined_pool_weighted,
          "mean",
          weights = area_rast,
          na.rm = TRUE
        )[1, 1],
        k_T_pool_weighted_mean = terra::global(
          k_T_pool_weighted,
          "mean",
          weights = area_rast,
          na.rm = TRUE
        )[1, 1],
        k_combined_mean_monthly = terra::global(
          k_combined_depth_mean,
          "mean",
          weights = area_rast,
          na.rm = TRUE
        )[1, 1]
      )
      
      mineralised_year_list[[length(mineralised_year_list) + 1]] <-
        mineralised_this_month
      
      bioavailable_year_list[[length(bioavailable_year_list) + 1]] <-
        bioavailable_this_month
      
      k_T_year_list[[length(k_T_year_list) + 1]] <-
        k_T_depth_mean
      
      k_combined_year_list[[length(k_combined_year_list) + 1]] <-
        k_combined_depth_mean
      
      rm(
        k_T_depths,
        k_T_depth_mean,
        k_combined_depths,
        k_combined_depth_mean,
        mineralisation_result,
        mineralised_this_month,
        bioavailable_this_month,
        inorg_rapid_this_month,
        area_rast
      )
      gc()
    }
    
    if (write_yearly_rasters) {
      
      r_min <- rast(mineralised_year_list)
      terra::time(r_min) <- seq(as.Date(paste0(yr, "-01-16")), by = "month", length.out = 12)
      names(r_min) <- month.abb
      
      terra::writeRaster(
        r_min,
        file.path(out_dir_yearly, paste0("region_mineralised_N_monthly_", ssp, "_", yr, "_w_temp_sm.nc")),
        overwrite = TRUE,
        filetype = "NetCDF"
      )
      
      r_bio <- rast(bioavailable_year_list)
      terra::time(r_bio) <- seq(as.Date(paste0(yr, "-01-16")), by = "month", length.out = 12)
      names(r_bio) <- month.abb
      
      terra::writeRaster(
        r_bio,
        file.path(out_dir_yearly, paste0("region_bioavailable_N_monthly_", ssp, "_", yr, "_w_temp_sm.nc")),
        overwrite = TRUE,
        filetype = "NetCDF"
      )
      
      r_kT <- rast(k_T_year_list)
      terra::time(r_kT) <- seq(as.Date(paste0(yr, "-01-16")), by = "month", length.out = 12)
      names(r_kT) <- month.abb
      
      terra::writeRaster(
        r_kT,
        file.path(out_dir_yearly, paste0("region_k_T_temperature_only_monthly_", ssp, "_", yr, ".nc")),
        overwrite = TRUE,
        filetype = "NetCDF"
      )
      
      r_kcombined <- rast(k_combined_year_list)
      terra::time(r_kcombined) <- seq(as.Date(paste0(yr, "-01-16")), by = "month", length.out = 12)
      names(r_kcombined) <- month.abb
      
      terra::writeRaster(
        r_kcombined,
        file.path(out_dir_yearly, paste0("region_k_combined_monthly_", ssp, "_", yr, "_w_temp_sm.nc")),
        overwrite = TRUE,
        filetype = "NetCDF"
      )
      
      rm(r_min, r_bio, r_kT, r_kcombined)
    }
    
    organic_pool_remaining <- app(organic_pool_depths, sum, na.rm = TRUE)
    organic_pool_remaining <- terra::mask(organic_pool_remaining, land_mask)
    
    yearly_mineralised_N <- app(rast(mineralised_year_list), sum, na.rm = TRUE)
    yearly_mineralised_N <- terra::mask(yearly_mineralised_N, land_mask)
    
    net_organic_pool_change <- new_thawed_organic_N - yearly_mineralised_N
    net_organic_pool_change <- terra::mask(net_organic_pool_change, land_mask)
    
    
    previous_max_ALD <- terra::ifel(
      current_ALD > previous_max_ALD,
      current_ALD,
      previous_max_ALD
    )
    previous_max_ALD <- terra::mask(previous_max_ALD, land_mask)
    
    previous_baseline_relative_thawed <- current_baseline_relative_thawed
    
    
    organic_pool_yearly_list[[i]] <- organic_pool_remaining
    thawed_increment_yearly_list[[i]] <- thawed_increment
    max_ALD_yearly_list[[i]] <- previous_max_ALD
    net_organic_pool_change_yearly_list[[i]] <- net_organic_pool_change
    
    rm(
      mineralised_year_list,
      bioavailable_year_list,
      k_T_year_list,
      k_combined_year_list,
      current_thawed,
      current_ALD,
      thawed_increment,
      net_organic_pool_change,
      new_thawed_total_N,
      new_thawed_organic_N,
      newly_thawed_org_depths,
      inorg_rapid_yearly,
      NContentLayer_year,
      organic_pool_remaining,
      theta_year,
      f_theta,
      k_moisture,
      current_baseline_relative_thawed
    )
    gc()
  }
  
  list(
    diagnostic_df = dplyr::bind_rows(diagnostic_rows),
    organic_pool_remaining_yearly = rast(organic_pool_yearly_list),
    net_organic_pool_change_yearly = rast(net_organic_pool_change_yearly_list),
    thawed_N_increment_yearly = rast(thawed_increment_yearly_list),
    new_thawed_organic_N_yearly = rast(new_thawed_organic_N_yearly_list),
    new_thawed_total_N_yearly = rast(new_thawed_total_N_yearly_list),
    rapid_bioavailable_pool_yearly = rast(rapid_bioavailable_pool_yearly_list),
    max_ALD_yearly = rast(max_ALD_yearly_list)
  )
}



# compute_flux_monthly_depth_pools <- function(
    #     permafrost_total_thawed,
#     ALD,
#     temp_full,
#     temp_month_dates,
#     k_base_mineralisation,
#     Ea,
#     Ed,
#     theta,
#     k_T_ref_scalar,
#     f_theta_ref_scalar,
#     start_year = 1850,
#     t_opt_C = 28,
#     f_inorg_rapid = 0.1136,
#     n_depths,
#     average_arctic_thaw_length = 0.3,
#     active_mineralisation_months = c(6, 7, 8, 9),
#     depth_bounds,
#     params,
#     profile_depth_max,
#     LC,
#     write_yearly_rasters = FALSE,
#     use_temperature_scaling = TRUE
# ) {
#   
#   f_org <- 1 - f_inorg_rapid
#   
#   n_years <- nlyr(permafrost_total_thawed)
#   years <- start_year:(start_year + n_years - 1)
#   temp_years <- as.integer(format(temp_month_dates, "%Y"))
#   
#   land_mask <- !is.na(permafrost_total_thawed[[n_years]])
#   
#   if (nlyr(ALD) != n_years) {
#     stop("ALD and permafrost_total_thawed do not have the same number of layers.")
#   }
#   
#   permafrost_total_thawed <- terra::mask(permafrost_total_thawed, land_mask)
#   ALD <- terra::mask(ALD, land_mask)
#   
#   baseline_idx <- which(years >= 1850 & years <= 1900)
#   
#   baseline_thawed_mean <- app(
#     permafrost_total_thawed[[baseline_idx]],
#     mean,
#     na.rm = TRUE
#   )
#   
#   baseline_thawed_mean <- terra::mask(baseline_thawed_mean, land_mask)
#   
#   diagnostic_rows <- list()
#   
#   organic_pool_yearly_list <- vector("list", n_years)
#   thawed_increment_yearly_list <- vector("list", n_years)
#   new_thawed_organic_N_yearly_list <- vector("list", n_years)
#   new_thawed_total_N_yearly_list <- vector("list", n_years)
#   max_ALD_yearly_list <- vector("list", n_years)
#   net_organic_pool_change_yearly_list <- vector("list", n_years)
#   
#   base_mineralisation_rate_monthly <-
#     k_base_mineralisation *
#     average_arctic_thaw_length /
#     length(active_mineralisation_months)
#   
#   cat("Active mineralisation months:",
#       paste(active_mineralisation_months, collapse = ", "), "\n")
#   cat("Base mineralisation rate per active month:",
#       base_mineralisation_rate_monthly, "\n")
#   
#   zero_layer <- terra::mask(permafrost_total_thawed[[1]] * 0, land_mask)
#   
#   previous_baseline_relative_thawed <- zero_layer
#   
#   previous_max_ALD <- terra::mask(ALD[[1]], land_mask)
#   
#   organic_pool_depths <- rast(
#     replicate(n_depths, zero_layer, simplify = FALSE)
#   )
#   names(organic_pool_depths) <- paste0("depth_", seq_len(n_depths))
#   organic_pool_depths <- terra::mask(organic_pool_depths, land_mask)
#   
#   for (i in seq_len(n_years)) {
#     
#     yr <- years[i]
#     cat("Processing year", yr, "(", i, "of", n_years, ")\n")
#     
#     mineralised_year_list <- list()
#     bioavailable_year_list <- list()
#     k_T_year_list <- list()
#     k_combined_year_list <- list()
#     
#     current_thawed <- terra::mask(permafrost_total_thawed[[i]], land_mask)
#     
#     current_baseline_relative_thawed <- current_thawed - baseline_thawed_mean
#     current_baseline_relative_thawed <- terra::ifel(
#       current_baseline_relative_thawed < 0,
#       0,
#       current_baseline_relative_thawed
#     )
#     current_baseline_relative_thawed <- terra::mask(
#       current_baseline_relative_thawed,
#       land_mask
#     )
#     
#     current_ALD <- terra::mask(ALD[[i]], land_mask)
#     
#     new_thawed_total_N <- current_baseline_relative_thawed -
#       previous_baseline_relative_thawed
#     
#     new_thawed_total_N <- terra::ifel(
#       new_thawed_total_N < 0,
#       0,
#       new_thawed_total_N
#     )
#     new_thawed_total_N <- terra::mask(new_thawed_total_N, land_mask)
#     
#     new_thawed_organic_N <- new_thawed_total_N * f_org
#     new_thawed_organic_N <- terra::mask(new_thawed_organic_N, land_mask)
#     
#     thawed_increment <- new_thawed_total_N
#     
#     if (i == 1) {
#       NContentLayer_year <- rast(
#         replicate(n_depths, zero_layer, simplify = FALSE)
#       )
#       names(NContentLayer_year) <- paste0("depth_", seq_len(n_depths))
#       NContentLayer_year <- terra::mask(NContentLayer_year, land_mask)
#     } else {
#       NContentLayer_year <- make_depth_weights_ALD_increment(
#         previous_max_ALD = previous_max_ALD,
#         current_ALD = current_ALD,
#         LC = LC,
#         depth_bounds = depth_bounds,
#         n_depths = n_depths,
#         params = params,
#         profile_depth_max = profile_depth_max,
#         land_mask = land_mask
#       )
#     }
#     
#     new_thawed_total_N_yearly_list[[i]] <- new_thawed_total_N
#     new_thawed_organic_N_yearly_list[[i]] <- new_thawed_organic_N
#     
#     newly_thawed_org_depths <- new_thawed_organic_N * NContentLayer_year
#     newly_thawed_org_depths <- terra::mask(newly_thawed_org_depths, land_mask)
#     
#     organic_pool_depths <- organic_pool_depths + newly_thawed_org_depths
#     organic_pool_depths <- terra::mask(organic_pool_depths, land_mask)
#     
#     inorg_rapid_yearly <- new_thawed_total_N * f_inorg_rapid
#     inorg_rapid_yearly <- terra::mask(inorg_rapid_yearly, land_mask)
#     
#     theta_year <- theta[[min(i, nlyr(theta))]]
#     theta_year <- terra::mask(theta_year, land_mask)
#     
#     f_theta <- terra::clamp(1 - theta_year, lower = 0, upper = 1)
#     k_moisture <- f_theta / f_theta_ref_scalar
#     k_moisture <- terra::ifel(
#       land_mask,
#       terra::ifel(is.na(k_moisture), 1, k_moisture),
#       NA
#     )
#     k_moisture <- terra::mask(k_moisture, land_mask)
#     
#     month_idx_year <- which(temp_years == yr)
#     
#     if (length(month_idx_year) == 0) {
#       stop(paste("No temperature months found for year", yr))
#     }
#     
#     for (j in month_idx_year) {
#       
#       month_number <- as.integer(format(temp_month_dates[j], "%m"))
#       mineralisation_active <- month_number %in% active_mineralisation_months
#       
#       k_T_depths <- compute_k_factor_depths(
#         temp_full = temp_full,
#         month_index_selected = j,
#         temp_idx_months = temp_idx_months,
#         n_depths_full = n_depths_full,
#         depth_keep = depth_keep,
#         Ea = Ea,
#         Ed = Ed,
#         t_opt_C = t_opt_C,
#         k_T_ref_scalar = k_T_ref_scalar,
#         land_mask = land_mask,
#         use_temperature_scaling = use_temperature_scaling
#       )
#       
#       k_T_depths <- terra::mask(k_T_depths, land_mask)
#       
#       k_T_depth_mean <- weighted_depth_mean(k_T_depths, depth_weights)
#       k_T_depth_mean <- terra::mask(k_T_depth_mean, land_mask)
#       
#       organic_pool_sum <- app(organic_pool_depths, sum, na.rm = TRUE)
#       
#       k_T_pool_weighted <- app(
#         k_T_depths * organic_pool_depths,
#         sum,
#         na.rm = TRUE
#       ) / organic_pool_sum
#       
#       k_T_pool_weighted <- terra::ifel(
#         organic_pool_sum > 0,
#         k_T_pool_weighted,
#         NA
#       )
#       k_T_pool_weighted <- terra::mask(k_T_pool_weighted, land_mask)
#       
#       k_combined_depths <- k_T_depths * k_moisture
#       k_combined_depths <- terra::mask(k_combined_depths, land_mask)
#       
#       k_combined_pool_weighted <- app(
#         k_combined_depths * organic_pool_depths,
#         sum,
#         na.rm = TRUE
#       ) / organic_pool_sum
#       
#       k_combined_pool_weighted <- terra::ifel(
#         organic_pool_sum > 0,
#         k_combined_pool_weighted,
#         NA
#       )
#       k_combined_pool_weighted <- terra::mask(
#         k_combined_pool_weighted,
#         land_mask
#       )
#       
#       k_combined_depth_mean <- weighted_depth_mean(
#         k_combined_depths,
#         depth_weights
#       )
#       k_combined_depth_mean <- terra::mask(k_combined_depth_mean, land_mask)
#       
#       if (mineralisation_active) {
#         
#         mineralisation_result <- mineralise_depth_pools(
#           organic_pool_depths = organic_pool_depths,
#           k_factor_depths = k_combined_depths,
#           base_mineralisation_rate_monthly = base_mineralisation_rate_monthly
#         )
#         
#         mineralised_this_month <- terra::mask(
#           mineralisation_result$mineralised_this_month,
#           land_mask
#         )
#         
#         organic_pool_depths <- terra::mask(
#           mineralisation_result$organic_pool_depths,
#           land_mask
#         )
#         
#       } else {
#         
#         mineralised_this_month <- zero_layer
#       }
#       
#       inorg_rapid_this_month <- zero_layer
#       
#       if (month_number %in% 6:8) {
#         inorg_rapid_this_month <- inorg_rapid_yearly / 3
#       }
#       
#       bioavailable_this_month <- terra::mask(
#         mineralised_this_month + inorg_rapid_this_month,
#         land_mask
#       )
#       
#       bioavailable_this_month <- terra::ifel(
#         bioavailable_this_month < 0,
#         0,
#         bioavailable_this_month
#       )
#       
#       area_rast <- terra::cellSize(k_T_depth_mean, unit = "m")
#       
#       diagnostic_rows[[length(diagnostic_rows) + 1]] <- data.frame(
#         Date = temp_month_dates[j],
#         Year = yr,
#         Month = month_number,
#         mineralisation_active = mineralisation_active,
#         mineralised_pg_monthly = terra::global(
#           mineralised_this_month * area_rast,
#           "sum",
#           na.rm = TRUE
#         )[1, 1] / 1e12,
#         bioavailable_pg_monthly = terra::global(
#           bioavailable_this_month * area_rast,
#           "sum",
#           na.rm = TRUE
#         )[1, 1] / 1e12,
#         k_T_mean_monthly = terra::global(
#           k_T_depth_mean,
#           "mean",
#           weights = area_rast,
#           na.rm = TRUE
#         )[1, 1],
#         k_moisture_mean_monthly = terra::global(
#           k_moisture,
#           "mean",
#           weights = area_rast,
#           na.rm = TRUE
#         )[1, 1],
#         k_combined_pool_weighted_mean = terra::global(
#           k_combined_pool_weighted,
#           "mean",
#           weights = area_rast,
#           na.rm = TRUE
#         )[1, 1],
#         k_T_pool_weighted_mean = terra::global(
#           k_T_pool_weighted,
#           "mean",
#           weights = area_rast,
#           na.rm = TRUE
#         )[1, 1],
#         k_combined_mean_monthly = terra::global(
#           k_combined_depth_mean,
#           "mean",
#           weights = area_rast,
#           na.rm = TRUE
#         )[1, 1]
#       )
#       
#       mineralised_year_list[[length(mineralised_year_list) + 1]] <-
#         mineralised_this_month
#       
#       bioavailable_year_list[[length(bioavailable_year_list) + 1]] <-
#         bioavailable_this_month
#       
#       k_T_year_list[[length(k_T_year_list) + 1]] <-
#         k_T_depth_mean
#       
#       k_combined_year_list[[length(k_combined_year_list) + 1]] <-
#         k_combined_depth_mean
#       
#       rm(
#         k_T_depths,
#         k_T_depth_mean,
#         k_combined_depths,
#         k_combined_depth_mean,
#         mineralised_this_month,
#         bioavailable_this_month,
#         inorg_rapid_this_month,
#         area_rast
#       )
#       
#       if (exists("mineralisation_result")) {
#         rm(mineralisation_result)
#       }
#       
#       gc()
#     }
#     
#     organic_pool_remaining <- app(organic_pool_depths, sum, na.rm = TRUE)
#     organic_pool_remaining <- terra::mask(organic_pool_remaining, land_mask)
#     
#     yearly_mineralised_N <- app(rast(mineralised_year_list), sum, na.rm = TRUE)
#     yearly_mineralised_N <- terra::mask(yearly_mineralised_N, land_mask)
#     
#     net_organic_pool_change <- new_thawed_organic_N - yearly_mineralised_N
#     net_organic_pool_change <- terra::mask(net_organic_pool_change, land_mask)
#     
#     previous_max_ALD <- terra::ifel(
#       current_ALD > previous_max_ALD,
#       current_ALD,
#       previous_max_ALD
#     )
#     previous_max_ALD <- terra::mask(previous_max_ALD, land_mask)
#     
#     previous_baseline_relative_thawed <- current_baseline_relative_thawed
#     
#     organic_pool_yearly_list[[i]] <- organic_pool_remaining
#     thawed_increment_yearly_list[[i]] <- thawed_increment
#     max_ALD_yearly_list[[i]] <- previous_max_ALD
#     net_organic_pool_change_yearly_list[[i]] <- net_organic_pool_change
#     
#     rm(
#       mineralised_year_list,
#       bioavailable_year_list,
#       k_T_year_list,
#       k_combined_year_list,
#       current_thawed,
#       current_ALD,
#       thawed_increment,
#       net_organic_pool_change,
#       new_thawed_total_N,
#       new_thawed_organic_N,
#       newly_thawed_org_depths,
#       inorg_rapid_yearly,
#       NContentLayer_year,
#       organic_pool_remaining,
#       theta_year,
#       f_theta,
#       k_moisture,
#       current_baseline_relative_thawed
#     )
#     gc()
#   }
#   
#   list(
#     diagnostic_df = dplyr::bind_rows(diagnostic_rows),
#     organic_pool_remaining_yearly = rast(organic_pool_yearly_list),
#     net_organic_pool_change_yearly = rast(net_organic_pool_change_yearly_list),
#     thawed_N_increment_yearly = rast(thawed_increment_yearly_list),
#     new_thawed_organic_N_yearly = rast(new_thawed_organic_N_yearly_list),
#     new_thawed_total_N_yearly = rast(new_thawed_total_N_yearly_list),
#     max_ALD_yearly = rast(max_ALD_yearly_list)
#   )
# }

flux_result <- compute_flux_monthly_depth_pools(
  permafrost_total_thawed = thawed_total,
  ALD = ALD,
  temp_full = temp_full,
  temp_month_dates = temp_month_dates,
  k_base_mineralisation = k_base,
  Ea = Ea_nitrification,
  Ed = Ed_nitrification,
  theta = theta,
  k_T_ref_scalar = k_T_ref_scalar,
  f_theta_ref_scalar = f_theta_ref_scalar,
  start_year = start_year,
  t_opt_C = t_opt_C,
  f_inorg_rapid = f_inorg_rapid,
  n_depths = n_depths,
  average_arctic_thaw_length = average_arctic_thaw_length,
  depth_bounds = depth_bounds,
  params = params,
  profile_depth_max = profile_depth_max,
  LC = LC,
  write_yearly_rasters = FALSE,
  use_temperature_scaling = use_temperature_scaling
)

# ============================================================
# 11. Extract outputs
# ============================================================
diagnostic_df <- flux_result$diagnostic_df

# diagnostic_df %>%
#   filter(Year >= 2000, Year <= 2014) %>%
#   summarise(mean_ref = mean(k_t_mean_monthly, na.rm = TRUE))

new_thawed_organic_N <- prepare_yearly_raster(
  flux_result$new_thawed_organic_N_yearly,
  start_year = start_year
)

net_organic_pool_change <- prepare_yearly_raster(
  flux_result$net_organic_pool_change_yearly,
  start_year = start_year
)

new_thawed_total_N <- prepare_yearly_raster(
  flux_result$new_thawed_total_N_yearly,
  start_year = start_year
)

organic_pool_remaining <- prepare_yearly_raster(
  flux_result$organic_pool_remaining_yearly,
  start_year = start_year
)

thawed_N_increment <- prepare_yearly_raster(
  flux_result$thawed_N_increment_yearly,
  start_year = start_year
)

rapid_bioavailable_pool <- prepare_yearly_raster(
  
  flux_result$rapid_bioavailable_pool_yearly,
  
  start_year = start_year
  
)

max_ALD_yearly <- prepare_yearly_raster(
  flux_result$max_ALD_yearly,
  start_year = start_year
)

# ============================================================
# 12. Diagnostics
# ============================================================

cat("Diagnostics...\n")

area_rast <- cellSize(organic_pool_remaining[[1]], unit = "m")

organic_pool_pg_yearly <- global(
  organic_pool_remaining * area_rast,
  "sum",
  na.rm = TRUE
)[, 1] / 1e12

thawed_increment_pg_yearly <- global(
  thawed_N_increment * area_rast,
  "sum",
  na.rm = TRUE
)[, 1] / 1e12

yearly_extra_df <- data.frame(
  Year = as.integer(format(time(organic_pool_remaining), "%Y")),
  organic_pool_remaining_pg = organic_pool_pg_yearly,
  thawed_increment_pg = thawed_increment_pg_yearly
)

diagnostic_df <- diagnostic_df %>%
  left_join(yearly_extra_df, by = "Year")

print(head(diagnostic_df))
print(summary(diagnostic_df))

write.csv(
  diagnostic_df,
  file.path(out_dir, paste0("region_monthly_diagnostics_w_temp", ssp, "_1850_2100.csv")),
  row.names = FALSE
)
# Monthly bioavailable N plot

diagnostic_df %>%
  filter(Year >= 2000, Year <= 2014) %>%
  summarise(
    mean_k_T = mean(k_T_mean_monthly, na.rm = TRUE),
    mean_k_moisture = mean(k_moisture_mean_monthly, na.rm = TRUE),
    mean_k_combined = mean(k_combined_mean_monthly, na.rm = TRUE)
  )



# ============================================================
# Check: cumulative new_thawed_total_N equals thawed permafrost N
# ============================================================
years <- as.integer(format(time(new_thawed_total_N), "%Y"))
cat("Checking cumulative new_thawed_total_N vs baseline-relative thawed N...\n")

# cumulative sum of yearly new thawed N
cumulative_new_thawed_total_N <- cumsum(new_thawed_total_N)

# baseline-relative thawed N from the model input
baseline_idx <- which(years >= 1850 & years <= 1900)

baseline_thawed_mean_check <- app(
  thawed_total[[baseline_idx]],
  mean,
  na.rm = TRUE
)

permafrost_thawed_N_check <- thawed_total - baseline_thawed_mean_check
permafrost_thawed_N_check <- terra::ifel(
  permafrost_thawed_N_check < 0,
  0,
  permafrost_thawed_N_check
)
permafrost_thawed_N_check <- mask(permafrost_thawed_N_check, land_mask)

# difference
cum_error <- cumulative_new_thawed_total_N - permafrost_thawed_N_check
cum_error <- mask(cum_error, land_mask)

area_rast <- cellSize(cum_error[[1]], unit = "m")

cum_check_df <- data.frame(
  Year = years,
  cumulative_new_thawed_total_pg = global(
    cumulative_new_thawed_total_N * area_rast,
    "sum",
    na.rm = TRUE
  )[, 1] / 1e12,
  permafrost_thawed_pg = global(
    permafrost_thawed_N_check * area_rast,
    "sum",
    na.rm = TRUE
  )[, 1] / 1e12,
  error_pg = global(
    cum_error * area_rast,
    "sum",
    na.rm = TRUE
  )[, 1] / 1e12
)

cum_check_df <- cum_check_df %>%
  mutate(abs_error_pg = abs(error_pg))

print(summary(cum_check_df$error_pg))

write.csv(
  cum_check_df,
  file.path(out_dir, paste0("region_cumulative_new_thawed_check_", ssp, ".csv")),
  row.names = FALSE
)

p_cum_check <- ggplot(cum_check_df, aes(x = Year)) +
  geom_line(aes(y = cumulative_new_thawed_total_pg, color = "Cumulative new thawed N"), linewidth = 1) +
  geom_line(aes(y = permafrost_thawed_pg, color = "Baseline-relative thawed N"), linewidth = 1, linetype = "dashed") +
  theme_bw() +
  labs(
    title = "Check: cumulative new thawed N equals baseline-relative thawed N",
    subtitle = paste0("SSP", ssp),
    x = "Year",
    y = "N pool [Pg N]",
    color = NULL
  )

ggsave(
  file.path(out_dir, paste0("region_cumulative_new_thawed_check_", ssp, ".png")),
  p_cum_check,
  width = 8,
  height = 5,
  dpi = 300
)

p_cum_error <- ggplot(cum_check_df, aes(x = Year, y = error_pg)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 1) +
  theme_bw() +
  labs(
    title = "Cumulative new thawed N conservation error",
    subtitle = "Cumulative new_thawed_total_N - baseline-relative thawed N",
    x = "Year",
    y = "Error [Pg N]"
  )

ggsave(
  file.path(out_dir, paste0("region_cumulative_new_thawed_error_", ssp, ".png")),
  p_cum_error,
  width = 8,
  height = 5,
  dpi = 300
)

# ============================================================
# 13. Write outputs
# ============================================================

cat("Writing outputs...\n")

#out_dir_yearly <- file.path(out_dir, paste0("_w_temp_yearly_nc_", ssp))
#dir.create(out_dir_yearly, recursive = TRUE, showWarnings = FALSE)

years <- 1850:2099
writeRaster(
  new_thawed_organic_N,
  file.path(
    out_dir,
    paste0(
      "region_new_thawed_organic_N_yearly_",
      ssp,
      "_mean_1850_2099_w_temp.nc"
    )
  ),
  overwrite = TRUE,
  filetype = "NetCDF"
)

writeRaster(
  net_organic_pool_change,
  file.path(
    out_dir,
    paste0(
      "region_net_organic_N_pool_change_yearly_",
      ssp,
      "_mean_1850_2099_w_temp.nc"
    )
  ),
  overwrite = TRUE,
  filetype = "NetCDF"
)

writeRaster(
  new_thawed_total_N,
  file.path(
    out_dir,
    paste0(
      "region_new_thawed_total_N_yearly_",
      ssp,
      "_mean_1850_2099_w_temp.nc"
    )
  ),
  overwrite = TRUE,
  filetype = "NetCDF"
)

writeRaster(
  organic_pool_remaining,
  file.path(
    out_dir,
    paste0(
      "region_organic_N_pool_remaining_yearly_",
      ssp,
      "_mean_1850_2099_w_temp.nc"
    )
  ),
  overwrite = TRUE,
  filetype = "NetCDF"
)

writeRaster(
  thawed_N_increment,
  file.path(
    out_dir,
    paste0(
      "region_thawed_N_yearly_increment_",
      ssp,
      "_mean_1850_2099_w_temp.nc"
    )
  ),
  overwrite = TRUE,
  filetype = "NetCDF"
)

writeRaster(
  rapid_bioavailable_pool,
  file.path(
    out_dir_yearly,
    paste0(
      "region_rapid_bioavailable_N_yearly_",
      ssp,
      "_mean_1850_2099_w_temp.nc"
    )
  ),
  overwrite = TRUE,
  filetype = "NetCDF"
)

writeRaster(
  max_ALD_yearly,
  file.path(
    out_dir,
    paste0(
      "region_max_ALD_yearly_",
      ssp,
      "_1850_2099_w_temp.nc"
    )
  ),
  overwrite = TRUE,
  filetype = "NetCDF"
)

cat("Finished writing raster outputs.\n")
cat("Finished SSP", ssp, "\n")
cat("Log file:", log_file, "\n")