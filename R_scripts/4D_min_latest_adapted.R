# ============================================================
# latest fixes: 
# k_baseline = 0.02 = 2% annually -- 0.16 % monthly
# f_inorg = 0.01 = 1% 
# removed rapid organic N ~11 %; only distinction between rapid inorganic fraction = 0.01 - 3 %, and 
# mineralised N and organic N pool; 
#
# rapid inorganic = fraction of what is immediately bioavailable following thawing, 
# between 0.01 - 3% of total soil N stocks (Strauss 2024)

# bioavailable = renamed as total inorganic; total inorg N =
# monthly mineralised organic N + rapid inorganic fraction of newly thawed N

# renamed k_combined to k_env (environment)= k_T * k_moisture

#Monthly mineralised and inorganic N
# Years: 1850–2099
# Dynamic depth-pool version, with soil moisture

# Main logic:
# Thawed permafrost N is calculated relative to the
# 1850–1900 mean baseline at each grid cell.
#
# baseline-relative thawed N =
# total thawed N in current year - mean total thawed N from 1850–1900
##
# Organic N is mineralised monthly using:
# - temperature modifier k_T
# - soil moisture modifier k_moisture
# - base mineralisation rate
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

 # if (is.na(ssp) || length(ssp) == 0) {
 #   ssp <- "585"
 #  }
run <- if (length(args) >= 2) args[2] else "mean"     # "mean", "plus" or "minus"
sd_factor <- c(mean = 0, plus = 1, minus = -1)[[run]]

out_dir <- file.path("monthly_mineralised", paste0(run, "_2perc_baseline"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_dir_yearly <- file.path(out_dir, paste0("yearly_nc_", ssp))
dir.create(out_dir_yearly, recursive = TRUE, showWarnings = FALSE)


log_dir <- "logs"
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(
  log_dir,
  paste0("arctic_mineralisation_w_temp_w_sm_", run, "_", ssp, "_",
         format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
)
 

log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, append = TRUE, type = "message")

start_time <- Sys.time()


cat("====================================================\n")
cat("STARTED SCRIPT\n")
cat("Time:", as.character(start_time), "\n")
cat("SSP:", ssp, "\n")
cat("Run:", run, "\n")          # add this line
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
# mean of 1-3% (Annually, about 1% to 3% of the N present in SOM is mineralized to nitrate (NO3).)

k_base <- 0.02
# this will be converted to a monthly rate in the mineralisation calculation
# 0ld value = 0.002 = 0.2%

Ea_nitrification <- 80000
Ed_nitrification <- 200000
t_opt_C <- 28

#f_inorg_rapid_OLD <- 0.1136
# f_inorg_rapid <- mean of 0.01 - 3 % (Shaver et al 1992, Buckeridge et al 2010)
f_inorg_rapid <-0.01
f_org <- 1 - f_inorg_rapid

#test_years <- 1850:2014
simulation_years <- 1850:2099
start_year <- min(simulation_years)

baseline_start <- 1850
baseline_end <- 1900

T_ref_C <- 25
theta_ref <- 0.65

#ref_start <- 2000
#ref_end <- 2014

profile_depth_max <- 5


# -----------------------------
# 2. File paths
# -----------------------------
# for 60 deg N: 

ald_mean_file    <- paste0("60N/mean_ALD_ssp", ssp, "_60deg.nc")               # adjust
ald_sd_file      <- paste0("60N/std_ALD_ssp",  ssp, "_60deg.nc")
thawed_mean_file <- paste0("60N/arctic_total_thawed_", ssp, "_60N_mean.nc")    # adjust
thawed_sd_file   <- paste0("60N/arctic_total_thawed_", ssp, "_60N_std.nc")

temp_file <- paste0("60N/mean_verttemp_ssp", ssp, "_60deg.nc")

# # soil moisture data
theta_file <- paste0("60N/sm_", ssp, "_60deg.nc")
common_extent <- ext(-179.95, 179.95, 60, 90)


fill_na_focal <- function(r, max_iter = 10, w = 3) {
  # w = window size (3 = 3x3 focal window); increase if gaps are larger
  r_filled <- r
  
  for (iter in seq_len(max_iter)) {
    
    na_before <- global(is.na(r_filled), "sum", na.rm = TRUE)[,1]
    
    if (all(na_before == 0)) {
      cat("No NAs left after", iter - 1, "iterations.\n")
      break
    }
    
    # focal mean of the 8 (or more) neighbors, ignoring NA
    filled_vals <- terra::focal(
      r_filled,
      w = w,
      fun = "mean",
      na.policy = "only",   # only fill cells that are currently NA
      na.rm = TRUE
    )
    
    r_filled <- terra::cover(r_filled, filled_vals)
    
    cat("Iteration", iter, "- NAs remaining per layer:",
        paste(global(is.na(r_filled), "sum", na.rm = TRUE)[,1], collapse = ", "),
        "\n")
  }
  
  r_filled
}
# -----------------------------
# 3. Load data
# -----------------------------

cat("Loading rasters...\n")

ALD          <- crop(rast(ald_mean_file),    common_extent)
thawed_total <- crop(rast(thawed_mean_file), common_extent)

if (sd_factor != 0) {
  
  ALD_sd    <- crop(rast(ald_sd_file),    common_extent)
  thawed_sd <- crop(rast(thawed_sd_file), common_extent)
  
  if (nlyr(ALD_sd) != nlyr(ALD) || nlyr(thawed_sd) != nlyr(thawed_total)) {
    stop("Mean and SD files have different numbers of layers.")
  }
  
  # missing SD = no uncertainty in that cell (keeps the land mask unchanged)
  ALD_sd    <- classify(ALD_sd,    cbind(NA, 0))
  thawed_sd <- classify(thawed_sd, cbind(NA, 0))
  
  ald_time    <- time(ALD)
  thawed_time <- time(thawed_total)
  
  # perturb both inputs in the same direction (no clamping)
  ALD          <- ALD          + sd_factor * ALD_sd
  thawed_total <- thawed_total + sd_factor * thawed_sd
  
  time(ALD)          <- ald_time
  time(thawed_total) <- thawed_time
  
  cat("Inputs perturbed by", sd_factor, "x SD\n")
}


land_mask <- !is.na(thawed_total[[nlyr(thawed_total)]])
plot(land_mask)

temp_full <- rast(temp_file)
plot(temp_full[[2]])

crs(temp_full) <- crs(thawed_total)
#plot(temp_full[[2]])
#temp_full<-mask(temp_full, land_mask)

theta <- rast(theta_file)
theta <- crop(theta, common_extent)

cat("Gap-filling theta before resampling...\n")
theta <- fill_na_focal(theta, max_iter = 20, w = 3)

theta <- resample(theta, thawed_total, method = "bilinear")
crs(theta) <- crs(thawed_total)
plot(theta[[2]])

LC <- rast("LC_remapnn_corr.nc")
LC <- crop(LC, common_extent)
LC <- resample(LC, thawed_total[[1]], method = "near")
plot(LC)
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


frac_capped <- sapply(seq_len(nlyr(ALD)), function(i) {
  global(ALD[[i]] >= profile_depth_max, "mean", na.rm = TRUE)[1,1]
})
plot(start_year:(start_year+length(frac_capped)-1), frac_capped, type = "l",
     xlab = "Year", ylab = "Fraction of land cells with ALD >= profile_depth_max")

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

test_mode <- FALSE

if (test_mode) {
  
  cat("RUNNING SMALL REGION TEST MODE\n")
  
  test_ext <- ext(90, 90.3, 75, 75.3)
  
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
### save workspace
#save.image(file='4D_min_after_readingfiles.RData')

### read workspace previously saved
#load(file='4D_min_after_readingfiles.RData')

# -----------------------------
# 4. Subset thawed_total and ALD to 1850–2014
# -----------------------------

thawed_years <- as.numeric(format(time(thawed_total), "%Y"))

years_to_keep <- which(
  thawed_years >= min(simulation_years) &
    thawed_years <= max(simulation_years)
)

if (length(years_to_keep) == 0) {
  stop("No thawed_total layers found for 1850–2014.")
}

thawed_total <- thawed_total[[years_to_keep]]

ALD_years <- as.numeric(format(time(ALD), "%Y"))

ALD_idx <- which(
  ALD_years >= min(simulation_years) &
    ALD_years <= max(simulation_years)
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

temp_idx_months <- which(temp_years_all %in% simulation_years)

temp_month_dates <- temp_dates_all[temp_idx_months]
n_months <- length(temp_month_dates)

cat("Using", n_depths, "depth layers and", n_months, "months.\n")
print(depth_vals)

#plot(temp_full[[2]])


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
  
  # always fix the bottom boundary at max_depth, regardless of
  # what the midpoint-extrapolation would have given
  depth_bounds[length(depth_bounds)] <- max_depth
  
  # safety: clip any interior bound that somehow exceeds max_depth
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

make_depth_weights_ALD_annual_change <- function(previous_ALD,
                                             current_ALD,
                                             LC,
                                             depth_bounds,
                                             n_depths,
                                             params,
                                             profile_depth_max,
                                             land_mask) {
  
  previous_ALD <- terra::ifel(
    previous_ALD > profile_depth_max,
    profile_depth_max,
    previous_ALD
  )
  
  current_ALD <- terra::ifel(
    current_ALD > profile_depth_max,
    profile_depth_max,
    current_ALD
  )
  
  # Full cumulative depth band between the fixed preindustrial ALD
  # and the current year's ALD.
  band_lo <- terra::ifel(previous_ALD < current_ALD, previous_ALD, current_ALD)  # min
  band_hi <- terra::ifel(previous_ALD > current_ALD, previous_ALD, current_ALD)  # max
  
  taiga_mask <- LC %in% taiga_classes
  tundra_mask <- LC %in% tundra_classes
  wetlands_mask <- LC %in% wetlands_classes
  barren_mask <- LC %in% barren_classes
  
  layer_mass_list <- vector("list", n_depths)
  
  for (d in seq_len(n_depths)) {
    
    z1_layer <- depth_bounds[d]
    z2_layer <- depth_bounds[d + 1]
    
    z1 <- terra::ifel(band_lo > z1_layer, band_lo, z1_layer)
    z2 <- terra::ifel(band_hi < z2_layer, band_hi, z2_layer)
    
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



unclassified <- !(LC %in% c(taiga_classes, tundra_classes, wetlands_classes, barren_classes))
plot(mask(unclassified, land_mask))


#============================================================
 # Moisture reference
#============================================================
moisture_response <- function(
    theta,
    theta_opt = theta_ref
) {
  
  if (inherits(theta, "SpatRaster")) {
    
    # to ensure meaningful values:
    theta <- terra::clamp(
      theta,
      lower = 0,
      upper = 1,
      values = TRUE
    )
    
    response <- terra::ifel(
      theta <= theta_opt,
      theta / theta_opt,
      (1 - theta) / (1 - theta_opt)
    )
    
    response <- terra::clamp(
      response,
      lower = 0,
      upper = 1,
      values = TRUE
    )
    
  } else {
    
    theta <- pmin(pmax(theta, 0), 1)
    
    response <- ifelse(
      theta <= theta_opt,
      theta / theta_opt,
      (1 - theta) / (1 - theta_opt)
    )
  }
  
  response
}

theta_present <- terra::app(
  theta,
  mean,
  na.rm = TRUE
)

theta_present <- terra::mask(
  theta_present,
  land_mask
)

k_moisture_present <- moisture_response(
  theta = theta_present,
  theta_opt = 0.65
)

k_moisture_present <- terra::mask(
  k_moisture_present,
  land_mask
)

k_moisture_present <- terra::ifel(
  is.na(k_moisture_present) & land_mask,
  1,
  k_moisture_present
)

plot(k_moisture_present)
# Before/after comparison
theta_raw <- rast(theta_file)
theta_raw <- crop(theta_raw, common_extent)

cat("NAs before filling:", global(is.na(theta_raw[[1]]), "sum", na.rm=TRUE)[1,1], "\n")
cat("NAs after filling:", global(is.na(theta[[1]]), "sum", na.rm=TRUE)[1,1], "\n")

# Check specifically over land_mask after the full pipeline
k_moisture_gaps <- is.na(k_moisture_present) & land_mask
plot(k_moisture_gaps, main = "Remaining k_moisture NAs over land")





#######



get_temp_layers_from_full <- function(month_index_selected,
                                      temp_idx_months,
                                      n_depths_full,
                                      depth_keep) {

  month_index_full <- temp_idx_months[month_index_selected]

  ((month_index_full - 1) * n_depths_full) + depth_keep
}



# ============================================================
# Old: 8. k_T reference scalar
# ============================================================
# compute_k_T_ref_scalar_surface_summer <- function(temp_full,
#                                                   temp_month_dates,
#                                                   temp_idx_months,
#                                                   n_depths_full,
#                                                   depth_keep,
#                                                   Ea,
#                                                   Ed,
#                                                   t_opt_C = 28,
#                                                   ref_start = 2000,
#                                                   ref_end = 2014,
#                                                   ref_months_use = 6:8,
#                                                   land_mask,
#                                                   use_fixed_T_ref = TRUE,
#                                                   fixed_T_ref_C = 8) {
# 
#   if (use_fixed_T_ref) {
# 
#     T_ref_C <- fixed_T_ref_C
# 
#   } else {
# 
#     temp_years  <- as.integer(format(temp_month_dates, "%Y"))
#     temp_months <- as.integer(format(temp_month_dates, "%m"))
# 
#     ref_months <- which(
#       temp_years >= ref_start &
#         temp_years <= ref_end &
#         temp_months %in% ref_months_use
#     )
# 
#     if (length(ref_months) == 0) {
#       stop("No reference temperature months found.")
#     }
# 
#     area_rast <- terra::cellSize(temp_full[[1]], unit = "m")
#     area_rast <- terra::mask(area_rast, land_mask)
# 
#     # get surface layer
#     surface_depth_keep <- depth_keep[1]
# 
#     ref_temps <- numeric(length(ref_months))
# 
#     for (ii in seq_along(ref_months)) {
# 
#       j <- ref_months[ii]
# 
#       cat(
#         "Surface summer reference month",
#         ii, "of", length(ref_months),
#         "date:", as.character(temp_month_dates[j]), "\n"
#       )
# 
#       surface_layer <- get_temp_layers_from_full(
#         month_index_selected = j,
#         temp_idx_months = temp_idx_months,
#         n_depths_full = n_depths_full,
#         depth_keep = surface_depth_keep
#       )
# 
#       temp_surface_C <- temp_full[[surface_layer]] - 273.15
#       temp_surface_C <- terra::mask(temp_surface_C, land_mask)
# 
#       ref_temps[ii] <- terra::global(
#         temp_surface_C,
#         "mean",
#         weights = area_rast,
#         na.rm = TRUE
#       )[1, 1]
# 
#       rm(temp_surface_C)
#       gc()
#     }
# 
#     T_ref_C <- mean(ref_temps, na.rm = TRUE)
#   }
# 
#   k_T_ref_scalar <- peaked_arrhenius(
#     temp_C = T_ref_C,
#     Ea = Ea,
#     Ed = Ed,
#     t_opt_C = t_opt_C
#   )
# 
#   cat("Surface summer T_ref_C:", T_ref_C, "\n")
#   cat("Surface summer k_T_ref_scalar:", k_T_ref_scalar, "\n")
# 
#   if (!is.finite(k_T_ref_scalar) || k_T_ref_scalar <= 0) {
#     stop("Surface-summer k_T_ref_scalar is invalid.")
#   }
# 
#   k_T_ref_scalar
# }


###############################################################################
###  helper functions:

use_temperature_scaling <- TRUE


compute_k_factor_depths <- function(
    temp_full,
    month_index_selected,
    temp_idx_months,
    n_depths_full,
    depth_keep,
    Ea,
    Ed,
    t_opt_C,
    T_ref_C = 25,
    land_mask,
    use_temperature_scaling = TRUE
) {
  
  depth_layers <- get_temp_layers_from_full(
    month_index_selected = month_index_selected,
    temp_idx_months = temp_idx_months,
    n_depths_full = n_depths_full,
    depth_keep = depth_keep
  )
  
  temp_month_depths <- temp_full[[depth_layers]] - 273.15
  
  # No temperature scaling: activity is 1 in thawed soil
  # and 0 in frozen soil.
  if (!use_temperature_scaling) {
    
    k_T_depths <- terra::ifel(
      is.na(temp_month_depths),
      NA,
      terra::ifel(temp_month_depths > 0, 1, 0)
    )
    
    k_T_depths <- terra::mask(
      k_T_depths,
      land_mask
    )
    
    names(k_T_depths) <- paste0(
      "depth_",
      seq_along(depth_keep)
    )
    
    return(k_T_depths)
  }
  
  # Temperature response at each grid cell and depth
  k_T_raw <- peaked_arrhenius(
    temp_C = temp_month_depths,
    Ea = Ea,
    Ed = Ed,
    t_opt_C = t_opt_C
  )
  
  # Temperature response under the conditions represented by k_base
  k_T_reference <- peaked_arrhenius(
    temp_C = T_ref_C,
    Ea = Ea,
    Ed = Ed,
    t_opt_C = t_opt_C
  )
  
  # Normalize so that k_T = 1 at 25°C
  k_T_depths <- k_T_raw / k_T_reference
  
  # No mineralisation in frozen soil
  k_T_depths <- terra::ifel(
    temp_month_depths > 0,
    k_T_depths,
    0
  )
  
  k_T_depths <- terra::ifel(
    is.na(temp_month_depths),
    NA,
    k_T_depths
  )
  
  k_T_depths <- terra::clamp(
    k_T_depths,
    lower = 0,
    upper = Inf,
    values = TRUE
  )
  
  k_T_depths <- terra::mask(
    k_T_depths,
    land_mask
  )
  
  names(k_T_depths) <- paste0(
    "depth_",
    seq_along(depth_keep)
  )
  
  k_T_depths
}


# with temperature old:
# compute_k_factor_depths <- function(
#     temp_full,
#     month_index_selected,
#     temp_idx_months,
#     n_depths_full,
#     depth_keep,
#     Ea,
#     Ed,
#     t_opt_C,
#     #k_T_ref_scalar,
#     land_mask,
#     use_temperature_scaling = TRUE
# ) {
#   
#   depth_layers <- get_temp_layers_from_full(
#     month_index_selected = month_index_selected,
#     temp_idx_months = temp_idx_months,
#     n_depths_full = n_depths_full,
#     depth_keep = depth_keep
#   )
#   
#   # Monthly soil temperature at every selected depth
#   temp_month_depths <- temp_full[[depth_layers]] - 273.15
#   
#   # ----------------------------------------------------------
#   # No temperature scaling:
#   # kT = 1 where soil is thawed
#   # kT = 0 where soil is frozen
#   # ----------------------------------------------------------
#   if (!use_temperature_scaling) {
#     
#     k_factor_depths <- terra::ifel(
#       is.na(temp_month_depths),
#       NA,
#       terra::ifel(temp_month_depths > 0, 1, 0)
#     )
#     
#     k_factor_depths <- terra::mask(
#       k_factor_depths,
#       land_mask
#     )
#     
#     names(k_factor_depths) <- paste0(
#       "depth_",
#       seq_along(depth_keep)
#     )
#     
#     return(k_factor_depths)
#   }
#   
#   # ----------------------------------------------------------
#   # Temperature-scaled version
#   # ----------------------------------------------------------
#   k_T_response_depths <- peaked_arrhenius(
#     temp_month_depths,
#     Ea = Ea,
#     Ed = Ed,
#     t_opt_C = t_opt_C
#   )
#   
#   # No activity in frozen soil
#   k_T_response_depths <- terra::ifel(
#     temp_month_depths > 0,
#     k_T_response_depths,
#     0
#   )
#   
#   k_T_response_depths <- terra::ifel(
#     is.na(temp_month_depths),
#     NA,
#     k_T_response_depths
#   )
#   
#   # Normalize against reference temperature
#   k_factor_depths <- k_T_response_depths / k_T_ref_scalar
#   
#   k_factor_depths <- terra::mask(
#     k_factor_depths,
#     land_mask
#   )
#   
#   names(k_factor_depths) <- paste0(
#     "depth_",
#     seq_along(depth_keep)
#   )
#   
#   k_factor_depths
# }


mineralise_depth_pools <- function(
    organic_pool_depths,
    k_factor_depths,
    base_mineralisation_rate_monthly
) {
  k_depths <- base_mineralisation_rate_monthly * k_factor_depths
  
  mineralised_depths_signed <- organic_pool_depths * k_depths
  
  organic_pool_depths <- organic_pool_depths - mineralised_depths_signed
  
  mineralised_this_month <- terra::app(
    mineralised_depths_signed,
    sum,
    na.rm = TRUE
  )
  
  list(
    mineralised_this_month = mineralised_this_month,
    organic_pool_depths = organic_pool_depths,
    # NEW: un-summed, depth-resolved mineralised N stack (n_depths layers)
    # kept separate from mineralised_this_month so existing callers are unaffected
    mineralised_depths_signed = mineralised_depths_signed
  )
}


# k_T_ref_scalar <- compute_k_T_ref_scalar_surface_summer(
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
#   ref_months_use = 6:8,
#   land_mask = land_mask,
#   use_fixed_T_ref = TRUE,
#   fixed_T_ref_C = 8
# )


# ============================================================
# 9. Main calculation


compute_flux_monthly_depth_pools <- function(
    permafrost_total_thawed,
    ALD,
    temp_full,
    temp_month_dates,
    k_base_mineralisation,
    Ea,
    Ed,
    k_moisture_present,
    T_ref_C,
    #k_T_ref_scalar,
    start_year = 1850,
    t_opt_C = 28,
    f_inorg_rapid = 0.01,
    n_depths,
    depth_bounds,
    params,
    profile_depth_max,
    LC,
    write_yearly_rasters = TRUE,
    use_temperature_scaling = TRUE) {
  
  f_org <- 1 - f_inorg_rapid
  n_years <- nlyr(permafrost_total_thawed)
  years <- start_year:(start_year + n_years - 1)
  temp_years <- as.integer(format(temp_month_dates, "%Y"))
  
  # Common land mask setup
  land_mask <- terra::ifel(
    !is.na(permafrost_total_thawed[[nlyr(permafrost_total_thawed)]]),
    1,
    NA
  )
  
  if (nlyr(ALD) != n_years) {
    stop("ALD and permafrost_total_thawed do not have the same number of layers.")
  }
  
  permafrost_total_thawed <- terra::mask(permafrost_total_thawed, land_mask)
  ALD <- terra::mask(ALD, land_mask)
  
  # Baseline calculations (1850 - 1900)
  baseline_idx <- which(years >= 1850 & years <= 1900)
  baseline_thawed_mean <- terra::mask(
    app(permafrost_total_thawed[[baseline_idx]], mean, na.rm = TRUE),
    land_mask
  )
  
  # Pre-allocate diagnostic and raster tracking lists
  diagnostic_rows <- list()
  depth_distribution_check_rows <- vector("list", n_years)
  
  organic_pool_yearly_list         <- vector("list", n_years)
  thawed_change_yearly_list        <- vector("list", n_years)
  thawed_permafrost_N_yearly_list  <- vector("list", n_years)
  new_thawed_organic_N_yearly_list <- vector("list", n_years)
  new_thawed_total_N_yearly_list   <- vector("list", n_years)
  rapid_inorg_pool_yearly_list <- vector("list", n_years)
  newly_thawed_depth_pg_rows <- vector("list", n_years)
  #mineralised_N_yearly_list  <- vector("list", n_years)
  #total_inorg_N_yearly_list  <- vector("list", n_years)
  
  base_mineralisation_rate_monthly <- k_base_mineralisation / 12
  print("base mineralization")
  print(base_mineralisation_rate_monthly)
  
  zero_layer <- terra::mask(permafrost_total_thawed[[1]] * 0, land_mask)
  na_layer   <- zero_layer * NA_real_
  
  # Year 1 Initialization
  first_thawed_total_N <- permafrost_total_thawed[[1]]
  first_baseline_relative_thawed_N <- first_thawed_total_N - baseline_thawed_mean
  previous_baseline_relative_thawed <- first_baseline_relative_thawed_N
  
  preindustrial_ALD <- app(ALD[[baseline_idx]], mean, na.rm = TRUE)
  
  # try: set organic pool year 1 to organic pool relative to preindustrial
  initial_ALD_lower <- ALD[[1]] * 0
  
  initial_depth_weights <- make_depth_weights_ALD_annual_change(
    previous_ALD = initial_ALD_lower,
    current_ALD = ALD[[1]],
    LC = LC,
    depth_bounds = depth_bounds,
    n_depths = n_depths,
    params = params,
    profile_depth_max = profile_depth_max,
    land_mask = land_mask
  )
  
  layer_sums <- app(initial_depth_weights, sum, na.rm = TRUE)
  initial_depth_weights <- ifel(layer_sums > 0,
                                initial_depth_weights / layer_sums,
                                0)
  
  organic_pool_depths <-
    (first_baseline_relative_thawed_N * f_org) * initial_depth_weights
  
  names(organic_pool_depths) <- paste0("depth_", seq_len(n_depths))
  
  
  #organic_pool_depths <- rast(replicate(n_depths, zero_layer))
  #names(organic_pool_depths) <- paste0("depth_", seq_len(n_depths))
  
  # n_years=50 --> only for debugging
  # Main Yearly Loop
  for (i in seq_len(n_years)) {
    yr <- years[i]
    cat("Processing year", yr, "(", i, "of", n_years, ")\n")
    
    mineralised_year_list  <- vector("list", 12)
    total_inorg_year_list <- vector("list", 12)
    k_T_year_list          <- vector("list", 12)
    k_env_year_list   <- vector("list", 12)
    mineralised_depths_year_list <- vector("list", 12)   # NEW: depth-resolved, n_depths layers per month
    
    current_thawed <- permafrost_total_thawed[[i]]
    current_baseline_relative_thawed <- current_thawed - baseline_thawed_mean
    current_ALD <- ALD[[i]]
    
    thawed_permafrost_N_yearly_list[[i]] <- current_baseline_relative_thawed

    if (i == 1L) {
      new_thawed_total_N     <- current_baseline_relative_thawed
      new_thawed_organic_N   <- current_baseline_relative_thawed * f_org
      #thawed_N_annual_change <- na_layer
      thawed_N_annual_change <- current_baseline_relative_thawed
      previous_ALD = ALD[[1]]
    }
     else {
      thawed_N_annual_change <- current_baseline_relative_thawed - previous_baseline_relative_thawed
      new_thawed_total_N     <- thawed_N_annual_change 
      new_thawed_organic_N   <- new_thawed_total_N * f_org
    }
    #print("new thawed total")
    #print(new_thawed_total_N[1,1,])
    
    # Weight normalisation across depths
    NContentLayer_year <- make_depth_weights_ALD_annual_change(
      previous_ALD = previous_ALD,
      current_ALD = current_ALD,
      LC = LC,
      depth_bounds = depth_bounds,
      n_depths = n_depths,
      params = params,
      profile_depth_max = profile_depth_max,
      land_mask = land_mask
    )
    
    layer_sums <- terra::app(NContentLayer_year, sum, na.rm = TRUE)
    #print("lay_sum")
    #print(layer_sums[1,1,])
    NContentLayer_year <- terra::ifel(layer_sums > 0, NContentLayer_year / layer_sums, 0)
    #print("NContentLayer_year")
    #print(NContentLayer_year[1,1,])
    
    new_thawed_total_N_yearly_list[[i]]   <- new_thawed_total_N
    new_thawed_organic_N_yearly_list[[i]] <- new_thawed_organic_N
    
    newly_thawed_org_depths <- new_thawed_organic_N * NContentLayer_year
    
    #### ------------------------------------------------------------------------ 
    ## Depth distribution diagnostics
    distributed_organic_N_sum <- terra::app(newly_thawed_org_depths, sum, na.rm = TRUE)
    distribution_difference    <- distributed_organic_N_sum - new_thawed_organic_N
    distribution_abs_difference <- abs(distribution_difference)
    
    distribution_relative_difference <- terra::ifel(
      abs(new_thawed_organic_N) > 1e-12,
      distribution_difference / new_thawed_organic_N,
      NA
    )
    
    distribution_area_rast <- terra::cellSize(new_thawed_organic_N, unit = "m")
    
    organic_before_distribution_pg <- terra::global(new_thawed_organic_N * distribution_area_rast, "sum", na.rm = TRUE)[1, 1] / 1e12
    organic_after_distribution_pg  <- terra::global(distributed_organic_N_sum * distribution_area_rast, "sum", na.rm = TRUE)[1, 1] / 1e12
    organic_distribution_difference_pg <- organic_after_distribution_pg - organic_before_distribution_pg
    
    organic_distribution_relative_error <- if (abs(organic_before_distribution_pg) > 1e-15) {
      organic_distribution_difference_pg / organic_before_distribution_pg
    } else {
      NA_real_
    }
    
    weight_sum_check <- terra::app(NContentLayer_year, sum, na.rm = TRUE)
    weight_sum_on_input <- terra::mask(weight_sum_check, abs(new_thawed_organic_N) > 1e-12)
    weight_sum_range <- terra::global(weight_sum_on_input, "range", na.rm = TRUE)
    
    depth_distribution_check_rows[[i]] <- data.frame(
      Year = yr,
      organic_before_distribution_pg = organic_before_distribution_pg,
      organic_after_distribution_pg  = organic_after_distribution_pg,
      difference_pg                  = organic_distribution_difference_pg,
      relative_error                 = organic_distribution_relative_error,
      mean_abs_cell_difference_kg_m2 = terra::global(distribution_abs_difference, "mean", na.rm = TRUE)[1, 1]
  
    )
    # check end
    ## ------------------------------------------------------------------------------------------
    
    
    newly_thawed_depth_pg_year <- terra::global(
      newly_thawed_org_depths * distribution_area_rast,
      "sum",
      na.rm = TRUE
    )[, 1] / 1e12
    
    newly_thawed_depth_pg_rows[[i]] <- data.frame(
      Year = yr,
      depth_index = seq_len(n_depths),
      depth_mid_m = depth_vals,
      depth_thickness_m = depth_thickness,
      newly_thawed_organic_N_pg = newly_thawed_depth_pg_year
    )
    
    # Pool addition
    organic_pool_depths <- organic_pool_depths + newly_thawed_org_depths
    inorg_rapid_yearly  <- new_thawed_total_N * f_inorg_rapid
    rapid_inorg_pool_yearly_list[[i]] <- inorg_rapid_yearly
    
    # Environmental factors
    k_moisture <- k_moisture_present
    #print(k_moisture[1,1,])
    
    month_idx_year <- which(temp_years == yr)
    if (length(month_idx_year) == 0) {
      stop(paste("No temperature months found for year", yr))
    }
    
    month_counter <- 1
    
    # Monthly Processing Loop
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
        T_ref_C = T_ref_C,
        #k_T_ref_scalar = k_T_ref_scalar,
        land_mask = land_mask,
        use_temperature_scaling = use_temperature_scaling
      )
      #print("Kt depths")
      #print(k_T_depths[1,1,])
      
      k_T_depth_mean   <- weighted_depth_mean(k_T_depths, depth_weights)
      #organic_pool_sum <- app(organic_pool_depths, sum, na.rm = TRUE)
      organic_pool_abs     <- abs(organic_pool_depths)
      organic_pool_abs_sum <- terra::app(organic_pool_abs, sum, na.rm = TRUE)
      
      #k_T_pool_weighted <- app(k_T_depths * organic_pool_depths, sum, na.rm = TRUE) / organic_pool_sum
      k_T_pool_weighted <- terra::app(k_T_depths * organic_pool_abs, sum, na.rm = TRUE) / organic_pool_abs_sum
      k_T_pool_weighted <- terra::ifel(organic_pool_abs_sum > 0, k_T_pool_weighted, NA)
      
      k_env_depths <- k_T_depths * k_moisture
      k_env_depths <- terra::clamp(
        k_env_depths,
        lower = 0,
        upper = Inf,
        values = TRUE
      )
      
      #k_env_pool_weighted <- app(k_env_depths * organic_pool_depths, sum, na.rm = TRUE) / organic_pool_sum
      k_env_pool_weighted <- terra::app(k_env_depths * organic_pool_abs, sum, na.rm = TRUE) / organic_pool_abs_sum
      k_env_pool_weighted <- terra::ifel(organic_pool_abs_sum > 0, k_env_pool_weighted, NA)

      k_env_depth_mean <- weighted_depth_mean(k_env_depths, depth_weights)
      
      mineralisation_result  <- mineralise_depth_pools(
        organic_pool_depths = organic_pool_depths,
        k_factor_depths = k_env_depths,
        base_mineralisation_rate_monthly = base_mineralisation_rate_monthly
      )
      #print("mineralisation")
      #print(mineralisation_result$mineralised_this_month[1,1,])
      
      mineralised_this_month <- mineralisation_result$mineralised_this_month
      organic_pool_depths    <- mineralisation_result$organic_pool_depths
      mineralised_depths_signed_this_month <- mineralisation_result$mineralised_depths_signed  # NEW
      
      inorg_rapid_this_month <- if (month_number %in% 6:8) inorg_rapid_yearly / 3 else zero_layer
      total_inorg_this_month <- mineralised_this_month + inorg_rapid_this_month
      
      area_rast <- terra::cellSize(k_T_depth_mean, unit = "m")
      # ------------------------------------------------------------------
      # Collapse k_env_depths across depth BEFORE calling global() on it.
      # Without this, global() returns one row per depth layer, and [1,1]
      # silently grabs only the shallowest layer (depth_1 ~ 0.01 m) instead
      # of a genuine depth-unweighted mean.
      # ------------------------------------------------------------------
      k_env_depth_unweighted <- terra::app(k_env_depths, mean, na.rm = TRUE)
      
      diagnostic_rows[[length(diagnostic_rows) + 1]] <- data.frame(
        Date = temp_month_dates[j],
        Year = yr,
        Month = month_number,
        mineralised_pg_monthly   = terra::global(mineralised_this_month * area_rast, "sum", na.rm = TRUE)[1, 1] / 1e12,
        total_inorg_pg_monthly   = terra::global(total_inorg_this_month * area_rast, "sum", na.rm = TRUE)[1, 1] / 1e12,
        inorg_rapid_available    = terra::global(inorg_rapid_this_month * area_rast, "sum", na.rm = TRUE)[1, 1] / 1e12,
        k_T_mean_monthly         = terra::global(k_T_depth_mean, "mean", weights = area_rast, na.rm = TRUE)[1, 1],
        k_moisture_mean_monthly  = terra::global(k_moisture, "mean", weights = area_rast, na.rm = TRUE)[1, 1],
        k_env_pool_weighted_mean = terra::global(k_env_pool_weighted, "mean", weights = area_rast, na.rm = TRUE)[1, 1],
        k_T_pool_weighted_mean   = terra::global(k_T_pool_weighted, "mean", weights = area_rast, na.rm = TRUE)[1, 1],
        k_env_depth_unweighted_mean = terra::global(k_env_depth_unweighted, "mean", weights = area_rast, na.rm = TRUE)[1, 1]
      )
      
      mineralised_year_list[[month_counter]]  <- mineralised_this_month
      total_inorg_year_list[[month_counter]] <- total_inorg_this_month
      k_T_year_list[[month_counter]]          <- k_T_depth_mean
      k_env_year_list[[month_counter]]   <- weighted_depth_mean(k_env_depths, depth_weights)
      mineralised_depths_year_list[[month_counter]] <- mineralised_depths_signed_this_month  # NEW
      
      month_counter <- month_counter + 1
    }
    
    # Optional Raster Exporting
    if (write_yearly_rasters) {
      write_nc_layer <- function(layer_list, name_suffix) {
        r_stack <- rast(layer_list)
        terra::time(r_stack) <- seq(as.Date(paste0(yr, "-01-16")), by = "month", length.out = 12)
        names(r_stack) <- month.abb
        terra::writeRaster(
          r_stack,
          file.path(out_dir_yearly, paste0("arctic_", name_suffix, "_", ssp, "_", yr, ".nc")),
          overwrite = TRUE,
          filetype = "NetCDF"
        )
      }
      
      # NEW: depth-resolved writer.
      # Each element of depth_month_list is a SpatRaster with n_depths layers
      # (one month's worth of per-depth mineralised N). Stacking all 12 elements
      # gives n_depths * 12 layers total, ordered month-major:
      # month1_depth1, month1_depth2, ..., month1_depthN, month2_depth1, ...
      # NetCDF has no native 4th (depth) dimension here, so depth and month are
      # both encoded in the layer name instead; recover them by parsing names().
      # write_nc_layer_depth_resolved <- function(depth_month_list, name_suffix) {
      #   r_stack <- rast(depth_month_list)
      #   layer_names <- unlist(lapply(seq_len(12), function(m) {
      #     paste0("month_", sprintf("%02d", m), "_depth_", seq_len(n_depths))
      #   }))
      #   names(r_stack) <- layer_names
      #   terra::writeRaster(
      #     r_stack,
      #     file.path(out_dir_yearly, paste0("arctic_", name_suffix, "_", ssp, "_", yr, ".nc")),
      #     overwrite = TRUE,
      #     filetype = "NetCDF",
      #     datatype = "FLT4S",
      #     gdal = c("FORMAT=NC4", "COMPRESS=DEFLATE", "ZLEVEL=6")
      #   )
      # }
      
      write_nc_layer(mineralised_year_list,  "mineralised_N_monthly_w_temp_sm")
      write_nc_layer(total_inorg_year_list, "total_inorg_N_monthly_w_temp_sm")
      write_nc_layer(k_T_year_list,          "k_T_temperature_only_monthly")
      write_nc_layer(k_env_year_list,   "k_env_monthly_w_temp_sm")
    }
    #   # NEW: write the depth-resolved mineralised N (n_depths x 12 layers)
    #   write_nc_layer_depth_resolved(
    #     mineralised_depths_year_list,
    #     "mineralised_N_by_depth_monthly_w_temp_sm"
    #   )
    # }
    
    # Annual Summary (now safely bounded to exactly 12 layers)
    organic_pool_remaining <- app(organic_pool_depths, sum, na.rm = TRUE)
    yearly_mineralised_N   <- app(rast(mineralised_year_list), sum, na.rm = TRUE)
    yearly_total_inorg_N <- app(rast(total_inorg_year_list), sum, na.rm = TRUE)
    
        #print("organic_pool_remaining")
    #print(organic_pool_remaining[1,1,])   
    #print("yearly_mineralised_N")
    #print(yearly_mineralised_N[1,1,])
    
    previous_baseline_relative_thawed <- current_baseline_relative_thawed
    previous_ALD = ALD[[i]]
    
    organic_pool_yearly_list[[i]]             <- organic_pool_remaining
    thawed_change_yearly_list[[i]]            <- thawed_N_annual_change
    #mineralised_N_yearly_list[[i]] <- yearly_mineralised_N   # kg N m-2 yr-1, depth-integrated
    #total_inorg_N_yearly_list[[i]] <- yearly_total_inorg_N   # kg N m-2 yr-1
    
    gc()
  }
  
  # Return Structured Outputs
  list(
    diagnostic_df                  = dplyr::bind_rows(diagnostic_rows),
    depth_distribution_check_df   = dplyr::bind_rows(depth_distribution_check_rows),
    thawed_permafrost_N_yearly    = rast(thawed_permafrost_N_yearly_list),
    organic_pool_remaining_yearly = rast(organic_pool_yearly_list),
    thawed_change_yearly          = rast(thawed_change_yearly_list),
    new_thawed_organic_N_yearly   = rast(new_thawed_organic_N_yearly_list),
    new_thawed_total_N_yearly     = rast(new_thawed_total_N_yearly_list),
    rapid_inorg_pool_yearly= rast(rapid_inorg_pool_yearly_list),
    newly_thawed_depth_pg_df = dplyr::bind_rows(newly_thawed_depth_pg_rows) 
    #mineralised_N_yearly = rast(mineralised_N_yearly_list),
    #total_inorg_N_yearly = rast(total_inorg_N_yearly_list)
  )
  
}


flux_result <- compute_flux_monthly_depth_pools(
  permafrost_total_thawed = thawed_total,
  ALD = ALD,
  temp_full = temp_full,
  temp_month_dates = temp_month_dates,
  k_base_mineralisation = k_base,
  Ea = Ea_nitrification,
  Ed = Ed_nitrification,
  k_moisture_present = k_moisture_present,
  T_ref_C = T_ref_C,
  #k_T_ref_scalar = k_T_ref_scalar,
  start_year = start_year,
  t_opt_C = t_opt_C,
  f_inorg_rapid = f_inorg_rapid,
  n_depths = n_depths,
  depth_bounds = depth_bounds,
  params = params,
  profile_depth_max = profile_depth_max,
  LC = LC,
  write_yearly_rasters = TRUE,
  use_temperature_scaling = use_temperature_scaling
)

# ============================================================
# 11. Extract outputs
# ============================================================
diagnostic_df <- flux_result$diagnostic_df

# mineralised_N_yearly <- prepare_yearly_raster(flux_result$mineralised_N_yearly, 
#                                               start_year = start_year)
# total_inorg_N_yearly <- prepare_yearly_raster(flux_result$total_inorg_N_yearly, 
#                                               start_year = start_year)

thawed_permafrost_N <- prepare_yearly_raster(
  flux_result$thawed_permafrost_N_yearly,
  start_year = start_year
)

new_thawed_organic_N <- prepare_yearly_raster(
  flux_result$new_thawed_organic_N_yearly,
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

thawed_N_annual_change <- prepare_yearly_raster(
  flux_result$thawed_change_yearly,
  start_year = start_year
)

rapid_inorg_pool <- prepare_yearly_raster(
  flux_result$rapid_inorg_pool_yearly,
  start_year = start_year
)


depth_distribution_check_df <-
  flux_result$depth_distribution_check_df

print(depth_distribution_check_df)

write.csv(
  depth_distribution_check_df,
  file.path(
    out_dir,
    paste0(
      "depth_distribution_conservation_check_",
      ssp,
      "_1850_2099.csv"
    )
  ),
  row.names = FALSE
)


newly_thawed_depth_org_df <-
  flux_result$newly_thawed_depth_pg_df

print(newly_thawed_depth_org_df)

write.csv(
  newly_thawed_depth_org_df,
  file.path(
    out_dir,
    paste0(
      "depth_distribution_organic_pool_",
      ssp,
      "_1850_2099.csv"
    )
  ),
  row.names = FALSE
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

thawed_change_pg_yearly <- global(
  thawed_N_annual_change * area_rast,
  "sum",
  na.rm = TRUE
)[, 1] / 1e12

yearly_extra_df <- data.frame(
  Year = as.integer(format(time(organic_pool_remaining), "%Y")),
  organic_pool_remaining_pg = organic_pool_pg_yearly,
  thawed_N_annual_change_pg = thawed_change_pg_yearly
)

diagnostic_df <- diagnostic_df %>%
  left_join(yearly_extra_df, by = "Year")

print(head(diagnostic_df))
print(summary(diagnostic_df))

write.csv(
  diagnostic_df,
  file.path(out_dir, paste0("arctic_monthly_diagnostics_w_temp", ssp, "_1850_2100.csv")),
  row.names = FALSE
)
# Monthly total inorg N plot

present_k_summary <- diagnostic_df %>%
  filter(Year >= 2000, Year <= 2014) %>%
  summarise(
    mean_k_T        = mean(k_T_mean_monthly, na.rm = TRUE),
    mean_k_moisture = mean(k_moisture_mean_monthly, na.rm = TRUE),
    mean_k_env      = mean(k_env_depth_unweighted_mean, na.rm = TRUE)
  )

print(present_k_summary)


# ============================================================
# 13. Write outputs
# ============================================================

cat("Writing outputs...\n")

out_dir_yearly <- file.path(out_dir, paste0("yearly_nc_", ssp))
dir.create(out_dir_yearly, recursive = TRUE, showWarnings = FALSE)

years <- 1850:2099

# for (nm in c("mineralised_N", "total_inorg_N")) {
#   writeRaster(
#     get(paste0(nm, "_yearly")),
#     file.path(out_dir, paste0("arctic_", nm, "_yearly_", ssp, "_mean_1850_2099_w_temp.nc")),
#     overwrite = TRUE, filetype = "NetCDF", datatype = "FLT4S",
#     gdal = c("FORMAT=NC4", "COMPRESS=DEFLATE", "ZLEVEL=6")
#   )
# }

writeRaster(
  thawed_permafrost_N,
  file.path(
    out_dir,
    paste0(
      "arctic_permafrost_thawed_total_N_",
      ssp,
      "_", run, "_1850_2099_w_temp.nc"
    )
  ),
  overwrite = TRUE,
  filetype = "NetCDF",
  datatype = "FLT4S",
  gdal = c("FORMAT=NC4", "COMPRESS=DEFLATE", "ZLEVEL=6")
)

writeRaster(
  new_thawed_organic_N,
  file.path(
    out_dir,
    paste0(
      "arctic_new_thawed_organic_N_yearly_",
      ssp,
      "_", run, "_1850_2099_w_temp.nc"
    )
  ),
  overwrite = TRUE,
  filetype = "NetCDF",
  datatype = "FLT4S",
  gdal = c("FORMAT=NC4", "COMPRESS=DEFLATE", "ZLEVEL=6")
)



writeRaster(
  new_thawed_total_N,
  file.path(
    out_dir,
    paste0(
      "arctic_new_thawed_total_N_yearly_",
      ssp,
      "_", run, "_1850_2099_w_temp.nc"
    )
  ),
  overwrite = TRUE,
  filetype = "NetCDF",
  datatype = "FLT4S",
  gdal = c("FORMAT=NC4", "COMPRESS=DEFLATE", "ZLEVEL=6")
)

writeRaster(
  organic_pool_remaining,
  file.path(
    out_dir,
    paste0(
      "arctic_organic_N_pool_remaining_yearly_",
      ssp,
      "_", run, "_1850_2099_w_temp.nc"
    )
  ),
  overwrite = TRUE,
  filetype = "NetCDF",
  datatype = "FLT4S",
  gdal = c("FORMAT=NC4", "COMPRESS=DEFLATE", "ZLEVEL=6"))
# 
# writeRaster(
#   thawed_N_annual_change,
#   file.path(
#     out_dir,
#     paste0(
#       "arctic_thawed_N_annual_change_yearly_",
#       ssp,
#       "_mean_1850_2099_w_temp.nc"
#     )
#   ),
#   overwrite = TRUE,
#   filetype = "NetCDF",
#   datatype = "FLT4S",
#   gdal = c("FORMAT=NC4", "COMPRESS=DEFLATE", "ZLEVEL=6")
# )

writeRaster(
  rapid_inorg_pool,
  file.path(
    out_dir,
    paste0(
      "arctic_rapid_inorg_N_yearly_",
      ssp,
      "_", run, "_1850_2099_w_temp.nc"
    )
  ),
  overwrite = TRUE,
  filetype = "NetCDF",
  datatype = "FLT4S",
  gdal = c("FORMAT=NC4", "COMPRESS=DEFLATE", "ZLEVEL=6")
)


cat("Finished writing raster outputs.\n")
cat("Finished SSP", ssp, "\n")
cat("Log file:", log_file, "\n")

