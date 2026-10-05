terra::gdal(drivers = TRUE)

# --- 0️⃣ Setup ---
library(terra)
library(ggplot2)
library(tidyr)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
ssp <- args[1]  # SSP passed from SLURM

# in case of just one ssp:
ssp <- "126"

# --- 1️⃣ Set parameters ---
k_base <- 0.002  # either 0.019 0.01 = 1 / 100; 1/500 = 0.002

Ea_nitrification <- 80000
Ed_nitrification <- 200000
t_opt_C <- 28
f_inorg_rapid <- 0.1136

# File paths
ald_file <- paste0("mean_ssp", ssp, "_fillnans.nc")
thawed_file <- paste0("total_thawed_extended/arctic_total_thawed_", ssp, "_no_lim_mean.nc")

# --- 2️⃣ Load data ---
if (!file.exists(ald_file)) stop(paste("ALD file not found:", ald_file))
ALD <- rast(ald_file)

common_extent <- ext(-179.95, 179.95, 45, 90)
ALD<-crop(ALD, common_extent)
#plot(ALD[[2]])
# surface Temperature
#temp_file <- paste0("mean_ST_", ssp, "_clean_final.nc")

# ALD Temperature
#temp_file <- paste0("ALD_temp/MRI_", ssp, "_clean.nc")
#if (!file.exists(temp_file)) stop(paste("Temperature file not found:", temp_file))

# deep soil temperature
temp_file <- paste0("deep_temp/mean_deep_temp_", ssp, ".nc")
if (!file.exists(temp_file)) stop(paste("Temperature file not found:", temp_file))
temp_data <- rast(temp_file)
temp_data <- crop(temp_data, common_extent)
crs(temp_data) <- crs(ALD)
temp_data <- temp_data[[1:250]] - 273.15  # Convert to °C
plot(temp_data[[2]])

# Fixation & deposition
if (ssp == "245") {
  n_fixation <- 0
  n_deposition <- 0
} else {
  n_fixation <- rast(paste0("fix_1850_2100_ssp", ssp, ".nc")) * 3600*24*365
  n_fixation <- crop(n_fixation, common_extent)
  n_deposition <- rast(paste0("dep_1850_2100_ssp", ssp, "_totalN.nc")) * 3600*24*365
  n_deposition <- crop(n_deposition, common_extent)
}
# Thawed N pool
if (!file.exists(thawed_file)) stop(paste("Thawed N file not found:", thawed_file))
thawed_total <- rast(thawed_file)

thawed_total <- crop(thawed_total, common_extent)
plot(thawed_total[[2]])
temp_data <- resample(temp_data, thawed_total, method = "bilinear")
plot(thawed_total[[250]])
# # soil moisture data
theta_file <- paste0("soil_moisture/sm_", ssp, "_present.nc")

theta <- rast(theta_file)
theta <- resample(theta, thawed_total, method = "bilinear")
plot(theta[[2]])

# --- 3️⃣ Peaked Arrhenius function ---
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

ext_60N <- common_extent

Ea<- Ea_nitrification
Ed<-Ed_nitrification
t_opt_C<-t_opt_C

years <- 1850:2015

ref_idx <- which(years >= 2000 & years <= 2015)
temp_present <- temp_data[[ref_idx]]
theta_present <- theta #already just present day
temp_present <- crop(temp_present, ext_60N)
theta_present <- crop(theta_present, ext_60N)
plot(theta_present[[2]])
plot(thawed_total[[2]])
k_T_present <- peaked_arrhenius(temp_present, Ea, Ed, t_opt_C)
#k_T_present <- mask(k_T_present, !is.na(thawed_total[[1]]))
#plot(k_T_present)
#k_T_ref_scalar <- mean(values(k_T_present), na.rm = TRUE)
k_T_ref_scalar <- mean(k_T_present,na.rm = TRUE) 
plot(k_T_ref_scalar)

# Moisture reference

f_theta_present <- clamp(1 - theta_present, 0, 1)
#f_theta_present <- mask(f_theta_present, !is.na(thawed_total[[1]]))
plot(f_theta_present[[2]])
f_theta_ref_scalar <- mean(values(f_theta_present), na.rm = TRUE)
#f_theta_ref_scalar <- mean(f_theta_present)


# # --- 4️⃣ Compute mineralisation flux ---
# compute_flux <- function(total_thawed,
#                          temp_data,
#                          n_deposition,
#                          n_fixation,
#                          k_base_mineralisation,
#                          Ea,
#                          Ed,
#                          theta,
#                          k_T_ref_scalar,
#                          f_theta_ref_scalar,
#                          t_opt_C = 28,
#                          f_inorg_rapid = 0.1136) {
# 
#   f_org <- 1 - f_inorg_rapid
#   land_mask <- !is.na(total_thawed)
#   # --- Temperature response ---
#   k_T_year <- peaked_arrhenius(temp_data, Ea, Ed, t_opt_C)
#   k_T_year <- terra::mask(k_T_year, land_mask)
#   k_factor <- k_T_year / k_T_ref_scalar
#   k_factor <- terra::ifel(land_mask,
#                           terra::ifel(is.na(k_factor), 1, k_factor),
#                           NA)
#   
#   #to switch temperature effect on / off:
#   #k_factor <- k_factor * 0 + 1
#   
#   # --- Moisture response ---
#   f_theta <- clamp(1 - theta, 0, 1)
#   k_moisture <- f_theta / f_theta_ref_scalar
#   k_moisture <- terra::ifel(land_mask,
#                             terra::ifel(is.na(k_moisture), 1, k_moisture),
#                             NA)
# 
#   # --- Pools ---
#   org_pool <- total_thawed * f_org
#   inorg_rapid_pool <- total_thawed * f_inorg_rapid
#   # --- Flux ---
#   k_t <- k_base_mineralisation * k_factor
#   
#   mineralised_N <- org_pool * k_t #* k_moisture
# 
#   bioavailable_N_total <- inorg_rapid_pool +
#     mineralised_N #+
#   #n_fixation +
#   #n_deposition
#   list(
#     mineralised_N = mineralised_N,
#     bioavailable_N_total = bioavailable_N_total,
#     k_t_values = k_factor,
#     k_moisture = k_moisture
#   )
# }

# version that calculates: organic_pool_year_t = organic_pool_previous_year - mineralised_N
compute_flux <- function(total_thawed,
                         temp_data,
                         n_deposition,
                         n_fixation,
                         k_base_mineralisation,
                         Ea,
                         Ed,
                         theta,
                         k_T_ref_scalar,
                         f_theta_ref_scalar,
                         t_opt_C = 28,
                         f_inorg_rapid = 0.1136) {
  
  f_org <- 1 - f_inorg_rapid
  n_years <- nlyr(total_thawed)
  
  land_mask <- !is.na(total_thawed[[n_years]])
  
  mineralised_list <- vector("list", n_years)
  bioavailable_list <- vector("list", n_years)
  k_t_list <- vector("list", n_years)
  k_moisture_list <- vector("list", n_years)
  organic_pool_list <- vector("list", n_years)
  
  # Initial organic pool in 1850
  organic_pool <- total_thawed[[1]] * f_org
  organic_pool <- terra::mask(organic_pool, land_mask)
  
  for (i in seq_len(n_years)) {
    
    cat("Processing year", i, "of", n_years, "\n")
    
    current_thawed <- terra::mask(total_thawed[[i]], land_mask)
    current_temp <- temp_data[[i]]
    
    # Temperature response
    k_T_year <- peaked_arrhenius(current_temp, Ea, Ed, t_opt_C)
    k_T_year <- terra::mask(k_T_year, land_mask)
    
    k_factor <- k_T_year / k_T_ref_scalar
    
    k_factor <- terra::ifel(
      land_mask,
      terra::ifel(is.na(k_factor), 1, k_factor),
      NA
    )
    
    # Moisture response
    f_theta <- clamp(1 - theta[[min(i, nlyr(theta))]], 0, 1)
    
    k_moisture <- f_theta / f_theta_ref_scalar
    
    k_moisture <- terra::ifel(
      land_mask,
      terra::ifel(is.na(k_moisture), 1, k_moisture),
      NA
    )
    
    # Mineralisation rate
    k_t <- k_base_mineralisation * k_factor * k_moisture
    
    # Mineralise from remaining organic pool
    mineralised_N <- organic_pool * k_t
    
    # Safety: cannot mineralise more than remaining pool
    mineralised_N <- terra::ifel(
      mineralised_N > organic_pool,
      organic_pool,
      mineralised_N
    )
    
    mineralised_N <- terra::ifel(
      mineralised_N < 0,
      0,
      mineralised_N
    )
    
    # Update organic pool after mineralisation loss
    organic_pool <- organic_pool - mineralised_N
    
    organic_pool <- terra::ifel(
      organic_pool < 0,
      0,
      organic_pool
    )
    
    organic_pool <- terra::mask(organic_pool, land_mask)
    
    # Rapid inorganic pool still from total thawed pool
    inorg_rapid_pool <- current_thawed * f_inorg_rapid
    
    bioavailable_N_total <- inorg_rapid_pool + mineralised_N
    
    mineralised_list[[i]] <- mineralised_N
    bioavailable_list[[i]] <- bioavailable_N_total
    k_t_list[[i]] <- k_factor
    k_moisture_list[[i]] <- k_moisture
    organic_pool_list[[i]] <- organic_pool
    
    rm(
      current_thawed,
      current_temp,
      k_T_year,
      k_factor,
      f_theta,
      k_moisture,
      k_t,
      mineralised_N,
      inorg_rapid_pool,
      bioavailable_N_total
    )
    gc()
  }
  
  list(
    mineralised_N = rast(mineralised_list),
    bioavailable_N_total = rast(bioavailable_list),
    k_t_values = rast(k_t_list),
    k_moisture = rast(k_moisture_list),
    organic_pool_remaining = rast(organic_pool_list)
  )
}



# # --- 5️⃣ Update total thawed pools ---
# update_total_thawed <- function(total_thawed, mineralised_N) {
#   terra::ifel(total_thawed - mineralised_N < 0, 0, total_thawed - mineralised_N)
# }

# --- 6️⃣ Run flux calculation ---
flux_result <- compute_flux(
  thawed_total,
  temp_data,
  n_deposition,
  n_fixation,
  k_base_mineralisation = k_base,
  Ea = Ea_nitrification,
  Ed = Ed_nitrification,
  theta = theta,
  k_T_ref_scalar = k_T_ref_scalar,
  f_theta_ref_scalar = f_theta_ref_scalar,
  t_opt_C = t_opt_C,
  f_inorg_rapid = f_inorg_rapid
)


# --- 1️⃣ Helper function ---
prepare_raster <- function(r, start_year = 1850) {
  n <- nlyr(r)
  years <- start_year:(start_year + n - 1)
  time_dates <- as.Date(paste0(years, "-07-01"))
  names(r) <- as.character(years)
  terra::time(r) <- time_dates
  return(r)
}

# --- 2️⃣ Prepare all rasters ---
bioavailable_N_total <- prepare_raster(flux_result$bioavailable_N_total)
mineralised_N        <- prepare_raster(flux_result$mineralised_N)
k_t_values           <- prepare_raster(flux_result$k_t_values)
k_moisture          <- prepare_raster(flux_result$k_moisture)
organic_pool_remaining <- prepare_raster(flux_result$organic_pool_remaining)

global(k_t_values, "mean", na.rm = TRUE)

# # Update pools first, THEN assign time
# total_thawed_after_mineralisation <- update_total_thawed(
#   thawed_total,
#   flux_result$mineralised_N
# )
# total_thawed_after_mineralisation <- prepare_raster(total_thawed_after_mineralisation)

#--- 3️⃣ Write NetCDF files ---
writeCDF(
  bioavailable_N_total,
  paste0("mineralised_results/no_ALD_limit/0.2perc/deep_temp/arctic_bioavailable_N_pool_", ssp, "_no_lim_no_sm_w_min_losses_mean.nc"),
  overwrite = TRUE,
  varname = "bio_available_N",
  longname = "Bioavailable soil N",
  unit = "kgN_m2"
)

writeCDF(
  mineralised_N,
  paste0("mineralised_results/no_ALD_limit/0.2perc/deep_temp/arctic_mineralised_N_", ssp, "_no_lim_no_sm_w_min_losses_mean.nc"),
  overwrite = TRUE,
  varname = "mineralised_N",
  longname = "Mineralised soil nitrogen",
  unit = "kgN_m2"
)

writeCDF(
  k_t_values,
  paste0("mineralised_results/no_ALD_limit/1perc/deep_temp/soil_moisure/arctic_k_t_", ssp, "_no_lim_w_min_losses_mean.nc"),
  overwrite = TRUE,
  varname = "k_t",
  longname = "Temperature response factor",
  unit = "-"
)

writeCDF(
  total_thawed_after_mineralisation,
  paste0("mineralised_results/no_ALD_limit/0.2perc/deep_temp/soil_moisure/arctic_total_thawed_after_mineralisation_", ssp, "_no_lim_w_min_losses_mean.nc"),
  overwrite = TRUE,
  varname = "total_thawed_N_after_mineralisation",
  longname = "Total thawed soil nitrogen after mineralisation",
  unit = "kgN_m2"
)

# --- 4️⃣ Done ---
#cat("Finished SSP", ssp, "\n")


# Area per grid cell
area_rast <- cellSize(organic_pool_remaining[[1]], unit = "m")

# Pg N per year
organic_N_pg <- global(
  organic_pool_remaining * area_rast,
  "sum",
  na.rm = TRUE
)[, 1] / 1e12

# Years
years <- as.integer(names(organic_N_pg))

organic_N_df <- data.frame(
  Year = years,
  organic_N_pg = organic_N_pg
)

print(head(bioavailable_df))
print(summary(bioavailable_df))

# Plot
p_bio <- ggplot(bioavailable_df, aes(x = Year, y = bioavailable_PgN)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  labs(
    title = "Bioavailable N through time",
    x = "Year",
    y = "Bioavailable N (Pg N)"
  )

print(p_bio)




# ============================================================

# Sum of total bioavailable N from 1900–2014

# ============================================================

# Select years

years <- as.integer(names(bioavailable_N_total))

idx <- which(years >= 1900 & years <= 2014)

# Subset raster

bioavailable_1900_2014 <- bioavailable_N_total[[idx]]

# Sum through time

bioavailable_sum <- app(
  
  bioavailable_1900_2014,
  
  sum,
  
  na.rm = TRUE
  
)

print(bioavailable_sum)

# ============================================================

# Pan-Arctic total in Pg N

# ============================================================

area_rast <- cellSize(
  
  bioavailable_sum,
  
  unit = "m"
  
)

total_pgN <- global(
  
  bioavailable_sum * area_rast,
  
  "sum",
  
  na.rm = TRUE
  
)[1,1] / 1e12

cat("Total cumulative bioavailable N 1900–2014 (Pg N):\n")

print(total_pgN)

years <- as.integer(names(bioavailable_N_total))

idx_2014 <- which(years == 2014)

bioavailable_2014 <- bioavailable_N_total[[idx_2014]]

print(bioavailable_2014)

area_rast <- cellSize(
  bioavailable_2014,
  unit = "m"
)

total_pgN <- global(
  bioavailable_2014 * area_rast,
  "sum",
  na.rm = TRUE
)[1,1] / 1e12

cat("Cumulative bioavailable N by 2014 (Pg N):\n")
print(total_pgN)
