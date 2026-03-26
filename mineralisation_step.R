# code to compute bioavailable n pool with temperature - dependent mineralisation flux + n deposition + BNF
library(terra)
#terra::gdal()
library(ggplot2)

# Load ALD files
ALD <- rast("mean_ssp126_corr_new.nc") # insert ssp 126, 126, 126, 126
#ALD <- clamp(ALD, upper = 3, values = TRUE) # for bioavailable N 
# Set common extent
common_extent <- ext(-179.95, 179.95, 30, 90)
# Load and prepare soil surface temperature data
temp_data <- rast("ALD_temp/MRI_126_clean.nc")
#temp_data <- rast("mean_ST_126_clean_final.nc")
ext(temp_data) <- common_extent
crs(temp_data)<-crs(ALD)
temp_data <- crop(temp_data, common_extent)
temp_data <- temp_data[[1:250]]
temp_data <- temp_data - 273.15

plot(temp_data[[2]])

# Define Nitrogen loss parameters
f_n_loss <- 0.4927715
f_denitrification_nitrate <- 0.4927401
f_runoff <- 3.137 * 10^(-05)

# Parameters
Ea_nitrification <- 80000   # J/mol
Ed_nitrification <- 200000  # J/mol
Ea_denitrification <- 47000   # J/mol
Ed_denitrification <- 200000  # J/mol
t_opt_C <- 28



# Define peaked_arrhenius function with your specific formulation
peaked_arrhenius <- function(temp_C, Ea, Ed, t_opt_C = 28) {
  # Gas constant
  R_gas <- 8.314462618
  
  # Convert to Kelvin
  temp_K <- temp_C + 273.15
  t_opt_K <- t_opt_C + 273.15
  
  # Helper terms
  hlp1 <- temp_K - t_opt_K
  hlp2 <- temp_K * t_opt_K * R_gas
  
  # Core Arrhenius terms
  numerator <- Ed * exp(Ea * hlp1 / hlp2)
  denominator <- Ed - Ea * (1 - exp(Ed * hlp1 / hlp2))
  
  # Calculate result
  result <- numerator / denominator
  return(result)
}

###
#load fixation, deposition, (but not for ssp 126, not available; comment lines)
#then convert from kg N / m2 * s to kg N / m2 * yr
# 
n_fixation<-rast("fix_1850_2100_ssp126.nc")
n_deposition<-rast("dep_1850_2100_ssp126_totalN.nc")

n_fixation<- n_fixation* 3600 * 24 * 365
plot(n_fixation[[2]])
ext(n_fixation) <- common_extent
n_deposition<- n_deposition* 3600 * 24 * 365
ext(n_deposition)<-common_extent
# ######


# --- 2 Compute mineralisation flux ---
#######

years<- c(1850:2099)
# 

compute_flux <- function(total_thawed,
                         temp_data,
                         n_deposition,
                         n_fixation,
                         ref_start = 2000,
                         ref_end   = 2020,
                         k_base_mineralisation,
                         Ea,
                         Ed,
                         t_opt_C = 28,
                         f_inorg_rapid = 0.1235) {
  
  f_org <- 1 - f_inorg_rapid
  land_mask <- !is.na(total_thawed)
  
  # temperature response
  k_T_year <- peaked_arrhenius(temp_data, Ea, Ed, t_opt_C)
  k_T_year <- terra::mask(k_T_year, land_mask, maskvalues = FALSE)
  
  # baseline
  ref_idx <- which(years >= ref_start & years <= ref_end)
  k_T_ref_raster <- terra::mean(k_T_year[[ref_idx]], na.rm = TRUE)
  
  # scaling factor
  k_factor <- k_T_year / k_T_ref_raster
  
  k_factor <- terra::ifel(
    land_mask,
    terra::ifel(is.na(k_factor), 1, k_factor),
    NA
  )
  
  k_t <- k_base_mineralisation * k_factor
  
  # pools
  org_pool         <- total_thawed * f_org
  inorg_rapid_pool <- total_thawed * f_inorg_rapid
  
  # mineralisation flux
  mineralised_N <- org_pool * k_t
  
  # --- ADDITION HERE ---
  bioavailable_N_total <- inorg_rapid_pool +
    mineralised_N +
   n_fixation +
   n_deposition
  
  list(
    mineralised_N     = mineralised_N,
    bioavailable_N_total  = bioavailable_N_total,
    inorg_rapid_pool  = inorg_rapid_pool,
    k_t_values        = k_factor,
    land_mask         = land_mask,
    k_T_ref_raster    = k_T_ref_raster
  )
}
update_total_thawed <- function(total_thawed, mineralised_N) {
  # Direct arithmetic with protection against negative values
  total_thawed_after_mineralisation <- terra::ifel(total_thawed - mineralised_N < 0, 0, total_thawed - mineralised_N)
  
  list(
    total_thawed_after_mineralisation = total_thawed_after_mineralisation
  )
}


# Step 2: Mineralisation flux
thawed_total<-rast("1.9%_min/arctic_total_thawed_126_no_lim.tif") # add mean 
plot(thawed_total[[2]])
#flux_result <- compute_flux(thawed_total, temp_data, ref_start=2000, ref_end=2020, k_base_mineralisation=0.019,Ea_nitrification, Ed_nitrification, t_opt_C= 28, f_inorg_rapid=0.1235)
flux_result <- compute_flux(
  thawed_total,
  temp_data,
  n_fixation,
  n_deposition,
  ref_start = 2000,
  ref_end   = 2020,
  k_base_mineralisation = 0.019,
  Ea_nitrification,
  Ed_nitrification,
  t_opt_C = 28,
  f_inorg_rapid = 0.1235
)

print(flux_result)
writeRaster(flux_result$bioavailable_N_total, "1.9%_min/ALD_temp/arctic_bioavailable_N_pool_126_no_lim_mean.tif")
writeRaster(flux_result$mineralised_N, "1.9%_min/ALD_temp/arctic_mineralised_N_126_no_lim_mean.tif")
#writeRaster(flux_result$inorg_rapid_pool,"1.9%_min/arctic_baseline_mineralised_N_126_no_lim_mean.tif")
writeRaster(flux_result$k_t_values, "1.9%_min/ALD_temp/arctic_k_t_126_no_lim_mean.tif")

rm(flux_result)

mineralised_N<-rast("1.9%_min/ALD_temp/arctic_mineralised_N_126_no_lim_mean.tif")

# Then update the pools using the computed flux
pools_result <- update_total_thawed(thawed_total, mineralised_N)

writeRaster(pools_result$total_thawed_after_mineralisation, "1.9%_min/ALD_temp/arctic_total_thawed_after_mineralisation_126_no_lim_mean.tif")
