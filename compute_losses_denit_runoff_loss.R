f_runoff <- 0.000051

# Parameters
Ea_nitrification <- 80000   # J/mol
Ed_nitrification <- 200000  # J/mol
Ea_denitrification <- 47000   # J/mol
Ed_denitrification <- 200000  # J/mol
k_base_nitrification <- 0.1235 # from QUINCY output
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

# Step 3: losses
library(here)
library(terra)
mineralised_N<-rast("1.9%_min/3m_lim/arctic_min_N_pool_245_3m_lim_mean.tif")
total_thawed_after_mineralisation<-rast("1.9%_min/3m_lim/arctic_total_thawed_after_mineralisation_245_3m_lim_mean.tif")
temp_data <- rast("mean_ST_245_clean_final.nc")

# Set common extent
common_extent <- ext(-179.95, 179.95, 30, 90)
ext(temp_data) <- common_extent
crs(temp_data)<-crs(mineralised_N)
# Test with subset of data

test_temp <- temp_data[[1:250]]
test_temp<- test_temp - 273.15

k_base_denitrification<-0.0308 # 

###
# #load fixation, deposition, then convert from kg N / m2 * s to kg N / m2 * yr
# n_fixation<-rast("fix_1850_2100_ssp245.nc")
# n_deposition<-rast("dep_1850_2100_ssp245_totalN.nc")
# 
# n_fixation<- n_fixation* 3600 * 24 * 365
# ext(n_fixation) <- common_extent
# n_deposition<- n_deposition* 3600 * 24 * 365
# ext(n_deposition)<-common_extent
######

# mask only land pixels
land_mask <- !is.na(mineralised_N[[1]])
#plot(land_mask)
test_temp<-test_temp[[1:250]]
mineralised_N<-mineralised_N[[1:250]]
total_thawed_after_mineralisation<-total_thawed_after_mineralisation[[1:250]]

years <- 1850:2099

ref_start <- 2000
ref_end   <- 2020
ref_idx <- which(years >= ref_start & years <= ref_end)

land_mask <- !is.na(mineralised_N[[1]])  # fixed footprint

# temperature response for every year (250 layers)
k_T_year <- peaked_arrhenius(test_temp, Ea_denitrification, Ed_denitrification, t_opt_C)

# enforce fixed land footprint across all years
k_T_year <- terra::mask(k_T_year, land_mask, maskvalues = FALSE)

# per-cell baseline (mean across 2000–2020)
k_T_ref_raster <- terra::mean(k_T_year[[ref_idx]], na.rm = TRUE)

# normalized multiplier (baseline ~ 1)
k_factor_stack <- k_T_year / k_T_ref_raster

# fill NA temps on land with 1, keep ocean NA
k_factor_stack <- terra::ifel(
  land_mask,
  terra::ifel(is.na(k_factor_stack), 1, k_factor_stack),
  NA
)

# optional: free big intermediates
rm(k_T_year, k_t_ref_raster); gc()

######

# --- INITIAL STATE (year 1850) ---
fast_pool <- mineralised_N[[1]]                       # state
slow_pool <- total_thawed_after_mineralisation[[1]]   # state


# preallocate lists
n_years <- length(years)
cumulative_total_loss_list      <- vector("list", n_years)
fast_pool_list       <- vector("list", n_years)
slow_pool_list       <- vector("list", n_years)

# Initialize lists for individual loss components
cumulative_denit_list <- vector("list", n_years)
cumulative_runoff_list <- vector("list", n_years)

# Initialize cumulative totals
cumulative_denit <- mineralised_N[[1]] * 0
cumulative_runoff <- mineralised_N[[1]] * 0


compute_losses <- function(
    mineralised_N,
    total_thawed,
    k_factor,
    f_runoff,
    k_base_denitrification
    #n_fixation,          # raster or scalar
    #n_deposition       # raster or scalar
) {
  
  k_t <- k_base_denitrification * k_factor
  
  # add external inputs to fast pool
  #non_permafrost <- n_fixation + n_deposition
  inorg_pool <- mineralised_N #+ non_permafrost
  org_pool<-total_thawed #+ non_permafrost
  # losses (apply to updated fast pool)
  fast_loss_denit  <- (inorg_pool*0.1)  * k_t
  fast_loss_runoff <- inorg_pool * f_runoff
  slow_loss_runoff <- org_pool * f_runoff
  
  fast_after_loss <- inorg_pool - fast_loss_denit - fast_loss_runoff
  slow_after_loss <- org_pool - slow_loss_runoff
  
  fast_after_loss <- terra::ifel(fast_after_loss < 0, 0, fast_after_loss)
  slow_after_loss <- terra::ifel(slow_after_loss < 0, 0, slow_after_loss)
  
  list(
    fast_after_loss = fast_after_loss,
    slow_after_loss = slow_after_loss,
    total_denitrification_loss = fast_loss_denit,
    total_runoff_loss = fast_loss_runoff + slow_loss_runoff,
    k_t = k_t
    #fast_inputs = non_permafrost
  )
}


library(terra)

n_years <- length(years)
# Reset progress bar
pb <- txtProgressBar(min = 0, max = n_years, style = 3)

for (i in seq_along(years)) {
  
  mineralised  <- mineralised_N[[i]]
  total_thawed <- total_thawed_after_mineralisation[[i]]
  kfac_layer   <- k_factor_stack[[i]]
  
  #nfix_i <- n_fixation[[i]]      # or just n_fixation if scalar
  #ndep_i <- n_deposition[[i]]    # or just n_deposition if scalar
  
  res <- compute_losses(
    mineralised_N = mineralised,
    total_thawed  = total_thawed,
    k_factor      = kfac_layer,
    f_runoff      = f_runoff,
    k_base_denitrification = k_base_denitrification
    #n_fixation    = nfix_i,
    #n_deposition  = ndep_i
  )
  
  cumulative_denit  <- cumulative_denit  + res$total_denitrification_loss[[1]]
  cumulative_runoff <- cumulative_runoff + res$total_runoff_loss[[1]]
  
  cumulative_denit_list[[i]]  <- cumulative_denit
  cumulative_runoff_list[[i]] <- cumulative_runoff
  fast_pool_list[[i]] <- res$fast_after_loss[[1]]
  slow_pool_list[[i]] <- res$slow_after_loss[[1]]
  
  setTxtProgressBar(pb, i)
}


close(pb)
# Save individual loss components as separate raster stacks
denitrification_cumulative_r <- rast(cumulative_denit_list)
writeRaster(denitrification_cumulative_r, "1.9%_min/3m_lim/denitrification_loss_cumulative_245_3m_lim_mean.tif")
rm(denitrification_cumulative_r, cumulative_denit_list)

runoff_cumulative_r <- rast(cumulative_runoff_list)
writeRaster(runoff_cumulative_r, "1.9%_min/3m_lim/runoff_loss_cumulative_245_3m_lim_mean.tif")
rm(runoff_cumulative_r, cumulative_runoff_list)

fast_pool_r       <- rast(fast_pool_list)
writeRaster(fast_pool_r, "1.9%_min/3m_lim/mineralised_N_pool_245_3m_lim_mean.tif")
rm(fast_pool_r, fast_pool_list)

slow_pool_r       <- rast(slow_pool_list)
writeRaster(slow_pool_r, "1.9%_min/3m_lim/thawed_organic_N_pool_245_3m_lim_mean.tif")
rm(slow_pool_r, slow_pool_list)
