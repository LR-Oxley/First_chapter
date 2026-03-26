# clean version of computing thawed N, the N losses through denitrification and runoff, 
# and the increased mineralisation due to higher soil temp;
install.packages("here", repos="https://cloud.r-project.org")
f_n_loss <- 0.4927715
f_denitrification_nitrate <- 0.4927401
f_runoff <- 3.137 * 10^(-05)

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
mineralised_N<-rast("arctic_mineralised_N_585_3m_lim.tif")
total_thawed_after_mineralisation<-rast("arctic_total_thawed_after_mineralisation_585_3m_lim.tif")
temp_data <- rast("mean_ST_585_clean_final.nc")
# Set common extent
common_extent <- ext(-179.95, 179.95, 30, 90)
ext(temp_data) <- common_extent
crs(temp_data)<-crs(mineralised_N)
# Test with subset of data

test_temp <- temp_data[[1:250]]  
test_temp<- test_temp - 273.15

k_base_denitrification<-0.091
######

# mask only land pixels
land_mask <- !is.na(mineralised_N[[1]])
#plot(land_mask)
test_temp<-test_temp[[1:250]]
mineralised_N<-mineralised_N[[1:250]]
total_thawed_after_mineralisation<-total_thawed_after_mineralisation[[1:250]]
rm()
years <- 1850:2099
#years <- as.character(years)         # if factor
#years_date <- as.Date(paste0(years, "-07-01"))
time(test_temp) <- years
time(mineralised_N) <- years
time(total_thawed_after_mineralisation) <- years


######

# --- INITIAL STATE (year 1850) ---
fast_pool <- mineralised_N[[1]]                       # state
slow_pool <- total_thawed_after_mineralisation[[1]]   # state

#years <- 2095:2099    # or whatever range you want

# preallocate lists
n_years <- length(years)
cumulative_total_loss_list      <- vector("list", n_years)
fast_pool_list       <- vector("list", n_years)
slow_pool_list       <- vector("list", n_years)
cumulative_total_loss <- mineralised_N[[1]] * 0
#years<-1949:2099

compute_losses <- function(
    mineralised_N,              # fast pool
    total_thawed,               # slow pool
    test_temp,
    f_denitrification_nitrate,
    f_runoff,
    Ea_denitrification,
    Ed_denitrification,
    k_base_denitrification,
    t_opt_C
) {
  # temperature multiplier (peaked Arrhenius)
  k_factor <- peaked_arrhenius(
    test_temp,
    Ea_denitrification,
    Ed_denitrification,
    t_opt_C
  )
  
  # initialise k_t with NA
  k_t <- mineralised_N * NA
  
  # pixels where test_temp is NA → use base rate
  k_t <- terra::ifel(land_mask & is.na(k_factor), 
                     k_base_denitrification, 
                     k_t)
  
  # pixels with valid temperature → k_factor * k_base
  k_t <- terra::ifel(land_mask & !is.na(k_factor),
                     k_factor * k_base_denitrification,
                     k_t)
  
  # --- compute losses ---
  
  # denitrification only acts on nitrate = 10% of mineralised N
  fast_loss_denit <- (mineralised_N * 0.10) *f_denitrification_nitrate * k_t
  
  fast_loss_runoff <- mineralised_N * f_runoff
  slow_loss_runoff <- total_thawed * f_runoff
  
  # updated pools
  fast_after_loss <- mineralised_N - fast_loss_denit - fast_loss_runoff
  slow_after_loss <- total_thawed - slow_loss_runoff
  
  # clamp negatives
  fast_after_loss <- terra::ifel(fast_after_loss < 0, 0, fast_after_loss)
  slow_after_loss <- terra::ifel(slow_after_loss < 0, 0, slow_after_loss)
  
  total_loss <- fast_loss_denit + fast_loss_runoff + slow_loss_runoff
  
  list(
    fast_after_loss = fast_after_loss,
    slow_after_loss = slow_after_loss,
    total_loss = total_loss,
    k_t = k_t
  )
}

# progress bar
pb <- txtProgressBar(min = 0, max = n_years, style = 3)
for (i in seq_along(years)) {
  yr <- years[i]
  idx <- which(time(test_temp) == yr)
  
  
  # --- background data for the current year (already increasing) ---
  mineralised  <- mineralised_N[[idx]]               # 1-layer
  total_thawed <- total_thawed_after_mineralisation[[idx]]  # 1-layer
  temp_layer   <- test_temp[[idx]]                         # 1-layer
  
  # --- compute losses ---
  res <- compute_losses(
    mineralised_N = mineralised,
    total_thawed  = total_thawed,
    test_temp     = temp_layer,
    f_denitrification_nitrate = f_denitrification_nitrate,
    f_runoff      = f_runoff,
    Ea_denitrification = Ea_denitrification,
    Ed_denitrification = Ed_denitrification,
    k_base_denitrification = k_base_denitrification,
    t_opt_C = t_opt_C
  )
  
  # --- save outputs (1-layer rasters per year) ---
  layer_loss <- res$total_loss[[1]]
  cumulative_total_loss <- cumulative_total_loss + layer_loss  # update
  
  # --- save outputs (1-layer rasters per year) ---
  cumulative_total_loss_list[[i]] <- cumulative_total_loss  # cumulative
  fast_pool_list[[i]]       <- res$fast_after_loss[[1]]
  slow_pool_list[[i]]       <- res$slow_after_loss[[1]]
  
  # --- update state for next iteration ---
  fast_pool <- res$fast_after_loss[[1]]
  slow_pool <- res$slow_after_loss[[1]]
  
  ## --- update progress bar ---
  setTxtProgressBar(pb, i)
}



total_loss_r      <- rast(cumulative_total_loss_list)
writeRaster(total_loss_r, "total_loss_r_585_mean.tif")
rm(total_loss_r, cumulative_total_loss_list)

fast_pool_r       <- rast(fast_pool_list)
writeRaster(fast_pool_r, "mineralised_N_pool_585.tif")
rm(fast_pool_r, fast_pool_list)

slow_pool_r       <- rast(slow_pool_list)
writeRaster(slow_pool_r, "thawed_organic_N_pool_585.tif")
rm(slow_pool_r, slow_pool_list)

