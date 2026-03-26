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
mineralised_N<-rast("arctic_mineralised_N_126_3m_lim_std.tif")
total_thawed_after_mineralisation<-rast("arctic_total_thawed_after_mineralisation_126_3m_lim_std.tif")
temp_data <- rast("mean_ST_126_clean_final.nc")
# Set common extent
common_extent <- ext(-179.95, 179.95, 30, 90)
ext(temp_data) <- common_extent
crs(temp_data)<-crs(mineralised_N)
# Test with subset of data

test_temp <- temp_data[[1:250]]  
test_temp<- test_temp - 273.15

k_base_denitrification<-0.091 
## k_base_denitrification:baseline denitrification rate constant under reference conditions.
# ~9.1% of the available nitrate pool can be denitrified per year before temperature scaling
######

# mask only land pixels
land_mask <- !is.na(mineralised_N[[1]])
#plot(land_mask)
test_temp<-test_temp[[1:250]]
mineralised_N<-mineralised_N[[1:250]]
total_thawed_after_mineralisation<-total_thawed_after_mineralisation[[1:250]]
rm()
years <- 1850:2099


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

# Initialize lists for individual loss components
cumulative_denit_list <- vector("list", n_years)
cumulative_runoff_list <- vector("list", n_years)

# Initialize cumulative totals
cumulative_denit <- mineralised_N[[1]] * 0
cumulative_runoff <- mineralised_N[[1]] * 0
#years<-1949:2099

# Modified compute_losses function to return individual loss components
compute_losses <- function(
    mineralised_N,              # fast pool (SpatRaster)
    total_thawed,               # slow pool (SpatRaster)
    test_temp,                  # temperature (SpatRaster, same layers as pools)
    years = 1850:2099,
    ref_start = 2000,
    ref_end   = 2020,
    f_denitrification_nitrate,
    f_runoff,
    Ea_denitrification,
    Ed_denitrification,
    k_base_denitrification,
    t_opt_C = 28
) {
  
  land_mask <- !is.na(total_thawed)
  
  # --- compute temperature response for every year ---
  k_T_year <- peaked_arrhenius(
    test_temp,
    Ea_denitrification,
    Ed_denitrification,
    t_opt_C
  )
  
  # enforce fixed land footprint across all years
  k_T_year <- terra::mask(k_T_year, land_mask, maskvalues = FALSE)
  
  # baseline indices (2000–2020)
  ref_idx <- which(years >= ref_start & years <= ref_end)
  if (length(ref_idx) == 0) stop("Reference years not found in 'years' vector.")
  
  # per-cell baseline
  k_T_ref_raster <- terra::mean(k_T_year[[ref_idx]], na.rm = TRUE)
  
  # normalized multiplier (same method as compute_flux)
  k_factor <- k_T_year / k_T_ref_raster
  
  # fill NA temps on land with 1, keep ocean NA
  k_factor <- terra::ifel(
    land_mask,
    terra::ifel(is.na(k_factor), 1, k_factor),
    NA
  )
  
  # effective denitrification rate constant
  k_t <- k_base_denitrification * k_factor
  
  # --- compute losses ---
  fast_loss_denit  <- (mineralised_N * 0.10) * f_denitrification_nitrate * k_t
  fast_loss_runoff <- mineralised_N * f_runoff
  slow_loss_runoff <- total_thawed  * f_runoff
  
  total_denitrification_loss <- fast_loss_denit
  total_runoff_loss <- fast_loss_runoff + slow_loss_runoff
  
  fast_after_loss <- mineralised_N - fast_loss_denit - fast_loss_runoff
  slow_after_loss <- total_thawed  - slow_loss_runoff
  
  # clamp negatives
  fast_after_loss <- terra::ifel(fast_after_loss < 0, 0, fast_after_loss)
  slow_after_loss <- terra::ifel(slow_after_loss < 0, 0, slow_after_loss)
  
  list(
    fast_after_loss = fast_after_loss,
    slow_after_loss = slow_after_loss,
    fast_loss_denit = fast_loss_denit,
    fast_loss_runoff = fast_loss_runoff,
    slow_loss_runoff = slow_loss_runoff,
    total_denitrification_loss = total_denitrification_loss,
    total_runoff_loss = total_runoff_loss,
    k_t = k_t,
    k_factor = k_factor,
    k_T_ref_raster = k_T_ref_raster
  )
}

#compute_losses <- function(
 #   mineralised_N,              # fast pool
  #  total_thawed,               # slow pool
   # test_temp,
    #f_denitrification_nitrate,
    #f_runoff,
    #Ea_denitrification,
    #Ed_denitrification,
    #k_base_denitrification,
    #t_opt_C
#) {
  # temperature multiplier (peaked Arrhenius)
 # k_factor <- peaked_arrhenius(
  #  test_temp,
   # Ea_denitrification,
    #Ed_denitrification,
    #t_opt_C
#  )
  
  # initialise k_t with NA
 # k_t <- mineralised_N * NA
  # k_t is the effective denitrification rate constant at the current soil temperature.
  # pixels where test_temp is NA → use base rate
  #k_t <- terra::ifel(land_mask & is.na(k_factor), 
   #                  k_base_denitrification, 
    #                 k_t)
  
  # pixels with valid temperature → k_factor * k_base
#  k_t <- terra::ifel(land_mask & !is.na(k_factor),
 #                    k_factor * k_base_denitrification,
  #                   k_t)
  
  # --- compute losses ---
  
  # denitrification only acts on nitrate = 10% of mineralised N
  #fast_loss_denit <- (mineralised_N * 0.10) * f_denitrification_nitrate * k_t
  
  #fast_loss_runoff <- mineralised_N * f_runoff
  #slow_loss_runoff <- total_thawed * f_runoff
  
  # total losses
  #total_denitrification_loss <- fast_loss_denit
  #total_runoff_loss <- fast_loss_runoff + slow_loss_runoff
  
  # updated pools
  #fast_after_loss <- mineralised_N - fast_loss_denit - fast_loss_runoff
  #slow_after_loss <- total_thawed - slow_loss_runoff
  
  # clamp negatives
  #fast_after_loss <- terra::ifel(fast_after_loss < 0, 0, fast_after_loss)
  #slow_after_loss <- terra::ifel(slow_after_loss < 0, 0, slow_after_loss)
  
  #list(
   # fast_after_loss = fast_after_loss,
    #slow_after_loss = slow_after_loss,
    #fast_loss_denit = fast_loss_denit,          # denitrification from fast pool
    #fast_loss_runoff = fast_loss_runoff,        # runoff from fast pool
    #slow_loss_runoff = slow_loss_runoff,        # runoff from slow pool
    #total_denitrification_loss = total_denitrification_loss,  # all denitrification
    #total_runoff_loss = total_runoff_loss,      # all runoff
    #k_t = k_t
#  )
#}

# Reset progress bar
pb <- txtProgressBar(min = 0, max = n_years, style = 3)

for (i in seq_along(years)) {
  # Just use i as the index (assuming same order for all rasters)
  mineralised  <- mineralised_N[[i]]               # 1-layer
  total_thawed <- total_thawed_after_mineralisation[[i]]  # 1-layer
  temp_layer   <- test_temp[[i]]                          # 1-layer
  
  # --- compute losses ---
  res <- compute_losses(
    mineralised_N = mineralised,
    total_thawed  = total_thawed,
    test_temp     = temp_layer,
    years = 1850:2099,
    ref_start = 2000,
    ref_end   = 2020,
    f_denitrification_nitrate = f_denitrification_nitrate,
    f_runoff      = f_runoff,
    Ea_denitrification = Ea_denitrification,
    Ed_denitrification = Ed_denitrification,
    k_base_denitrification = k_base_denitrification,
    t_opt_C = t_opt_C
  )
  
  # --- update cumulative totals ---
  cumulative_denit <- cumulative_denit + res$total_denitrification_loss[[1]]
  cumulative_runoff <- cumulative_runoff + res$total_runoff_loss[[1]]
  
  # --- save outputs ---
  # Save cumulative losses
  cumulative_denit_list[[i]] <- cumulative_denit
  cumulative_runoff_list[[i]] <- cumulative_runoff
  
  # Save pool states (keep your existing code)
  fast_pool_list[[i]] <- res$fast_after_loss[[1]]
  slow_pool_list[[i]] <- res$slow_after_loss[[1]]
  
  # --- update state for next iteration ---
  fast_pool <- res$fast_after_loss[[1]]
  slow_pool <- res$slow_after_loss[[1]]
  
  # --- update progress bar ---
  setTxtProgressBar(pb, i)
}

close(pb)

# Save individual loss components as separate raster stacks
denitrification_cumulative_r <- rast(cumulative_denit_list)
writeRaster(denitrification_cumulative_r, "denitrification_loss_cumulative_126_std.tif")
rm(denitrification_cumulative_r, cumulative_denit_list)

runoff_cumulative_r <- rast(cumulative_runoff_list)
writeRaster(runoff_cumulative_r, "runoff_loss_cumulative_126_std.tif")
rm(runoff_cumulative_r, cumulative_runoff_list)

fast_pool_r       <- rast(fast_pool_list)
writeRaster(fast_pool_r, "mineralised_N_pool_126_std.tif")
rm(fast_pool_r, fast_pool_list)

slow_pool_r       <- rast(slow_pool_list)
writeRaster(slow_pool_r, "thawed_organic_N_pool_126_std.tif")
rm(slow_pool_r, slow_pool_list)


library(terra)
library(dplyr)

library(tidyr)
library(ggplot2)

# ---- load stacks ----
denit  <- rast("denitrification_loss_cumulative_126_no_lim_mean.tif")
runoff <- rast("runoff_loss_cumulative_126_no_lim_mean.tif")
fast   <- rast("mineralised_N_pool_126_no_lim_mean.tif")
slow   <- rast("thawed_organic_N_pool_126_no_lim_mean.tif")

# years (adjust if your stacks don't start at 1850)
years <- 1850:(1850 + nlyr(denit) - 1)

# ---- land/ocean mask from any stack (layer 1 is fine) ----
land_mask <- !is.na(denit[[1]])

# ---- cell area weights, masked to land only ----
area_m2 <- cellSize(denit[[1]], unit = "m")
area_m2 <- mask(area_m2, land_mask, maskvalues = FALSE)

# helper: weighted mean time series for a SpatRaster stack
weighted_arctic_mean_ts <- function(x, area_m2) {
  # ensure same footprint
  x <- mask(x, area_m2, maskvalues = NA)
  terra::global(x, fun = "mean", weights = area_m2, na.rm = TRUE)[, 1]
}

ts_denit  <- weighted_arctic_mean_ts(denit,  area_m2)
ts_runoff <- weighted_arctic_mean_ts(runoff, area_m2)
ts_fast   <- weighted_arctic_mean_ts(fast,   area_m2)
ts_slow   <- weighted_arctic_mean_ts(slow,   area_m2)

# ---- build long df for plotting ----
df <- tibble(
  Year = years,
  Denitrification_cum = ts_denit,
  Runoff_cum          = ts_runoff,
  Fast_pool           = ts_fast,
  Slow_pool           = ts_slow
) %>%
  pivot_longer(-Year, names_to = "Variable", values_to = "Mean")

# ---- plot ----
ggplot(df, aes(x = Year, y = Mean, color = Variable)) +
  geom_line(linewidth = 0.8) +
  theme_minimal() +
  labs(
    x = "Year",
    y = "Area-weighted Arctic mean (masked to land)",
    color = NULL
  ) +
  theme(legend.position = "bottom")
