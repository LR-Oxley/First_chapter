####################################################################################################
# Compute cumulative N losses (denitrification + runoff)
# SSP-aware version (for Slurm job arrays)
# Cumulative mineralised N since 1850
####################################################################################################

library(terra)
library(here)

# -----------------------------
# 0) READ SSP ARGUMENT
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("No SSP provided")


ssp <- args[1]
message("Running losses for SSP", ssp)

# -----------------------------
# 1) CONSTANTS & PARAMETERS
# -----------------------------
f_runoff <- 0.000051

Ea_denitrification <- 47000
Ed_denitrification <- 200000
k_base_denitrification <- 0.0308
t_opt_C <- 28

years <- 1850:2099
n_years <- length(years)

common_extent <- ext(-179.95, 179.95, 30, 90)

# -----------------------------
# 2) FUNCTIONS
# -----------------------------
peaked_arrhenius <- function(temp_C, Ea, Ed, t_opt_C = 28) {
  R <- 8.314462618
  T  <- temp_C + 273.15
  To <- t_opt_C + 273.15
  
  hlp1 <- T - To
  hlp2 <- T * To * R
  
  num <- Ed * exp(Ea * hlp1 / hlp2)
  den <- Ed - Ea * (1 - exp(Ed * hlp1 / hlp2))
  
  num / den
}

compute_losses <- function(
    mineralised_N,
    total_thawed,
    k_factor,
    f_runoff,
    k_base_denitrification,
    n_fixation,          # raster or scalar
    n_deposition       # raster or scalar
) {
  
  k_t <- k_base_denitrification * k_factor
  
  # add external inputs to fast pool
  non_permafrost <- n_fixation + n_deposition
  inorg_pool <- mineralised_N + non_permafrost
  org_pool<-total_thawed + non_permafrost
  # losses (apply to updated fast pool)
  fast_loss_denit  <- (inorg_pool*0.1) * k_t
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
    k_t = k_t, 
    fast_inputs = non_permafrost
  )
}

# -----------------------------
# 3) INPUT DATA
# -----------------------------
mineralised_N <- rast(
  paste0("1.9%_min/3m_lim/arctic_min_N_pool_", ssp, "_3m_lim_std.tif")
)
total_thawed <- rast(
  paste0("1.9%_min/3m_lim/arctic_total_thawed_after_mineralisation_", ssp, "_3m_lim_std.tif")
)

temp_K <- rast(paste0("mean_ST_", ssp, "_clean_final.nc"))
ext(temp_K) <- common_extent
crs(temp_K) <- crs(mineralised_N)
temp_C <- temp_K - 273.15

n_fixation <- rast(paste0("fix_1850_2100_ssp", ssp, ".nc")) * 3600 * 24 * 365
n_deposition <- rast(paste0("dep_1850_2100_ssp", ssp, "_totalN.nc")) * 3600 * 24 * 365
ext(n_fixation) <- ext(n_deposition) <- common_extent

# -----------------------------
# 4) TEMPERATURE MULTIPLIER
# -----------------------------
land_mask <- !is.na(mineralised_N[[1]])

k_T <- peaked_arrhenius(temp_C, Ea_denitrification, Ed_denitrification, t_opt_C)
k_T <- mask(k_T, land_mask, maskvalues = FALSE)

ref_idx <- which(years >= 2000 & years <= 2020)
k_T_ref <- mean(k_T[[ref_idx]], na.rm = TRUE)

k_factor_stack <- k_T / k_T_ref
k_factor_stack <- ifel(land_mask, ifel(is.na(k_factor_stack), 1, k_factor_stack), NA)

rm(k_T); gc()

# -----------------------------
# 5) INITIAL STATE
# -----------------------------
cumulative_denit  <- mineralised_N[[1]] * 0
cumulative_runoff <- mineralised_N[[1]] * 0

cumulative_denit_list  <- vector("list", n_years)
cumulative_runoff_list <- vector("list", n_years)
fast_pool_list <- vector("list", n_years)
slow_pool_list <- vector("list", n_years)

pb <- txtProgressBar(min = 0, max = n_years, style = 3)

# -----------------------------
# 6) TIME LOOP
# -----------------------------
for (i in seq_along(years)) {
  
  mineralised_i  <- mineralised_N[[i]]
  total_thawed_i <- total_thawed[[i]]
  kfac_layer     <- k_factor_stack[[i]]
  
  nfix_i <- n_fixation[[i]]
  ndep_i <- n_deposition[[i]]
  
  res <- compute_losses(
    mineralised_N = mineralised_i,
    total_thawed  = total_thawed_i,
    k_factor      = kfac_layer,
    f_runoff      = f_runoff,
    k_base_denitrification = k_base_denitrification,
    n_fixation    = nfix_i,
    n_deposition  = ndep_i
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

# -----------------------------
# 7) OUTPUT
# -----------------------------
writeRaster(
  rast(cumulative_denit_list),
  paste0("1.9%_min/3m_lim/denitrification_loss_cumulative_", ssp, "_3m_lim_std.tif"),
  overwrite = TRUE
)

writeRaster(
  rast(cumulative_runoff_list),
  paste0("1.9%_min/3m_lim/runoff_loss_cumulative_", ssp, "_3m_lim_std.tif"),
  overwrite = TRUE
)

writeRaster(
  rast(fast_pool_list),
  paste0("1.9%_min/3m_lim/mineralised_N_pool_", ssp, "_3m_lim_std.tif"),
  overwrite = TRUE
)

writeRaster(
  rast(slow_pool_list),
  paste0("1.9%_min/3m_lim/thawed_organic_N_pool_", ssp, "_3m_lim_std.tif"),
  overwrite = TRUE
)

message("Finished SSP", ssp)