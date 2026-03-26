library(terra)


# -----------------------------
# 0) ARRAY TASK ID → SSP
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)
task_id <- if (length(args) > 0) as.integer(args[1]) else as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "0"))

ssp_levels <- c("126", "245", "370", "585")
if (!(task_id %in% 0:3)) stop("SLURM_ARRAY_TASK_ID must be 0–3")

ssp <- ssp_levels[task_id + 1]
message("Running SSP", ssp)

# -----------------------------
# 1) INPUTS (scenario-dependent)
# -----------------------------
arctic_path <- paste0("1.9%_min/3m_lim/mineralised_N_pool_", ssp, "_3m_lim_std.tif")
arctic <- rast(arctic_path)

LC <- rast("LC_remapnn_corr.nc")
common_extent <- ext(-179.95, 179.95, 30, 90)
ext(LC) <- common_extent
LC <- resample(LC, arctic, method = "near")

temp_path <- paste0("mean_ST_", ssp, "_clean_final.nc")
temp_data <- rast(temp_path)

ext(temp_data) <- common_extent
crs(temp_data) <- crs(arctic)
temp_C <- temp_data - 273.15
rm(temp_data)

# -----------------------------
# 2) YEARS + LAND MASK
# -----------------------------
years <- 1850:2099
names(arctic) <- paste0("Y", years)
land_mask <- !is.na(arctic[[1]])

# -----------------------------
# 3) LAND-COVER MASKS
# -----------------------------
taiga_mask    <- LC %in% c(1,2,3,4,5,8,9)
tundra_mask   <- LC %in% c(6,7,10)
wetlands_mask <- LC == 11
barren_mask   <- LC %in% c(15,16)

taiga_min    <- mask(arctic, taiga_mask,    maskvalue = 0)
tundra_min   <- mask(arctic, tundra_mask,   maskvalue = 0)
wetlands_min <- mask(arctic, wetlands_mask, maskvalue = 0)
barren_min   <- mask(arctic, barren_mask,   maskvalue = 0)
rm(arctic, LC)

# -----------------------------
# 4) TEMPERATURE MULTIPLIER (peaked Arrhenius)
# -----------------------------
Ea_denitrification <- 47000
Ed_denitrification <- 200000
t_opt_C <- 28

peaked_arrhenius <- function(temp_C, Ea, Ed, t_opt_C = 28) {
  R_gas <- 8.314462618
  temp_K  <- temp_C + 273.15
  t_opt_K <- t_opt_C + 273.15
  hlp1 <- temp_K - t_opt_K
  hlp2 <- temp_K * t_opt_K * R_gas
  numerator   <- Ed * exp(Ea * hlp1 / hlp2)
  denominator <- Ed - Ea * (1 - exp(Ed * hlp1 / hlp2))
  numerator / denominator
}

k_T_year <- peaked_arrhenius(temp_C, Ea_denitrification, Ed_denitrification, t_opt_C)
k_T_year <- mask(k_T_year, land_mask, maskvalues = FALSE)
names(k_T_year) <- paste0("Y", years)

ref_idx <- which(years >= 2000 & years <= 2020)
k_T_ref <- mean(k_T_year[[ref_idx]], na.rm = TRUE)
k_factor_stack <- k_T_year / k_T_ref

k_factor_stack <- ifel(
  land_mask,
  ifel(is.na(k_factor_stack), 1, k_factor_stack),
  NA
)

rm(temp_C, k_T_year, k_T_ref)

# -----------------------------
# 5) UNIT CONVERSION (kg → mg)
# -----------------------------
tundra_mg   <- tundra_min   * 1e6
taiga_mg    <- taiga_min    * 1e6
wetlands_mg <- wetlands_min * 1e6
barren_mg   <- barren_min   * 1e6
rm(tundra_min, wetlands_min, barren_min)

# -----------------------------
# 6) EMISSION FACTORS (Voigt)
# -----------------------------
season_days <- 100

ef_dry_veg  <- 0.0001113896 * season_days
ef_wet_veg  <- 0.0001314326 * season_days
ef_dry_bare <- 0.00369822   * season_days


# only affects nitrate part of inorganic N!
n2o_tundra <- (tundra_mg  *0.1) * ef_dry_veg  * k_factor_stack
n2o_taiga  <- (taiga_mg  *0.1)  * ef_dry_veg  * k_factor_stack
n2o_wet    <- (wetlands_mg *0.1)* ef_wet_veg  * k_factor_stack
n2o_bare   <- (barren_mg *0.1)  * ef_dry_bare * k_factor_stack

rm(tundra_mg, wetlands_mg, barren_mg)

# -----------------------------
# 7) TOTAL + OPTIONAL WINTER SCALING
# -----------------------------
n2o_tundra[is.na(n2o_tundra)] <- 0
n2o_taiga[is.na(n2o_taiga)] <- 0
n2o_wet[is.na(n2o_wet)] <- 0
n2o_bare[is.na(n2o_bare)] <- 0

n2o_total <- n2o_taiga+ n2o_tundra  + n2o_wet + n2o_bare 
n2o_total <- ifel(n2o_total == 0, NA, n2o_total)

winter_factor <- 2
n2o_total_annual <- n2o_total * winter_factor

# -----------------------------
# 8) OUTPUT (scenario-tagged)
# -----------------------------
dir.create("n2o", showWarnings = FALSE, recursive = TRUE)

writeRaster(
  n2o_total,
  paste0("n2o/arctic_total_n2o_ssp", ssp, "_std_new.tif"),
  overwrite = TRUE
)

writeRaster(
  n2o_total_annual,
  paste0("n2o/arctic_total_n2o_ssp", ssp, "_std_annual_new.tif"),
  overwrite = TRUE
)

message("Finished SSP", ssp)