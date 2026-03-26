library(terra)

# -----------------------------
# 1) INPUTS
# -----------------------------
arctic <- rast("1.9%_min/3m_lim/mineralised_N_pool_370_3m_lim_std.tif")  # cumulative since 1850 (units as in your file)

LC <- rast("LC_remapnn_corr.nc")
common_extent <- ext(-179.95, 179.95, 30, 90)
ext(LC) <- common_extent
LC <- resample(LC, arctic, method = "near")

# Temperature (K -> °C)
temp_data <- rast("mean_ST_370_clean_final.nc")
ext(temp_data) <- common_extent
crs(temp_data) <- crs(arctic)
temp_C <- temp_data - 273.15
rm(temp_data)
# -----------------------------
# 2) YEARS + CONVERT CUMULATIVE N -> ANNUAL INCREMENT (Option 1)
# -----------------------------

years <- 1850:2099
names(arctic) <- paste0("Y", years)

# Land footprint (fixed): use first increment layer
land_mask <- !is.na(arctic[[1]])

# -----------------------------
# 3) LAND-COVER MASKS
# -----------------------------
taiga_mask    <- LC %in% c(1,2,3,4,5,8,9)
tundra_mask   <- LC %in% c(6,7,10)
wetlands_mask <- LC == 11
barren_mask   <- LC %in% c(15,16)

# mask() keeps values where mask is TRUE, sets others to NA
taiga_min    <- mask(arctic, taiga_mask,    maskvalue = 0)
tundra_min   <- mask(arctic, tundra_mask,   maskvalue = 0)
wetlands_min <- mask(arctic, wetlands_mask, maskvalue = 0)
barren_min   <- mask(arctic, barren_mask,   maskvalue = 0)
rm(arctic)

# -----------------------------
# 4) TEMPERATURE MULTIPLIER (peaked Arrhenius), BASELINED TO 2000–2020
#    IMPORTANT: align to 1851–2099 to match dN years
# -----------------------------
Ea_denitrification <- 47000   # J/mol
Ed_denitrification <- 200000  # J/mol
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

k_T_year_full <- peaked_arrhenius(temp_C, Ea_denitrification, Ed_denitrification, t_opt_C = t_opt_C)
k_T_year_full <- mask(k_T_year_full, land_mask, maskvalues = FALSE)

# Use same years as dN (1851–2099)
k_T_year <- k_T_year_full[[1:nlyr(k_T_year_full)]]
names(k_T_year) <- paste0("Y", years)

# baseline 2000–2020 within 1851–2099 range
ref_idx <- which(years >= 2000 & years <= 2020)
k_T_ref <- mean(k_T_year[[ref_idx]], na.rm = TRUE)

k_factor_stack <- k_T_year / k_T_ref

# make NA temps on land become 1, keep ocean NA

rm(k_T_year_full, k_T_year, k_T_ref)
# -----------------------------
# 5) UNIT CONVERSION: thawed/mineralised N -> mg
#    NOTE: This assumes your arctic raster is in kg (of N) per unit area.
#    You currently use *1e6. Keep it as in your script, but double-check your units.
# -----------------------------
thawed_tundra_bioavailable_mg   <- tundra_min   * 1e6
thawed_taiga_bioavailable_mg    <- taiga_min    * 1e6
thawed_wetlands_bioavailable_mg <- wetlands_min * 1e6
thawed_barren_bioavailable_mg   <- barren_min   * 1e6

rm(tundra_min, taiga_min, wetlands_min, barren_min)

# -----------------------------
# 6) EMISSION FACTORS (Voigt) + ANNUAL (100 season days) + TEMP MULTIPLIER
# -----------------------------
dry_vegetated_rate_1m_mean <- 0.0001113896
wet_vegetated_rate_1m_mean <- 0.0001314326
dry_bare_rate_1m_mean      <- 0.00369822

emission_factor_year_dry_vegetated <- dry_vegetated_rate_1m_mean * 100
emission_factor_year_wet_vegetated <- wet_vegetated_rate_1m_mean * 100
emission_factor_year_dry_bare      <- dry_bare_rate_1m_mean      * 100

# mg N2O m-2 yr-1 (per year emissions driven by annual increment pool)
n2o_tundra_per_m2 <- thawed_tundra_bioavailable_mg   * emission_factor_year_dry_vegetated * k_factor_stack
n2o_taiga_per_m2  <- thawed_taiga_bioavailable_mg    * emission_factor_year_dry_vegetated * k_factor_stack
n2o_wet_per_m2    <- thawed_wetlands_bioavailable_mg * emission_factor_year_wet_vegetated * k_factor_stack
n2o_bare_per_m2   <- thawed_barren_bioavailable_mg   * emission_factor_year_dry_bare      * k_factor_stack

# add tundra and taiga together as both are dry vegetated
rm(thawed_tundra_bioavailable_mg,
   thawed_wetlands_bioavailable_mg, thawed_barren_bioavailable_mg)
n2o_tundra_per_m2[is.na(n2o_tundra_per_m2)] <- 0
n2o_taiga_per_m2[is.na(n2o_taiga_per_m2)] <- 0
n2o_wet_per_m2[is.na(n2o_wet_per_m2)] <- 0
n2o_bare_per_m2[is.na(n2o_bare_per_m2)] <- 0
plot(n2o_bare_per_m2[[2]])
n2o_total <- n2o_tundra_per_m2 + n2o_wet_per_m2+ n2o_bare_per_m2
n2o_total <- ifel(n2o_total == 0, NA, n2o_total)
plot(n2o_total[[2]])
n2o_total_annual <- n2o_total * 2
writeRaster(n2o_total_annual, "n2o/arctic_total_n2o_370_std_annual_new.tif", overwrite = TRUE)
writeRaster(n2o_total, "n2o/arctic_total_n2o_370_std_new.tif", overwrite = TRUE)


# rm(tundra_thawed_total,barren_thawed_total,wetlands_thawed_total, taiga_thawed_total)
# rm(n2o_tundra_per_m2, n2o_taiga_per_m2, n2o_wet_per_m2, n2o_bare_per_m2)
# 
# # -----------------------------
# # 7) CONVERT TO PER-CELL TOTALS (mg/cell/yr) — MULTIPLY BY AREA ONCE
# # -----------------------------
# area_m2<-cellSize(n2o_total, mask = TRUE, unit = "m")
# plot(n2o_total[[2]])
# 
# n2o_total_area<-n2o_total*area_m2
# 
# # Optional winter doubling (as in your original script)
# n2o_cell_total_annual <- n2o_total_area * 1
# plot(n2o_cell_total_annual[[2]])
# rm(n2o_cell_total)
# 
# names(n2o_cell_total_annual) <- paste0("Y", years)
# 
# # -----------------------------
# # 8) ARCTIC TOTALS TIME SERIES (Tg N / yr) + ANOMALY RELATIVE TO 2000–2020
# # -----------------------------
# total_area_emissions_mg_annual <- global(n2o_cell_total_annual, "sum", na.rm = TRUE)[,1]
# 
# # Convert mg N2O -> Tg N (N2O-N) using 0.636 and mg->Tg (1e15 mg = 1 Tg)
# total_area_emissions_Tg_Neq <- total_area_emissions_mg_annual * 0.636 / 1e15
# 
# df_out <- data.frame(
#   years = years,
#   TgN = total_area_emissions_Tg_Neq
# )
# 
# baseline_mean <- mean(df_out$TgN[df_out$years %in% 2000:2020], na.rm = TRUE)
# df_out$anomaly <- df_out$TgN - baseline_mean
# 
# write.csv(df_out, "n2o-n_emissions_Tg_mean_370.csv", row.names = FALSE)
