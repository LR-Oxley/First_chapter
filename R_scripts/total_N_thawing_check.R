
# ============================================================
# Total active-layer N stock:
#   1) Present day, CMIP6 (2000-2020 mean, NOT bias-corrected)
#   2) Present day, ESA observational reference (for comparison)
#   3) Present day, CMIP6 bias-corrected against ESA
#   4) Hypothetical: ALD = 3 m everywhere (also the profile's
#      calibration depth -- a useful internal sanity check)
#   5) Hypothetical: ALD = 5 m everywhere
#
# thawed N (state) = f(ALD, TN_profile)
# ============================================================

library(terra)

# ------------------------------------------------------------
# 1. Settings
# ------------------------------------------------------------

common_extent <- ext(-179.95, 179.95, 60, 90)
f_inorg <- 0.01

params <- list(
  taiga  = c(a = 0.007, b = 0.097, k = 2.7),
  tundra = c(a = 0.01,  b = 0.017, k = 1.9),
  barren = c(a = 0,     b = 0.0161, k = 1.6)
)


# ------------------------------------------------------------
# 2. Load rasters
# ------------------------------------------------------------

ALD_raw <- rast("60N/mean_ALD_ssp585_60deg.nc")
#ALD_raw <- rast("mean_ssp585_corr_new.nc") old values
ALD_raw <- crop(ALD_raw, common_extent)
plot(ALD_raw[[2]])
n_layers <- nlyr(ALD_raw)
years <- seq(1850, by = 1, length.out = n_layers)
present_idx <- which(years >= 2000 & years <= 2020)

N_data <- rast("TN_30deg_corr.nc", lyr = 1)
N_data <- crop(N_data, common_extent)

LC <- rast("LC_remapnn_corr.nc")
LC <- crop(LC, common_extent)


# ------------------------------------------------------------
# 3. Present-day ALD, two sources, NO bias correction
# ------------------------------------------------------------
# Source 1: CMIP6 model climatology, 2000-2020 mean, as-is
ALD_present_cmip6 <- mean(ALD_raw[[present_idx]], na.rm = TRUE)

# Source 2: ESA observational reference (independent of the model)
ALD_present_esa <- rast("ESA_ALD/ESA_ALD_present_day.nc", lyr = 1)   # lyr 1 = mean, lyr 2 = std
ALD_present_esa <- crop(ALD_present_esa, common_extent)

# Source 3: CMIP6, bias-corrected against the ESA reference
# (shift the FULL model time series by a constant offset, calibrated so
# the model's own 2000-2020 mean matches the ESA present-day value, then
# take the present-day mean of that corrected series)
#
# NOTE: because the correction is calibrated using the present-day period
# itself, ALD_present_cmip6_corr's PRESENT-DAY MEAN would be mathematically
# forced to exactly equal ALD_present_esa if it weren't for the negative-value
# clipping applied per-year (ALD_raw_corr[ALD_raw_corr < 0] <- 0) BEFORE
# averaging -- any individual year within 2000-2020 that the correction would
# push below zero gets floored first, which nudges the resulting mean slightly
# above the pure algebraic ESA match. Expect this row's present-day total to
# come out very close to, but not perfectly identical to, Source 2 (ESA).
# It only becomes a genuinely different estimate at OTHER time periods (e.g.
# future decades), where the additive bias correction carries real information.
if (!compareGeom(ALD_present_esa, ALD_present_cmip6, stopOnError = FALSE)) {
  ALD_present_esa_for_bias <- resample(ALD_present_esa, ALD_present_cmip6, method = "bilinear")
} else {
  ALD_present_esa_for_bias <- ALD_present_esa
}

bias <- ALD_present_cmip6 - ALD_present_esa_for_bias

ALD_raw_corr <- ALD_raw - bias
ALD_raw_corr[ALD_raw_corr < 0] <- 0

ALD_present_cmip6_corr <- mean(ALD_raw_corr[[present_idx]], na.rm = TRUE)


# ------------------------------------------------------------
# 4. Align grids (resample onto N_data's grid, don't just force
#    extent metadata -- ext<- alone does not reproject/regrid)
# ------------------------------------------------------------

if (!compareGeom(LC, N_data, stopOnError = FALSE)) {
  LC <- resample(LC, N_data, method = "near")
}

if (!compareGeom(ALD_present_cmip6, N_data, stopOnError = FALSE)) {
  ALD_present_cmip6 <- resample(ALD_present_cmip6, N_data, method = "bilinear")
}

if (!compareGeom(ALD_present_esa, N_data, stopOnError = FALSE)) {
  ALD_present_esa <- resample(ALD_present_esa, N_data, method = "bilinear")
}

if (!compareGeom(ALD_present_cmip6_corr, N_data, stopOnError = FALSE)) {
  ALD_present_cmip6_corr <- resample(ALD_present_cmip6_corr, N_data, method = "bilinear")
}


# ------------------------------------------------------------
# 5. Land-cover masks
# ------------------------------------------------------------

taiga_mask    <- LC %in% c(1, 2, 3, 4, 5, 8, 9)
tundra_mask   <- LC %in% c(6, 7, 10, 12, 14)
wetlands_mask <- LC == 11
barren_mask   <- LC %in% c(13, 15, 16)

land_mask <- taiga_mask | tundra_mask | wetlands_mask | barren_mask


# ------------------------------------------------------------
# 6. Normalise the measured N stock to a per-metre profile density
#    (dividing out the known 3 m reference-depth integral)
# ------------------------------------------------------------

normalize_N <- function(N, a, b, k) {
  A_3m <- 3 * a + (b / k) * (1 - exp(-3 * k))
  N / A_3m
}

taiga_N    <- normalize_N(N_data * taiga_mask,  params$taiga["a"],  params$taiga["b"],  params$taiga["k"])
tundra_N   <- normalize_N(N_data * tundra_mask, params$tundra["a"], params$tundra["b"], params$tundra["k"])
barren_N   <- normalize_N(N_data * barren_mask, params$barren["a"], params$barren["b"], params$barren["k"])
wetlands_N <- (N_data * wetlands_mask) / 3


# ------------------------------------------------------------
# 7. Thaw functions
# ------------------------------------------------------------

compute_thawed_N <- function(ALD, N, a, b, k, f_inorg) {
  A_ALD <- a * ALD + (b / k) * (1 - exp(-k * ALD))
  total_thawed <- ifel(N == 0, NA, N * A_ALD)
  total_thawed / (1 - f_inorg)
}

compute_thawed_N_wetlands <- function(ALD, N, f_inorg) {
  total_thawed <- ifel(N == 0, NA, N * ALD)
  total_thawed / (1 - f_inorg)
}


# ------------------------------------------------------------
# 8. Wrapper: total N stock [Pg] for any given ALD raster/scenario
# ------------------------------------------------------------

compute_total_N_pg <- function(ALD, label) {
  
  thawed_taiga <- compute_thawed_N(
    ALD, taiga_N,
    params$taiga["a"], params$taiga["b"], params$taiga["k"],
    f_inorg
  )
  
  thawed_tundra <- compute_thawed_N(
    ALD, tundra_N,
    params$tundra["a"], params$tundra["b"], params$tundra["k"],
    f_inorg
  )
  
  thawed_barren <- compute_thawed_N(
    ALD, barren_N,
    params$barren["a"], params$barren["b"], params$barren["k"],
    f_inorg
  )
  
  thawed_wetlands <- compute_thawed_N_wetlands(ALD, wetlands_N, f_inorg)
  
  thawed_taiga[is.na(thawed_taiga)]       <- 0
  thawed_tundra[is.na(thawed_tundra)]     <- 0
  thawed_barren[is.na(thawed_barren)]     <- 0
  thawed_wetlands[is.na(thawed_wetlands)] <- 0
  
  combined <- thawed_taiga + thawed_tundra + thawed_barren + thawed_wetlands
  combined <- ifel(combined == 0, NA, combined)
  combined <- mask(combined, land_mask, maskvalues = 0)
  
  weights <- cellSize(combined, unit = "m")
  
  total_mass <- global(combined * weights, "sum", na.rm = TRUE)[1, 1]
  
  total_pg <- total_mass / 1e12
  
  cat(label, ":", round(total_pg, 4), "Pg N\n")
  
  list(raster = combined, total_pg = total_pg)
}


# ------------------------------------------------------------
# 9. Run the scenarios
# ------------------------------------------------------------

present_day_cmip6 <- compute_total_N_pg(
  ALD_present_cmip6,
  "Present-day, CMIP6 (2000-2020 mean, not bias-corrected)"
)
plot(present_day_cmip6$raster)

present_day_esa <- compute_total_N_pg(
  ALD_present_esa,
  "Present-day, ESA observational reference"
)
plot(present_day_esa$raster)

present_day_cmip6_corr <- compute_total_N_pg(
  ALD_present_cmip6_corr,
  "Present-day, CMIP6 bias-corrected against ESA"
)
plot(present_day_cmip6_corr$raster)

# Sanity check: ALD = 3 m matches the profile's calibration depth, so this
# run should closely reproduce the original N_data inventory. If it doesn't,
# something upstream (grid alignment, masking, units) needs another look
# before trusting the 5 m scenario.
ALD_3m <- ifel(land_mask, 3, NA)
hyp_3m <- compute_total_N_pg(ALD_3m, "Hypothetical: ALD = 3 m everywhere")

ALD_5m <- ifel(land_mask, 5, NA)
hyp_5m <- compute_total_N_pg(ALD_5m, "Hypothetical: ALD = 5 m everywhere")


# ------------------------------------------------------------
# 10. Summary
# ------------------------------------------------------------

summary_df <- data.frame(
  Scenario = c(
    "Present-day, CMIP6 (not bias-corrected)",
    "Present-day, ESA reference",
    "Present-day, CMIP6 (bias-corrected against ESA)",
    "Uniform ALD = 3 m",
    "Uniform ALD = 5 m"
  ),
  Total_N_Pg = c(
    present_day_cmip6$total_pg,
    present_day_esa$total_pg,
    present_day_cmip6_corr$total_pg,
    hyp_3m$total_pg,
    hyp_5m$total_pg
  )
)

print(summary_df)




library(dplyr)
series <- readr::read_csv(file.path("monthly_mineralised/",
                                    "diagnostic_csv_analysis_60N",
                                    "pan_arctic_annual_from_csv.csv"))

series %>%
  filter(run == "mean", Year >= 2000, Year <= 2020,
         var %in% c("total_inorg_pg_monthly", "organic_pool_remaining_pg")) %>%
  group_by(SSP, var) %>%
  summarise(value = mean(value), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = var, values_from = value) %>%
  mutate(effective_rate_pct = 100 * total_inorg_pg_monthly / organic_pool_remaining_pg)
