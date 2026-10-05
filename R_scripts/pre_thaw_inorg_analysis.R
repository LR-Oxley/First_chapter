
# ============================================================
# Pre-thaw inorganic N release (inorg_rapid_available):
# present day (2000-2020) and 2080-2099, absolute and per m2
# uncertainty = |plus - minus| / 2
# ============================================================

library(dplyr)
library(tidyr)
library(readr)

# annual table written by plot_diagnostics_csv.R (Step 3)
annual_file <- file.path("monthly_mineralised/",
                         "diagnostic_csv_analysis_60N",
                         "pan_arctic_annual_from_csv.csv")

# domain area for the per-m2 values: same land mask as the model
domain_file <- file.path("total_thawed_extended", "arctic_total_thawed_126_60N_mean.nc")

r_dom <- terra::crop(terra::rast(domain_file), terra::ext(-179.95, 179.95, 60, 90))
land  <- !is.na(r_dom[[terra::nlyr(r_dom)]])
domain_area_m2 <- terra::global(terra::mask(terra::cellSize(r_dom[[1]], unit = "m"), land, maskvalues = 0),
                                "sum", na.rm = TRUE)[1, 1]
cat("Domain area [10^6 km2]:", round(domain_area_m2 / 1e12, 2), "\n")

annual <- read_csv(annual_file, show_col_types = FALSE)

in_p <- function(y, a, b) y >= a & y <= b

prethaw <- annual %>%
  filter(var == "inorg_rapid_available") %>%          # Pg N yr-1 (sum of 12 months)
  group_by(SSP, run) %>%
  summarise(
    present = mean(value[in_p(Year, 2000, 2020)]),
    future  = mean(value[in_p(Year, 2080, 2099)]),
    .groups = "drop"
  ) %>%
  mutate(anomaly = future - present) %>%
  pivot_longer(c(present, future, anomaly), names_to = "period", values_to = "PgN_yr") %>%
  pivot_wider(names_from = run, values_from = PgN_yr) %>%
  mutate(
    sd      = abs(plus - minus) / 2,
    Tg_yr   = mean * 1000,
    Tg_sd   = sd   * 1000,
    g_m2_yr = mean * 1e15 / domain_area_m2,
    g_m2_sd = sd   * 1e15 / domain_area_m2,
    text_Tg = sprintf("%.2f \u00b1 %.2f Tg N/yr", Tg_yr, Tg_sd),
    text_g  = sprintf("%.3f \u00b1 %.3f g N/m2/yr", g_m2_yr, g_m2_sd),
    period  = factor(period, levels = c("present", "future", "anomaly"))
  ) %>%
  arrange(period, SSP)

print(prethaw %>% select(period, SSP, text_Tg, text_g), n = Inf, width = Inf)
write_csv(prethaw, "prethaw_inorganic_N_present_future.csv")

