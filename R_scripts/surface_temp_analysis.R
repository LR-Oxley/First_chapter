library(terra)
library(ggplot2)
library(dplyr)
library(purrr)
library(zoo)
library(readr)


# ============================================================
# Pan-Arctic (60-90 N) surface temperature from CMIP6
#   - AREA-WEIGHTED annual means (grid cells shrink towards the pole)
#   - anomaly relative to 2000-2020, 20-year running mean, plot
#   - historical warming (1850-1900 vs 2000-2020) and 2080-2099 anomaly
#   - comparison with Biskaborn et al. (2019): trend over the reference
#     decade 2007-2016, anomalies relative to 1981-2010
# ============================================================


# ============================================================
# Setup
# ============================================================

ext_60N <- ext(-180, 180, 60, 90)

ssp_files <- list(
  "SSP1-2.6" = "60N/mean_ST_126_clean_final.nc",
  "SSP2-4.5" = "60N/mean_ST_245_clean_final.nc",
  "SSP3-7.0" = "60N/mean_ST_370_clean_final.nc",
  "SSP5-8.5" = "60N/mean_ST_585_clean_final.nc"
)

ssp_colors <- c("SSP1-2.6" = "blue", "SSP2-4.5" = "orange",
                "SSP3-7.0" = "#D73027", "SSP5-8.5" = "#7B3294")

# periods
base_period   <- c(2000, 2020)   # anomaly reference for the plot
preind_period <- c(1850, 1900)   # preindustrial
future_period <- c(2080, 2099)   # end of century

# Biskaborn et al. (2019): 0.86 +/- 0.84 degC per decade, continuous
# permafrost zone, reference decade 2007-2016, anomalies vs 1981-2010
bisk_ref_period   <- c(1981, 2010)
bisk_trend_period <- c(2007, 2016)
biskaborn         <- c(mean = 0.86, sd = 0.84)

# optional: raster of the continuous permafrost zone (1 = continuous),
# used ONLY for the Biskaborn comparison. NULL = all cells 60-90 N.
permafrost_mask_file <- NULL


# ============================================================
# Helpers
# ============================================================

in_p <- function(y, p) y >= p[1] & y <= p[2]

# linear trend per year (NA-safe), used per grid cell
slope_fun <- function(y, t) {
  ok <- is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  t <- t[ok]; y <- y[ok]
  sum((t - mean(t)) * (y - mean(y))) / sum((t - mean(t))^2)
}

wmean <- function(x, w) { ok <- is.finite(x) & is.finite(w); sum(x[ok] * w[ok]) / sum(w[ok]) }
wsd   <- function(x, w) { ok <- is.finite(x) & is.finite(w); m <- wmean(x, w)
sqrt(sum(w[ok] * (x[ok] - m)^2) / sum(w[ok])) }


# ============================================================
# Function: process one SSP file
#   -> annual area-weighted pan-Arctic mean (degC) + Biskaborn trends
# ============================================================

process_ssp <- function(file_path, scenario_name) {
  
  cat("Processing", scenario_name, "\n")
  
  r <- rast(file_path)
  r <- crop(r, ext_60N)
  r <- r - 273.15                                  # Kelvin -> Celsius
  
  years  <- as.numeric(format(time(r), "%Y"))
  annual <- tapp(r, index = years, fun = "mean", na.rm = TRUE)
  yrs    <- sort(unique(years))
  
  # cell areas as weights
  area <- cellSize(annual[[1]], unit = "km")
  
  # --- area-weighted annual pan-Arctic mean ---
  ts <- data.frame(
    Year       = yrs,
    MeanTemp_C = global(annual, fun = "mean", weights = area, na.rm = TRUE)[, 1],
    Scenario   = scenario_name
  )
  
  # --- Biskaborn comparison (optionally continuous permafrost only) ---
  annual_b <- annual
  if (!is.null(permafrost_mask_file)) {
    pf <- resample(crop(rast(permafrost_mask_file), ext_60N), annual[[1]], method = "near")
    annual_b <- mask(annual_b, pf, maskvalues = c(0, NA))
  }
  
  ref_idx   <- which(in_p(yrs, bisk_ref_period))
  trend_idx <- which(in_p(yrs, bisk_trend_period))
  t         <- yrs[trend_idx]
  anom_b    <- annual_b - mean(annual_b[[ref_idx]], na.rm = TRUE)
  
  # pan-Arctic trend of the area-weighted mean anomaly
  pan <- global(anom_b[[trend_idx]], "mean", weights = area, na.rm = TRUE)[, 1]
  fit <- lm(pan ~ t)
  
  # trend in every grid cell -> area-weighted mean and SD across cells
  cell_slope <- app(anom_b[[trend_idx]], function(y) slope_fun(y, t)) * 10
  vals <- as.data.frame(c(cell_slope, area), na.rm = TRUE)
  names(vals) <- c("slope", "area")
  
  trend <- data.frame(
    Scenario            = scenario_name,
    pan_trend_C_decade  = unname(coef(fit)[2]) * 10,
    pan_trend_se        = summary(fit)$coefficients[2, 2] * 10,
    cell_trend_mean     = wmean(vals$slope, vals$area),
    cell_trend_sd       = wsd(vals$slope, vals$area),
    anomaly_2007_2016_C = mean(pan)
  )
  
  list(ts = ts, trend = trend)
}


# ============================================================
# Apply to all four SSP files and combine
# ============================================================

results <- imap(ssp_files, process_ssp)

df_all <- map(results, "ts")    |> bind_rows()
trends <- map(results, "trend") |> bind_rows()


# ============================================================
# Anomaly relative to 2000-2020 (per scenario), period, 20-yr running mean
# ============================================================

long_data <- df_all %>%
  group_by(Scenario) %>%
  arrange(Year, .by_group = TRUE) %>%
  mutate(
    baseline_mean     = mean(MeanTemp_C[in_p(Year, base_period)], na.rm = TRUE),
    Anomaly_C         = MeanTemp_C - baseline_mean,
    Rolling_Anomaly_C = rollmean(Anomaly_C, k = 20, fill = NA, align = "center"),
    Period            = ifelse(Year <= 2014, "Before 2015", "After 2015")
  ) %>%
  ungroup() %>%
  mutate(Scenario = factor(Scenario, levels = names(ssp_files)))

write_csv(long_data, "surface_temp_pan_arctic_area_weighted.csv")


# ============================================================
# Historical warming (1850-1900 vs 2000-2020) and future anomaly
# ============================================================

historical_increase <- df_all %>%
  group_by(Scenario) %>%
  summarise(
    mean_preindustrial = mean(MeanTemp_C[in_p(Year, preind_period)], na.rm = TRUE),
    mean_present       = mean(MeanTemp_C[in_p(Year, base_period)],   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(increase_C = mean_present - mean_preindustrial,
         label      = sprintf("%.2f \u00B0C", increase_C))

summary_future_period <- long_data %>%
  filter(in_p(Year, future_period)) %>%
  group_by(Scenario) %>%
  summarise(period_mean_anomaly = mean(Anomaly_C, na.rm = TRUE), .groups = "drop") %>%
  mutate(label = sprintf("%.2f \u00B0C", period_mean_anomaly))

cat("\n=== Historical warming, 1850-1900 to 2000-2020 (area-weighted) ===\n")
print(historical_increase)
cat("\n=== Anomaly 2080-2099 relative to 2000-2020 (area-weighted) ===\n")
print(summary_future_period)


# ============================================================
# Comparison with Biskaborn et al. (2019)
# ============================================================

comparison <- trends %>%
  mutate(
    this_study = sprintf("pan-Arctic %.2f \u00b1 %.2f; cells %.2f \u00b1 %.2f \u00B0C/decade",
                         pan_trend_C_decade, pan_trend_se, cell_trend_mean, cell_trend_sd),
    biskaborn_text = sprintf("%.2f \u00b1 %.2f \u00B0C/decade",
                             .env$biskaborn["mean"], .env$biskaborn["sd"]),
    within_1sd = abs(cell_trend_mean - .env$biskaborn["mean"]) <= .env$biskaborn["sd"]
  )

cat("\n=== Trend 2007-2016 (anomalies vs 1981-2010) vs Biskaborn et al. (2019) ===\n")
print(tibble::as_tibble(comparison) %>%
        select(Scenario, this_study, biskaborn_text, within_1sd),
      width = Inf)
write_csv(trends, "temperature_trend_2007_2016_vs_biskaborn.csv")


# ============================================================
# Plot
# ============================================================

temp_plot <- ggplot(long_data, aes(x = Year, y = Anomaly_C, group = Scenario)) +
  
  # Before 2015: black line
  geom_line(data = filter(long_data, Period == "Before 2015"),
            color = "black", linewidth = 0.2) +
  
  # After 2015: coloured lines
  geom_line(data = filter(long_data, Period == "After 2015"),
            aes(color = Scenario), linewidth = 0.2) +
  
  # Rolling mean before 2015 (black)
  geom_line(data = filter(long_data, Year < 2015),
            aes(y = Rolling_Anomaly_C), color = "black", linewidth = 0.2, na.rm = TRUE) +
  
  # Rolling mean after 2015 (coloured)
  geom_line(data = filter(long_data, Year >= 2015),
            aes(y = Rolling_Anomaly_C, color = Scenario), linewidth = 0.2, na.rm = TRUE) +
  
  geom_vline(xintercept = 2015, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  
  labs(x = "Year", y = "Surface temperature anomaly [\u00B0C]", title = "", color = "") +
  scale_color_manual(values = ssp_colors) +
  
  theme_minimal(base_size = 15) +
  theme(
    legend.position   = "bottom",
    legend.margin     = margin(t = 0),
    legend.box.margin = margin(1, 1, 1, 1),
    plot.margin       = margin(3, 3, 2, 3),
    legend.key.size   = unit(0.3, "cm"),
    legend.text       = element_text(size = 5),
    legend.title      = element_text(size = 5),
    axis.title        = element_text(size = 6),
    axis.text         = element_text(size = 6),
    axis.title.x      = element_text(margin = margin(t = 2)),
    axis.title.y      = element_text(margin = margin(r = 2)),
    axis.text.x       = element_text(margin = margin(t = 1)),
    axis.text.y       = element_text(margin = margin(r = 1))
  )

print(temp_plot)

ggsave(
  "plot_surface_temp.pdf",
  temp_plot,
  width = 8.5,
  height = 8,
  units = "cm",
  dpi = 500
)