# ============================================================
# Site comparison with Beermann et al. (2017):
# pre-thaw inorganic N release and ALD increase at three sites in
# eastern Siberia (LEN, IND, KOL), extracted from the model output.
#
# - values at the grid cell containing each site, and as the mean over
#   a buffer around it (delta / river sites often fall on water or
#   coastline cells at 0.1 deg resolution)
# - N release converted to mg N m-2 yr-1 (unit used by Beermann et al.)
# - period means and ALD change, plus time series plots
# ============================================================

library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(zoo)
library(readr)


# ------------------------------------------------------------
# 1. Settings
# ------------------------------------------------------------
ssp <- "585"

# pre-thaw inorganic N release, one layer per year, kg N m-2 yr-1
rapid_file <- file.path("monthly_mineralised/mean_2perc_baseline",
                        sprintf("arctic_rapid_inorg_N_yearly_%s_mean_1850_2099_w_temp.nc", ssp))
# bias-corrected ALD (m), one layer per year; use the file that drove the model
ald_file <- file.path("60N", sprintf("mean_ALD_ssp%s_60deg.nc", ssp))
# if you don't have the corrected file, the raw multi-model mean:
# ald_file <- sprintf("60N/mean_ALD_ssp%s_60deg.nc", ssp)

start_year <- 1850         # used if the files carry no time information
buffer_m   <- 25000        # radius of the buffer around each site (m)

present <- c(2000, 2020)
future  <- c(2080, 2099)
proj    <- c(2015, 2099)   # projection period, for a mean annual release

sites <- tibble(
  site = c("LEN", "IND", "KOL"),
  name = c("Samoylov Island, Lena Delta", "Indigirka Lowlands (Kytalyk)", "Kolyma Delta (Pokhodsk)"),
  lat  = c(72.3676, 70.8296, 69.0790),
  lon  = c(126.4838, 147.4895, 160.9634)
)

out_dir <- "site_comparison_beermann"
dir.create(out_dir, showWarnings = FALSE)


# ------------------------------------------------------------
# 2. Helpers
# ------------------------------------------------------------
layer_years <- function(r) {
  tt <- time(r)
  yrs <- if (inherits(tt, "Date") || inherits(tt, "POSIXt")) as.integer(format(tt, "%Y")) else as.integer(tt)
  if (length(yrs) != nlyr(r) || any(is.na(yrs))) yrs <- start_year + seq_len(nlyr(r)) - 1
  yrs
}

# extract a yearly stack at the site cells and as buffer means -> long table
extract_sites <- function(r, value_name) {
  yrs <- layer_years(r)
  pts <- vect(sites, geom = c("lon", "lat"), crs = "EPSG:4326")
  
  at_cell <- terra::extract(r, pts, ID = FALSE)
  in_buf  <- terra::extract(r, buffer(pts, width = buffer_m), fun = mean, na.rm = TRUE, ID = FALSE)
  
  to_long <- function(m, type) {
    as_tibble(m) %>%
      setNames(as.character(yrs)) %>%
      mutate(site = sites$site, type = type) %>%
      pivot_longer(-c(site, type), names_to = "Year", values_to = value_name) %>%
      mutate(Year = as.integer(Year))
  }
  bind_rows(to_long(at_cell, "cell"), to_long(in_buf, "buffer"))
}


# ------------------------------------------------------------
# 3. Extract N release and ALD
# ------------------------------------------------------------
rapid <- rast(rapid_file)
ald   <- rast(ald_file)

site_data <- extract_sites(rapid, "rapid_kg_m2") %>%
  full_join(extract_sites(ald, "ALD_m"), by = c("site", "type", "Year")) %>%
  mutate(rapid_mg_m2 = rapid_kg_m2 * 1e6)          # kg N m-2 yr-1 -> mg N m-2 yr-1

# warn if a site cell has no data (then use the buffer values)
missing_cells <- site_data %>%
  filter(type == "cell") %>%
  group_by(site) %>%
  summarise(all_na_N = all(is.na(rapid_mg_m2)), all_na_ALD = all(is.na(ALD_m)))
print(missing_cells)

write_csv(site_data, file.path(out_dir, sprintf("site_timeseries_%s.csv", ssp)))


# ------------------------------------------------------------
# 4. Period means
# ------------------------------------------------------------
in_p <- function(y, p) y >= p[1] & y <= p[2]

site_summary <- site_data %>%
  group_by(site, type) %>%
  summarise(
    N_present_mg   = mean(rapid_mg_m2[in_p(Year, present)], na.rm = TRUE),
    N_future_mg    = mean(rapid_mg_m2[in_p(Year, future)],  na.rm = TRUE),
    N_proj_mean_mg = mean(rapid_mg_m2[in_p(Year, proj)],    na.rm = TRUE),
    ALD_present_m  = mean(ALD_m[in_p(Year, present)], na.rm = TRUE),
    ALD_future_m   = mean(ALD_m[in_p(Year, future)],  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    ALD_change_m     = ALD_future_m - ALD_present_m,
    ALD_rate_cm_yr   = 100 * ALD_change_m / (mean(future) - mean(present))
  ) %>%
  left_join(sites %>% select(site, name), by = "site") %>%
  arrange(type, site)

cat("\n=== Site summary, SSP", ssp, "(N release in mg N m-2 yr-1) ===\n")
print(site_summary %>% mutate(across(where(is.numeric), ~ round(.x, 2))), width = Inf)

write_csv(site_summary, file.path(out_dir, sprintf("site_summary_%s.csv", ssp)))


# ------------------------------------------------------------
# 5. Time series plots (buffer means, 20-year running mean)
# ------------------------------------------------------------
plot_data <- site_data %>%
  filter(type == "buffer") %>%
  group_by(site) %>%
  arrange(Year, .by_group = TRUE) %>%
  mutate(N_20   = rollapply(rapid_mg_m2, 20, mean, align = "center", fill = NA, na.rm = TRUE),
         ALD_20 = rollapply(ALD_m,       20, mean, align = "center", fill = NA, na.rm = TRUE)) %>%
  ungroup()

p_N <- ggplot(plot_data, aes(Year)) +
  # Beermann et al. (2017): 8-81 mg N m-2 yr-1 under RCP8.5 (range across sites)
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 8, ymax = 81, alpha = 0.15, fill = "darkgreen") +
  geom_line(aes(y = rapid_mg_m2), colour = "grey60", linewidth = 0.2, na.rm = TRUE) +
  geom_line(aes(y = N_20), colour = "black", linewidth = 0.5, na.rm = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.2) +
  facet_wrap(~site, ncol = 3) +
  labs(x = "Year", y = expression("Pre-thaw inorganic N release [mg N m"^-2 * " yr"^-1 * "]"),
       caption = "Green band: range estimated by Beermann et al. (2017), RCP8.5") +
  theme_minimal(base_size = 8)

p_ALD <- ggplot(plot_data, aes(Year)) +
  geom_line(aes(y = ALD_m), colour = "grey60", linewidth = 0.2, na.rm = TRUE) +
  geom_line(aes(y = ALD_20), colour = "steelblue", linewidth = 0.5, na.rm = TRUE) +
  scale_y_reverse() +
  facet_wrap(~site, ncol = 3) +
  labs(x = "Year", y = "ALD [m]") +
  theme_minimal(base_size = 8)

print(p_N)
print(p_ALD)
ggsave(file.path(out_dir, sprintf("site_N_release_%s.pdf", ssp)), p_N,  width = 17, height = 7, units = "cm")
ggsave(file.path(out_dir, sprintf("site_ALD_%s.pdf", ssp)),       p_ALD, width = 17, height = 7, units = "cm")

cat("\nDone. Outputs in", out_dir, "\n")