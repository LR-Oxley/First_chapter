# ============================================================
# Arctic N fluxes: pan-Arctic and by biome
#
# Compares three N inputs to Arctic soils:
#   - Total inorganic N following permafrost thaw (own model, monthly files)
#   - N fixation      (CMIP6, kg N m-2 s-1)
#   - N deposition    (CMIP6, kg N m-2 s-1)
#
# for present day (2000-2020) and the future anomaly (2080-2099 minus
# 2000-2020), pan-Arctic and split into four biomes
# (Taiga, Tundra, Wetlands, Barren).
#
# Uncertainty
#   thaw N:               three model runs (mean, mean + 1 SD, mean - 1 SD);
#                         every quantity is computed in each run, and
#                         uncertainty = |plus - minus| / 2
#   fixation, deposition: year-to-year SD (present day); for the anomaly
#                         sqrt(SD_future^2 + SD_present^2)
#
# Present-day fixation and deposition = mean of observations and CMIP6
# (setting pool_obs_present); future anomalies are CMIP6 future minus
# CMIP6 present day.
#
# Outputs (in out_dir):
#   Figure 1: pan-Arctic comparison          N_budget_pan_arctic_<unit>.png/.pdf
#   Figure 2: comparison split by biome, two versions (choose in Step 1):
#     "contribution": biome contribution to the pan-Arctic flux
#                                            N_budget_by_biome_<unit>.png/.pdf
#     "intensity":    biome flux per unit area relative to the pan-Arctic
#                     mean (1 = average; >1 = more per m2 than its area
#                     share would suggest)   N_budget_by_biome_intensity.png/.pdf
#   Tables:   annual values, present-day, anomalies, biome shares,
#             relative contributions (.csv)
#
# The script runs top to bottom in 13 steps. Change settings in Step 1 only.
# ============================================================


# ------------------------------------------------------------
# Step 0: packages
# ------------------------------------------------------------
library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(patchwork)   # stacks the two panels of Figure 2; install.packages("patchwork")


# ------------------------------------------------------------
# Step 1: settings
# ------------------------------------------------------------

# folders
out_dir <- file.path("coding_first_chapter", "diagnostic_csv_analysis_biome")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# study region (60-90 N)
ext_sub <- ext(-179.95, 179.95, 60, 90)

# scenarios (no CMIP6 fixation / deposition for SSP2-4.5)
ssps        <- c("126", "370", "585")
ssp_labels  <- c("126" = "SSP1-2.6", "370" = "SSP3-7.0", "585" = "SSP5-8.5")
ssp_colours <- c("SSP1-2.6" = "blue", "SSP3-7.0" = "#D73027", "SSP5-8.5" = "#7B3294")

# periods
baseline_start <- 2000
baseline_end   <- 2020
future_start   <- 2080
future_end     <- 2099
years_needed   <- c(baseline_start:baseline_end, future_start:future_end)

# flux names (order = bottom to top on the y axis)
min_name <- "Total inorg. due to permafrost thaw"
fix_name <- "N Fixation"
dep_name <- "N Deposition"
variable_order <- c(min_name, fix_name, dep_name)

# thaw-N input: "mineralised_N" or "total_inorg_N"
min_flux_var <- "total_inorg_N"

# model runs: mean, mean + 1 SD, mean - 1 SD
# (set plus and minus to NULL to plot without thaw-N uncertainty)
min_dirs <- list(
  mean  = file.path("monthly_mineralised", "mean_2perc_baseline"),
  plus  = file.path("monthly_mineralised", "plus_2perc_baseline"),
  minus = file.path("monthly_mineralised", "minus_2perc_baseline")
)
# file name: variable, SSP, year
min_file_fmt <- "arctic_%s_monthly_w_temp_sm_%s_%d.nc"

# deposition and fixation files (%s = SSP code)
dep_file_fmt <- "dep_1850_2100_ssp%s_totalN.nc"
fix_file_fmt <- "fix_1850_2100_ssp%s.nc"
seconds_per_year <- 31556926

# land mask (model domain) and land cover
thawed_file <- file.path("total_thawed_extended", "arctic_total_thawed_370_60N_mean.nc")

lc_file     <- "LC_remapnn_corr.nc"

# land-cover classes per biome; all other codes (water, NA, ...) are EXCLUDED
biome_classes <- list(
  Taiga    = c(1, 2, 3, 4, 5, 8, 9),
  Tundra   = c(6, 7, 10, 12, 14),
  Wetlands = 11,
  Barren   = c(13, 15, 16)
)
biome_levels  <- names(biome_classes)
biome_colours <- c(Taiga = "#1B7837", Tundra = "#A6DBA0",
                   Wetlands = "#2166AC", Barren = "#BF812D")

# observations for present day (pan-Arctic only), g N m-2 yr-1
obs <- tibble(
  Variable = c(fix_name, dep_name),
  obs_mean = c(0.275, 0.140),    
  obs_sd   = c(0.125, 0.110)
)

# TRUE = present-day N fixation and deposition = mean of observations and
#        CMIP6 (in both figures and all tables); the biome split follows the
#        CMIP6 pattern. Future anomalies stay CMIP6-only (future minus CMIP6
#        present day). FALSE = CMIP6 only.
pool_obs_present <- TRUE

# figures
fig_width_cm  <- 15
fig_height_cm <- 12
fig_dpi       <- 500
base_pt  <- 7       # axis titles, panel titles
axis_pt  <- 6       # tick labels, legend
label_pt <- 5       # numbers next to bars in Figure 1
biome_label_pt   <- 4      # numbers next to bars in Figure 2
label_biome_bars <- TRUE   # TRUE = also write values next to the biome bars in Figure 2
bar_alpha_total  <- 0.7    # transparency of the pan-Arctic / SSP bars in Figure 2 (0-1)
bar_alpha_biome  <- 0.7    # transparency of the biome bars in Figure 2 (0-1)

# Figure 2 versions to draw (one or both):
#   "contribution" = biome contribution to the pan-Arctic flux (bars add up
#                    to the pan-Arctic bar); drawn once per unit
#   "intensity"    = biome flux per m2 / pan-Arctic flux per m2
#                    (= share of flux / share of area); unit-free, drawn once
fig2_modes <- c("contribution", "intensity")

# intensity is left empty (n/a) when the pan-Arctic flux or change is
# smaller than this (g N m-2 yr-1): dividing by ~0 gives meaningless ratios
min_abs_for_ratio <- 0.02

units <- c("g_m2_yr", "Tg_yr")
unit_axis_labels <- list(
  g_m2_yr = expression("N [g N m"^{-2} * " yr"^{-1} * "]"),
  Tg_yr   = expression("N [Tg N yr"^{-1} * "]")
)
intensity_axis_label <- "Flux per unit area relative to the pan-Arctic mean (pan-Arctic = 1)"


# ------------------------------------------------------------
# Step 2: study domain and biome map
#
# biome_r: raster with value 1-4 (= position in biome_levels) in every cell
#          of the study domain, NA everywhere else.
# ------------------------------------------------------------

# 2a. model land mask
thawed    <- crop(rast(thawed_file), ext_sub)
model_land <- !is.na(thawed[[nlyr(thawed)]])          # TRUE on model land cells

# 2b. land cover on the same grid
LC <- crop(rast(lc_file)[[1]], ext_sub)
LC <- resample(LC, thawed[[1]], method = "near")

# 2c. reclassification table: land-cover code -> biome number
#     (one row per code, e.g. code 11 -> 3 = Wetlands)
rcl <- NULL
for (i in seq_along(biome_classes)) {
  rcl <- rbind(rcl, cbind(biome_classes[[i]], i))
}

# 2d. biome map; codes not in the table become NA (= excluded)
biome_r <- classify(round(LC), rcl, others = NA)
biome_r <- mask(biome_r, model_land, maskvalues = 0)
names(biome_r) <- "zone"

# 2e. areas
cell_area   <- cellSize(thawed[[1]], unit = "m")
area_model  <- global(mask(cell_area, model_land, maskvalues = 0), "sum", na.rm = TRUE)[1, 1]
area_domain <- global(mask(cell_area, biome_r), "sum", na.rm = TRUE)[1, 1]

cat("Model land area      [10^6 km2]:", round(area_model / 1e12, 3), "\n")
cat("Study domain (biomes)[10^6 km2]:", round(area_domain / 1e12, 3), "\n")
cat("Excluded (no biome)  [10^6 km2]:", round((area_model - area_domain) / 1e12, 3),
    sprintf("(%.2f %%)", 100 * (area_model - area_domain) / area_model), "\n")

biome_area <- zonal(cell_area, biome_r, "sum", na.rm = TRUE)
biome_area$Biome <- biome_levels[biome_area$zone]
biome_area$area_Mkm2 <- biome_area$area / 1e12
print(biome_area[, c("Biome", "area_Mkm2")])


# ------------------------------------------------------------
# Step 3: helper -- sum a stack of annual flux maps per biome
#
# Input:  r_gm2  = raster, one layer per year, in g N m-2 yr-1
#         years  = the year of each layer
# Output: table with one row per biome and year, plus "Pan-Arctic":
#         Biome, Year, g_N_yr (total mass per year), area_m2 (area with data)
# ------------------------------------------------------------
biome_totals <- function(r_gm2, years) {
  
  names(r_gm2) <- paste0("y", years)
  
  # biome map on the grid of this raster (only resampled if grids differ)
  zones <- biome_r
  if (!compareGeom(r_gm2, biome_r, stopOnError = FALSE)) {
    zones <- resample(biome_r, r_gm2, method = "near")
  }
  area <- cellSize(r_gm2[[1]], unit = "m")
  
  # mass per cell = flux x cell area; area counted only where there is data
  mass_per_cell <- r_gm2 * area
  area_with_data <- (!is.na(r_gm2)) * area
  
  # sum per biome (one column per year) and reshape to long format
  mass_tab <- zonal(mass_per_cell,  zones, "sum", na.rm = TRUE) %>%
    pivot_longer(-zone, names_to = "Year", values_to = "g_N_yr")
  area_tab <- zonal(area_with_data, zones, "sum", na.rm = TRUE) %>%
    pivot_longer(-zone, names_to = "Year", values_to = "area_m2")
  
  by_biome <- left_join(mass_tab, area_tab, by = c("zone", "Year")) %>%
    mutate(Year  = as.integer(sub("y", "", Year)),
           Biome = biome_levels[zone]) %>%
    select(Biome, Year, g_N_yr, area_m2)
  
  # pan-Arctic = sum of the biomes
  pan <- by_biome %>%
    group_by(Year) %>%
    summarise(g_N_yr = sum(g_N_yr), area_m2 = sum(area_m2), .groups = "drop") %>%
    mutate(Biome = "Pan-Arctic")
  
  bind_rows(by_biome, pan)
}


# ------------------------------------------------------------
# Step 4: thaw N (mean, plus and minus run)
#
# Each yearly file has 12 layers (Jan-Dec) in kg N m-2 month-1.
# Sum the 12 months -> kg N m-2 yr-1, x 1000 -> g N m-2 yr-1.
# ------------------------------------------------------------
min_tables <- list()

for (s in ssps) {
  for (run in names(min_dirs)) {
    
    if (is.null(min_dirs[[run]])) next
    
    folder <- file.path(min_dirs[[run]], paste0("yearly_nc_", s))
    annual_maps <- list()
    
    for (yr in years_needed) {
      f <- file.path(folder, sprintf(min_file_fmt, min_flux_var, s, yr))
      if (!file.exists(f)) {
        warning("Missing file: ", f)
        next
      }
      monthly <- rast(f)
      if (xmax(monthly) > 180) monthly <- rotate(monthly)
      monthly <- crop(monthly, ext_sub)
      annual_maps[[as.character(yr)]] <- sum(monthly, na.rm = TRUE) * 1000
    }
    
    if (length(annual_maps) == 0) next
    cat("Thaw N", run, "run, SSP", s, ":", length(annual_maps), "years loaded\n")
    
    tab <- biome_totals(rast(annual_maps), as.integer(names(annual_maps)))
    tab$SSP      <- ssp_labels[[s]]
    tab$Variable <- min_name
    tab$run      <- run
    min_tables[[length(min_tables) + 1]] <- tab
  }
}
min_table <- bind_rows(min_tables)


# ------------------------------------------------------------
# Step 5: N deposition and N fixation (CMIP6, one run)
#
# Files hold one layer per year from 1850 on, in kg N m-2 s-1.
# x seconds per year x 1000 -> g N m-2 yr-1.
# ------------------------------------------------------------
depfix_tables <- list()

for (s in ssps) {
  for (v in c(dep_name, fix_name)) {
    
    f <- if (v == dep_name) sprintf(dep_file_fmt, s) else sprintf(fix_file_fmt, s)
    if (!file.exists(f)) {
      warning("Missing file: ", f)
      next
    }
    
    r <- crop(rast(f), ext_sub)
    file_years <- 1850 + seq_len(nlyr(r)) - 1
    keep <- which(file_years %in% years_needed)
    r_gm2 <- r[[keep]] * seconds_per_year * 1000
    
    tab <- biome_totals(r_gm2, file_years[keep])
    tab$SSP      <- ssp_labels[[s]]
    tab$Variable <- v
    tab$run      <- "mean"
    depfix_tables[[length(depfix_tables) + 1]] <- tab
  }
}
depfix_table <- bind_rows(depfix_tables)


# ------------------------------------------------------------
# Step 6: one annual table with all fluxes and runs, in both units
#
# value   = flux of this biome (or pan-Arctic)
#           g_m2_yr: mass / own area;  Tg_yr: mass / 1e12
# contrib = contribution to the pan-Arctic flux
#           g_m2_yr: mass / pan-Arctic area (biomes add up to pan-Arctic)
#           Tg_yr:   same as value
# run     = mean / plus / minus (fixation and deposition: mean only)
# ------------------------------------------------------------
all_runs <- bind_rows(min_table, depfix_table)

# pan-Arctic area per run / SSP / flux / year (needed for the contributions)
pan_area <- all_runs %>%
  filter(Biome == "Pan-Arctic") %>%
  select(run, SSP, Variable, Year, pan_area_m2 = area_m2)

all_runs <- all_runs %>%
  left_join(pan_area, by = c("run", "SSP", "Variable", "Year"))

# both units below each other (column "unit")
annual <- bind_rows(
  all_runs %>% mutate(unit = "g_m2_yr",
                      value   = g_N_yr / area_m2,
                      contrib = g_N_yr / pan_area_m2),
  all_runs %>% mutate(unit = "Tg_yr",
                      value   = g_N_yr / 1e12,
                      contrib = g_N_yr / 1e12)
) %>%
  select(unit, run, SSP, Variable, Biome, Year, value, contrib)

write_csv(annual, file.path(out_dir, "annual_N_fluxes_by_biome.csv"))

# helper: |plus - minus| / 2 of a column, per group (NA if runs are missing)
spread_pm <- function(d, keys, col) {
  w <- d %>%
    filter(run %in% c("plus", "minus")) %>%
    select(all_of(c(keys, "run", col))) %>%
    pivot_wider(names_from = run, values_from = all_of(col))
  for (k in c("plus", "minus")) if (!k %in% names(w)) w[[k]] <- NA_real_
  w %>%
    mutate(sd_runs = abs(plus - minus) / 2) %>%
    select(all_of(c(keys, "sd_runs")))
}


# ------------------------------------------------------------
# Step 7: present-day means (2000-2020, all SSPs pooled), per run
#
# SD: thaw N                -> |plus - minus| / 2
#     deposition / fixation -> year-to-year SD
# ------------------------------------------------------------
present_runs <- annual %>%
  filter(Year >= baseline_start, Year <= baseline_end) %>%
  group_by(unit, run, Variable, Biome) %>%
  summarise(mean     = mean(value, na.rm = TRUE),
            contrib  = mean(contrib, na.rm = TRUE),
            sd_years = sd(value, na.rm = TRUE),
            .groups = "drop")

present <- present_runs %>%
  filter(run == "mean") %>%
  select(-run) %>%
  left_join(spread_pm(present_runs, c("unit", "Variable", "Biome"), "mean"),
            by = c("unit", "Variable", "Biome")) %>%
  mutate(sd = ifelse(is.finite(sd_runs), sd_runs, sd_years))

present_model <- present      # CMIP6 / model values only, kept for reference
write_csv(present_model, file.path(out_dir, "present_day_2000_2020_model_only.csv"))


# ------------------------------------------------------------
# Step 7b: present-day fixation and deposition = mean of observations
#          and CMIP6 (if pool_obs_present = TRUE)
#
# Observations exist only for the pan-Arctic, so every biome value is
# scaled by the same factor:
#   factor = ((CMIP6 pan-Arctic + observation) / 2) / CMIP6 pan-Arctic
# -> the pan-Arctic value becomes the mean of both, and the biome split
#    (shares, intensities) follows CMIP6.
# SD of the pan-Arctic value: (SD_CMIP6 + SD_obs) / 2
# Thaw N is not changed.
# ------------------------------------------------------------
if (pool_obs_present) {
  
  obs_units <- bind_rows(
    obs %>% mutate(unit = "g_m2_yr"),
    obs %>% mutate(unit = "Tg_yr",
                   obs_mean = obs_mean * area_domain / 1e12,
                   obs_sd   = obs_sd   * area_domain / 1e12)
  )
  
  pool_factors <- present_model %>%
    filter(Biome == "Pan-Arctic") %>%
    inner_join(obs_units, by = c("unit", "Variable")) %>%
    transmute(unit, Variable,
              model_pan = mean,
              obs_mean,
              pooled_pan = (mean + obs_mean) / 2,
              factor     = pooled_pan / model_pan,
              sd_pooled  = (sd + obs_sd) / 2)
  
  cat("\nPresent-day fixation and deposition, CMIP6 vs observations vs mean of both:\n")
  print(pool_factors, width = Inf)
  
  present <- present_model %>%
    left_join(pool_factors %>% select(unit, Variable, factor, sd_pooled),
              by = c("unit", "Variable")) %>%
    mutate(mean    = ifelse(is.na(factor), mean,    mean    * factor),
           contrib = ifelse(is.na(factor), contrib, contrib * factor),
           sd      = ifelse(!is.na(sd_pooled) & Biome == "Pan-Arctic", sd_pooled, sd)) %>%
    select(-factor, -sd_pooled)
}

write_csv(present, file.path(out_dir, "present_day_2000_2020.csv"))


# ------------------------------------------------------------
# Step 8: future anomaly (2080-2099 mean minus present-day mean), per run
#
# SD: thaw N                -> |anomaly_plus - anomaly_minus| / 2
#     deposition / fixation -> sqrt(SD_future^2 + SD_present^2)
# ------------------------------------------------------------
future_runs <- annual %>%
  filter(Year >= future_start, Year <= future_end) %>%
  group_by(unit, run, SSP, Variable, Biome) %>%
  summarise(mean_future    = mean(value, na.rm = TRUE),
            contrib_future = mean(contrib, na.rm = TRUE),
            sd_future      = sd(value, na.rm = TRUE),
            .groups = "drop")

anomaly_runs <- future_runs %>%
  left_join(present_runs %>% select(unit, run, Variable, Biome,
                                    mean_present = mean,
                                    contrib_present = contrib,
                                    sd_present = sd_years),
            by = c("unit", "run", "Variable", "Biome")) %>%
  mutate(anomaly         = mean_future - mean_present,
         contrib_anomaly = contrib_future - contrib_present)

anomaly <- anomaly_runs %>%
  filter(run == "mean") %>%
  select(-run) %>%
  left_join(spread_pm(anomaly_runs, c("unit", "SSP", "Variable", "Biome"), "anomaly"),
            by = c("unit", "SSP", "Variable", "Biome")) %>%
  mutate(anomaly_sd = ifelse(is.finite(sd_runs), sd_runs,
                             sqrt(sd_future^2 + sd_present^2)))

future <- future_runs %>% filter(run == "mean") %>% select(-run)

write_csv(anomaly, file.path(out_dir, "future_anomaly_2080_2099.csv"))


# ------------------------------------------------------------
# Step 9: biome shares (% of the pan-Arctic mass flux), for the text
# ------------------------------------------------------------
shares_present <- present %>%
  filter(unit == "Tg_yr", Biome != "Pan-Arctic") %>%
  group_by(Variable) %>%
  mutate(period = "Present day", share_pct = 100 * mean / sum(mean)) %>%
  ungroup() %>%
  select(period, Variable, Biome, Tg_N_yr = mean, share_pct)

shares_future <- future %>%
  filter(unit == "Tg_yr", Biome != "Pan-Arctic") %>%
  group_by(SSP, Variable) %>%
  mutate(period = paste(SSP, "2080-2099"), share_pct = 100 * mean_future / sum(mean_future)) %>%
  ungroup() %>%
  select(period, Variable, Biome, Tg_N_yr = mean_future, share_pct)

write_csv(bind_rows(shares_present, shares_future),
          file.path(out_dir, "biome_shares.csv"))


# ------------------------------------------------------------
# Step 10: relative contributions -- flux share vs area share
#
# For every biome (present day, and future anomaly per SSP):
#   area_pct     = share of the study-domain area
#   flux_pct     = share of the pan-Arctic flux (or change)
#   density      = flux per m2 of the biome itself (g N m-2 yr-1)
#   intensity    = density / pan-Arctic density
#                  (~ flux_pct / area_pct; 1 = average, >1 = above average)
# Ratios are set to NA when the pan-Arctic value is below min_abs_for_ratio.
# ------------------------------------------------------------
area_share <- biome_area %>%
  transmute(Biome, area_Mkm2, area_pct = 100 * area / sum(area))

relative_to_pan <- function(d, panel_name) {
  d %>%
    group_by(group, Variable) %>%
    mutate(pan_density = density[Biome == "Pan-Arctic"],
           pan_contrib = contrib[Biome == "Pan-Arctic"]) %>%
    ungroup() %>%
    filter(Biome != "Pan-Arctic") %>%
    mutate(panel     = panel_name,
           ok        = abs(pan_density) >= min_abs_for_ratio,
           flux_pct  = ifelse(ok, 100 * contrib / pan_contrib, NA_real_),
           intensity = ifelse(ok, density / pan_density, NA_real_)) %>%
    select(-ok, -pan_contrib)
}

rel <- bind_rows(
  present %>%
    filter(unit == "g_m2_yr") %>%
    transmute(group = "Present day", Variable, Biome, density = mean, contrib) %>%
    relative_to_pan("Present day"),
  anomaly %>%
    filter(unit == "g_m2_yr") %>%
    transmute(group = SSP, Variable, Biome, density = anomaly, contrib = contrib_anomaly) %>%
    relative_to_pan("Future anomaly")
) %>%
  left_join(area_share, by = "Biome") %>%
  select(panel, group, Variable, Biome, area_Mkm2, area_pct, flux_pct,
         density, pan_density, intensity)

cat("\nRelative contributions (flux share vs area share):\n")
print(rel %>% mutate(across(where(is.numeric), ~ round(.x, 3))), n = Inf)
write_csv(rel, file.path(out_dir, "relative_contributions.csv"))


# ------------------------------------------------------------
# Step 11: Figure 1 -- pan-Arctic comparison
#
# Present-day bars = `present` (fixation and deposition pooled with the
# observations in Step 7b if pool_obs_present = TRUE). Thaw N = model only.
# ------------------------------------------------------------
for (u in units) {
  
  # present-day bars
  fig1_present <- present %>%
    filter(unit == u, Biome == "Pan-Arctic") %>%
    mutate(value = mean,
           Type  = "Present day",
           fill  = "Present day")
  
  # future bars: one per SSP
  fig1_future <- anomaly %>%
    filter(unit == u, Biome == "Pan-Arctic") %>%
    mutate(value = anomaly,
           sd    = anomaly_sd,
           Type  = "Future anomaly",
           fill  = SSP)
  
  fig1_data <- bind_rows(fig1_present, fig1_future) %>%
    select(Type, Variable, fill, value, sd) %>%
    mutate(Type     = factor(Type, levels = c("Present day", "Future anomaly")),
           Variable = factor(Variable, levels = variable_order),
           fill     = factor(fill, levels = c("Present day", unname(ssp_labels))),
           label_x  = value + ifelse(is.na(sd), 0, sd) * sign(value),
           hjust    = ifelse(value >= 0, -0.3, 1.3))
  
  write_csv(fig1_data, file.path(out_dir, paste0("figure1_data_", u, ".csv")))
  
  dodge <- position_dodge(width = 0.8)
  
  fig1 <- ggplot(fig1_data, aes(x = value, y = Variable, fill = fill, group = fill)) +
    geom_col(position = dodge, width = 0.75) +
    geom_errorbar(aes(xmin = value - sd, xmax = value + sd),
                  position = dodge, width = 0.2, linewidth = 0.3) +
    geom_text(aes(x = label_x, label = sprintf("%.2f", value), hjust = hjust),
              position = dodge, size = label_pt / .pt) +
    geom_vline(xintercept = 0, linewidth = 0.4) +
    facet_wrap(~Type, nrow = 2, scales = "free_y") +
    scale_x_continuous(expand = expansion(mult = c(0.08, 0.15))) +
    scale_fill_manual(values = c("Present day" = "grey70", ssp_colours),
                      breaks = unname(ssp_labels), name = "") +
    labs(x = unit_axis_labels[[u]], y = "") +
    theme_minimal(base_size = base_pt) +
    theme(panel.grid.major.y = element_blank(),
          strip.text   = element_text(face = "bold", size = base_pt),
          axis.text    = element_text(size = axis_pt),
          legend.text  = element_text(size = axis_pt),
          legend.position = "bottom")
  
  ggsave(file.path(out_dir, paste0("N_budget_pan_arctic_", u, ".png")), fig1,
         width = fig_width_cm, height = fig_height_cm, units = "cm", dpi = fig_dpi)
  ggsave(file.path(out_dir, paste0("N_budget_pan_arctic_", u, ".pdf")), fig1,
         width = fig_width_cm, height = fig_height_cm, units = "cm", dpi = fig_dpi)
  print(fig1)
}


# ------------------------------------------------------------
# Step 12: Figure 2 -- the same comparison split by biome
#
# Per flux: pan-Arctic bar (grey / SSP colour), and below it one thin bar
# per biome. Present-day fixation and deposition as in Step 7b (mean of
# observations and CMIP6 if pool_obs_present = TRUE); everything else is
# model output. Two versions (see fig2_modes in Step 1):
#
#   "contribution": biome bars = contribution to the pan-Arctic flux, so
#                   they add up to the pan-Arctic bar; error bars on totals.
#                   Drawn once per unit.
#   "intensity":    biome bars = flux per m2 relative to the pan-Arctic
#                   mean (from Step 10); pan-Arctic bar = 1 (dashed line).
#                   Unit-free, drawn once. n/a = pan-Arctic value ~0.
# ------------------------------------------------------------
gap_rows <- 1.5                                   # empty space between flux blocks
group_levels <- c("Present day", rev(unname(ssp_labels)))   # SSP5-8.5 on top
item_levels  <- c("Pan-Arctic", biome_levels)

# list of figures to draw: (mode, unit)
fig2_jobs <- list()
if ("contribution" %in% fig2_modes) {
  for (u in units) fig2_jobs[[length(fig2_jobs) + 1]] <- list(mode = "contribution", unit = u)
}
if ("intensity" %in% fig2_modes) {
  fig2_jobs[[length(fig2_jobs) + 1]] <- list(mode = "intensity", unit = "g_m2_yr")
}

for (job in fig2_jobs) {
  
  mode <- job$mode
  u    <- job$unit
  tag  <- if (mode == "contribution") u else "intensity"   # used in file names
  
  # 12a. one row per bar
  if (mode == "contribution") {
    
    bars_present <- present %>%
      filter(unit == u) %>%
      transmute(panel = "Present day", Variable, group = "Present day",
                item = Biome, value = contrib, sd = sd)
    
    bars_future <- anomaly %>%
      filter(unit == u) %>%
      transmute(panel = "Future anomaly", Variable, group = SSP,
                item = Biome, value = contrib_anomaly, sd = anomaly_sd)
    
    bars <- bind_rows(bars_present, bars_future)
    
  } else {
    
    bars_biomes <- rel %>%
      transmute(panel, Variable, group, item = Biome, value = intensity, sd = NA_real_)
    
    # reference bar = 1 (left empty if all biome ratios are n/a)
    bars_totals <- bars_biomes %>%
      group_by(panel, Variable, group) %>%
      summarise(value = ifelse(all(is.na(value)), NA_real_, 1), .groups = "drop") %>%
      mutate(item = "Pan-Arctic", sd = NA_real_)
    
    bars <- bind_rows(bars_totals, bars_biomes)
  }
  
  bars <- bars %>%
    mutate(
      kind     = ifelse(item == "Pan-Arctic", "total", "biome"),
      sd       = ifelse(kind == "total", sd, NA),            # error bars on totals only
      fill_key = ifelse(kind == "total", group, item),
      name     = ifelse(kind == "total",
                        ifelse(group == "Present day", "Pan-Arctic", group),
                        item),
      label    = paste0(name, "  ", ifelse(is.na(value), "n/a", sprintf("%.2f", value))),
      label_x  = pmax(value, value + ifelse(is.na(sd), 0, sd), 0, na.rm = TRUE)
    )
  
  # 12b. check: biome bars add up to the pan-Arctic bar (contribution only)
  if (mode == "contribution") {
    check <- bars %>%
      group_by(panel, Variable, group) %>%
      summarise(pan_arctic = sum(value[kind == "total"]),
                sum_biomes = sum(value[kind == "biome"]), .groups = "drop") %>%
      mutate(diff = sum_biomes - pan_arctic)
    cat("\nFigure 2 check (", u, "): biome sum - pan-Arctic\n")
    print(check, n = Inf)
  }
  
  # 12c. y position of every bar
  #      blocks of rows per flux (bottom to top = variable_order),
  #      inside a block: group (SSP) by group, pan-Arctic bar first, then biomes
  bars <- bars %>%
    mutate(block       = match(Variable, variable_order),
           group_order = match(group, group_levels),
           item_order  = match(item, item_levels)) %>%
    arrange(panel, block, group_order, item_order) %>%
    group_by(panel, Variable) %>%
    mutate(row = row_number(), n_rows = n()) %>%
    ungroup() %>%
    mutate(y = (block - 1) * (n_rows + gap_rows) + (n_rows - row + 1))
  
  write_csv(bars, file.path(out_dir, paste0("figure2_data_", tag, ".csv")))
  
  # 12d. where the flux names go on the y axis (middle of each block)
  y_ticks <- bars %>%
    group_by(panel, Variable) %>%
    summarise(centre = mean(range(y)), .groups = "drop")
  ticks_present <- y_ticks %>% filter(panel == "Present day")
  ticks_future  <- y_ticks %>% filter(panel == "Future anomaly")
  
  # 12e. same x range in both panels, so the zero lines line up
  x_lo <- min(0, bars$value - ifelse(is.na(bars$sd), 0, bars$sd), na.rm = TRUE)
  x_hi <- max(0, bars$value + ifelse(is.na(bars$sd), 0, bars$sd), na.rm = TRUE)
  
  x_title  <- if (mode == "contribution") unit_axis_labels[[u]] else intensity_axis_label
  ref_line <- if (mode == "intensity") {
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3, colour = "grey40")
  } else NULL
  
  # 12f. plot layers used by both panels
  shared_layers <- list(
    geom_col(data = function(d) filter(d, kind == "total"),
             aes(x = value, y = y, fill = fill_key), orientation = "y", width = 0.85,
             alpha = bar_alpha_total, na.rm = TRUE),
    geom_col(data = function(d) filter(d, kind == "biome"),
             aes(x = value, y = y, fill = fill_key), orientation = "y", width = 0.6,
             alpha = bar_alpha_biome, na.rm = TRUE),
    geom_errorbar(data = function(d) filter(d, !is.na(sd)),
                  aes(xmin = value - sd, xmax = value + sd, y = y),
                  orientation = "y", width = 0.5, linewidth = 0.3),
    geom_text(data = function(d) filter(d, kind == "total" | label_biome_bars),
              aes(x = label_x, y = y, label = label,
                  fontface = ifelse(kind == "total", "bold", "plain")),
              hjust = -0.05, size = biome_label_pt / .pt),
    geom_vline(xintercept = 0, linewidth = 0.3),
    ref_line,
    scale_x_continuous(limits = c(x_lo, x_hi * 1.35), expand = expansion(mult = c(0.02, 0))),
    scale_fill_manual(values = c("Present day" = "grey60", ssp_colours, biome_colours),
                      limits = c("Present day", unname(ssp_labels), biome_levels),
                      breaks = c(unname(ssp_labels), biome_levels),   # "Present day" not in legend
                      guide  = guide_legend(nrow = 2, byrow = TRUE),
                      name = ""),
    labs(x = x_title, y = ""),
    theme_minimal(base_size = base_pt),
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          plot.subtitle   = element_text(face = "bold", size = base_pt, hjust = 0.5),
          axis.text       = element_text(size = axis_pt),
          legend.text     = element_text(size = axis_pt),
          legend.key.size = unit(0.3, "cm"),
          legend.position = "bottom")
  )
  
  # 12g. the two panels
  p_present <- ggplot(filter(bars, panel == "Present day")) +
    shared_layers +
    scale_y_continuous(breaks = ticks_present$centre, labels = ticks_present$Variable,
                       expand = expansion(add = 0.8)) +
    labs(subtitle = "Present day (2000-2020)") +
    theme(axis.title.x = element_blank())
  
  p_future <- ggplot(filter(bars, panel == "Future anomaly")) +
    shared_layers +
    scale_y_continuous(breaks = ticks_future$centre, labels = ticks_future$Variable,
                       expand = expansion(add = 0.8)) +
    labs(subtitle = "Future anomaly (2080-2099 vs 2000-2020)")
  
  # 12h. stack them; panel heights follow the number of bars
  fig2 <- p_present / p_future +
    plot_layout(heights = c(max(p_present$data$y), max(p_future$data$y)),
                guides = "collect") &
    theme(legend.position = "none")
  
  ggsave(file.path(out_dir, paste0("N_budget_by_biome_", tag, ".png")), fig2,
         width = fig_width_cm, height = fig_height_cm, units = "cm", dpi = fig_dpi)
  ggsave(file.path(out_dir, paste0("N_budget_by_biome_", tag, ".pdf")), fig2,
         width = fig_width_cm, height = fig_height_cm, units = "cm", dpi = fig_dpi)
  print(fig2)
}


# ------------------------------------------------------------
# Step 13: done
# ------------------------------------------------------------
cat("\nDone. Outputs in", out_dir, "\n")
