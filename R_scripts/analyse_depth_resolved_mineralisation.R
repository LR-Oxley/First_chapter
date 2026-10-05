# g N / m3 * yr analysis, with top-soil (0.1 m) comparison to Ramm et al. (2022)

# ============================================================
# Depth-resolved mineralised N analysis, ONE SSP per run
#
# This script is designed to run once per SSP (e.g. as one task
# of a SLURM array job). It saves per-SSP output files only.
# Run combine_ssp_results.R afterward to merge all SSPs into
# comparison tables and plots.
#
# Produces (per SSP):
# 1. Annual + cumulative mineralised N (Pg N), whole-column
#    profile-mean areal (g N m-2 yr-1) and volumetric (g N m-3 yr-1)
#    rate, PLUS top-soil-only (0.1 m) areal and volumetric annual
#    rate - full available record
# 2. Period comparison (present-day / mid-century / end-century)
#    of mean rates, whole-column and top-soil
# 3. Mean monthly seasonal cycle of areal and volumetric flux,
#    whole-column and top-soil, by period, for years within the
#    three periods
# 4. Top-soil growing-season-only total (matching Ramm et al. 2022's
#    100-day growing-season basis), by period
# ============================================================

library(terra)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)

# -----------------------------
# 0. Setup
# -----------------------------

args <- commandArgs(trailingOnly = TRUE)
ssp <- args[1]

if (is.na(ssp) || length(ssp) == 0 || nchar(ssp) == 0) {
  stop("No SSP provided. Usage: Rscript arctic_mineralisation_by_depth.R <ssp>")
}

base_dir <- "monthly_mineralised/mean_2perc_baseline_new"
out_dir <- "monthly_mineralised/mineralisation_by_depth_analysis"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_dir <- "logs"
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(
  log_dir,
  paste0(
    "arctic_mineralisation_by_depth_",
    ssp,
    "_",
    format(Sys.time(), "%Y%m%d_%H%M%S"),
    ".log"
  )
)

log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, append = TRUE, type = "message")

start_time <- Sys.time()

cat("====================================================\n")
cat("STARTED SCRIPT\n")
cat("Time:", as.character(start_time), "\n")
cat("SSP:", ssp, "\n")
cat("====================================================\n")

on.exit({
  end_time <- Sys.time()
  cat("====================================================\n")
  cat("SCRIPT FINISHED\n")
  cat("End:", as.character(end_time), "\n")
  cat("Runtime:", as.character(end_time - start_time), "\n")
  cat("Log file:", log_file, "\n")
  cat("====================================================\n")
  sink(type = "message")
  sink()
  close(log_con)
}, add = TRUE)

ssp_labels <- c(
  "126" = "SSP1-2.6",
  "245" = "SSP2-4.5",
  "370" = "SSP3-7.0",
  "585" = "SSP5-8.5"
)

if (!ssp %in% names(ssp_labels)) {
  stop("Unrecognised SSP: '", ssp, "'. Expected one of: ", paste(names(ssp_labels), collapse = ", "))
}

profile_depth_max <- 5

periods <- list(
  "Present-day (2000-2020)"    = 2000:2020,
  "Mid-century (2040-2060)"    = 2040:2060,
  "End-of-century (2080-2099)" = 2080:2099
)

# Top-soil target depth for comparison with Ramm et al. (2022)
target_depth_topsoil <- 0.1   # metres

# Growing-season months used for the direct Ramm et al. (2022) comparison
# (their 552 g N m-3 figure is a 100-day growing-season TOTAL, not annualised).
# Adjust this if you want a different definition of the growing season.
growing_season_months <- 6:8   # June, July, August

# ============================================================
# 1. Recover the physical depth structure
#    (assumes depth levels are identical across SSPs;
#    uses this run's own SSP's temperature file)
# ============================================================

get_depth_bounds <- function(depth_mid, max_depth) {
  depth_bounds <- numeric(length(depth_mid) + 1)
  depth_bounds[1] <- 0
  for (i in 2:length(depth_mid)) {
    depth_bounds[i] <- (depth_mid[i - 1] + depth_mid[i]) / 2
  }
  depth_bounds[length(depth_bounds)] <- max_depth
  depth_bounds[depth_bounds > max_depth] <- max_depth
  depth_bounds
}

get_depth_structure <- function(ssp_for_depth) {
  temp_file <- paste0("60N/mean_verttemp_ssp", ssp_for_depth, "_60deg.nc")
  temp_full <- rast(temp_file)
  n_depths_full <- length(unique(as.numeric(depth(temp_full))))
  all_depth_vals <- as.numeric(depth(temp_full))[1:n_depths_full]
  depth_keep <- which(all_depth_vals <= profile_depth_max)
  depth_vals <- all_depth_vals[depth_keep]
  n_depths <- length(depth_vals)
  depth_bounds <- get_depth_bounds(depth_vals, max_depth = profile_depth_max)
  depth_thickness <- diff(depth_bounds)
  depth_weights <- depth_thickness / sum(depth_thickness)
  rm(temp_full); gc()
  list(
    n_depths = n_depths,
    depth_thickness = depth_thickness,
    depth_weights = depth_weights,
    depth_bounds = depth_bounds
  )
}

depth_struct <- get_depth_structure(ssp)
n_depths <- depth_struct$n_depths
depth_thickness <- depth_struct$depth_thickness
depth_weights <- depth_struct$depth_weights
depth_bounds <- depth_struct$depth_bounds

cat("Recovered", n_depths, "depth layers. Thickness (m):\n")
print(round(depth_thickness, 3))

# ============================================================
# 2. Top-soil overlap fractions
#    Each depth layer's fractional overlap with the top
#    target_depth_topsoil band, so a layer straddling the
#    boundary is pro-rated rather than crudely in/excluded.
#    Assumes N is uniformly distributed within each layer.
# ============================================================

get_topsoil_overlap_fractions <- function(depth_bounds, target_depth) {
  
  n_depths_local <- length(depth_bounds) - 1
  overlap_fraction <- numeric(n_depths_local)
  
  for (d in seq_len(n_depths_local)) {
    z1 <- depth_bounds[d]
    z2 <- depth_bounds[d + 1]
    layer_thickness <- z2 - z1
    
    overlap_z1 <- max(z1, 0)
    overlap_z2 <- min(z2, target_depth)
    overlap_thickness <- max(overlap_z2 - overlap_z1, 0)
    
    overlap_fraction[d] <- if (layer_thickness > 0) overlap_thickness / layer_thickness else 0
  }
  
  overlap_fraction
}

topsoil_overlap_fraction <- get_topsoil_overlap_fractions(depth_bounds, target_depth_topsoil)

# Actual depth represented by the overlap-weighted top-soil sum, in case
# layer boundaries don't align exactly with target_depth_topsoil
actual_topsoil_depth <- sum(topsoil_overlap_fraction * depth_thickness)

cat("Top-soil overlap fractions per depth layer:\n")
print(round(topsoil_overlap_fraction, 3))
cat("Actual top-soil depth represented:", round(actual_topsoil_depth, 4), "m",
    "(target was", target_depth_topsoil, "m)\n")

# ============================================================
# 3. File-listing helper
# ============================================================

list_depth_resolved_files <- function(ssp, years_filter = NULL) {
  
  yearly_dir <- file.path(base_dir, paste0("yearly_nc_", ssp))
  
  file_pattern <- paste0(
    "^arctic_mineralised_N_by_depth_monthly_w_temp_sm_", ssp, "_[0-9]{4}\\.nc$"
  )
  
  nc_files <- list.files(yearly_dir, pattern = file_pattern, full.names = TRUE)
  
  if (length(nc_files) == 0) {
    warning("No depth-resolved files found for SSP ", ssp, " in: ", yearly_dir)
    return(NULL)
  }
  
  file_years <- as.integer(str_extract(basename(nc_files), "(?<=_)[0-9]{4}(?=\\.nc$)"))
  file_df <- data.frame(file = nc_files, year = file_years) %>% arrange(year)
  
  if (!is.null(years_filter)) {
    file_df <- file_df %>% filter(year %in% years_filter)
  }
  
  file_df
}

# ============================================================
# 4a. Per-year processing: ANNUAL summary
#     (Pg N, whole-column areal/volumetric, top-soil areal/volumetric)
# ============================================================

process_one_year_annual <- function(nc_path, yr, n_depths, depth_thickness, depth_weights,
                                    topsoil_overlap_fraction, actual_topsoil_depth) {
  
  r <- rast(nc_path)
  
  expected_nlyr <- n_depths * 12
  if (nlyr(r) != expected_nlyr) {
    stop(
      "Layer count mismatch in ", nc_path, ": found ", nlyr(r),
      ", expected ", expected_nlyr, " (n_depths=", n_depths, " x 12 months)."
    )
  }
  
  layer_position <- seq_len(nlyr(r))
  depth_of_layer <- ((layer_position - 1) %% n_depths) + 1
  
  area_rast <- cellSize(r[[1]], unit = "m")
  
  # Annual areal density per depth layer (kg N m-2 yr-1): sum 12 months per depth
  annual_by_depth <- terra::tapp(r, index = depth_of_layer, fun = "sum", na.rm = TRUE)
  names(annual_by_depth) <- paste0("depth_", seq_len(n_depths))
  
  # ---- Whole-column quantities ----
  
  annual_total_areal <- sum(annual_by_depth)   # kg N m-2 yr-1, all depths summed
  
  # Spatial sum -> Pg N (pan-Arctic total)
  annual_total_pg <- global(
    annual_total_areal * area_rast, "sum", na.rm = TRUE
  )[1, 1] / 1e12
  
  # Pan-Arctic spatial MEAN areal rate, whole column, g N m-2 yr-1
  profile_mean_areal_gN_m2_yr <- global(
    annual_total_areal * 1000, "mean", weights = area_rast, na.rm = TRUE
  )[1, 1]
  
  # Volumetric rate per depth layer, then thickness-weighted collapse across depth
  # (this collapses to whole-column-areal / whole-column-depth - see derivation notes;
  # kept here for continuity with earlier analysis versions)
  volumetric_by_depth <- (annual_by_depth * 1000) / depth_thickness
  weighted_layers <- volumetric_by_depth * depth_weights
  profile_volumetric <- sum(weighted_layers)
  
  profile_volumetric_mean <- global(
    profile_volumetric, "mean", weights = area_rast, na.rm = TRUE
  )[1, 1]
  
  # ---- Top-soil-only quantities (0.1 m, pro-rated across straddling layers) ----
  
  topsoil_areal_by_depth <- annual_by_depth * topsoil_overlap_fraction
  topsoil_areal_total <- sum(topsoil_areal_by_depth)   # kg N m-2 yr-1, top 0.1 m only
  
  topsoil_areal_mean_gN_m2_yr <- global(
    topsoil_areal_total * 1000, "mean", weights = area_rast, na.rm = TRUE
  )[1, 1]
  
  topsoil_volumetric <- (topsoil_areal_total * 1000) / actual_topsoil_depth   # g N m-3 yr-1
  
  topsoil_volumetric_mean_gN_m3_yr <- global(
    topsoil_volumetric, "mean", weights = area_rast, na.rm = TRUE
  )[1, 1]
  
  rm(r, annual_by_depth, annual_total_areal, volumetric_by_depth,
     weighted_layers, profile_volumetric, topsoil_areal_by_depth,
     topsoil_areal_total, topsoil_volumetric)
  gc()
  
  data.frame(
    Year = yr,
    annual_total_pg = annual_total_pg,
    profile_mean_areal_gN_m2_yr = profile_mean_areal_gN_m2_yr,
    profile_mean_volumetric_gN_m3_yr = profile_volumetric_mean,
    topsoil_areal_gN_m2_yr = topsoil_areal_mean_gN_m2_yr,
    topsoil_volumetric_gN_m3_yr = topsoil_volumetric_mean_gN_m3_yr
  )
}

# ============================================================
# 4b. Per-year processing: MONTHLY seasonal profile-mean flux
#     (whole-column and top-soil, areal and volumetric)
# ============================================================

process_one_year_monthly <- function(nc_path, yr, n_depths, depth_thickness, depth_weights,
                                     topsoil_overlap_fraction, actual_topsoil_depth) {
  
  r <- rast(nc_path)
  
  expected_nlyr <- n_depths * 12
  if (nlyr(r) != expected_nlyr) {
    stop(
      "Layer count mismatch in ", nc_path, ": found ", nlyr(r),
      ", expected ", expected_nlyr, " (n_depths=", n_depths, " x 12 months)."
    )
  }
  
  layer_position <- seq_len(nlyr(r))
  month_of_layer <- ceiling(layer_position / n_depths)
  depth_of_layer <- ((layer_position - 1) %% n_depths) + 1
  
  area_rast <- cellSize(r[[1]], unit = "m")
  
  thickness_vec <- depth_thickness[depth_of_layer]
  weight_vec    <- depth_weights[depth_of_layer]
  topsoil_frac_vec <- topsoil_overlap_fraction[depth_of_layer]
  
  # ---- Whole-column volumetric pass ----
  volumetric_r <- (r * 1000) / thickness_vec   # kg N m-2 -> g N m-3, per layer
  weighted_r   <- volumetric_r * weight_vec    # apply depth weight, per layer
  
  monthly_profile_volumetric_r <- terra::tapp(weighted_r, index = month_of_layer, fun = "sum", na.rm = TRUE)
  names(monthly_profile_volumetric_r) <- paste0("month_", sprintf("%02d", seq_len(12)))
  
  monthly_means_volumetric <- terra::global(
    monthly_profile_volumetric_r, "mean", weights = area_rast, na.rm = TRUE
  )[, 1]
  
  # ---- Whole-column areal pass ----
  areal_r <- r * 1000   # kg N m-2 -> g N m-2, per layer
  
  monthly_profile_areal_r <- terra::tapp(areal_r, index = month_of_layer, fun = "sum", na.rm = TRUE)
  names(monthly_profile_areal_r) <- paste0("month_", sprintf("%02d", seq_len(12)))
  
  monthly_means_areal <- terra::global(
    monthly_profile_areal_r, "mean", weights = area_rast, na.rm = TRUE
  )[, 1]
  
  # ---- Top-soil areal pass (pro-rated by overlap fraction) ----
  topsoil_areal_r <- areal_r * topsoil_frac_vec   # g N m-2, top-soil share, per layer
  
  monthly_topsoil_areal_r <- terra::tapp(topsoil_areal_r, index = month_of_layer, fun = "sum", na.rm = TRUE)
  names(monthly_topsoil_areal_r) <- paste0("month_", sprintf("%02d", seq_len(12)))
  
  monthly_means_topsoil_areal <- terra::global(
    monthly_topsoil_areal_r, "mean", weights = area_rast, na.rm = TRUE
  )[, 1]
  
  # ---- Top-soil volumetric pass ----
  monthly_topsoil_volumetric_r <- monthly_topsoil_areal_r / actual_topsoil_depth   # g N m-3
  
  monthly_means_topsoil_volumetric <- terra::global(
    monthly_topsoil_volumetric_r, "mean", weights = area_rast, na.rm = TRUE
  )[, 1]
  
  rm(r, volumetric_r, weighted_r, monthly_profile_volumetric_r,
     areal_r, monthly_profile_areal_r, topsoil_areal_r,
     monthly_topsoil_areal_r, monthly_topsoil_volumetric_r)
  gc()
  
  data.frame(
    Year = yr,
    Month = seq_len(12),
    profile_mean_areal_gN_m2_month = monthly_means_areal,
    profile_mean_volumetric_gN_m3_month = monthly_means_volumetric,
    topsoil_areal_gN_m2_month = monthly_means_topsoil_areal,
    topsoil_volumetric_gN_m3_month = monthly_means_topsoil_volumetric
  )
}

# ============================================================
# 5. Period-assignment helper
# ============================================================

assign_period <- function(yr) {
  for (p in names(periods)) if (yr %in% periods[[p]]) return(p)
  NA_character_
}

# ============================================================
# 6. PART 1: Annual + cumulative + period comparison, THIS SSP,
#    full available record
# ============================================================

cat("\n=== PART 1: Annual / cumulative / period comparison, SSP", ssp, "===\n")

process_one_ssp_annual <- function(ssp) {
  
  file_df <- list_depth_resolved_files(ssp)
  if (is.null(file_df)) return(NULL)
  
  cat("SSP", ssp, "- found", nrow(file_df), "yearly files, years",
      min(file_df$year), "-", max(file_df$year), "\n")
  
  results_list <- lapply(seq_len(nrow(file_df)), function(i) {
    cat("  Processing SSP", ssp, "year", file_df$year[i], "(annual)\n")
    process_one_year_annual(
      file_df$file[i], file_df$year[i],
      n_depths = n_depths, depth_thickness = depth_thickness, depth_weights = depth_weights,
      topsoil_overlap_fraction = topsoil_overlap_fraction,
      actual_topsoil_depth = actual_topsoil_depth
    )
  })
  
  bind_rows(results_list) %>%
    arrange(Year) %>%
    mutate(
      cumulative_pg = cumsum(annual_total_pg),
      SSP_raw = ssp,
      SSP = unname(ssp_labels[ssp])
    )
}

annual_df_ssp <- process_one_ssp_annual(ssp)

if (is.null(annual_df_ssp)) {
  stop("No annual results produced for SSP ", ssp, " - check input files in ",
       file.path(base_dir, paste0("yearly_nc_", ssp)))
}

write.csv(
  annual_df_ssp,
  file.path(out_dir, paste0("mineralised_N_annual_cumulative_", ssp, ".csv")),
  row.names = FALSE
)

annual_df_ssp <- annual_df_ssp %>%
  mutate(Period = vapply(Year, assign_period, character(1)))

period_summary_ssp <- annual_df_ssp %>%
  filter(!is.na(Period)) %>%
  mutate(Period = factor(Period, levels = names(periods))) %>%
  group_by(SSP, Period) %>%
  summarise(
    n_years = n(),
    mean_volumetric_gN_m3_yr = mean(profile_mean_volumetric_gN_m3_yr, na.rm = TRUE),
    sd_volumetric_gN_m3_yr   = sd(profile_mean_volumetric_gN_m3_yr, na.rm = TRUE),
    mean_areal_gN_m2_yr      = mean(profile_mean_areal_gN_m2_yr, na.rm = TRUE),
    sd_areal_gN_m2_yr        = sd(profile_mean_areal_gN_m2_yr, na.rm = TRUE),
    mean_topsoil_areal_gN_m2_yr      = mean(topsoil_areal_gN_m2_yr, na.rm = TRUE),
    sd_topsoil_areal_gN_m2_yr        = sd(topsoil_areal_gN_m2_yr, na.rm = TRUE),
    mean_topsoil_volumetric_gN_m3_yr = mean(topsoil_volumetric_gN_m3_yr, na.rm = TRUE),
    sd_topsoil_volumetric_gN_m3_yr   = sd(topsoil_volumetric_gN_m3_yr, na.rm = TRUE),
    mean_annual_pg           = mean(annual_total_pg, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Period)

write.csv(
  period_summary_ssp,
  file.path(out_dir, paste0("volumetric_rate_period_summary_", ssp, ".csv")),
  row.names = FALSE
)

cat("SSP", ssp, "period summary (annual-basis rates):\n")
print(period_summary_ssp)

# ============================================================
# 7. PART 2: Monthly seasonal cycle by period, THIS SSP
#    (only processes years within the three periods)
#    Also derives the growing-season-only top-soil total for
#    direct comparison with Ramm et al. (2022).
# ============================================================

cat("\n=== PART 2: Monthly seasonal cycle + growing-season top-soil total, SSP", ssp, "===\n")

years_needed <- unlist(periods, use.names = FALSE)

process_one_ssp_monthly <- function(ssp) {
  
  file_df <- list_depth_resolved_files(ssp, years_filter = years_needed)
  if (is.null(file_df) || nrow(file_df) == 0) {
    warning("No files found within the requested periods for SSP ", ssp)
    return(NULL)
  }
  
  cat("SSP", ssp, "- processing", nrow(file_df), "years within requested periods\n")
  
  results_list <- lapply(seq_len(nrow(file_df)), function(i) {
    cat("  SSP", ssp, "year", file_df$year[i], "(monthly)\n")
    process_one_year_monthly(
      file_df$file[i], file_df$year[i],
      n_depths = n_depths, depth_thickness = depth_thickness, depth_weights = depth_weights,
      topsoil_overlap_fraction = topsoil_overlap_fraction,
      actual_topsoil_depth = actual_topsoil_depth
    )
  })
  
  bind_rows(results_list) %>%
    mutate(SSP_raw = ssp, SSP = unname(ssp_labels[ssp]))
}

monthly_df_ssp <- process_one_ssp_monthly(ssp)

if (is.null(monthly_df_ssp)) {
  stop("No monthly results produced for SSP ", ssp, " within the requested periods.")
}

monthly_df_ssp <- monthly_df_ssp %>%
  mutate(Period = vapply(Year, assign_period, character(1)))

# ---- 7a. Seasonal cycle: mean across years, per calendar month ----

seasonal_cycle_ssp <- monthly_df_ssp %>%
  filter(!is.na(Period)) %>%
  mutate(Period = factor(Period, levels = names(periods))) %>%
  group_by(SSP, Period, Month) %>%
  summarise(
    mean_volumetric_gN_m3_month = mean(profile_mean_volumetric_gN_m3_month, na.rm = TRUE),
    sd_volumetric_gN_m3_month   = sd(profile_mean_volumetric_gN_m3_month, na.rm = TRUE),
    mean_areal_gN_m2_month      = mean(profile_mean_areal_gN_m2_month, na.rm = TRUE),
    sd_areal_gN_m2_month        = sd(profile_mean_areal_gN_m2_month, na.rm = TRUE),
    mean_topsoil_areal_gN_m2_month      = mean(topsoil_areal_gN_m2_month, na.rm = TRUE),
    sd_topsoil_areal_gN_m2_month        = sd(topsoil_areal_gN_m2_month, na.rm = TRUE),
    mean_topsoil_volumetric_gN_m3_month = mean(topsoil_volumetric_gN_m3_month, na.rm = TRUE),
    sd_topsoil_volumetric_gN_m3_month   = sd(topsoil_volumetric_gN_m3_month, na.rm = TRUE),
    n_years = n(),
    .groups = "drop"
  )

write.csv(
  seasonal_cycle_ssp,
  file.path(out_dir, paste0("monthly_volumetric_flux_seasonal_cycle_", ssp, ".csv")),
  row.names = FALSE
)

cat("SSP", ssp, "seasonal cycle summary:\n")
print(seasonal_cycle_ssp)

# ---- 7b. Growing-season-only top-soil total, per YEAR first,
#          then averaged across years within each period.
#          This is computed at the per-year level (summing across
#          growing_season_months within each year) BEFORE averaging
#          across years, so the reported SD reflects genuine
#          interannual variability of the seasonal total itself. ----

growing_season_yearly_ssp <- monthly_df_ssp %>%
  filter(Month %in% growing_season_months) %>%
  group_by(Year, SSP, Period) %>%
  summarise(
    topsoil_areal_gN_m2_season      = sum(topsoil_areal_gN_m2_month, na.rm = TRUE),
    topsoil_volumetric_gN_m3_season = sum(topsoil_volumetric_gN_m3_month, na.rm = TRUE),
    .groups = "drop"
  )

growing_season_summary_ssp <- growing_season_yearly_ssp %>%
  filter(!is.na(Period)) %>%
  mutate(Period = factor(Period, levels = names(periods))) %>%
  group_by(SSP, Period) %>%
  summarise(
    n_years = n(),
    mean_topsoil_areal_gN_m2_season      = mean(topsoil_areal_gN_m2_season, na.rm = TRUE),
    sd_topsoil_areal_gN_m2_season        = sd(topsoil_areal_gN_m2_season, na.rm = TRUE),
    mean_topsoil_volumetric_gN_m3_season = mean(topsoil_volumetric_gN_m3_season, na.rm = TRUE),
    sd_topsoil_volumetric_gN_m3_season   = sd(topsoil_volumetric_gN_m3_season, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Period)

write.csv(
  growing_season_summary_ssp,
  file.path(out_dir, paste0("topsoil_growing_season_summary_", ssp, ".csv")),
  row.names = FALSE
)

cat("SSP", ssp, "growing-season (months", paste(growing_season_months, collapse = ","),
    ") top-soil total, for comparison with Ramm et al. (2022):\n")
print(growing_season_summary_ssp)

cat("\nPer-SSP outputs for SSP", ssp, "saved in:\n", out_dir, "\n")
cat("Run combine_ssp_results.R after all SSPs have finished to produce\n")
cat("the multi-SSP comparison tables and plots.\n")