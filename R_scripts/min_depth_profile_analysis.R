# ============================================================
# Depth-profile analysis: volumetric mineralisation rate vs. depth
#
# Unlike profile_mean_volumetric_gN_m3_yr in the main analysis script,
# this does NOT collapse across depth with thickness-weighting (which
# mathematically cancels the depth information out - see derivation
# notes). Instead it keeps each depth layer's rate separate, so the
# resulting profile actually shows where mineralisation is concentrated.
#
# Runs once per SSP (array job). Produces one CSV with mean volumetric
# rate per depth layer per period. Run combine_depth_profiles.R
# afterward to make the faceted depth-profile plot across all SSPs.
# ============================================================

library(terra)
library(dplyr)
library(ggplot2)
library(stringr)

# -----------------------------
# 0. Setup
# -----------------------------

args <- commandArgs(trailingOnly = TRUE)
ssp <- args[1]

if (is.na(ssp) || length(ssp) == 0 || nchar(ssp) == 0) {
  stop("No SSP provided. Usage: Rscript depth_profile_analysis.R <ssp>")
}

base_dir <- "monthly_mineralised/mean_2perc_baseline_new"
out_dir <- "monthly_mineralised/depth_profile_analysis"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_dir <- "logs"
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(
  log_dir,
  paste0("depth_profile_", ssp, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
)

log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, append = TRUE, type = "message")

start_time <- Sys.time()
cat("====================================================\n")
cat("STARTED depth-profile analysis, SSP", ssp, "\n")
cat("Time:", as.character(start_time), "\n")
cat("====================================================\n")

on.exit({
  end_time <- Sys.time()
  cat("Runtime:", as.character(end_time - start_time), "\n")
  cat("Log file:", log_file, "\n")
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

# ============================================================
# 1. Recover the physical depth structure
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

temp_file <- paste0("60N/mean_verttemp_ssp", ssp, "_60deg.nc")
temp_full <- rast(temp_file)
n_depths_full <- length(unique(as.numeric(depth(temp_full))))
all_depth_vals <- as.numeric(depth(temp_full))[1:n_depths_full]
depth_keep <- which(all_depth_vals <= profile_depth_max)
depth_vals <- all_depth_vals[depth_keep]     # nominal midpoint depth of each layer (m)
n_depths <- length(depth_vals)
depth_bounds <- get_depth_bounds(depth_vals, max_depth = profile_depth_max)
depth_thickness <- diff(depth_bounds)         # thickness of each layer (m)
rm(temp_full); gc()

cat("Recovered", n_depths, "depth layers.\n")
cat("Depth midpoints (m):\n"); print(round(depth_vals, 3))
cat("Depth thickness (m):\n"); print(round(depth_thickness, 3))

# ============================================================
# 2. File-listing helper
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
# 3. Per-year processing: volumetric rate PER DEPTH LAYER
#    (no thickness-weighted collapsing across depth - each
#    layer's rate is kept separate and spatially averaged)
# ============================================================

process_one_year_depth_profile <- function(nc_path, yr, n_depths, depth_thickness, depth_vals) {
  
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
  
  # Volumetric rate per depth layer: kg N m-2 yr-1 -> g N m-3 yr-1
  # (each layer normalised by its OWN thickness only - no cross-depth
  # weighting or summation, so the depth signal is preserved)
  volumetric_by_depth <- (annual_by_depth * 1000) / depth_thickness
  
  # Spatial (area-weighted) mean, per depth layer - one number per layer
  spatial_mean_by_depth <- terra::global(
    volumetric_by_depth, "mean", weights = area_rast, na.rm = TRUE
  )[, 1]
  
  rm(r, annual_by_depth, volumetric_by_depth)
  gc()
  
  data.frame(
    Year = yr,
    depth_index = seq_len(n_depths),
    depth_mid_m = depth_vals,
    volumetric_gN_m3_yr = spatial_mean_by_depth
  )
}

# ============================================================
# 4. Process all years within the requested periods, this SSP
# ============================================================

assign_period <- function(yr) {
  for (p in names(periods)) if (yr %in% periods[[p]]) return(p)
  NA_character_
}

years_needed <- unlist(periods, use.names = FALSE)

file_df <- list_depth_resolved_files(ssp, years_filter = years_needed)

if (is.null(file_df) || nrow(file_df) == 0) {
  stop("No files found within the requested periods for SSP ", ssp)
}

cat("SSP", ssp, "- processing", nrow(file_df), "years within requested periods\n")

results_list <- lapply(seq_len(nrow(file_df)), function(i) {
  cat("  SSP", ssp, "year", file_df$year[i], "(depth profile)\n")
  process_one_year_depth_profile(
    file_df$file[i], file_df$year[i],
    n_depths = n_depths, depth_thickness = depth_thickness, depth_vals = depth_vals
  )
})

depth_profile_yearly <- bind_rows(results_list) %>%
  mutate(
    Period = vapply(Year, assign_period, character(1)),
    SSP_raw = ssp,
    SSP = unname(ssp_labels[ssp])
  )

# ============================================================
# 5. Average across years within each period, per depth layer
# ============================================================

depth_profile_by_period <- depth_profile_yearly %>%
  filter(!is.na(Period)) %>%
  mutate(Period = factor(Period, levels = names(periods))) %>%
  group_by(SSP, Period, depth_index, depth_mid_m) %>%
  summarise(
    mean_volumetric_gN_m3_yr = mean(volumetric_gN_m3_yr, na.rm = TRUE),
    sd_volumetric_gN_m3_yr   = sd(volumetric_gN_m3_yr, na.rm = TRUE),
    n_years = n(),
    .groups = "drop"
  ) %>%
  arrange(Period, depth_index)

write.csv(
  depth_profile_by_period,
  file.path(out_dir, paste0("depth_profile_by_period_", ssp, ".csv")),
  row.names = FALSE
)

cat("SSP", ssp, "depth profile by period:\n")
print(depth_profile_by_period)

cat("\nPer-SSP depth profile saved for SSP", ssp, "in:\n", out_dir, "\n")
cat("Run combine_depth_profiles.R after all SSPs finish to make the faceted plot.\n")



# ============================================================
# Combine per-SSP depth-profile results into a single faceted
# depth-profile plot (rate vs. depth, coloured by period,
# faceted by SSP). Run after all SSP runs of
# depth_profile_analysis.R have completed.
# ============================================================

library(dplyr)
library(ggplot2)
library(readr)

out_dir <- "monthly_mineralised/depth_profile_analysis"

ssps <- c("126", "245", "370", "585")

ssp_labels <- c(
  "126" = "SSP1-2.6",
  "245" = "SSP2-4.5",
  "370" = "SSP3-7.0",
  "585" = "SSP5-8.5"
)

period_colors <- c(
  "Present-day (2000-2020)"    = "steelblue",
  "Mid-century (2040-2060)"    = "darkorange",
  "End-of-century (2080-2099)" = "firebrick"
)

period_levels <- names(period_colors)

# ============================================================
# 1. Combine per-SSP CSVs
# ============================================================

profile_files <- file.path(out_dir, paste0("depth_profile_by_period_", ssps, ".csv"))
found <- file.exists(profile_files)

if (!all(found)) {
  warning("Missing depth profile files for SSP(s): ", paste(ssps[!found], collapse = ", "))
}

if (!any(found)) {
  stop("No depth profile files found. Nothing to combine.")
}

depth_profile_all <- bind_rows(
  lapply(profile_files[found], read_csv, show_col_types = FALSE)
) %>%
  mutate(
    SSP = factor(SSP, levels = unname(ssp_labels[ssps])),
    Period = factor(Period, levels = period_levels)
  )

write_csv(depth_profile_all, file.path(out_dir, "depth_profile_by_period_all_SSPs.csv"))

cat("Combined depth profile data:", nrow(depth_profile_all), "rows across",
    length(unique(depth_profile_all$SSP)), "SSPs\n")

# ============================================================
# 2. Plot: classic soil-profile style
#    x = rate, y = depth (INVERTED so surface is at top, like a
#    real soil profile), colour = period, facet = SSP
# ============================================================

p_depth_profile <- ggplot(
  depth_profile_all,
  aes(x = mean_volumetric_gN_m3_yr, y = depth_mid_m, colour = Period)
) +
  geom_path(linewidth = 1, orientation = "y") +
  geom_point(size = 1.8) +
  geom_ribbon(
    aes(
      xmin = mean_volumetric_gN_m3_yr - sd_volumetric_gN_m3_yr,
      xmax = mean_volumetric_gN_m3_yr + sd_volumetric_gN_m3_yr,
      fill = Period
    ),
    alpha = 0.12, colour = NA
  ) +
  scale_y_reverse(name = "Depth (m)") +   # surface (0 m) at top, like a real soil profile
  scale_colour_manual(values = period_colors, name = "Time period") +
  scale_fill_manual(values = period_colors, guide = "none") +
  facet_wrap(~ SSP, ncol = 2) +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom") +
  labs(
    title = "Depth profile of mineralisation rate by time period and scenario",
    subtitle = "Pan-Arctic spatial mean, mean ± SD across years within each period",
    x = expression("Volumetric mineralisation rate (g N m"^{-3}*" yr"^{-1}*")")
  )

print(p_depth_profile)

ggsave(
  file.path(out_dir, "depth_profile_by_period_all_SSPs.png"),
  p_depth_profile, width = 10, height = 8, dpi = 300
)

cat("Combined depth-profile plot saved in:\n", out_dir, "\n")