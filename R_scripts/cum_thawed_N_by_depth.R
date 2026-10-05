# ============================================================
# Cumulative newly thawed N by soil depth
# All SSPs
# Updated to match the depth-weighting / baseline logic used in
# the latest monthly mineralisation script:
#   - uses 60N input files
#   - depth weights based on previous YEAR's ALD (not a running
#     maximum), via make_depth_weights_ALD_annual_change()
#   - baseline-relative thawed N is NOT clamped to >= 0
# ============================================================

library(terra)
library(dplyr)
library(ggplot2)

terraOptions(memfrac = 0.4, progress = 1)

# ------------------------------------------------------------
# 1. Setup
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Provide an SSP code: 126, 245, 370, or 585")
}

ssp <- args[1]

if (!ssp %in% c("126", "245", "370", "585")) {
  stop("Invalid SSP: ", ssp)
}

ssp_labels <- c(
  "126" = "SSP1-2.6",
  "245" = "SSP2-4.5",
  "370" = "SSP3-7.0",
  "585" = "SSP5-8.5"
)

ssp_colors <- c(
  "SSP1-2.6" = "blue",
  "SSP2-4.5" = "orange",
  "SSP3-7.0" = "#D73027",
  "SSP5-8.5" = "#7B3294"
)

start_year <- 1850
end_year <- 2099

profile_depth_max <- 5


f_inorg_rapid <- 0.01
f_org <- 1 - f_inorg_rapid

out_dir <- file.path(
  "monthly_mineralised",
  "cumulative_newly_thawed_N_depth_all_ssps"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Input files
# ------------------------------------------------------------
# Updated to the 60N inputs used by the latest working script

ald_files <- setNames(
  paste0("60N/mean_ALD_ssp", ssp, "_60deg.nc"),
  ssp
)

thawed_files <- setNames(
  paste0(
    "60N/arctic_total_thawed_",
    ssp,
    "_60N_mean.nc"
  ),
  ssp
)

LC_file <- "LC_remapnn_corr.nc"

common_extent <- ext(-179.95, 179.95, 60, 90)


# ------------------------------------------------------------
# 3. Depth structure from vertical temperature file
# ------------------------------------------------------------

temp_file <- paste0(
  "60N/mean_verttemp_ssp",
  ssp,
  "_60deg.nc"
)

if (!file.exists(temp_file)) {
  stop("Missing vertical temperature file: ", temp_file)
}

temp_full <- rast(temp_file)

# Extract depth midpoints from the NetCDF
depth_vals <- as.numeric(depth(temp_full))

cat("Depth midpoints read from:", temp_file, "\n")
print(depth_vals)

# Keep only layers whose midpoint is <= 5 m
depth_vals <- depth_vals[
  depth_vals <= profile_depth_max
]

n_depths <- length(depth_vals)

cat("Number of depth layers:", n_depths, "\n")
cat("Depth midpoints (m):\n")
print(depth_vals)


# ------------------------------------------------------------
# Calculate layer boundaries
# ------------------------------------------------------------

get_depth_bounds <- function(depth_mid, max_depth) {
  
  depth_bounds <- numeric(length(depth_mid) + 1)
  
  depth_bounds[1] <- 0
  
  for (i in 2:length(depth_mid)) {
    depth_bounds[i] <-
      (depth_mid[i - 1] + depth_mid[i]) / 2
  }
  
  depth_bounds[length(depth_bounds)] <-
    depth_mid[length(depth_mid)] +
    (
      depth_mid[length(depth_mid)] -
        depth_bounds[length(depth_bounds) - 1]
    )
  
  depth_bounds[depth_bounds > max_depth] <- max_depth
  
  depth_bounds
}

depth_bounds <- get_depth_bounds(
  depth_mid = depth_vals,
  max_depth = profile_depth_max
)


# ------------------------------------------------------------
# Calculate layer thickness
# ------------------------------------------------------------

depth_thickness <- diff(depth_bounds)

cat("Depth layer thicknesses (m):\n")
print(depth_thickness)

if (length(depth_thickness) != n_depths) {
  stop(
    "Number of depth thicknesses does not match ",
    "number of depth layers."
  )
}

# ------------------------------------------------------------
# 4. Land-cover classes and Palmtag parameters
# ------------------------------------------------------------

taiga_classes <- c(1, 2, 3, 4, 5, 8, 9)
tundra_classes <- c(6, 7, 10, 12, 14)
wetlands_classes <- 11
barren_classes <- c(13, 15, 16)

params <- list(
  taiga = c(
    a = 0.007,
    b = 0.097,
    k = 2.7
  ),
  tundra = c(
    a = 0.010,
    b = 0.017,
    k = 1.9
  ),
  barren = c(
    a = 0.000,
    b = 0.0161,
    k = 1.6
  )
)

profile_integral_raster <- function(z1, z2, a, b, k) {
  
  a * (z2 - z1) +
    (b / k) * (
      exp(-k * z1) -
        exp(-k * z2)
    )
}

# ------------------------------------------------------------
# 5. Function to make depth weights
# ------------------------------------------------------------

# This version matches the latest working script: it uses the
# PREVIOUS YEAR's ALD (not a running max) and distributes mass
# across the full band between last year's and this year's ALD,
# whichever direction it moved.

make_depth_weights_ALD_annual_change <- function(previous_ALD,
                                                 current_ALD,
                                                 LC,
                                                 depth_bounds,
                                                 n_depths,
                                                 params,
                                                 profile_depth_max,
                                                 land_mask) {
  
  previous_ALD <- terra::ifel(
    previous_ALD > profile_depth_max,
    profile_depth_max,
    previous_ALD
  )
  
  current_ALD <- terra::ifel(
    current_ALD > profile_depth_max,
    profile_depth_max,
    current_ALD
  )
  
  # Full cumulative depth band between last year's ALD and the
  # current year's ALD (works whether the active layer deepened
  # or shallowed).
  band_lo <- terra::ifel(previous_ALD < current_ALD, previous_ALD, current_ALD)  # min
  band_hi <- terra::ifel(previous_ALD > current_ALD, previous_ALD, current_ALD)  # max
  
  taiga_mask <- LC %in% taiga_classes
  tundra_mask <- LC %in% tundra_classes
  wetlands_mask <- LC %in% wetlands_classes
  barren_mask <- LC %in% barren_classes
  
  layer_mass_list <- vector("list", n_depths)
  
  for (d in seq_len(n_depths)) {
    
    z1_layer <- depth_bounds[d]
    z2_layer <- depth_bounds[d + 1]
    
    z1 <- terra::ifel(band_lo > z1_layer, band_lo, z1_layer)
    z2 <- terra::ifel(band_hi < z2_layer, band_hi, z2_layer)
    
    valid_overlap <- z2 > z1
    
    taiga_mass <- profile_integral_raster(
      z1, z2,
      params$taiga["a"],
      params$taiga["b"],
      params$taiga["k"]
    )
    
    tundra_mass <- profile_integral_raster(
      z1, z2,
      params$tundra["a"],
      params$tundra["b"],
      params$tundra["k"]
    )
    
    barren_mass <- profile_integral_raster(
      z1, z2,
      params$barren["a"],
      params$barren["b"],
      params$barren["k"]
    )
    
    wetland_mass <- z2 - z1
    
    layer_mass <- terra::ifel(
      taiga_mask, taiga_mass,
      terra::ifel(
        tundra_mask, tundra_mass,
        terra::ifel(
          wetlands_mask, wetland_mass,
          terra::ifel(
            barren_mask, barren_mass,
            NA
          )
        )
      )
    )
    
    layer_mass <- terra::ifel(valid_overlap, layer_mass, 0)
    layer_mass <- terra::mask(layer_mass, land_mask)
    
    layer_mass_list[[d]] <- layer_mass
  }
  
  layer_mass_stack <- rast(layer_mass_list)
  names(layer_mass_stack) <- paste0("depth_", seq_len(n_depths))
  
  total_mass <- app(layer_mass_stack, sum, na.rm = TRUE)
  
  weights <- layer_mass_stack / total_mass
  
  weights <- terra::ifel(
    total_mass > 0,
    weights,
    0
  )
  
  weights <- terra::mask(weights, land_mask)
  names(weights) <- paste0("depth_", seq_len(n_depths))
  
  weights
}

# ------------------------------------------------------------
# 6. Function for one SSP
# ------------------------------------------------------------

calculate_cumulative_depth_ssp <- function(ssp) {
  
  cat("\n====================================================\n")
  cat("Processing", ssp_labels[[ssp]], "\n")
  cat("====================================================\n")
  
  ald_file <- ald_files[[ssp]]
  thawed_file <- thawed_files[[ssp]]
  
  if (!file.exists(ald_file)) {
    stop("Missing ALD file: ", ald_file)
  }
  
  if (!file.exists(thawed_file)) {
    stop("Missing thawed-N file: ", thawed_file)
  }
  
  # ----------------------------------------------------------
  # Read rasters
  # ----------------------------------------------------------
  
  ALD <- rast(ald_file)
  thawed_total <- rast(thawed_file)
  LC <- rast(LC_file)
  
  ALD <- crop(ALD, common_extent)
  thawed_total <- crop(thawed_total, common_extent)
  LC <- crop(LC, common_extent)
  
  LC <- resample(
    LC,
    thawed_total[[1]],
    method = "near"
  )
  
  # ----------------------------------------------------------
  # Select years
  # ----------------------------------------------------------
  
  thawed_years <- as.integer(
    format(time(thawed_total), "%Y")
  )
  
  ald_years <- as.integer(
    format(time(ALD), "%Y")
  )
  
  thawed_idx <- which(
    thawed_years >= start_year &
      thawed_years <= end_year
  )
  
  ald_idx <- which(
    ald_years >= start_year &
      ald_years <= end_year
  )
  
  if (length(thawed_idx) == 0) {
    stop("No thawed-N years found for SSP ", ssp)
  }
  
  if (length(ald_idx) == 0) {
    stop("No ALD years found for SSP ", ssp)
  }
  
  thawed_total <- thawed_total[[thawed_idx]]
  ALD <- ALD[[ald_idx]]
  
  years <- thawed_years[thawed_idx]
  
  if (nlyr(ALD) != nlyr(thawed_total)) {
    stop(
      "ALD and thawed N have different layer counts for SSP ",
      ssp
    )
  }
  
  # ----------------------------------------------------------
  # Land mask and baseline correction
  # ----------------------------------------------------------
  
  land_mask <- !is.na(
    thawed_total[[nlyr(thawed_total)]]
  )
  
  land_mask <- terra::ifel(
    land_mask,
    1,
    NA
  )
  
  thawed_total <- terra::mask(
    thawed_total,
    land_mask
  )
  
  ALD <- terra::mask(
    ALD,
    land_mask
  )
  
  LC <- terra::mask(
    LC,
    land_mask
  )
  
  baseline_idx <- which(
    years >= 1850 &
      years <= 1900
  )
  
  if (length(baseline_idx) == 0) {
    stop("No 1850-1900 baseline found for SSP ", ssp)
  }
  
  baseline_thawed_N <- app(
    thawed_total[[baseline_idx]],
    mean,
    na.rm = TRUE
  )
  
  baseline_thawed_N <- terra::mask(
    baseline_thawed_N,
    land_mask
  )
  
  # Baseline-relative thawed N. NOT clamped to >= 0, to match
  # the latest working script (keeps the signed anomaly so
  # negative year-to-year changes are preserved rather than
  # floored at zero).
  thawed_permafrost_N <-
    thawed_total - baseline_thawed_N
  
  thawed_permafrost_N <- terra::mask(
    thawed_permafrost_N,
    land_mask
  )
  
  # ----------------------------------------------------------
  # Area raster
  # ----------------------------------------------------------
  
  area_rast <- cellSize(
    thawed_permafrost_N[[1]],
    unit = "m"
  )
  
  area_rast <- terra::mask(
    area_rast,
    land_mask
  )
  
  # ----------------------------------------------------------
  # Initialise cumulative totals
  # ----------------------------------------------------------
  
  # ----------------------------------------------------------
  # Year-1 contribution.
  # Matches script 2's initial organic pool: the anomaly already
  # present in year 1 (thawed_permafrost_N[[1]] relative to the
  # 1850-1900 baseline) is distributed by depth using weights
  # built from an ALD band of 0 -> ALD[[1]], and its organic
  # fraction (f_org) is added into the cumulative total up front.
  # Without this step the loop below (which starts at i = 2)
  # telescopes to thawed_permafrost_N[[last]] - thawed_permafrost_N[[1]],
  # silently dropping the year-1 anomaly.
  # ----------------------------------------------------------
  
  first_thawed_organic <- thawed_permafrost_N[[1]] * f_org
  
  zero_ALD <- ALD[[1]] * 0
  
  initial_depth_weights <- make_depth_weights_ALD_annual_change(
    previous_ALD = zero_ALD,
    current_ALD = ALD[[1]],
    LC = LC,
    depth_bounds = depth_bounds,
    n_depths = n_depths,
    params = params,
    profile_depth_max = profile_depth_max,
    land_mask = land_mask
  )
  
  first_thawed_by_depth <- first_thawed_organic * initial_depth_weights
  
  first_depth_pg <- terra::global(
    first_thawed_by_depth * area_rast,
    "sum",
    na.rm = TRUE
  )[, 1] / 1e12
  
  first_expected_pg <- terra::global(
    first_thawed_organic * area_rast,
    "sum",
    na.rm = TRUE
  )[1, 1] / 1e12
  
  cumulative_depth_pg <- first_depth_pg
  
  total_expected_pg <- first_expected_pg
  total_distributed_pg <- sum(first_depth_pg, na.rm = TRUE)
  
  # annual_depth_list now spans every year (including year 1),
  # not just years 2..n
  annual_depth_list <- vector(
    "list",
    nlyr(thawed_permafrost_N)
  )
  
  annual_depth_list[[1]] <- data.frame(
    Year = years[1],
    SSP = unname(ssp_labels[[ssp]]),
    Depth_layer = seq_len(n_depths),
    Depth_top_m = depth_bounds[1:n_depths],
    Depth_bottom_m = depth_bounds[2:(n_depths + 1)],
    Depth_mid_m = (
      depth_bounds[1:n_depths] +
        depth_bounds[2:(n_depths + 1)]
    ) / 2,
    Annual_newly_thawed_N_Pg = first_depth_pg
  )
  
  previous_thawed <- thawed_permafrost_N[[1]]
  
  # Previous YEAR's ALD (not a running maximum) -- matches the
  # latest working script's depth-weighting logic.
  previous_ALD <- ALD[[1]]
  
  # ----------------------------------------------------------
  # Annual loop
  # ----------------------------------------------------------
  
  for (i in 2:nlyr(thawed_permafrost_N)) {
    
    yr <- years[i]
    
    cat(
      ssp_labels[[ssp]],
      "- year",
      yr,
      "\n"
    )
    
    current_thawed <- thawed_permafrost_N[[i]]
    current_ALD <- ALD[[i]]
    
    # Year-to-year N increment (signed; not forced positive)
    newly_thawed_N <-
      (current_thawed - previous_thawed) * f_org
    
    newly_thawed_N <- terra::mask(
      newly_thawed_N,
      land_mask
    )
    
    # Depth distribution based on the ALD band between last
    # year and this year
    depth_weights_year <-
      make_depth_weights_ALD_annual_change(
        previous_ALD = previous_ALD,
        current_ALD = current_ALD,
        LC = LC,
        depth_bounds = depth_bounds,
        n_depths = n_depths,
        params = params,
        profile_depth_max = profile_depth_max,
        land_mask = land_mask
      )
    
    newly_thawed_by_depth <-
      newly_thawed_N * depth_weights_year
    
    annual_depth_pg <- terra::global(
      newly_thawed_by_depth * area_rast,
      "sum",
      na.rm = TRUE
    )[, 1] / 1e12
    
    cumulative_depth_pg <-
      cumulative_depth_pg + annual_depth_pg
    
    expected_pg <- terra::global(
      newly_thawed_N * area_rast,
      "sum",
      na.rm = TRUE
    )[1, 1] / 1e12
    
    distributed_pg <- sum(
      annual_depth_pg,
      na.rm = TRUE
    )
    
    total_expected_pg <-
      total_expected_pg + expected_pg
    
    total_distributed_pg <-
      total_distributed_pg + distributed_pg
    
    annual_depth_list[[1]] <- data.frame(
      Year = yr,
      SSP = unname(ssp_labels[[ssp]]),
      Depth_layer = seq_len(n_depths),
      Depth_top_m = depth_bounds[1:n_depths],
      Depth_bottom_m = depth_bounds[2:(n_depths + 1)],
      Depth_mid_m = (
        depth_bounds[1:n_depths] +
          depth_bounds[2:(n_depths + 1)]
      ) / 2,
      Annual_newly_thawed_N_Pg = annual_depth_pg
    )
    
    # Update previous thaw state
    previous_thawed <- current_thawed
    
    # Update previous ALD to THIS year's ALD (year-to-year,
    # not a running maximum)
    previous_ALD <- current_ALD
    
    rm(
      current_thawed,
      current_ALD,
      newly_thawed_N,
      depth_weights_year,
      newly_thawed_by_depth
    )
    
    gc()
  }
  
  # ----------------------------------------------------------
  # Cumulative table
  # ----------------------------------------------------------
  
  # ----------------------------------------------------------
  # Cumulative table
  # ----------------------------------------------------------
  
  cumulative_df <- data.frame(
    SSP_code = ssp,
    SSP = unname(ssp_labels[[ssp]]),
    Depth_layer = seq_len(n_depths),
    Depth_top_m = depth_bounds[1:n_depths],
    Depth_bottom_m = depth_bounds[2:(n_depths + 1)],
    Depth_mid_m = (
      depth_bounds[1:n_depths] +
        depth_bounds[2:(n_depths + 1)]
    ) / 2,
    
    # Actual cumulative N in each layer
    Cumulative_newly_thawed_N_Pg =
      cumulative_depth_pg,
    
    # Thickness of each layer
    Layer_thickness_m =
      depth_thickness,
    
    # N normalized by soil thickness
    Cumulative_newly_thawed_N_Pg_per_m =
      cumulative_depth_pg / depth_thickness
    
  ) %>%
    mutate(
      Share_percent =
        100 *
        Cumulative_newly_thawed_N_Pg /
        sum(
          Cumulative_newly_thawed_N_Pg,
          na.rm = TRUE
        )
    )
  
  annual_df <- bind_rows(annual_depth_list)
  
  conservation_df <- data.frame(
    SSP_code = ssp,
    SSP = unname(ssp_labels[[ssp]]),
    Expected_cumulative_Pg = total_expected_pg,
    Distributed_cumulative_Pg = total_distributed_pg,
    Difference_Pg =
      total_distributed_pg - total_expected_pg,
    Distributed_fraction =
      total_distributed_pg / total_expected_pg
  )
  
  cat("\nConservation check for", ssp_labels[[ssp]], "\n")
  print(conservation_df)
  
  # Save SSP-specific outputs
  write.csv(
    cumulative_df,
    file.path(
      out_dir,
      paste0(
        "cumulative_newly_thawed_N_by_depth_",
        ssp,
        "_",
        start_year,
        "_",
        end_year,
        ".csv"
      )
    ),
    row.names = FALSE
  )
  
  write.csv(
    annual_df,
    file.path(
      out_dir,
      paste0(
        "annual_newly_thawed_N_by_depth_",
        ssp,
        "_",
        start_year,
        "_",
        end_year,
        ".csv"
      )
    ),
    row.names = FALSE
  )
  
  rm(
    ALD,
    thawed_total,
    thawed_permafrost_N,
    LC,
    land_mask,
    area_rast
  )
  
  gc()
  
  list(
    cumulative = cumulative_df,
    annual = annual_df,
    conservation = conservation_df
  )
}
# 
# # ------------------------------------------------------------
# # 7. Run one SSP (this script is called once per SSP, e.g. as
# #    an array job with `ssp` passed on the command line)
# # ------------------------------------------------------------
# 
# result <- calculate_cumulative_depth_ssp(ssp)
# 
# cumulative_df <- result$cumulative
# annual_df <- result$annual
# conservation_df <- result$conservation
# 
# cat("\nFinished", ssp_labels[[ssp]], "\n")
# cat("Outputs saved in:\n", out_dir, "\n")
# 
# print(conservation_df)


# 
# # ============================================================
# # Plot cumulative newly thawed N by depth
# # Combine outputs from all SSP array jobs
# #
# # Run this AFTER cumulative_newly_thawed_N_by_depth_updated.R
# # has been run once for each SSP (126, 245, 370, 585), since it
# # reads their CSV outputs from disk.
# # ============================================================
# 
# library(dplyr)
# library(ggplot2)
# 
# # ------------------------------------------------------------
# # 1. Setup
# # ------------------------------------------------------------
# 
# ssps <- c("126", "245", "370", "585")
# 
# ssp_labels <- c(
#   "126" = "SSP1-2.6",
#   "245" = "SSP2-4.5",
#   "370" = "SSP3-7.0",
#   "585" = "SSP5-8.5"
# )
# 
# ssp_colors <- c(
#   "SSP1-2.6" = "blue",
#   "SSP2-4.5" = "orange",
#   "SSP3-7.0" = "#D73027",
#   "SSP5-8.5" = "#7B3294"
# )
# 
# start_year <- 1850
# end_year <- 2099
# profile_depth_max <- 5
# 
# out_dir <- file.path(
#   "monthly_mineralised",
#   "cumulative_newly_thawed_N_depth_all_ssps"
# )
# 
# # ------------------------------------------------------------
# # 2. Read cumulative depth CSVs
# # ------------------------------------------------------------
# 
# cumulative_list <- lapply(ssps, function(ssp) {
# 
#   file <- file.path(
#     out_dir,
#     paste0(
#       "cumulative_newly_thawed_N_by_depth_",
#       ssp,
#       "_",
#       start_year,
#       "_",
#       end_year,
#       ".csv"
#     )
#   )
# 
#   if (!file.exists(file)) {
#     stop("Missing file: ", file)
#   }
# 
#   df <- read.csv(file)
# 
#   df$SSP_code <- ssp
#   df$SSP <- unname(ssp_labels[ssp])
# 
#   df
# })
# 
# cumulative_all <- bind_rows(cumulative_list)
# 
# cumulative_all$SSP <- factor(
#   cumulative_all$SSP,
#   levels = unname(ssp_labels)
# )
# 
# 
# # Check combined data
# print(cumulative_all)
# print(
#   cumulative_all %>%
#     group_by(SSP) %>%
#     summarise(
#       Total_Pg = sum(Cumulative_newly_thawed_N_Pg, na.rm = TRUE),
#       .groups = "drop"
#     )
# )
# 
# # ------------------------------------------------------------
# # 3. Combined depth-by-SSP bar plot
# # ------------------------------------------------------------
# 
# p_depth_bars <- ggplot(
#   cumulative_all,
#   aes(
#     x = factor(Plot_mid_m),
#     y = Cumulative_newly_thawed_N_Pg,
#     fill = SSP
#   )
# ) +
#   geom_col(
#     position = position_dodge(width = 0.7),
#     width = 0.6
#   ) +
#   coord_flip() +
#   scale_x_discrete(
#     labels = function(x) sub("^-", "", x)
#   ) +
#   scale_fill_manual(values = ssp_colors) +
#   theme_bw(base_size = 13) +
#   labs(
#     title = "Cumulative newly thawed N by soil depth",
#     subtitle = paste0(start_year, "-", end_year),
#     x = "Soil depth (m)",
#     y = "Cumulative newly thawed N (Pg N)",
#     fill = "Scenario"
#   )
# print(p_depth_bars)
# 
# ggsave(
#   filename = file.path(
#     out_dir,
#     paste0(
#       "cumulative_newly_thawed_N_depth_bars_all_SSPs_",
#       start_year,
#       "_",
#       end_year,
#       ".png"
#     )
#   ),
#   plot = p_depth_bars,
#   width = 8,
#   height = 7,
#   dpi = 300
# )
# 
# # ------------------------------------------------------------
# # 3b. Combined depth-profile line plot
# # ------------------------------------------------------------
# 
# p_depth_lines <- ggplot(
#   cumulative_all,
#   aes(
#     x = Cumulative_newly_thawed_N_Pg,
#     y = Plot_mid_m,
#     colour = SSP,
#     group = SSP
#   )
# ) +
#   geom_path(linewidth = 1.2) +
#   geom_point(size = 2.5) +
#   scale_colour_manual(values = ssp_colors) +
#   scale_y_continuous(
#     limits = c(-profile_depth_max, 0),
#     breaks = seq(-profile_depth_max, 0, by = 0.5),
#     expand = expansion(mult = c(0, 0))
#   ) +
#   scale_x_continuous(
#     expand = expansion(mult = c(0.02, 0.08))
#   ) +
#   theme_bw(base_size = 13) +
#   theme(
#     legend.position = "bottom",
#     panel.grid.minor = element_blank()
#   ) +
#   labs(
#     title = "Cumulative newly thawed N by soil depth",
#     subtitle = paste0(start_year, "-", end_year),
#     x = "Cumulative newly thawed N (Pg N)",
#     y = "Soil depth (m)",
#     colour = "Scenario"
#   )
# 
# print(p_depth_lines)
# 
# ggsave(
#   filename = file.path(
#     out_dir,
#     paste0(
#       "cumulative_newly_thawed_N_depth_lines_all_SSPs_",
#       start_year,
#       "_",
#       end_year,
#       ".png"
#     )
#   ),
#   plot = p_depth_lines,
#   width = 8,
#   height = 7,
#   dpi = 300
# )
# 
# # ------------------------------------------------------------
# # 4. Faceted soil-layer plot
# # ------------------------------------------------------------
# 
# p_depth_facets <- ggplot(cumulative_all) +
#   geom_rect(
#     aes(
#       xmin = 0,
#       xmax = Cumulative_newly_thawed_N_Pg,
#       ymin = Plot_bottom_m,
#       ymax = Plot_top_m,
#       fill = SSP
#     ),
#     colour = "black",
#     linewidth = 0.2
#   ) +
#   facet_wrap(
#     ~SSP,
#     ncol = 2,
#     scales = "free_x"
#   ) +
#   scale_fill_manual(values = ssp_colors) +
#   scale_y_continuous(
#     limits = c(-profile_depth_max, 0),
#     breaks = seq(-profile_depth_max, 0, by = 0.5),
#     expand = expansion(mult = c(0, 0))
#   ) +
#   scale_x_continuous(
#     expand = expansion(mult = c(0, 0.08))
#   ) +
#   theme_bw(base_size = 13) +
#   theme(
#     legend.position = "none",
#     panel.grid.minor = element_blank()
#   ) +
#   labs(
#     title = "Cumulative newly thawed N by soil layer",
#     subtitle = paste0(start_year, "-", end_year),
#     x = "Cumulative newly thawed N in layer (Pg N)",
#     y = "Soil depth (m)"
#   )
# 
# print(p_depth_facets)
# 
# ggsave(
#   filename = file.path(
#     out_dir,
#     paste0(
#       "cumulative_newly_thawed_N_depth_facets_all_SSPs_",
#       start_year,
#       "_",
#       end_year,
#       ".png"
#     )
#   ),
#   plot = p_depth_facets,
#   width = 10,
#   height = 8,
#   dpi = 300
# )
# 
# # ------------------------------------------------------------
# # 5. Save combined CSV
# # ------------------------------------------------------------
# 
# write.csv(
#   cumulative_all,
#   file.path(
#     out_dir,
#     paste0(
#       "cumulative_newly_thawed_N_by_depth_all_SSPs_",
#       start_year,
#       "_",
#       end_year,
#       ".csv"
#     )
#   ),
#   row.names = FALSE
# )
# 
# cat("Combined plots and CSV saved in:\n", out_dir, "\n")