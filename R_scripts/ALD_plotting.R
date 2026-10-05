
library(terra)
library(here)
library(ggplot2)
library(dplyr)
library(tidyr)
ext_60N <- ext(-179.95, 179.95, 60, 90)


# Load present-day satellite ALD
ALD_ESA_REF_mean <- rast("ESA_ALD/ESA_ALD_regrid.nc", lyr = 1)
ALD_ESA_REF_std <- rast("ESA_ALD/ESA_ALD_regrid.nc", lyr = 2)

plot(ALD_ESA_REF_mean)
# ============================================================
# FUNCTION: Bias-correct SSP scenario rasters
# ============================================================
bias_correct_ssp <- function(mean_rast, std_rast, esa_mean, esa_std, 
                             start_year = 2000, end_year = 2020) {
  
  # Get time indices for reference period
  rast_time <- time(mean_rast)
  present_idx <- which(rast_time >= start_year & rast_time <= end_year)
  
  # Calculate mean of model during reference period
  model_present_mean <- mean(mean_rast[[present_idx]], na.rm = TRUE)
  
  # Calculate bias (model - observation)
  bias_mean <- model_present_mean - esa_mean
  
  # Calculate standard deviation bias (if needed)
  model_present_std <- mean(std_rast[[present_idx]], na.rm = TRUE)
  bias_std <- model_present_std - esa_std
  
  # Apply bias correction
  mean_corrected <- mean_rast - bias_mean
  std_corrected <-  std_rast - bias_std
  
  # Return list with corrected rasters and bias
  return(list(mean_corrected = mean_corrected, 
              std_corrected = std_corrected,
              bias_mean = bias_mean,
              bias_std = bias_std))
}

# ============================================================
# APPLY BIAS CORRECTION TO EACH SCENARIO
# ============================================================

# SSP585
rast_mean_585 <- rast("mean_ssp585_corr_new.nc")
rast_std_585 <- rast("std_ssp585_corr_new.nc")
rast_mean_585 <- crop(rast("mean_ssp585_corr_new.nc"), ext_60N)
rast_std_585  <- crop(rast("std_ssp585_corr_new.nc"), ext_60N)

corr_585 <- bias_correct_ssp(rast_mean_585, rast_std_585, 
                             ALD_ESA_REF_mean, ALD_ESA_REF_std)

# mean ALD all the cell areas (total area of the grid)
mean_ssp585 <- global(corr_585$mean_corrected, fun = "mean", na.rm = TRUE)
std_ssp585<- global(corr_585$std_corrected, fun = "mean", na.rm = TRUE)

mean_ssp585_old <- global(rast_mean_585, fun = "mean", na.rm = TRUE)
std_ssp585_old <- global(rast_std_585, fun = "mean", na.rm = TRUE)



# SSP370
rast_mean_370 <- rast("mean_ssp370_corr_new.nc")
rast_std_370 <- rast("std_ssp370_corr_new.nc")
rast_mean_370 <- crop(rast("mean_ssp370_corr_new.nc"), ext_60N)
rast_std_370  <- crop(rast("std_ssp370_corr_new.nc"), ext_60N)

corr_370 <- bias_correct_ssp(rast_mean_370, rast_std_370, 
                             ALD_ESA_REF_mean, ALD_ESA_REF_std)

# mean ALD all the cell areas (total area of the grid)
mean_ssp370 <- global(corr_370$mean_corrected, fun = "mean", na.rm = TRUE)
std_ssp370<- global(corr_370$std_corrected, fun = "mean", na.rm = TRUE)

# SSP245
rast_mean_245 <- rast("mean_ssp245_corr_new.nc")
rast_std_245 <- rast("std_ssp245_corr_new.nc")
rast_mean_245 <- crop(rast("mean_ssp245_corr_new.nc"), ext_60N)
rast_std_245  <- crop(rast("std_ssp245_corr_new.nc"), ext_60N)

corr_245 <- bias_correct_ssp(rast_mean_245, rast_std_245, 
                             ALD_ESA_REF_mean, ALD_ESA_REF_std)

# mean ALD all the cell areas (total area of the grid)
mean_ssp245 <- global(corr_245$mean_corrected, fun = "mean", na.rm = TRUE)
std_ssp245 <- global(corr_245$std_corrected, fun = "mean", na.rm = TRUE)

# SSP126
rast_mean_126 <- rast("mean_ssp126_corr_new.nc")
rast_std_126 <- rast("std_ssp126_corr_new.nc")
rast_mean_126 <- crop(rast("mean_ssp126_corr_new.nc"), ext_60N)
rast_std_126  <- crop(rast("std_ssp126_corr_new.nc"), ext_60N)

corr_126 <- bias_correct_ssp(rast_mean_126, rast_std_126, 
                             ALD_ESA_REF_mean, ALD_ESA_REF_std)
# mean ALD all the cell areas (total area of the grid)
mean_ssp126 <- global(corr_126$mean_corrected, fun = "mean", na.rm = TRUE)
std_ssp126 <- global(corr_126$std_corrected, fun = "mean", na.rm = TRUE)


df<-data.frame(mean_ssp585, mean_ssp370, mean_ssp245, mean_ssp126, std_ssp585, std_ssp370, std_ssp245, std_ssp126)
df$Year<-rep(1850:2099)
names(df)<-c("mean_585", "mean_370", "mean_245", "mean_126", "std_585", "std_370", "std_245","std_126", "Year")
write.csv(df, "ALD_60N_bias_corrected.csv")

corr<-read.csv("ALD_no_lim_bias_corr.csv")

tail(corr2)
ald_old<-read.csv("ALD_final.csv")
tail(ald_old)
# ============================================================
# VISUALIZE BIAS MAPS
# ============================================================


library(terra)
library(sf)
library(ggplot2)
library(dplyr)
library(rnaturalearth)

# ============================================================
# Define target polar projection
# ============================================================
polar_crs <- "+proj=stere +lat_0=90 +lat_ts=71 +lon_0=0 +datum=WGS84 +units=m"

# ============================================================
# Reproject each bias raster to polar stereographic
# ============================================================
project_bias <- function(rast, res_m = 25000) {
  project(rast, polar_crs, res = res_m, method = "bilinear")
}

bias_585_p <- project_bias(corr_585$bias_mean)
bias_370_p <- project_bias(corr_370$bias_mean)
bias_245_p <- project_bias(corr_245$bias_mean)
bias_126_p <- project_bias(corr_126$bias_mean)

# ============================================================
# Convert projected rasters to dataframes
# ============================================================
raster_to_df <- function(rast, scenario_name) {
  df <- as.data.frame(rast, xy = TRUE, na.rm = TRUE)
  names(df)[3] <- "bias"
  df$scenario <- scenario_name
  return(df)
}

bias_all <- bind_rows(
  raster_to_df(bias_585_p, "SSP5-8.5"),
  raster_to_df(bias_370_p, "SSP3-7.0"),
  raster_to_df(bias_245_p, "SSP2-4.5"),
  raster_to_df(bias_126_p, "SSP1-2.6")
)

bias_all$scenario <- factor(bias_all$scenario,
                            levels = c("SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5"))

max_abs_bias <- max(abs(bias_all$bias), na.rm = TRUE)

# ============================================================
# Get coastlines, transform to same polar projection
# ============================================================
coastlines <- ne_countries(scale = "medium", returnclass = "sf")
coastlines_polar <- st_transform(coastlines, crs = polar_crs)

# Determine plot extent from the raster bounds (so coastlines don't extend
# beyond the 60-90N crop)
plot_extent <- st_bbox(c(
  xmin = min(bias_all$x), xmax = max(bias_all$x),
  ymin = min(bias_all$y), ymax = max(bias_all$y)
), crs = st_crs(polar_crs))

# ============================================================
# Plot
# ============================================================
p <- ggplot() +
  geom_raster(data = bias_all, aes(x = x, y = y, fill = bias)) +
  geom_sf(data = coastlines_polar, fill = NA, color = "grey20", linewidth = 0.2) +
  coord_sf(
    xlim = c(plot_extent["xmin"], plot_extent["xmax"]),
    ylim = c(plot_extent["ymin"], plot_extent["ymax"]),
    crs = polar_crs,
    expand = FALSE
  ) +
  scale_fill_gradient2(
    low = "#2166ac", mid = "white", high = "#b2182b",
    midpoint = 0,
    limits = c(-max_abs_bias, max_abs_bias),
    name = "ALD bias (m)"
  ) +
  facet_wrap(~ scenario, ncol = 2) +
  labs(
    title = "Model bias in Active Layer Depth (2000\u20132020)",
    subtitle = "Model mean minus ESA Permafrost_cci observed mean",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 6) +
  theme(
    panel.grid = element_line(color = "grey90", linewidth = 0.2),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.key.width = unit(2, "cm")
  )

print(p)

p <- p +
  theme(
    legend.margin = margin(t = 0),
    legend.box.margin = margin(1, 1, 1, 1),
    plot.margin = margin(3, 3, 2, 3),
    legend.position = "bottom", 
    legend.key.size = unit(0.3, "cm"),      # size of each colour swatch
    legend.text = element_text(size =5),     # label text
    legend.title = element_text(size = 5),    # legend title text
    
    # --- shrink the gap between axis titles and the panel ---
    axis.title.x = element_text(margin = margin(t = 2)),   # was default ~half a line
    axis.title.y = element_text(margin = margin(r = 2)),
    
    # --- pull tick labels closer to the axis (less padding) ---
    axis.text.x = element_text(margin = margin(t = 1)),
    axis.text.y = element_text(margin = margin(r = 1)),
    
    # --- optionally shrink the fonts a touch, which shrinks their box too ---
    axis.title = element_text(size = 6),
    axis.text  = element_text(size = 6)
  )

ggsave(
  "plot_ALD_bias.pdf",
  p,
  width = 15,
  height = 12,
  units = "cm", 
  dpi = 500
)

# CALCULATE GLOBAL AREA-WEIGHTED MEANS (CORRECTED)
# ============================================================

# Calculate area weights once (using any corrected raster)
area_weights <- cellSize(corr_585$mean_corrected[[1]], mask = TRUE, unit = "km")
total_area <- global(area_weights, "sum", na.rm = TRUE)[1,1]
area_weights_norm <- area_weights / total_area

# Function to calculate weighted mean time series
calc_weighted_ts <- function(mean_corrected, area_weights_norm) {
  n_years <- nlyr(mean_corrected)
  weighted_mean <- numeric(n_years)
  
  for(i in 1:n_years) {
    weighted_grid <- mean_corrected[[i]] * area_weights_norm
    weighted_mean[i] <- global(weighted_grid, "sum", na.rm = TRUE)[1,1]
  }
  
  return(weighted_mean)
}

# Calculate time series for all scenarios
years <- time(corr_585$mean_corrected)  # Extract years from raster

ts_585 <- data.frame(
  year = years,
  ssp585 = calc_weighted_ts(corr_585$mean_corrected, area_weights_norm)
)

ts_370 <- data.frame(
  year = years,
  ssp370 = calc_weighted_ts(corr_370$mean_corrected, area_weights_norm)
)

ts_245 <- data.frame(
  year = years,
  ssp245 = calc_weighted_ts(corr_245$mean_corrected, area_weights_norm)
)

ts_126 <- data.frame(
  year = years,
  ssp126 = calc_weighted_ts(corr_126$mean_corrected, area_weights_norm)
)

# Combine all scenarios
all_ts <- ts_126 %>%
  left_join(ts_245, by = "year") %>%
  left_join(ts_370, by = "year") %>%
  left_join(ts_585, by = "year")

# ============================================================
# PLOT TIME SERIES
# ============================================================

# Reshape for ggplot
ts_long <- all_ts %>%
  pivot_longer(cols = -year, names_to = "scenario", values_to = "ald_change")

# Plot all scenarios
ggplot(ts_long, aes(x = year, y = ald_change, color = scenario)) +
  geom_line(size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  scale_color_manual(values = c("ssp126" = "blue", 
                                "ssp245" = "green", 
                                "ssp370" = "orange", 
                                "ssp585" = "red")) +
  labs(
    title = "Area-Weighted Global Mean ALD Change (Bias-Corrected)",
    subtitle = "Relative to ESA present-day reference (2000-2020 baseline)",
    x = "Year",
    y = "ALD Change (units)",
    color = "SSP Scenario"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "bottom"
  )

# ============================================================
# COMPARE RAW VS CORRECTED FOR ONE SCENARIO
# ============================================================

# Calculate raw (uncorrected) weighted mean
raw_585 <- calc_weighted_ts(rast_mean_585, area_weights_norm)

comparison_df <- data.frame(
  year = years,
  raw = raw_585,
  corrected = ts_585$ssp585
) %>%
  pivot_longer(cols = c(raw, corrected), names_to = "type", values_to = "ald")

ggplot(comparison_df, aes(x = year, y = ald, color = type)) +
  geom_line(size = 1) +
  labs(
    title = "SSP585: Raw vs ESA Bias-Corrected ALD",
    x = "Year",
    y = "ALD [m]"
  ) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("raw" = "gray", "corrected" = "red"))


ALD_bias_corrected<-read.csv("ALD_no_lim_bias_corr.csv")
ALD_not_corrected<-read.csv("ALD_final.csv")

ALD_bias_corrected$Type <- "Bias corrected"
ALD_not_corrected$Type  <- "Not corrected"

combined_ALD <- dplyr::bind_rows(ALD_bias_corrected, ALD_not_corrected)

long_ALD <- combined_ALD %>%
  pivot_longer(
    cols = starts_with("mean_"),
    names_to = "SSP",
    values_to = "Mean_ALD"
  ) %>%
  dplyr::mutate(
    SSP = gsub("mean_", "SSP", SSP)
  )

ggplot(long_ALD, aes(x = Year, y = Mean_ALD, color = Type)) +
  geom_line() +
  facet_wrap(~SSP) +
  theme_minimal() +
  labs(y = "Active layer depth [m]", color = "")

compute_anomaly_present <- function(df) {
  
  ref <- df %>% dplyr::filter(Year >= 2000 & Year <= 2020)
  
  ref_mean <- sapply(df[,-1], function(x) mean(x[ref$Year >= 2000 & ref$Year <= 2020], na.rm = TRUE))
  
  df_anom <- df
  
  for (col in names(df)[-1]) {
    df_anom[[col]] <- df[[col]] - ref_mean[col]
  }
  
  return(df_anom)
}

convert_numeric <- function(df) {
  df$Year <- as.integer(df$Year)
  
  df[,-1] <- lapply(df[,-1], function(x) as.numeric(as.character(x)))
  
  return(df)
}

ALD_bias_corrected <- convert_numeric(ALD_bias_corrected)
ALD_not_corrected  <- convert_numeric(ALD_not_corrected)

ALD_bias_anom <- compute_anomaly_present(ALD_bias_corrected)
ALD_no_anom   <- compute_anomaly_present(ALD_not_corrected)

ALD_bias_anom$Type <- "Bias corrected"
ALD_no_anom$Type   <- "Not corrected"

combined <- dplyr::bind_rows(ALD_bias_anom, ALD_no_anom)

library(tidyr)

long_anom <- combined %>%
  pivot_longer(
    cols = starts_with("mean_"),
    names_to = "SSP",
    values_to = "Mean_ALD"
  ) %>%
  dplyr::mutate(SSP = gsub("mean_", "SSP", SSP))

library(ggplot2)

ggplot(long_anom, aes(x = Year, y = Mean_ALD, color = Type)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line() +
  facet_wrap(~SSP) +
  theme_minimal() +
  labs(
    y = "ALD anomaly relative to 2000–2020 [m]",
    color = ""
  )







df<-read.csv("ALD_no_lim_bias_corr.csv")

df<-read.csv("ALD_60N_bias_corrected.csv")

library(dplyr)

df <- read.csv("ALD_60N_bias_corrected.csv")


# ============================================================
# Compute period means AND period stds for each scenario
# ============================================================

summarise_period <- function(data, y_start, y_end) {
  data %>%
    filter(Year >= y_start, Year <= y_end) %>%
    summarise(
      across(starts_with("mean_"), \(x) mean(x, na.rm = TRUE)),
      #across(starts_with("std_"),  \(x) sqrt(mean(x^2, na.rm = TRUE)))  # RMS of variances -> pooled std
      across(starts_with("std_"),  \(x) (mean(x, na.rm = TRUE)))
     )
}

baseline_stats      <- summarise_period(df, 2000, 2020)
future_stats        <- summarise_period(df, 2080, 2100)
preindustrial_stats <- summarise_period(df, 1880, 1900)

# ============================================================
# Compute deltas (means) and propagated uncertainty (stds)
# ============================================================

scenarios <- c("585", "370", "245", "126")

df_change <- data.frame(
  setNames(
    lapply(scenarios, function(s) {
      future_stats[[paste0("mean_", s)]] - baseline_stats[[paste0("mean_", s)]]
    }),
    paste0("delta_", scenarios)
  ),
  setNames(
    lapply(scenarios, function(s) {
      sqrt(future_stats[[paste0("std_", s)]]^2 + baseline_stats[[paste0("std_", s)]]^2)
    }),
    paste0("delta_std_", scenarios)
  )
)

df_change_preindustrial <- data.frame(
  setNames(
    lapply(scenarios, function(s) {
      baseline_stats[[paste0("mean_", s)]] - preindustrial_stats[[paste0("mean_", s)]]
    }),
    paste0("delta_", scenarios)
  ),
  setNames(
    lapply(scenarios, function(s) {
      sqrt(baseline_stats[[paste0("std_", s)]]^2 + preindustrial_stats[[paste0("std_", s)]]^2)
    }),
    paste0("delta_std_", scenarios)
  )
)

# ============================================================
# Convert period-mean deltas (and their std) to a rate per decade
# ============================================================

midpoint_baseline      <- mean(2000:2020)
midpoint_future        <- mean(2080:2100)
midpoint_preindustrial <- mean(1880:1900)

decades_baseline_to_future        <- (midpoint_future - midpoint_baseline) / 10
decades_preindustrial_to_baseline <- (midpoint_baseline - midpoint_preindustrial) / 10



df_change <- df_change %>%
  mutate(
    rate_per_decade_585 = delta_585 / decades_baseline_to_future,
    rate_per_decade_370 = delta_370 / decades_baseline_to_future,
    rate_per_decade_245 = delta_245 / decades_baseline_to_future,
    rate_per_decade_126 = delta_126 / decades_baseline_to_future,
    # std scales linearly with a division by a constant
    rate_per_decade_std_585 = delta_std_585 / decades_baseline_to_future,
    rate_per_decade_std_370 = delta_std_370 / decades_baseline_to_future,
    rate_per_decade_std_245 = delta_std_245 / decades_baseline_to_future,
    rate_per_decade_std_126 = delta_std_126 / decades_baseline_to_future
  )

df_change_preindustrial <- df_change_preindustrial %>%
  mutate(
    rate_per_decade_585 = delta_585 / decades_preindustrial_to_baseline,
    rate_per_decade_370 = delta_370 / decades_preindustrial_to_baseline,
    rate_per_decade_245 = delta_245 / decades_preindustrial_to_baseline,
    rate_per_decade_126 = delta_126 / decades_preindustrial_to_baseline,
    rate_per_decade_std_585 = delta_std_585 / decades_preindustrial_to_baseline,
    rate_per_decade_std_370 = delta_std_370 / decades_preindustrial_to_baseline,
    rate_per_decade_std_245 = delta_std_245 / decades_preindustrial_to_baseline,
    rate_per_decade_std_126 = delta_std_126 / decades_preindustrial_to_baseline
  )

print(df_change %>% select(starts_with("delta"), starts_with("rate_per_decade")))
print(df_change_preindustrial %>% select(starts_with("delta"), starts_with("rate_per_decade")))

# Convert to cm
df_change <- df_change %>%
  mutate(
    rate_per_decade_585_cm = rate_per_decade_585 * 100,
    rate_per_decade_370_cm = rate_per_decade_370 * 100,
    rate_per_decade_245_cm = rate_per_decade_245 * 100,
    rate_per_decade_126_cm = rate_per_decade_126 * 100,
    rate_per_decade_std_585_cm = rate_per_decade_std_585 * 100,
    rate_per_decade_std_370_cm = rate_per_decade_std_370 * 100,
    rate_per_decade_std_245_cm = rate_per_decade_std_245 * 100,
    rate_per_decade_std_126_cm = rate_per_decade_std_126 * 100
  )





library(dplyr)

df <- read.csv("ALD_60N_bias_corrected.csv")

scenarios <- c("585", "370", "245", "126")

# ============================================================
# Get values at 2000 and 2024 for each scenario
# ============================================================

val_2000 <- df %>% filter(Year == 2000)
val_2024 <- df %>% filter(Year == 2024)

# ============================================================
# Compute endpoint difference, annual rate, and propagated std
# ============================================================

endpoint_trend <- data.frame(
  Scenario = paste0("SSP", scenarios),
  
  delta_m = sapply(scenarios, function(s) {
    val_2024[[paste0("mean_", s)]] - val_2000[[paste0("mean_", s)]]
  }),
  
  rate_m_per_yr = sapply(scenarios, function(s) {
    (val_2024[[paste0("mean_", s)]] - val_2000[[paste0("mean_", s)]]) / (2024 - 2000)
  }),
  
  delta_std_m = sapply(scenarios, function(s) {
    val_2024[[paste0("std_", s)]] - val_2000[[paste0("std_", s)]]
  }),
  
  rate_std_m_per_yr = sapply(scenarios, function(s) {
    (val_2024[[paste0("std_", s)]] - val_2000[[paste0("std_", s)]]) / (2024 - 2000)
  })
)

# Convert to cm for direct comparison with Streletskiy et al.
endpoint_trend <- endpoint_trend %>%
  mutate(
    rate_cm_per_yr = rate_m_per_yr * 100,
    rate_std_cm_per_yr = rate_std_m_per_yr * 100
  )

print(endpoint_trend)



library(dplyr)

df <- read.csv("ALD_60N_bias_corrected.csv")

scenarios <- c("585", "370", "245", "126")
# Function to compute anomaly for a given SSP scenario
compute_anomaly <- function(df, mean_col, std_col) {
  # Select relevant columns
  ssp_data <- df[, c("Year", mean_col, std_col)]
  
  # Compute reference period mean and standard deviation (1990-2010)
  ref_data <- ssp_data[ssp_data$Year >= 2000 & ssp_data$Year <= 2020, ]
  
  ref_mean <- mean(ref_data[[mean_col]], na.rm = TRUE)
  ref_std  <- mean(ref_data[[std_col]], na.rm = TRUE)
  
  # Compute anomalies
  ssp_data[[mean_col]] <- ssp_data[[mean_col]] - ref_mean
  ssp_data[[std_col]]  <- ssp_data[[std_col]] - ref_std  # Adjust standard deviation as well
  
  return(ssp_data)
}


# Apply the function to all SSPs
ssp585_a <- compute_anomaly(df, "mean_585", "std_585")
ssp370_a <- compute_anomaly(df, "mean_370", "std_370")
ssp245_a <- compute_anomaly(df, "mean_245", "std_245")
ssp126_a <- compute_anomaly(df, "mean_126", "std_126")

# Combine into a single data frame
anomaly_total <- Reduce(function(x, y) merge(x, y, by = "Year"), list(ssp585_a, ssp370_a, ssp245_a, ssp126_a))

# View the final dataset
head(anomaly_total)
anomaly_total<-anomaly_total[,c("Year", "mean_585", "std_585","mean_370","std_370","mean_245", "std_245", "mean_126","std_126")]


# Create a new column for color based on the Year
anomaly_total$color_period <- ifelse(anomaly_total$Year <= 2014, "black", "colored")

# Create separate dataframes for before and after 2015
before_2015 <- anomaly_total[anomaly_total$Year <= 2014, ]
after_2015 <- anomaly_total[anomaly_total$Year > 2014, ]

# Add period column
before_2015$Period <- "Before 2015"
after_2015$Period <- "After 2015"

# Combine datasets
combined_data <- bind_rows(before_2015, after_2015)

# Convert to long format
long_data <- combined_data %>%
  pivot_longer(cols = starts_with("mean_"), names_to = "SSP", values_to = "Mean_ALD") %>%
  pivot_longer(cols = starts_with("std_"), names_to = "SSP_std", values_to = "STD_ALD") %>%
  filter(gsub("mean_", "", SSP) == gsub("std_", "", SSP_std)) %>%
  dplyr::select(-SSP_std) %>%
  mutate(SSP = gsub("mean_", "SSP", SSP))  # Rename scenarios



# Define SSP colors for after 2015
ssp_colors <- c("SSP5-8.5" = "#7B3294", "SSP3-7.0" = "#D73027", "SSP2-4.5" = "orange", "SSP1-2.6" = "blue")


library(zoo)
library(zoo)

# Compute 20-year rolling mean and RMS-pooled standard deviation for each SSP
long_data <- long_data %>%
  group_by(SSP) %>%
  mutate(
    Rolling_Mean_ALD = rollmean(Mean_ALD, k = 20, fill = NA, align = "center"),
    Rolling_STD_ALD  = sqrt(rollapply(STD_ALD^2, width = 20, FUN = mean, fill = NA, align = "center"))
  ) %>%
  ungroup()


long_data$SSP <- factor(
  long_data$SSP,
  levels = c("SSP126", "SSP245", "SSP370", "SSP585"),
  labels = c("SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5")
)

long_data <- long_data %>%
  mutate(
    Mean_ALD = -Mean_ALD,
    Rolling_Mean_ALD = -Rolling_Mean_ALD
  )


# ============================================================
# Period-average version (e.g., mean over 2080-2099, matching
# your "future" period definition elsewhere in the analysis)
# ============================================================

summary_future_period <- long_data %>%
  filter(Year >= 2080 & Year <= 2099, Period == "After 2015") %>%
  group_by(SSP) %>%
  summarise(
    period_mean_ALD = mean(Mean_ALD, na.rm = TRUE),
    period_mean_SD  = mean(STD_ALD, na.rm = TRUE),   # average of yearly inter-model SDs
    .groups = "drop"
  ) %>%
  mutate(
    label = sprintf("%.2f ± %.2f m", period_mean_ALD, period_mean_SD)
  )

print(summary_future_period)


# Plot
ALD<-ggplot(long_data, aes(x = Year, y = Mean_ALD, group = SSP)) +
  
  # Before 2015: Black lines & grey ribbons
  geom_line(data = filter(long_data, Period == "Before 2015"), color = "black", linewidth = 0.8) +
  #geom_ribbon(data = filter(long_data, Period == "Before 2015"),
  #            aes(ymin = Mean_ALD - STD_ALD, ymax = Mean_ALD + STD_ALD),
  #            fill = "grey", alpha = 0.2) +
  
  # After 2015: Colored lines & ribbons
  geom_line(data = filter(long_data, Period == "After 2015"),
            aes(color = SSP), linewidth = 0.6) +
  #geom_ribbon(data = filter(long_data, Period == "After 2015"),
  #            aes(ymin = Mean_ALD - STD_ALD, ymax = Mean_ALD + STD_ALD, fill = SSP),
  #            alpha = 0.2) +
  
  # Rolling Mean before 2015 (Black line)
  geom_line(data = filter(long_data, Year < 2015),
            aes(y = Rolling_Mean_ALD), color = "black", linewidth = 0.6, linetype = "solid") +
  
  # Rolling Mean after 2015 (Colored lines)
  geom_line(data = filter(long_data, Year >= 2015),
            aes(y = Rolling_Mean_ALD, color = SSP), linewidth = 0.6, linetype = "solid") +
  
  # Rolling STD before 2015 (Grey ribbon)
  geom_ribbon(data = filter(long_data, Year < 2015),
              aes(ymin = Rolling_Mean_ALD - Rolling_STD_ALD, 
                  ymax = Rolling_Mean_ALD + Rolling_STD_ALD),
              fill = "grey", alpha = 0.1) +
  
  ylim(-0.75,3.0)+
  # Rolling STD after 2015 (Colored ribbons)
  geom_ribbon(data = filter(long_data, Year >= 2015),
  aes(ymin = Rolling_Mean_ALD - Rolling_STD_ALD, 
  ymax = Rolling_Mean_ALD + Rolling_STD_ALD, 
  fill = SSP),
  alpha = 0.1) +
  
  # Vertical line for 2015
  geom_vline(xintercept = 2015, linetype = "dashed", color = "black") +
  
  # Vertical line for 2015
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  
  # Labels and title
  labs(x = "Year", y = "Active layer depth [m]", title = "", 
       color = "SSP Scenario", fill = "SSP Scenario") +
  
  # Color scales
  scale_color_manual(values = ssp_colors) +
  scale_fill_manual(values = ssp_colors) +
  
   scale_y_continuous(
     breaks = c(-3, -2, -1, 0),
     labels = c("3", "2", "1", "0")
   )+
  # Theme settings
  theme_minimal(base_size = 15) +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 10),
    title = element_text(size = 10)
  )

ALD

ggsave(
  "monthly_mineralised/results_graphs/Fig1_ALD.png",
  plot = ALD,
  width = 12,
  height = 12,
  units = "cm",
  dpi = 900
)

##




# Load required libraries
library(terra)
library(ggplot2)
library(dplyr)
library(tidyr)


# Load NetCDF thawed nitrogen rasters
rast_mean_585 <- rast("60N/mean_ALD_ssp585_60deg.nc")
rast_std_585 <- rast("60N/std_ALD_ssp585_60deg.nc")

ext_60N <- ext(-179.95, 179.95, 60, 90)
# Load present-day satellite ALD
ALD_ESA_REF_mean <- rast("ESA_ALD_present_day.nc", lyr = 1)
ALD_ESA_REF_std <- rast("ESA_ALD_present_day.nc", lyr = 2)
ALD_ESA_REF_mean <- crop(ALD_ESA_REF_mean, ext_60N)
ALD_ESA_REF_std  <- crop(ALD_ESA_REF_std,  ext_60N)

present_idx <- which(time(rast_mean_585) >= 2000 & time(rast_mean_585) <= 2020)
ALD_present <- rast_mean_585[[present_idx]]
ALD_present_mean <- mean(ALD_present, na.rm = TRUE)
bias <- ALD_present_mean - ALD_ESA_REF_mean
plot(bias)

# bias correction of present-day ALD
rast_mean_585_corrected <- rast_mean_585 - bias
rast_mean_585_corrected[rast_mean_585_corrected < 0] <- 0
# Calculate area weights (grid cell area in km² or m²)
# This accounts for decreasing grid cell size toward poles
area_weights <- cellSize(rast_mean_585_corrected[[1]], mask = TRUE, unit = "km")
total_area <- global(area_weights, "sum", na.rm = TRUE)[1,1]
area_weights_normalized <- area_weights / total_area

# Calculate area-weighted global mean for each time step
n_years <- nlyr(rast_mean_585_corrected)
weighted_mean_ts <- numeric(n_years)

for(i in 1:n_years) {
  # Multiply each grid cell value by its normalized area weight
  weighted_grid <- rast_mean_585_corrected[[i]] * area_weights_normalized
  # Sum across all grid cells (now area-weighted)
  weighted_mean_ts[i] <- global(weighted_grid, "sum", na.rm = TRUE)[1,1]
}

# Also calculate weighted standard deviation (more complex)
# For uncertainty propagation with area weighting
weighted_std_ts <- numeric(n_years)
for(i in 1:n_years) {
  # Variance of weighted mean = sum(w_i^2 * sigma_i^2)
  var_weighted <- global((rast_std_585_corrected[[i]]^2) * (area_weights_normalized^2), 
                         "sum", na.rm = TRUE)[1,1]
  weighted_std_ts[i] <- sqrt(var_weighted)
}

# Create data frame for plotting
years <- seq(2020, by = 1, length.out = n_years)  # Adjust as needed
ts_data <- data.frame(
  year = years,
  mean_ald = weighted_mean_ts,
  std_ald = weighted_std_ts
)

# Plot area-weighted time series with uncertainty
ggplot(ts_data, aes(x = year, y = mean_ald)) +
  geom_ribbon(aes(ymin = mean_ald - std_ald, 
                  ymax = mean_ald + std_ald), 
              alpha = 0.3, fill = "steelblue") +
  geom_line(color = "steelblue", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.5) +
  labs(
    title = "Area-Weighted Global Mean Active Layer Depth (SSP585)",
    subtitle = "Corrected with present-day ESA ALD reference | Area-weighted for latitude",
    x = "Year",
    y = "Mean ALD Change (units)",
    caption = "Positive values = deeper active layer than present day"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

# For comparison, also calculate simple (unweighted) mean
simple_mean <- global(rast_mean_585_corrected, "mean", na.rm = TRUE)
ts_data$simple_mean <- simple_mean[,1]

# Compare weighted vs unweighted
comparison_data <- ts_data %>%
  select(year, mean_ald, simple_mean) %>%
  pivot_longer(cols = c(mean_ald, simple_mean), 
               names_to = "method", 
               values_to = "ald_value")

ggplot(comparison_data, aes(x = year, y = ald_value, color = method)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("mean_ald" = "steelblue", "simple_mean" = "orange"),
                     labels = c("mean_ald" = "Area-weighted", "simple_mean" = "Unweighted")) +
  labs(
    title = "Global Mean ALD: Area-Weighted vs Unweighted",
    subtitle = "Area-weighting is critical due to latitudinal variation in grid cell size",
    x = "Year",
    y = "Mean ALD Change (units)",
    color = "Method"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

# Optional: Calculate and plot zonal means (e.g., by latitude bands)
# Create latitude bands
lat <- init(rast_mean_585_corrected[[1]], "y")
lat_bands <- cut(lat, breaks = seq(-90, 90, by = 30), 
                 labels = c("90-60°S", "60-30°S", "30-0°S", 
                            "0-30°N", "30-60°N", "60-90°N"))

# Calculate area-weighted mean for each latitude band
zonal_weighted <- matrix(NA, nrow = n_years, ncol = 6)
colnames(zonal_weighted) <- levels(lat_bands)

for(i in 1:n_years) {
  for(j in 1:6) {
    # Mask to latitude band
    band_mask <- ifel(lat_bands == levels(lat_bands)[j], 1, NA)
    weighted_band <- (rast_mean_585_corrected[[i]] * area_weights_normalized) * band_mask
    zonal_weighted[i, j] <- global(weighted_band, "sum", na.rm = TRUE)[1,1]
  }
}




# ============================================================
# Direct comparison: SSP585 bias-corrected vs. not corrected
# Standalone script, pulled out of the main ALD pipeline so it
# can be run on its own.
#
# Two ways to get the comparison, pick whichever fits what you
# already have on disk:
#   A) From the raw rasters (recomputes everything for ssp585)
#   B) From the CSVs you already exported from the full pipeline
# ============================================================

library(terra)
library(ggplot2)
library(dplyr)
library(tidyr)

ext_60N <- ext(-179.95, 179.95, 60, 90)

# ------------------------------------------------------------
# A) FROM RASTERS
# ------------------------------------------------------------

# Present-day satellite reference
ALD_ESA_REF_mean <- rast("ESA_ALD_present_day.nc", lyr = 1)
ALD_ESA_REF_std  <- rast("ESA_ALD_present_day.nc", lyr = 2)
ALD_ESA_REF_mean <- crop(ALD_ESA_REF_mean, ext_60N)
ALD_ESA_REF_std  <- crop(ALD_ESA_REF_std,  ext_60N)

# SSP585 raw
rast_mean_585 <- rast("60N/mean_ALD_ssp585_60deg.nc")
rast_std_585  <- rast("60N/std_ALD_ssp585_60deg.nc")
plot(rast_mean_585[[2]])
bias_correct_ssp <- function(mean_rast, std_rast, esa_mean, esa_std,
                             start_year = 2000, end_year = 2020) {
  rast_time <- time(mean_rast)
  present_idx <- which(rast_time >= start_year & rast_time <= end_year)
  
  model_present_mean <- mean(mean_rast[[present_idx]], na.rm = TRUE)
  bias_mean <- model_present_mean - esa_mean
  
  model_present_std <- mean(std_rast[[present_idx]], na.rm = TRUE)
  bias_std <- model_present_std - esa_std
  
  list(
    mean_corrected = mean_rast - bias_mean,
    std_corrected  = std_rast - bias_std,
    bias_mean = bias_mean,
    bias_std  = bias_std
  )
}

corr_585 <- bias_correct_ssp(rast_mean_585, rast_std_585,
                             ALD_ESA_REF_mean, ALD_ESA_REF_std)

# Area weights, from the corrected raster's grid
area_weights <- cellSize(corr_585$mean_corrected[[1]], mask = TRUE, unit = "km")
total_area <- global(area_weights, "sum", na.rm = TRUE)[1, 1]
area_weights_norm <- area_weights / total_area

calc_weighted_ts <- function(mean_rast, area_weights_norm) {
  n_years <- nlyr(mean_rast)
  weighted_mean <- numeric(n_years)
  for (i in 1:n_years) {
    weighted_grid <- mean_rast[[i]] * area_weights_norm
    weighted_mean[i] <- global(weighted_grid, "sum", na.rm = TRUE)[1, 1]
  }
  weighted_mean
}

years <- time(rast_mean_585)

comparison_df <- data.frame(
  year      = years,
  raw       = calc_weighted_ts(rast_mean_585,          area_weights_norm),
  corrected = calc_weighted_ts(corr_585$mean_corrected, area_weights_norm)
) %>%
  pivot_longer(cols = c(raw, corrected), names_to = "type", values_to = "ald")

p_raster <- ggplot(comparison_df, aes(x = year, y = ald, color = type)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +
  scale_color_manual(
    values = c("raw" = "gray40", "corrected" = "red"),
    labels = c("raw" = "Not corrected", "corrected" = "Bias corrected")
  ) +
  labs(
    title = "SSP585: Raw vs. ESA Bias-Corrected ALD",
    subtitle = "Area-weighted global mean, 60-90°N",
    x = "Year", y = "ALD [m]", color = ""
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

print(p_raster)













# ============================================================
# ALD hotspots and bias analysis (60-90 N)
#
# Answers two questions:
#   1. Where are the hotspots of active-layer deepening?
#      Delta ALD = mean(2080-2099) - mean(2000-2020), bias-corrected,
#      per SSP; summarised by region and latitude band, top 5 % cells.
#   2. How large is the model bias, and where is it largest?
#      Bias = CMIP6 multi-model mean ALD (2000-2020) - ESA CCI ALD;
#      area-weighted mean, mean absolute bias, RMSE, relative bias,
#      summarised by region and latitude band.
#
# The bias correction is done exactly as in the thawed-N script:
#   ALD_corr = ALD - bias, negative values set to 0
#
# All spatial means are AREA-WEIGHTED (grid cells shrink towards the pole).
#
# Outputs (in out_dir):
#   bias_summary.csv, bias_by_region.csv, bias_by_latband.csv
#   dALD_summary.csv, dALD_by_region.csv, dALD_by_latband.csv,
#   dALD_hotspot_cells_by_region.csv
#   map_ALD_bias.pdf, map_dALD.pdf
# ============================================================

library(terra)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(sf)
library(rnaturalearth)


# ------------------------------------------------------------
# 1. Settings
# ------------------------------------------------------------
ssps <- c("126", "245", "370", "585")
ssp_labels <- c("126" = "SSP1-2.6", "245" = "SSP2-4.5",
                "370" = "SSP3-7.0", "585" = "SSP5-8.5")

ald_fmt  <- "mean_ssp%s_corr_new.nc"      # multi-model mean ALD (same as the thaw script)
esa_file <- "ESA_ALD/ESA_ALD_present_day.nc"           # layer 1 = mean (same as the thaw script)

ext_60N <- ext(-179.95, 179.95, 60, 90)

present_period <- c(2000, 2020)
future_period  <- c(2080, 2099)

hotspot_quantile  <- 0.95    # top 5 % of cells = hotspots
profile_depth_max <- 5       # depth limit of the mineralisation model (m)

lat_breaks <- c(60, 66.5, 75, 90)
lat_labels <- c("60-66.5 N (subarctic)", "66.5-75 N", "75-90 N")

out_dir <- "ALD_hotspots_bias"
dir.create(out_dir, showWarnings = FALSE)


# ------------------------------------------------------------
# 2. Helpers
# ------------------------------------------------------------
layer_years <- function(r) {
  tt <- time(r)
  if (inherits(tt, "Date") || inherits(tt, "POSIXt")) as.integer(format(tt, "%Y")) else as.integer(tt)
}

wmean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

# regions: countries, with Russia and Canada split by longitude
region_of <- function(admin, lon) {
  case_when(
    admin == "United States of America"         ~ "Alaska",
    admin == "Canada" & lon < -95               ~ "Western Canada",
    admin == "Canada"                           ~ "Eastern Canada / Hudson Bay",
    admin == "Greenland"                        ~ "Greenland",
    admin %in% c("Norway", "Sweden", "Finland") ~ "Fennoscandia",
    admin == "Russia" & lon >= 0  & lon < 60    ~ "European Russia",
    admin == "Russia" & lon >= 60 & lon < 90    ~ "Western Siberia",
    admin == "Russia" & lon >= 90 & lon < 120   ~ "Central Siberia",
    admin == "Russia"                           ~ "Eastern Siberia",
    TRUE                                        ~ "Other / unassigned"
  )
}


# ------------------------------------------------------------
# 3. Reference data: ESA ALD, cell area, country of each cell
# ------------------------------------------------------------
r_template <- crop(rast(sprintf(ald_fmt, ssps[1])), ext_60N)

esa <- crop(rast(esa_file, lyr = 1), ext_60N)
if (!compareGeom(esa, r_template[[1]], stopOnError = FALSE)) {
  esa <- resample(esa, r_template[[1]], method = "bilinear")
}

cell_area <- cellSize(r_template[[1]], unit = "km")

countries <- ne_countries(scale = "medium", returnclass = "sf")
country_r <- rasterize(vect(countries), r_template[[1]], field = "admin", touches = TRUE)
names(country_r) <- "admin"


# ------------------------------------------------------------
# 4. Per SSP: bias, bias-corrected ALD, delta ALD -> one table of cells
# ------------------------------------------------------------
cell_tables  <- list()
dALD_rasters <- list()     # NEW
bias_rasters <- list()     # NEW

for (s in ssps) {
  
  cat("Processing SSP", s, "\n")
  r <- crop(rast(sprintf(ald_fmt, s)), ext_60N)
  yrs <- layer_years(r)
  
  idx_present <- which(yrs >= present_period[1] & yrs <= present_period[2])
  idx_future  <- which(yrs >= future_period[1]  & yrs <= future_period[2])
  
  model_present <- mean(r[[idx_present]], na.rm = TRUE)
  bias <- model_present - esa
  
  # bias correction as in the thawed-N script
  r_corr <- r - bias
  r_corr[r_corr < 0] <- 0
  
  present_corr <- mean(r_corr[[idx_present]], na.rm = TRUE)
  future_corr  <- mean(r_corr[[idx_future]],  na.rm = TRUE)
  dALD <- future_corr - present_corr
  dALD_rasters[[ssp_labels[[s]]]] <- dALD
  bias_rasters[[ssp_labels[[s]]]] <- bias
  stack_s <- c(esa, model_present, bias, present_corr, future_corr, dALD, cell_area)
  names(stack_s) <- c("esa", "model_present", "bias", "present_corr", "future_corr", "dALD", "area_km2")
  
  df <- as.data.frame(c(stack_s, country_r), xy = TRUE, na.rm = FALSE) %>%
    filter(is.finite(esa), is.finite(model_present), is.finite(dALD)) %>%
    mutate(SSP      = ssp_labels[[s]],
           lon      = x,
           lat      = y,
           region   = region_of(as.character(admin), lon),
           lat_band = cut(lat, breaks = lat_breaks, labels = lat_labels, include.lowest = TRUE))
  
  cell_tables[[s]] <- df
}

cells <- bind_rows(cell_tables) %>%
  mutate(SSP = factor(SSP, levels = unname(ssp_labels)))


# ------------------------------------------------------------
# 5. Bias statistics
#    The bias only differs between SSPs through 2015-2020, so it is
#    reported per SSP but should be nearly identical.
# ------------------------------------------------------------
bias_stats <- function(d) {
  d %>%
    summarise(
      area_Mkm2          = sum(area_km2) / 1e6,
      mean_ESA_m         = wmean(esa, area_km2),
      mean_model_m       = wmean(model_present, area_km2),
      mean_bias_m        = wmean(bias, area_km2),
      mean_abs_bias_m    = wmean(abs(bias), area_km2),
      rmse_m             = sqrt(wmean(bias^2, area_km2)),
      relative_bias_pct  = 100 * mean_bias_m / mean_ESA_m,
      share_overest_pct  = 100 * sum(area_km2[bias > 0]) / sum(area_km2),
      p05_bias_m         = quantile(bias, 0.05),
      p95_bias_m         = quantile(bias, 0.95),
      .groups = "drop"
    )
}

bias_summary     <- cells %>% group_by(SSP) %>% bias_stats()
bias_by_region   <- cells %>% group_by(SSP, region) %>% bias_stats() %>% arrange(SSP, desc(mean_abs_bias_m))
bias_by_latband  <- cells %>% group_by(SSP, lat_band) %>% bias_stats()

# where are the largest biases? share of the top 5 % |bias| cells per region
big_bias <- cells %>%
  group_by(SSP) %>%
  mutate(top = abs(bias) >= quantile(abs(bias), hotspot_quantile)) %>%
  filter(top) %>%
  group_by(SSP, region) %>%
  summarise(area_Mkm2 = sum(area_km2) / 1e6,
            mean_bias_m = wmean(bias, area_km2), .groups = "drop") %>%
  group_by(SSP) %>%
  mutate(share_of_top_pct = 100 * area_Mkm2 / sum(area_Mkm2)) %>%
  arrange(SSP, desc(share_of_top_pct))

write_csv(bias_summary,    file.path(out_dir, "bias_summary.csv"))
write_csv(bias_by_region,  file.path(out_dir, "bias_by_region.csv"))
write_csv(bias_by_latband, file.path(out_dir, "bias_by_latband.csv"))
write_csv(big_bias,        file.path(out_dir, "bias_largest_cells_by_region.csv"))

cat("\n=== Bias (model 2000-2020 minus ESA), area-weighted ===\n")
print(bias_summary, width = Inf)
cat("\n=== Bias by region (sorted by mean absolute bias) ===\n")
print(bias_by_region %>% filter(SSP == "SSP5-8.5"), width = Inf)
cat("\n=== Regions holding the largest |bias| cells (top 5 %) ===\n")
print(big_bias %>% filter(SSP == "SSP5-8.5"), width = Inf)


# ------------------------------------------------------------
# 6. Delta ALD hotspots
# ------------------------------------------------------------
dALD_stats <- function(d) {
  d %>%
    summarise(
      area_Mkm2              = sum(area_km2) / 1e6,
      present_ALD_m          = wmean(present_corr, area_km2),
      future_ALD_m           = wmean(future_corr, area_km2),
      mean_dALD_m            = wmean(dALD, area_km2),
      p95_dALD_m             = quantile(dALD, 0.95),
      max_dALD_m             = max(dALD),
      share_dALD_gt_1m_pct   = 100 * sum(area_km2[dALD > 1]) / sum(area_km2),
      share_future_ge_5m_pct = 100 * sum(area_km2[future_corr >= profile_depth_max]) / sum(area_km2),
      .groups = "drop"
    )
}

dALD_summary    <- cells %>% group_by(SSP) %>% dALD_stats()
dALD_by_region  <- cells %>% group_by(SSP, region) %>% dALD_stats() %>% arrange(SSP, desc(mean_dALD_m))
dALD_by_latband <- cells %>% group_by(SSP, lat_band) %>% dALD_stats()

# hotspot cells: top 5 % of delta ALD per SSP, and where they are
hotspots <- cells %>%
  group_by(SSP) %>%
  mutate(threshold_m = quantile(dALD, hotspot_quantile),
         hotspot = dALD >= threshold_m) %>%
  ungroup()

hotspot_by_region <- hotspots %>%
  filter(hotspot) %>%
  group_by(SSP, region) %>%
  summarise(threshold_m     = first(threshold_m),
            area_Mkm2       = sum(area_km2) / 1e6,
            mean_dALD_m     = wmean(dALD, area_km2),
            mean_lat        = wmean(lat, area_km2),
            .groups = "drop") %>%
  group_by(SSP) %>%
  mutate(share_of_hotspots_pct = 100 * area_Mkm2 / sum(area_Mkm2)) %>%
  arrange(SSP, desc(share_of_hotspots_pct))

write_csv(dALD_summary,      file.path(out_dir, "dALD_summary.csv"))
write_csv(dALD_by_region,    file.path(out_dir, "dALD_by_region.csv"))
write_csv(dALD_by_latband,   file.path(out_dir, "dALD_by_latband.csv"))
write_csv(hotspot_by_region, file.path(out_dir, "dALD_hotspot_cells_by_region.csv"))

cat("\n=== Delta ALD (2080-2099 minus 2000-2020), area-weighted ===\n")
print(dALD_summary, width = Inf)
cat("\n=== Delta ALD by region (sorted by mean increase) ===\n")
print(dALD_by_region %>% filter(SSP == "SSP5-8.5"), width = Inf)
cat("\n=== Where the hotspot cells (top 5 %) are ===\n")
print(hotspot_by_region, n = Inf, width = Inf)


# ------------------------------------------------------------
# 7. Maps (polar stereographic, geom_raster)
# ------------------------------------------------------------
polar_crs <- "+proj=stere +lat_0=90 +lat_ts=71 +lon_0=0 +datum=WGS84 +units=m"
map_res_m <- 25000      # resolution of the projected maps (m); smaller = finer, slower

coast <- st_transform(ne_countries(scale = "medium", returnclass = "sf"), crs = polar_crs)

# project a raster to polar stereographic and convert it to a data frame
raster_to_polar_df <- function(r, value_name, label) {
  rp <- project(r, polar_crs, res = map_res_m, method = "bilinear")
  df <- as.data.frame(rp, xy = TRUE, na.rm = TRUE)
  names(df)[3] <- value_name
  df$SSP <- label
  df
}

map_theme <- theme_minimal(base_size = 6) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_line(colour = "grey90", linewidth = 0.2),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom", legend.key.width = unit(1.5, "cm"),
        legend.key.height = unit(0.25, "cm"))

# --- bias map (SSP5-8.5; nearly identical for all SSPs) ---
bias_df <- raster_to_polar_df(bias_rasters[["SSP5-8.5"]], "bias", "SSP5-8.5")
lim_x <- range(bias_df$x)
lim_y <- range(bias_df$y)
max_b <- quantile(abs(bias_df$bias), 0.99)       # clip colour scale at the 99th percentile

b_lim    <- round(max_b, 1)                                   # rounded end value
b_breaks <- c(-b_lim, pretty(c(-b_lim, b_lim), n = 4), b_lim)
b_breaks <- sort(unique(b_breaks[abs(b_breaks) <= b_lim]))
b_breaks <- b_breaks[!(abs(b_breaks) < b_lim & abs(b_breaks) > b_lim - 0.4)]  # drop breaks too close to the ends
b_labels <- ifelse(b_breaks == -b_lim, paste0("\u2264 ", -b_lim),
                   ifelse(b_breaks ==  b_lim, paste0("\u2265 ", b_lim),
                          as.character(b_breaks)))

p_bias <- ggplot() +
  geom_raster(data = bias_df, aes(x, y, fill = pmax(pmin(bias, max_b), -max_b))) +
  geom_sf(data = coast, fill = NA, colour = "grey20", linewidth = 0.05) +
  coord_sf(xlim = lim_x, ylim = lim_y, crs = polar_crs, expand = FALSE) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                       midpoint = 0, limits = c(-b_lim, b_lim),
                       breaks = b_breaks, labels = b_labels,
                       name = "ALD bias, model - ESA (m)",
                       guide = guide_colourbar(title.position = "top",
                                               title.hjust = 0.5,
                                               barwidth  = unit(4, "cm"),
                                               barheight = unit(0.2, "cm"))) +
  labs(x = NULL, y = NULL) +
  map_theme

ggsave(file.path(out_dir, "map_ALD_bias.pdf"), p_bias, width = 8.5, height = 6, units = "cm")

# --- delta ALD map, all SSPs ---
dALD_df <- bind_rows(lapply(names(dALD_rasters), function(lab) {
  raster_to_polar_df(dALD_rasters[[lab]], "dALD", lab)
})) %>%
  mutate(SSP = factor(SSP, levels = unname(ssp_labels)))

max_d <- quantile(dALD_df$dALD, 0.99)

p_dALD <- ggplot() +
  geom_raster(data = dALD_df, aes(x, y, fill = pmin(dALD, max_d))) +
  geom_sf(data = coast, fill = NA, colour = "grey20", linewidth = 0.05) +
  coord_sf(xlim = lim_x, ylim = lim_y, crs = polar_crs, expand = FALSE) +
  scale_fill_viridis_c(option = "magma", direction = -1,
                       name = "Change in ALD, 2080-2099 vs 2000-2020 (m)") +
  facet_wrap(~SSP, ncol = 2) +
  labs(x = NULL, y = NULL) +
  map_theme

ggsave(file.path(out_dir, "map_dALD.pdf"), p_dALD, width = 15, height = 12, units = "cm")

print(p_bias)
print(p_dALD)

cat("\nDone. Tables and maps in", out_dir, "\n")