library(terra)
library(here)
library(ggplot2)
library(dplyr)
library(tidyr)
getwd()
# Load NetCDF thawed nitrogen rasters
rast_mean_585 <- rast("mean_ssp585_corr_new.nc")
rast_std_585<- rast("std_ssp585_corr_new.nc")
plot(rast_std_585[[200]])

time(rast_mean_585)
# baseline and future periods
ald_2000_2020 <- rast_mean_585[[time(rast_mean_585) >= as.Date("2000-07-01") &
                                  time(rast_mean_585) <= as.Date("2020-07-01")]]

ald_2080_2099 <- rast_mean_585[[time(rast_mean_585) >= as.Date("2080-07-01") &
                                  time(rast_mean_585) <= as.Date("2099-07-01")]]

# mean ALD for each period
mean_2000_2020 <- mean(ald_2000_2020, na.rm = TRUE)
mean_2080_2099 <- mean(ald_2080_2099, na.rm = TRUE)

# delta ALD
delta_ald <- mean_2080_2099 - mean_2000_2020
names(delta_ald) <- "delta_ALD"
delta_df <- as.data.frame(delta_ald, xy = TRUE, na.rm = TRUE)

ggplot() +
  geom_raster(data = delta_df, aes(x = x, y = y, fill = delta_ALD)) +
  geom_sf(data = coastlines_sf, fill = NA, color = "black", linewidth = 0.2) +
  coord_sf(xlim = c(-180, 180), ylim = c(30, 90), expand = FALSE) +
  scale_fill_viridis_c(name = expression(Delta*"ALD [m]")) +
  labs(
    title = expression(Delta*"ALD for SSP5-8.5"),
    subtitle = "Mean 2080-2099 minus mean 2000-2020",
    x = "",
    y = ""
  ) +
  theme_minimal()
# mean ALD all the cell areas (total area of the grid)
mean_ssp585 <- global(rast_mean_585, fun = "mean", na.rm = TRUE)
std_ssp585<- global(rast_std_585, fun = "sd", na.rm = TRUE)

plot(rast_mean_585[[1]])

# Load NetCDF thawed nitrogen rasters
rast_mean_370 <- rast("mean_ssp370_corr_new.nc")
rast_std_370<- rast("std_ssp370_corr_new.nc")
# mean ALD all the cell areas (total area of the grid)
mean_ssp370 <- global(rast_mean_370, fun = "mean", na.rm = TRUE)
std_ssp370<- global(rast_std_370, fun = "std", na.rm = TRUE)


# Load NetCDF thawed nitrogen rasters
rast_mean_245 <- rast("mean_ssp245_corr_new.nc")
rast_std_245 <- rast("std_ssp245_corr_new.nc")
# mean ALD all the cell areas (total area of the grid)
mean_ssp245 <- global(rast_mean_245, fun = "mean", na.rm = TRUE)
std_ssp245<- global(rast_std_245, fun = "std", na.rm = TRUE)
plot(mean_ssp245)

# Load NetCDF thawed nitrogen rasters
rast_mean_126 <- rast("mean_ssp126_corr_new.nc")
rast_std_126 <- rast("std_ssp126_corr_new.nc")
# mean ALD all the cell areas (total area of the grid)
mean_ssp126 <- global(rast_mean_126, fun = "mean", na.rm = TRUE)
std_ssp126<- global(rast_std_126, fun = "std", na.rm = TRUE)
plot(mean_ssp245)


# compare ALD 1880-1900, 2000 - 2020, 2080-2100; 

library(terra)
year_idx <- function(years, y0, y1) which(years >= y0 & years <= y1)
yrs <- time(rast_mean_585)
yrs <- as.numeric(format(yrs, "%Y"))

area_weighted_mean <- function(r){
  a <- cellSize(r, unit = "m")
  a <- mask(a, r)  # same footprint
  num <- global(r * a, "sum", na.rm = TRUE)[[1]]
  den <- global(a, "sum", na.rm = TRUE)[[1]]
  num / den
}

period_stats_aw <- function(r, years, y0, y1){
  idx <- year_idx(years, y0, y1)
  r_mean_period <- mean(r[[idx]], na.rm = TRUE)
  m <- area_weighted_mean(r_mean_period)
  list(mean = m, raster = r_mean_period)
}

aw_1880_1900 <- period_stats_aw(rast_mean_585, yrs, 1880, 1900)
aw_2000_2020 <- period_stats_aw(rast_mean_585, yrs, 2000, 2020)
aw_2080_2100 <- period_stats_aw(rast_mean_585, yrs, 2080, 2100)

aw_1880_1900$mean; aw_2000_2020$mean; aw_2080_2100$mean

#





df<-data.frame(mean_ssp585, mean_ssp370, mean_ssp245, mean_ssp126, std_ssp585, std_ssp370, std_ssp245, std_ssp126)
df$Year<-rep(1850:2099)
names(df)<-c("mean_585", "mean_370", "mean_245", "mean_126", "std_585", "std_370", "std_245","std_126", "Year")



write.csv(df, "final_results/ALD_final.csv")



df<-read.csv("ALD_final.csv")

df_period <- df %>%
  mutate(
    Period = case_when(
      Year >= 2000 & Year <= 2020 ~ "baseline",
      Year >= 2080 & Year <= 2100 ~ "future",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Period))
df_mean <- df_period %>%
  group_by(Period) %>%
  summarise(
    mean_585 = mean(mean_585, na.rm = TRUE),
    mean_370 = mean(mean_370, na.rm = TRUE),
    mean_245 = mean(mean_245, na.rm = TRUE),
    mean_126 = mean(mean_126, na.rm = TRUE)
  )
df_change <- df_mean %>%
  pivot_wider(names_from = Period, values_from = starts_with("mean")) %>%
  mutate(
    delta_585 = mean_585_future - mean_585_baseline,
    delta_370 = mean_370_future - mean_370_baseline,
    delta_245 = mean_245_future - mean_245_baseline,
    delta_126 = mean_126_future - mean_126_baseline
  )

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
# Compute 20-year rolling mean and standard deviation for each SSP
long_data <- long_data %>%
  group_by(SSP) %>%
  mutate(
    Rolling_Mean_ALD = rollmean(Mean_ALD, k = 20, fill = NA, align = "center"),
    Rolling_STD_ALD  = rollapply(STD_ALD, width = 20, FUN = mean, fill = NA, align = "center")
  )
long_data$SSP <- factor(
  long_data$SSP,
  levels = c("SSP126", "SSP245", "SSP370", "SSP585"),
  labels = c("SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5")
)

# Plot
ggplot(long_data, aes(x = Year, y = Mean_ALD, group = SSP)) +
  
  # Before 2015: Black lines & grey ribbons
  geom_line(data = filter(long_data, Period == "Before 2015"), color = "black", linewidth = 0.8) +
  #geom_ribbon(data = filter(long_data, Period == "Before 2015"),
  #            aes(ymin = Mean_ALD - STD_ALD, ymax = Mean_ALD + STD_ALD),
  #            fill = "grey", alpha = 0.2) +
  
  # After 2015: Colored lines & ribbons
  geom_line(data = filter(long_data, Period == "After 2015"),
            aes(color = SSP), linewidth = 0.8) +
  #geom_ribbon(data = filter(long_data, Period == "After 2015"),
  #            aes(ymin = Mean_ALD - STD_ALD, ymax = Mean_ALD + STD_ALD, fill = SSP),
  #            alpha = 0.2) +
  
  # Rolling Mean before 2015 (Black line)
  geom_line(data = filter(long_data, Year < 2015),
            aes(y = Rolling_Mean_ALD), color = "black", linewidth = 1.2, linetype = "solid") +
  
  # Rolling Mean after 2015 (Colored lines)
  geom_line(data = filter(long_data, Year >= 2015),
            aes(y = Rolling_Mean_ALD, color = SSP), linewidth = 1.2, linetype = "solid") +
  
  # Rolling STD before 2015 (Grey ribbon)
  geom_ribbon(data = filter(long_data, Year < 2015),
              aes(ymin = Rolling_Mean_ALD - Rolling_STD_ALD, 
                  ymax = Rolling_Mean_ALD + Rolling_STD_ALD),
              fill = "grey", alpha = 0.15) +
  
  ylim(-0.75,3.0)+
  # Rolling STD after 2015 (Colored ribbons)
  geom_ribbon(data = filter(long_data, Year >= 2015),
  aes(ymin = Rolling_Mean_ALD - Rolling_STD_ALD, 
  ymax = Rolling_Mean_ALD + Rolling_STD_ALD, 
  fill = SSP),
  alpha = 0.15) +
  
  # Vertical line for 2015
  geom_vline(xintercept = 2015, linetype = "dashed", color = "black") +
  
  # Labels and title
  labs(x = "Year", y = "Active layer depth [m]", title = "", 
       color = "SSP Scenario", fill = "SSP Scenario") +
  
  # Color scales
  scale_color_manual(values = ssp_colors) +
  scale_fill_manual(values = ssp_colors) +
  
  # Theme settings
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14)
  )



