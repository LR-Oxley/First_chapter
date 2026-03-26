# R script for plotting graphs

###### ------------------ figure 1 ---------------------- ######

# plot ALD increase
df<-read.csv("ALD_final.csv")

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
plot_ALD<-ggplot(long_data, aes(x = Year, y = Mean_ALD, group = SSP)) +
  
  # Before 2015: Black lines & grey ribbons
  geom_line(data = filter(long_data, Period == "Before 2015"), color = "black", linewidth = 0.5) +
  #geom_ribbon(data = filter(long_data, Period == "Before 2015"),
  #            aes(ymin = Mean_ALD - STD_ALD, ymax = Mean_ALD + STD_ALD),
  #            fill = "grey", alpha = 0.2) +
  
  # After 2015: Colored lines & ribbons
  geom_line(data = filter(long_data, Period == "After 2015"),
            aes(color = SSP), linewidth = 0.5) +
  #geom_ribbon(data = filter(long_data, Period == "After 2015"),
  #            aes(ymin = Mean_ALD - STD_ALD, ymax = Mean_ALD + STD_ALD, fill = SSP),
  #            alpha = 0.2) +
  
  # Rolling Mean before 2015 (Black line)
  geom_line(data = filter(long_data, Year < 2015),
            aes(y = Rolling_Mean_ALD), color = "black", linewidth = 0.5, linetype = "solid") +
  
  # Rolling Mean after 2015 (Colored lines)
  geom_line(data = filter(long_data, Year >= 2015),
            aes(y = Rolling_Mean_ALD, color = SSP), linewidth = 0.5, linetype = "solid") +
  
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
  
  scale_color_manual(
    values = ssp_colors,
    guide = guide_legend(
      keyheight = unit(0.35, "cm"),
      keywidth  = unit(0.35, "cm"),
      override.aes = list(linewidth = 0.7)
    )
  ) +
  
  scale_fill_manual(
    values = ssp_colors,
    guide = guide_legend(
      keyheight = unit(0.35, "cm"),
      keywidth  = unit(0.35, "cm"),
      override.aes = list(alpha = 0.4)
    )
  ) +
  
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.spacing.y = unit(0.08, "cm"),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8)
  )




#### plotting total thawed N 

#plot total N release in Pg
library(terra)
total_thawed_585<-rast("n_thaw_no_lim/arctic_total_thawed_585_no_lim.tif")
total_thawed_370<-rast("n_thaw_no_lim/arctic_total_thawed_370_no_lim.tif")
total_thawed_245<-rast("n_thaw_no_lim/arctic_total_thawed_245_no_lim.tif")
total_thawed_126<-rast("n_thaw_no_lim/arctic_total_thawed_126_no_lim.tif")
total_thawed_585_std<-rast("n_thaw_no_lim/arctic_total_thawed_585_no_lim_std.tif")
total_thawed_370_std<-rast("n_thaw_no_lim/arctic_total_thawed_370_no_lim_std.tif")
total_thawed_245_std<-rast("n_thaw_no_lim/arctic_total_thawed_245_no_lim_std.tif")
total_thawed_126_std<-rast("n_thaw_no_lim/arctic_total_thawed_126_no_lim_std.tif")

library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(zoo)

# Helper function to extract time series for total Pg 
extract_ts <- function(mean_rast, sd_rast, ssp, years, area_rast) {
  total_mean <- global(mean_rast * area_rast, sum, na.rm = TRUE)[, 1]
  total_sd <- global(sd_rast * area_rast, sum, na.rm = TRUE)[, 1]
  total_mean <- total_mean / 1e12  # kg → Pg
  total_sd <- total_sd / 1e12
  
  tibble(
    Year = years,
    Mean = total_mean,
    SD = total_sd,
    SSP = ssp
  )
}

extract_total_pg <- function(mean_rast, sd_rast, ssp, years, area_rast) {
  
  # weights (m2), NA wherever mean_rast is NA (ocean)
  w <- terra::mask(area_rast, mean_rast)
  
  total_mean_kg <- terra::global(mean_rast * w, "sum", na.rm = TRUE)[[1]]
  total_sd_kg   <- terra::global(sd_rast   * w, "sum", na.rm = TRUE)[[1]]
  
  tibble::tibble(
    Year = years,
    Mean = total_mean_kg / 1e12,  # kg -> Pg
    SD   = total_sd_kg   / 1e12,
    SSP  = ssp
  )
}

# # in kg / m2
# extract_ts <- function(mean_rast, sd_rast, ssp, years, area_rast) {
#   
#   # total area (m2)
#   total_area <- global(area_rast, sum, na.rm = TRUE)[, 1]
#   
#   # area-weighted mean (kg / m2)
#   mean_weighted <- global(mean_rast * area_rast, sum, na.rm = TRUE)[, 1] / total_area
#   
#   # area-weighted SD (still in kg / m2)
#   sd_weighted <- global(sd_rast * area_rast, sum, na.rm = TRUE)[, 1] / total_area
#   
#   tibble(
#     Year = years,
#     Mean = mean_weighted,   # kg / m2
#     SD   = sd_weighted,     # kg / m2
#     SSP  = ssp
#   )
# }

# Define years (adjust if needed)
years <- 1850:(1850 + nlyr(total_thawed_370) - 1)

# Calculate cell area
cell_area <- cellSize(total_thawed_585, unit = "m")

# Extract time series for each scenario
df_585 <- extract_total_pg(total_thawed_585, total_thawed_585_std, "SSP585", years, cell_area)
rm(total_thawed_585,total_thawed_585_std)
df_370 <- extract_total_pg(total_thawed_370, total_thawed_370_std, "SSP370", years, cell_area)
rm(total_thawed_370, total_thawed_370_std)
df_245 <- extract_total_pg(total_thawed_245, total_thawed_245_std, "SSP245", years, cell_area)
rm(total_thawed_245, total_thawed_245_std)

df_126 <- extract_total_pg(total_thawed_126, total_thawed_126_std, "SSP126", years, cell_area)

# Combine all data
df_all <- bind_rows(df_585, df_370, df_245, df_126)
tail(df_585)

df_periods <- df_all %>%
  mutate(
    Period = case_when(
      Year >= 2000 & Year <= 2020 ~ "baseline",
      Year >= 2080 & Year <= 2099 ~ "future",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Period)) %>%
  group_by(SSP, Period) %>%
  summarise(
    Mean_ALD = mean(Mean, na.rm = TRUE),
    SD_ALD = mean(SD, na.rm = TRUE),
    .groups = "drop"
  )


df_change <- df_periods %>%
  dplyr::select(SSP, Period, Mean_ALD) %>%
  pivot_wider(names_from = Period, values_from = Mean_ALD) %>%
  dplyr::mutate(
    delta_ALD = future - baseline
  )

df_change





# Function to compute anomalies relative to 1880-1900
compute_anomaly <- function(df) {
  # Calculate baseline (1880-1900) - changed from 1860-1900 to match your description
  ref <- df %>%
    filter(Year >= 2000, Year <= 2020)  # or 1860-1900 if you prefer
  
  # Calculate reference period statistics
  ref_mean <- mean(ref$Mean, na.rm = TRUE)
  ref_std <- mean(ref$SD, na.rm = TRUE)
  
  # Compute anomalies
  df %>%
    mutate(
      Mean_anom = Mean - ref_mean,
      # Subtract reference period standard deviation
      SD_anom = SD - ref_std
    )
}

# Apply anomaly calculation
df_all_anom <- df_all %>%
  group_by(SSP) %>%
  group_modify(~ compute_anomaly(.)) %>%
  ungroup()

# Add period column (before/after 2015)
df_all_anom <- df_all_anom %>%
  mutate(Period = ifelse(Year < 2015, "Before 2015", "After 2015"))

df_all_anom <- df_all_anom %>%
  arrange(SSP, Year) %>%
  group_by(SSP) %>%
  mutate(
    Rolling_Mean_anom = rollapply(
      Mean_anom, width = 20, FUN = mean,
      align = "center", fill = NA, na.rm = TRUE
    ),
    Rolling_SD_anom = rollapply(
      SD_anom, width = 20, FUN = mean,
      align = "center", fill = NA, na.rm = TRUE
    )
  ) %>%
  ungroup() 


df_all_anom <- df_all_anom %>%
  mutate(
    SSP = recode(
      SSP,
      "SSP126" = "SSP1-2.6",
      "SSP245" = "SSP2-4.5",
      "SSP370" = "SSP3-7.0",
      "SSP585" = "SSP5-8.5"
    )
  )

ssp_colors <- c(
  "SSP3-7.0" = "#D73027",
  "SSP5-8.5" = "#7B3294",
  "SSP2-4.5" = "orange",
  "SSP1-2.6" = "blue"
)


# Create the plot
plot_total_thawed<-ggplot(df_all_anom, aes(x = Year, y = Mean_anom, group = SSP)) +
  
  # Before 2015: Black lines & grey ribbons
  geom_line(data = filter(df_all_anom, Period == "Before 2015"), 
            color = "black", linewidth = 0.5) +
  #geom_ribbon(data = filter(df_all_anom, Period == "Before 2015"),
  #           aes(ymin = Mean_anom - SD_anom, ymax = Mean_anom + SD_anom),
  #          fill = "grey", alpha = 0.2) +
  
  # After 2015: Colored lines & ribbons
  geom_line(data = filter(df_all_anom, Period == "After 2015"),
            aes(color = SSP), linewidth = 0.5) +
  #geom_ribbon(data = filter(df_all_anom, Period == "After 2015"),
  #           aes(ymin = Mean_anom - SD_anom, ymax = Mean_anom + SD_anom, fill = SSP),
  #          alpha = 0.15) +
  
  # Rolling Mean before 2015 (Black line)
  geom_line(data = filter(df_all_anom, Year < 2015),
            aes(y = Rolling_Mean_anom), color = "black", 
            linewidth = 0.5, linetype = "solid") +
  
  # Rolling Mean after 2015 (Colored lines)
  geom_line(data = filter(df_all_anom, Year >= 2015),
            aes(y = Rolling_Mean_anom, color = SSP), 
            linewidth = 0.5, linetype = "solid") +
  
  # Rolling STD before 2015 (Grey ribbon)
  geom_ribbon(data = filter(df_all_anom, Year < 2015),
              aes(ymin = Rolling_Mean_anom - Rolling_SD_anom, 
                  ymax = Rolling_Mean_anom + Rolling_SD_anom),
              fill = "grey", alpha = 0.15) +
  
  # Rolling STD after 2015 (Colored ribbons) - COMMENTED OUT as in your original
  geom_ribbon(data = filter(df_all_anom, Year >= 2015),
              aes(ymin = Rolling_Mean_anom - Rolling_SD_anom, 
                  ymax = Rolling_Mean_anom + Rolling_SD_anom, 
                  fill = SSP),
              alpha = 0.15) +
  
  # Vertical line for 2015
  geom_vline(xintercept = 2015, linetype = "dashed", color = "black") +
  
  # Labels and title - UPDATED for your nitrogen data
  labs(x = "Year", 
       y = "Total thawed N [Pg]", 
       color = "SSP Scenario", 
       fill = "SSP Scenario") +
  
  scale_color_manual(
    values = ssp_colors,
    guide = guide_legend(
      keyheight = unit(0.35, "cm"),
      keywidth  = unit(0.35, "cm"),
      override.aes = list(linewidth = 0.7)
    )
  ) +
  
  scale_fill_manual(
    values = ssp_colors,
    guide = guide_legend(
      keyheight = unit(0.35, "cm"),
      keywidth  = unit(0.35, "cm"),
      override.aes = list(alpha = 0.4)
    )
  ) +
  
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.spacing.y = unit(0.08, "cm"),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8)
  )

library(dplyr)

####

# to plot ALD and total thawed next to each other: 
library(patchwork)


combined_plot<-plot_ALD / plot_total_thawed + plot_annotation(tag_levels = "a", tag_suffix = ") ")
ggsave("Fig1.png", combined_plot,
       width = 12,
       height = 15,
       units = "cm", 
       #scale = 1.4,
       dpi = 600,  #70
       #device = cairo_pdf
)
  
######### ------------------------ figure 2 ---------------------- ##########
  
#### plotting min N 

#plot min N release in Pg
library(here)
here
setwd("coding_first_chapter")
getwd()
library(terra)
mineralised_585<-rast("1.9%_min/arctic_bioavailable_N_pool_585_no_lim_mean.tif")[[1:250]]
mineralised_370<-rast("1.9%_min/arctic_bioavailable_N_pool_370_no_lim_mean.tif")[[1:250]]
mineralised_245<-rast("1.9%_min/arctic_bioavailable_N_pool_245_no_lim_mean.tif")[[1:250]]
mineralised_126<-rast("1.9%_min/arctic_bioavailable_N_pool_126_no_lim_mean.tif")[[1:250]]
mineralised_585_std<-rast("1.9%_min/arctic_bioavailable_N_pool_585_no_lim_std.tif")[[1:250]]
mineralised_370_std<-rast("1.9%_min/arctic_bioavailable_N_pool_370_no_lim_std.tif")[[1:250]]
mineralised_245_std<-rast("1.9%_min/arctic_bioavailable_N_pool_245_no_lim_std.tif")[[1:250]]
mineralised_126_std<-rast("1.9%_min/arctic_bioavailable_N_pool_126_no_lim_std.tif")[[1:250]]

library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(zoo)
# Helper function to extract time series for total Pg 
extract_ts <- function(mean_rast,sd_rast, ssp, years, area_rast) {
  total_mean <- global(mean_rast * area_rast, sum, na.rm = TRUE)[, 1]
  total_sd <- global(sd_rast * area_rast, sum, na.rm = TRUE)[, 1]
  total_mean <- total_mean / 1e12  # kg → Pg
  total_sd <- total_sd / 1e12
  
  tibble(
    Year = years,
    Mean = total_mean,
    SD = total_sd,
    SSP = ssp
  )
}

# # in kg / m2
# extract_ts <- function(mean_rast, sd_rast, ssp, years, area_rast) {
#   
#   # total area (m2)
#   total_area <- global(area_rast, sum, na.rm = TRUE)[, 1]
#   
#   # area-weighted mean (kg / m2)
#   mean_weighted <- global(mean_rast * area_rast, sum, na.rm = TRUE)[, 1] / total_area
#   
#   # area-weighted SD (still in kg / m2)
#   sd_weighted <- global(sd_rast * area_rast, sum, na.rm = TRUE)[, 1] / total_area
#   
#   tibble(
#     Year = years,
#     Mean = mean_weighted,   # kg / m2
#     SD   = sd_weighted,     # kg / m2
#     SSP  = ssp
#   )
# }

# Define years (adjust if needed)
years <- 1850:(1850 + nlyr(mineralised_245) - 1)

# Calculate cell area
cell_area <- cellSize(mineralised_245, unit = "m")



# Extract time series for each scenario
df_585 <- extract_ts(mineralised_585, mineralised_585_std, "SSP585", years, cell_area)
df_370 <- extract_ts(mineralised_370, mineralised_370_std, "SSP370", years, cell_area)
df_245 <- extract_ts(mineralised_245, mineralised_245_std, "SSP245", years, cell_area)
df_126 <- extract_ts(mineralised_126, mineralised_126_std, "SSP126", years, cell_area)
df_all<-df_245
df_all <- bind_rows(df_585, df_370, df_245, df_126)
tail(df_all)

df_periods <- df_all %>%
  mutate(
    Period = case_when(
      Year >= 2000 & Year <= 2020 ~ "baseline",
      Year >= 2080 & Year <= 2099 ~ "future",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Period)) %>%
  group_by(SSP, Period) %>%
  summarise(
    Mean_ALD = mean(Mean, na.rm = TRUE),
    SD_ALD = mean(SD, na.rm = TRUE),
    .groups = "drop"
  )


df_change <- df_periods %>%
  dplyr::select(SSP, Period, Mean_ALD) %>%
  pivot_wider(names_from = Period, values_from = Mean_ALD) %>%
  dplyr::mutate(
    delta_ALD = future - baseline
  )

df_change




##############
# convert in g/m2
# df_all <- df_all %>%
#   mutate(
#     Mean = Mean * 1000,
#     SD   = SD * 1000
#   )

# Function to compute anomalies relative to 1880-1900
compute_anomaly <- function(df_all) {
  # Calculate baseline (1880-1900) - changed from 1860-1900 to match your description
  ref <- df_all %>%
    filter(Year >= 2000, Year <= 2020)  # or 1860-1900 if you prefer
  
  # Calculate reference period statistics
  ref_mean <- mean(ref$Mean, na.rm = TRUE)
  ref_std <- mean(ref$SD, na.rm = TRUE)
  
  # Compute anomalies
  df_all %>%
    mutate(
      Mean_anom = Mean - ref_mean,
      # Subtract reference period standard deviation
      SD_anom = SD - ref_std
    )
}

# Apply anomaly calculation
df_all_anom <- df_all %>%
  group_by(SSP) %>%
  group_modify(~ compute_anomaly(.)) %>%
  ungroup()

# Add period column (before/after 2015)
df_all_anom <- df_all_anom %>%
  mutate(Period = ifelse(Year < 2015, "Before 2015", "After 2015"))

df_all_anom <- df_all_anom %>%
  arrange(SSP, Year) %>%
  group_by(SSP) %>%
  mutate(
    Rolling_Mean_anom = rollapply(
      Mean_anom, width = 20, FUN = mean,
      align = "center", fill = NA, na.rm = TRUE
    ),
    Rolling_SD_anom = rollapply(
      SD_anom, width = 20, FUN = mean,
      align = "center", fill = NA, na.rm = TRUE
    )
  ) %>%
  ungroup() 


df_all_anom <- df_all_anom %>%
  mutate(
    SSP = recode(
      SSP,
      "SSP126" = "SSP1-2.6",
      "SSP245" = "SSP2-4.5",
      "SSP370" = "SSP3-7.0",
      "SSP585" = "SSP5-8.5"
    )
  )

ssp_colors <- c(
  "SSP3-7.0" = "#D73027",
  "SSP5-8.5" = "#7B3294",
  "SSP2-4.5" = "orange",
  "SSP1-2.6" = "blue"
)


# Create the plot
plot_min<-ggplot(df_all_anom, aes(x = Year, y = Mean_anom, group = SSP)) +
  
  # Before 2015: Black lines & grey ribbons
  geom_line(data = filter(df_all_anom, Period == "Before 2015"), 
            color = "black", linewidth = 0.8) +
  #geom_ribbon(data = filter(df_all_anom, Period == "Before 2015"),
   #          aes(ymin = Mean_anom - SD_anom, ymax = Mean_anom + SD_anom),
    #        fill = "grey", alpha = 0.2) +
  
  # After 2015: Colored lines & ribbons
  geom_line(data = filter(df_all_anom, Period == "After 2015"),
            aes(color = SSP), linewidth = 0.8) +
  #geom_ribbon(data = filter(df_all_anom, Period == "After 2015"),
   #          aes(ymin = Mean_anom - SD_anom, ymax = Mean_anom + SD_anom, fill = SSP),
   #         alpha = 0.15) +
  
  # Rolling Mean before 2015 (Black line)
  geom_line(data = filter(df_all_anom, Year < 2015),
            aes(y = Rolling_Mean_anom), color = "black", 
            linewidth = 1.2, linetype = "solid") +
  
  # Rolling Mean after 2015 (Colored lines)
  geom_line(data = filter(df_all_anom, Year >= 2015),
            aes(y = Rolling_Mean_anom, color = SSP), 
            linewidth = 1.2, linetype = "solid") +
  
  # Rolling STD before 2015 (Grey ribbon)
  geom_ribbon(data = filter(df_all_anom, Year < 2015),
              aes(ymin = Rolling_Mean_anom - Rolling_SD_anom, 
                  ymax = Rolling_Mean_anom + Rolling_SD_anom),
              fill = "grey", alpha = 0.15) +
  
  # Rolling STD after 2015 (Colored ribbons) - COMMENTED OUT as in your original
  geom_ribbon(data = filter(df_all_anom, Year >= 2015),
              aes(ymin = Rolling_Mean_anom - Rolling_SD_anom, 
                 ymax = Rolling_Mean_anom + Rolling_SD_anom, 
                  fill = SSP),
            alpha = 0.15) +
  
  # Vertical line for 2015
  geom_vline(xintercept = 2015, linetype = "dashed", color = "black") +
  
  # Labels and title - UPDATED for your nitrogen data
  labs(x = "Year", 
       y = "Mineralised N [Pg]", 
       color = "SSP Scenario", 
       fill = "SSP Scenario") +
  scale_y_continuous(
    breaks = seq(0, 8, by = 2)
  )+
  # Color scales
  scale_color_manual(values = ssp_colors) +
  scale_fill_manual(values = ssp_colors) +
  
  # Theme settings
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8)
  )



# make legend smaller: 

plot_min <- ggplot(df_all_anom, aes(x = Year, y = Mean_anom, group = SSP)) +
  
  # Before 2015: black mean line
  geom_line(
    data = filter(df_all_anom, Period == "Before 2015"),
    color = "black",
    linewidth = 0.5
  ) +
  
  # After 2015: colored mean lines
  geom_line(
    data = filter(df_all_anom, Period == "After 2015"),
    aes(color = SSP),
    linewidth = 0.5
  ) +
  
  # Rolling mean before 2015
  geom_line(
    data = filter(df_all_anom, Year < 2015),
    aes(y = Rolling_Mean_anom),
    color = "black",
    linewidth = 0.5
  ) +
  
  # Rolling mean after 2015
  geom_line(
    data = filter(df_all_anom, Year >= 2015),
    aes(y = Rolling_Mean_anom, color = SSP),
    linewidth = 0.5
  ) +
  
  # Rolling SD ribbon before 2015
  geom_ribbon(
    data = filter(df_all_anom, Year < 2015),
    aes(
      ymin = Rolling_Mean_anom - Rolling_SD_anom,
      ymax = Rolling_Mean_anom + Rolling_SD_anom
    ),
    fill = "grey",
    alpha = 0.15
  ) +
  
  # Rolling SD ribbons after 2015
  geom_ribbon(
    data = filter(df_all_anom, Year >= 2015),
    aes(
      ymin = Rolling_Mean_anom - Rolling_SD_anom,
      ymax = Rolling_Mean_anom + Rolling_SD_anom,
      fill = SSP
    ),
    alpha = 0.15
  ) +
  
  # 2015 divider
  geom_vline(
    xintercept = 2015,
    linetype = "dashed",
    color = "black"
  ) +
  
  labs(
    x = "Year",
    y = "Bio-available N [Pg]",
    color = "SSP Scenario",
    fill = "SSP Scenario"
  ) +
  
  scale_y_continuous(
    breaks = seq(0, 8, by = 2)
  ) +
  
  scale_color_manual(
    values = ssp_colors,
    guide = guide_legend(
      keyheight = unit(0.35, "cm"),
      keywidth  = unit(0.35, "cm"),
      override.aes = list(linewidth = 0.7)
    )
  ) +
  
  scale_fill_manual(
    values = ssp_colors,
    guide = guide_legend(
      keyheight = unit(0.35, "cm"),
      keywidth  = unit(0.35, "cm"),
      override.aes = list(alpha = 0.4)
    )
  ) +
  
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    legend.spacing.y = unit(0.08, "cm"),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8)
  )


# plot kt_values
kt_126<-rast("1%_min/arctic_k_t_126_no_lim_mean.tif")
kt_245<-rast("1%_min/arctic_k_t_245_no_lim_mean.tif")
kt_370<-rast("1%_min/arctic_k_t_370_no_lim_mean.tif")
kt_585<-rast("1%_min/arctic_k_t_585_no_lim_mean.tif")

kt_n2o<-rast("n2o/arctic_total_n2o_370_mean.tif")


cell_area <- cellSize(kt_n2o, unit = "m")

arctic_mean <- function(r, cell_area) {
  global(r, fun = "mean", weights = cell_area, na.rm = TRUE)[,1]
}

kt_n2o <- arctic_mean(kt_n2o, cell_area)

kt_126 <- arctic_mean(kt_126, cell_area)
kt_245 <- arctic_mean(kt_245, cell_area)
kt_370 <- arctic_mean(kt_370, cell_area)
kt_585 <- arctic_mean(kt_585, cell_area)
# ---- Wide -> long ----
df <- data.frame(
  Year   = years,
  SSP126 = kt_126,
  SSP245 = kt_245,
  SSP370 = kt_370,
  SSP585 = kt_585
)

df_long <- df %>%
  pivot_longer(starts_with("SSP"),
               names_to = "SSP",
               values_to = "Mean") %>%
  mutate(
    Period = ifelse(Year < 2015, "Before 2015", "After 2015"),
    SSP = recode(
      SSP,
      "SSP126" = "SSP1-2.6",
      "SSP245" = "SSP2-4.5",
      "SSP370" = "SSP3-7.0",
      "SSP585" = "SSP5-8.5"
    )
  )



ggplot(df_long, aes(x = Year, y = Mean, group = SSP)) +
  
  # Raw values
  geom_line(data = filter(df_long, Period == "Before 2015"),
            color = "black", linewidth = 0.6, alpha = 0.8) +
  
  geom_line(data = filter(df_long, Period == "After 2015"),
            aes(color = SSP), linewidth = 0.6, alpha = 0.8) +
  
  geom_vline(xintercept = 2015, linetype = "dashed") +
  
  labs(
    x = "Year",
    y = expression(paste("Arctic mean ", k[t])),
    color = "SSP Scenario"
  ) +
  
  scale_color_manual(values = ssp_colors) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text  = element_text(size = 10),
    axis.title = element_text(size = 10)
  )

# ---- Anomaly function (means only) ----
compute_anomaly_mean <- function(d, ref_start = 2000, ref_end = 2020) {
  ref_mean <- d %>%
    filter(Year >= ref_start, Year <= ref_end) %>%
    summarise(m = mean(Mean, na.rm = TRUE)) %>%
    pull(m)
  
  d %>%
    mutate(
      Mean_anom = Mean - ref_mean,
      Rolling_Mean_anom = zoo::rollapply(
        Mean_anom, width = 20, FUN = mean,
        align = "center", fill = NA, na.rm = TRUE
      )
    )
}

df_all_anom <- df_long %>%
  group_by(SSP) %>%
  group_modify(~ compute_anomaly_mean(.x, ref_start = 2000, ref_end = 2020)) %>%
  ungroup() %>%
  mutate(
    Period = ifelse(Year < 2015, "Before 2015", "After 2015"),
    SSP = recode(
      SSP,
      "SSP126" = "SSP1-2.6",
      "SSP245" = "SSP2-4.5",
      "SSP370" = "SSP3-7.0",
      "SSP585" = "SSP5-8.5"
    )
  )

# ---- Colors ----
ssp_colors <- c(
  "SSP3-7.0" = "#D73027",
  "SSP5-8.5" = "#7B3294",
  "SSP2-4.5" = "orange",
  "SSP1-2.6" = "blue"
)

# ---- Plot (no SD, no ribbons) ----
ggplot(df_all_anom, aes(x = Year, y = Mean_anom, group = SSP)) +
  
  # Raw anomalies: black before 2015
  geom_line(data = filter(df_all_anom, Period == "Before 2015"),
            color = "black", linewidth = 0.6, alpha = 0.8) +
  
  # Raw anomalies: colored after 2015
  geom_line(data = filter(df_all_anom, Period == "After 2015"),
            aes(color = SSP), linewidth = 0.6, alpha = 0.8) +
  
  # Rolling mean anomalies: black before 2015
  geom_line(data = filter(df_all_anom, Year < 2015),
            aes(y = Rolling_Mean_anom),
            color = "black", linewidth = 1.2) +
  
  # Rolling mean anomalies: colored after 2015
  geom_line(data = filter(df_all_anom, Year >= 2015),
            aes(y = Rolling_Mean_anom, color = SSP),
            linewidth = 1.2) +
  
  geom_vline(xintercept = 2015, linetype = "dashed", color = "black") +
  
  labs(
    x = "Year",
    y = expression(paste("Mean ", k[t], " anomaly (relative to 2000–2020)")),
    color = "SSP Scenario"
  ) +
  
  scale_color_manual(values = ssp_colors) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 10)
  )





## spatial plotting

## plot thawed nitrogen maps
# figures 1 a and b

library(terra)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(here)
library(sf)
library(viridis)
library(rnaturalearth)
library(raster)
library(rnaturalearthdata)
library(ggspatial)
library(tidyr)
library(viridis)


gc()
# Fig 1a and b)

# Load the NetCDF file
thawedN <- rast("1.9%_min/arctic_min_N_pool_585_no_lim_mean.tif")  # replace with your actual file name

# Extract years from names if necessary
years <- 1850:2099
names(thawedN) <- paste0("Y", years)

y0 <- 2000
y1 <- 2099   # last year in your stack

i0 <- which(years == y0)
i1 <- which(years == y1)

rate_per_year <- (thawedN[[i1]] - thawedN[[i0]]) / (y1 - y0)  # units: (your thawedN units) per year
plot(rate_per_year)

df_plot <- function(r, period) {
  df <- as.data.frame(r, xy = TRUE)
  names(df)[3] <- "thawed_N"
  df$Period <- period
  return(df)
}


df_min_n_future <- df_plot(rate_per_year, "2100")
# convert in g/m2
df_min_n_future <- df_min_n_future %>%
  dplyr::mutate(
    thawed_N = thawed_N * 1000
  )

## for a flat circular projection
# Step 1: Filter the data to include only the Arctic region
arctic_df <- subset(df_min_n_future)

# Step 2: Convert the filtered data to an sf object
arctic_sf <- st_as_sf(arctic_df, coords = c("x", "y"), crs = 4326)  # WGS84

# Step 3: Define the LAEA projection centered on the North Pole
laea_crs <- "+proj=laea +lat_0=90 +lon_0=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# Step 4: Reproject the data into the LAEA projection
arctic_sf_laea <- st_transform(arctic_sf, crs = laea_crs)
bbox <- st_bbox(arctic_sf_laea)
print(bbox)
# Define custom limits for the plot
xlim <- c(-4233688, 4592082)  # Adjust these values based on your data
ylim <- c(-2868318, 3264652)  # Adjust these values based on your data
coastlines <- ne_coastline(scale = "medium", returnclass = "sf")
summary(arctic_sf_laea$thawed_N)
# Step 5: Reproject the coastlines to match the LAEA projection
coastlines_laea <- st_transform(coastlines, crs = laea_crs)
arctic_sf_laea <- arctic_sf_laea[!is.na(arctic_sf_laea$thawed_N), ]
# Step 6: Create the plot with the LAEA projection

spatial_plot <- ggplot() +
  geom_sf(data = arctic_sf_laea, aes(color = thawed_N), size = 0.01) +
  
  scale_color_gradientn(
    colours = c("dodgerblue3", "#abd9e9", "orange"),
    values = scales::rescale(c(0, 2, 4, 6, 8, 10)),
    limits = c(0, 10),
    breaks = c(0, 2, 4, 6, 8, 10),
    labels = c("0", "2", "4", "6", "8", "10"),
    name = expression("N (g m"^{-2}*" yr"^{-1}*")"),
    oob = scales::squish,
    guide = guide_colorbar(
      barheight = unit(2.5, "cm"),
      barwidth  = unit(0.25, "cm"),
      ticks = FALSE
    )
  ) +
  geom_sf(data = coastlines_laea, color = "black", size = 0.05) +  
  coord_sf(crs = laea_crs, xlim = xlim, ylim = ylim) +  
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    axis.text = element_text(size = 8)
  ) +
  labs(x = "", y = "")


# figure 2: 
library(patchwork)


combined_plot
combined_plot<-spatial_plot / plot_min + plot_annotation(tag_levels = "a", tag_suffix = ") ")
ggsave("Fig2.png", combined_plot,
       width = 15,
       height = 15,
       units = "cm", 
       #scale = 1.4,
       dpi = 500,  #70
       #device = cairo_pdf
)


ggsave("Fig2_a.pdf", spatial_plot,
       width = 35,
       height = 30,
       units = "cm", 
       #scale = 1.4,
       #dpi = 500,  #70
       #device = cairo_pdf
)

#############################
library(terra)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(here)
library(sf)
library(viridis)
library(rnaturalearth)
library(raster)
library(rnaturalearthdata)
library(ggspatial)
library(tidyr)
library(viridis)
##### supplementary 

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


## for a flat circular projection
# Step 1: Filter the data to include only the Arctic region
arctic_df <- subset(delta_df)

# Step 2: Convert the filtered data to an sf object
arctic_sf <- st_as_sf(arctic_df, coords = c("x", "y"), crs = 4326)  # WGS84

# Step 3: Define the LAEA projection centered on the North Pole
laea_crs <- "+proj=laea +lat_0=90 +lon_0=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# Step 4: Reproject the data into the LAEA projection
arctic_sf_laea <- st_transform(arctic_sf, crs = laea_crs)
bbox <- st_bbox(arctic_sf_laea)
print(bbox)
# Define custom limits for the plot
xlim <- c(-4233688, 4592082)  # Adjust these values based on your data
ylim <- c(-2868318, 3264652)  # Adjust these values based on your data
coastlines <- ne_coastline(scale = "medium", returnclass = "sf")
summary(arctic_sf_laea$delta_ALD)
# Step 5: Reproject the coastlines to match the LAEA projection
coastlines_laea <- st_transform(coastlines, crs = laea_crs)
arctic_sf_laea <- arctic_sf_laea[!is.na(arctic_sf_laea$delta_ALD), ]
# Step 6: Create the plot with the LAEA projection
ggplot() +
  geom_sf(data = arctic_sf_laea, aes(color = delta_ALD), size = 0.1) +
  
  scale_color_gradientn(
    colours = c("darkblue", "gold", "darkred"),
    values = scales::rescale(c(-1, 0, 1,2,3, 4,5, 6)),
    limits = c(0, 6),
    breaks = c(-1,0, 1, 2, 3, 4, 5,6),
    labels = c("-1","0", "1", "2", "3", "4", "5", "6"),
    name = expression("ALD"),
    oob = scales::squish,
    guide = guide_colorbar(
      barheight = unit(2.5, "cm"),
      barwidth  = unit(0.25, "cm"),
      ticks = FALSE
    )
  ) +
  geom_sf(data = coastlines_laea, color = "black", size = 0.5) +  
  coord_sf(crs = laea_crs, xlim = xlim, ylim = ylim) +  
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    axis.text = element_text(size = 8)
  ) +
  labs(x = "", y = "")

ALD<-ggplot() +
  geom_sf(data = arctic_sf_laea, aes(color = delta_ALD), size = 0.1) +
  
  scale_color_gradient2(
    low = "#4575b4",   # muted blue
    mid = "white",
    high = "#d73027",  # muted red
    midpoint = 0,
    limits = c(0, 6),
    name = expression(Delta * "ALD [m]"),
    guide = guide_colorbar(
      barheight = unit(2.5, "cm"),
      barwidth  = unit(0.25, "cm"),
      ticks = TRUE
    )
  )+
  geom_sf(data = coastlines_laea, color = "black", size = 0.1) +  
  coord_sf(crs = laea_crs, xlim = xlim, ylim = ylim) +  
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    axis.text = element_text(size = 6)
  ) +
  labs(x = "", y = "")

ggsave("Fig_ALD.png", ALD,
       width = 8,
       height = 8,
       units = "cm", 
       #scale = 1.4,
       dpi = 600,  #70
       #device = cairo_pdf
)


# to plot all ssps together: 
library(terra)
library(sf)
library(ggplot2)
library(patchwork)
rast_mean_585 <- rast("mean_ssp585_corr_new.nc")
rast_mean_370 <- rast("mean_ssp370_corr_new.nc")
rast_mean_245 <- rast("mean_ssp245_corr_new.nc")
rast_mean_126 <- rast("mean_ssp126_corr_new.nc")

# ---- helper function: calculate delta ALD ----
calc_delta_ald <- function(r) {
  baseline <- r[[time(r) >= as.Date("2000-07-01") &
                   time(r) <= as.Date("2020-07-01")]]
  
  future <- r[[time(r) >= as.Date("2080-07-01") &
                 time(r) <= as.Date("2099-07-01")]]
  
  delta <- mean(future, na.rm = TRUE) - mean(baseline, na.rm = TRUE)
  names(delta) <- "delta_ALD"
  delta
}

# ---- calculate delta rasters ----
delta_126 <- calc_delta_ald(rast_mean_126)
delta_245 <- calc_delta_ald(rast_mean_245)
delta_370 <- calc_delta_ald(rast_mean_370)
delta_585 <- calc_delta_ald(rast_mean_585)

# ---- convert to sf/dataframe for plotting ----
make_sf_df <- function(r) {
  as.points(r, na.rm = TRUE) |>
    st_as_sf()
}

sf_126 <- make_sf_df(delta_126)
sf_245 <- make_sf_df(delta_245)
sf_370 <- make_sf_df(delta_370)
sf_585 <- make_sf_df(delta_585)

# ---- common plotting function ----
plot_delta_ald <- function(sf_obj, title_text) {
  ggplot() +
    geom_sf(data = sf_obj, aes(color = delta_ALD), size = 0.08) +
    geom_sf(data = coastlines_laea, color = "black", size = 0.1) +
    coord_sf(crs = laea_crs, xlim = xlim, ylim = ylim) +
    scale_color_viridis_c(
      option = "viridis",
      limits = c(-1, 6),
      breaks = c(-1,0, 1, 2, 3, 4, 5, 6),
      name = expression(Delta*"ALD [m]"),
      guide = guide_colorbar(
        barheight = unit(2.5, "cm"),
        barwidth  = unit(0.25, "cm"),
        ticks = FALSE
      )
    ) +
    labs(title = title_text, x = "", y = "") +
    theme_minimal() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 7),
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(size = 7, face = "bold")
    )
}

# ---- build plots ----
p126 <- p126 + labs(title = "a), SSP1-2.6")
p245 <- p245 + labs(title = "b), SSP2-4.5")
p370 <- p370 + labs(title = "c), SSP3-7.0")
p585 <- p585 + labs(title = "d), SSP5-8.5")

combined_plot <- (p126 + p245 + p370 + p585) +
  plot_layout(ncol = 2, guides = "collect")



combined_plot

ggsave("Fig9.png", combined_plot,
       width = 15,
       height = 15,
       units = "cm", 
       #scale = 1.4,
       dpi = 250,  #70
       #device = cairo_pdf
)










# plot Palmtag dataset: 

palmtag<-rast("TN_30deg_corr.nc")

df_plot <- function(r, period) {
  df <- as.data.frame(r, xy = TRUE)
  names(df)[3] <- "TN"
  df$Period <- period
  return(df)
}

df_thawedN<-df_plot(palmtag, "Total")

## for a flat circular projection
# Step 1: Filter the data to include only the Arctic region
arctic_df <- subset(df_thawedN)

# Step 2: Convert the filtered data to an sf object
arctic_sf <- st_as_sf(arctic_df, coords = c("x", "y"), crs = 4326)  # WGS84

# Step 3: Define the LAEA projection centered on the North Pole
laea_crs <- "+proj=laea +lat_0=90 +lon_0=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# Step 4: Reproject the data into the LAEA projection
arctic_sf_laea <- st_transform(arctic_sf, crs = laea_crs)
bbox <- st_bbox(arctic_sf_laea)
print(bbox)
# Define custom limits for the plot
xlim <- c(-4000000, 4000000)  # Adjust these values based on your data
ylim <- c(-4000000, 4000000)  # Adjust these values based on your data
coastlines <- ne_coastline(scale = "medium", returnclass = "sf")
summary(arctic_sf_laea$TN)
# Step 5: Reproject the coastlines to match the LAEA projection
coastlines_laea <- st_transform(coastlines, crs = laea_crs)
arctic_sf_laea <- arctic_sf_laea[!is.na(arctic_sf_laea$TN), ]
# Step 6: Create the plot with the LAEA projection
ggplot() +
  geom_sf(data = arctic_sf_laea, aes(color = TN), size = 0.1) +
  scale_color_gradientn(
    colours = c("lightblue", "slateblue3", "orange"),
    values = scales::rescale(
      c(min(arctic_sf_laea$TN, na.rm = TRUE),
        median(arctic_sf_laea$TN, na.rm = TRUE),
        max(arctic_sf_laea$TN, na.rm = TRUE))
    ),
    name = expression("N [kg m"^{-2}*"]")  )+
  geom_sf(data = coastlines_laea, color = "black", size = 0.1) +  
  coord_sf(crs = laea_crs, xlim = xlim, ylim = ylim) +  
  theme_minimal() +
  theme(legend.position = "right", 
        #panel.grid=element_blank()
  ) +  
  labs(
    x = "", y = "",
    #title="Yearly increase in mineralised soil nitrogen, 1880 - 2100",
    #subtitle="SSP 5-8.5"
  )







### plot surface temperature 
st_585<-rast("mean_ST_585_clean_final.nc")
st_370<-rast("mean_ST_370_clean_final.nc")
st_245<-rast("mean_ST_245_clean_final.nc")
st_126<-rast("mean_ST_126_clean_final.nc")

plot(st_585[[2]])
ext_sub <- ext(-179.95, 179.95, 60, 90)

st_585 <- crop(st_585, ext_sub)
st_370 <- crop(st_370, ext_sub)
st_245 <- crop(st_245, ext_sub)
st_126 <- crop(st_126, ext_sub)

st_585 <- st_585 - 273.15
st_126 <- st_126 - 273.15
st_245 <- st_245 - 273.15
st_370 <- st_370 - 273.15


cell_area <- cellSize(st_585, unit = "m")

arctic_mean <- function(r, cell_area) {
  global(r, fun = "mean", weights = cell_area, na.rm = TRUE)[,1]
}

ts_585 <- arctic_mean(st_585, cell_area)
ts_126 <- arctic_mean(st_126, cell_area)
ts_245 <- arctic_mean(st_245, cell_area)
ts_370 <- arctic_mean(st_370, cell_area)

years <- 1850:2099

df <- data.frame(
  Year = years,
  SSP126 = ts_126,
  SSP245 = ts_245,
  SSP370 = ts_370,
  SSP585 = ts_585
)

ggplot(df, aes(x = Year)) +
  geom_line(aes(y = SSP126, color = "SSP126"), linewidth = 1) +
  geom_line(aes(y = SSP245, color = "SSP245"), linewidth = 1) +
  geom_line(aes(y = SSP370, color = "SSP370"), linewidth = 1) +
  geom_line(aes(y = SSP585, color = "SSP585"), linewidth = 1) +
  labs(
    y = "Arctic mean surface temperature (°C)",
    x = "Year",
    color = "Scenario",
    title = "Arctic mean surface temperature (60–90°N)"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    text = element_text(size = 12)
  )

df_long <- df %>%
  pivot_longer(starts_with("SSP"),
               names_to = "SSP",
               values_to = "Mean") %>%
  mutate(
    Period = ifelse(Year < 2015, "Before 2015", "After 2015"),
    SSP = recode(
      SSP,
      "SSP126" = "SSP1-2.6",
      "SSP245" = "SSP2-4.5",
      "SSP370" = "SSP3-7.0",
      "SSP585" = "SSP5-8.5"
    )
  )


# ---- Anomaly function (means only) ----
compute_anomaly_mean <- function(d, ref_start = 2000, ref_end = 2020) {
  ref_mean <- d %>%
    filter(Year >= ref_start, Year <= ref_end) %>%
    summarise(m = mean(Mean, na.rm = TRUE)) %>%
    pull(m)
  
  d %>%
    mutate(
      Mean_anom = Mean - ref_mean,
      Rolling_Mean_anom = zoo::rollapply(
        Mean_anom, width = 20, FUN = mean,
        align = "center", fill = NA, na.rm = TRUE
      )
    )
}

df_all_anom <- df_long %>%
  group_by(SSP) %>%
  group_modify(~ compute_anomaly_mean(.x, ref_start = 2000, ref_end = 2020)) %>%
  ungroup() %>%
  mutate(
    Period = ifelse(Year < 2015, "Before 2015", "After 2015"),
    SSP = recode(
      SSP,
      "SSP126" = "SSP1-2.6",
      "SSP245" = "SSP2-4.5",
      "SSP370" = "SSP3-7.0",
      "SSP585" = "SSP5-8.5"
    )
  )

# ---- Colors ----
ssp_colors <- c(
  "SSP3-7.0" = "#D73027",
  "SSP5-8.5" = "#7B3294",
  "SSP2-4.5" = "orange",
  "SSP1-2.6" = "blue"
)

# ---- Plot (no SD, no ribbons) ----
ggplot(df_all_anom, aes(x = Year, y = Mean_anom, group = SSP)) +
  
  # Raw anomalies: black before 2015
  geom_line(data = filter(df_all_anom, Period == "Before 2015"),
            color = "black", linewidth = 0.6, alpha = 0.8) +
  
  # Raw anomalies: colored after 2015
  geom_line(data = filter(df_all_anom, Period == "After 2015"),
            aes(color = SSP), linewidth = 0.6, alpha = 0.8) +
  
  # Rolling mean anomalies: black before 2015
  geom_line(data = filter(df_all_anom, Year < 2015),
            aes(y = Rolling_Mean_anom),
            color = "black", linewidth = 1.2) +
  
  # Rolling mean anomalies: colored after 2015
  geom_line(data = filter(df_all_anom, Year >= 2015),
            aes(y = Rolling_Mean_anom, color = SSP),
            linewidth = 1.2) +
  
  geom_vline(xintercept = 2015, linetype = "dashed", color = "black") +
  
  labs(
    x = "Year",
    y = expression(paste("Surface temperature anomaly [°C]")),
    color = "SSP Scenario"
  ) +
  
  scale_color_manual(values = ssp_colors) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 10)
  )




### n2o results:

library(terra)
n2o_126_mean<-rast("n2o/arctic_total_n2o_ssp126_mean_new.tif")
n2o_245_mean<-rast("n2o/arctic_total_n2o_ssp245_mean_new.tif")
n2o_370_mean<-rast("n2o/arctic_total_n2o_ssp370_mean_new.tif")
n2o_585_mean<-rast("n2o/arctic_total_n2o_ssp585_mean_new.tif")

n2o_126_std<-rast("n2o/arctic_total_n2o_ssp126_std_new.tif")
n2o_245_std<-rast("n2o/arctic_total_n2o_ssp245_std_new.tif")
n2o_370_std<-rast("n2o/arctic_total_n2o_ssp370_std_new.tif")
n2o_585_std<-rast("n2o/arctic_total_n2o_ssp585_std_new.tif")

#convert to Tg N equivalent
# conversion factor: 
# M (N2O) = 44.013 g/mol
# M (N) = 14.007 g / mol
# Atomic weight of N (2 N atoms in N₂O) = 2 × 14.007 = 28.014 g/mol
# So, the fraction of N in N₂O is: 28.014 / 44.013 = 0.636

#total_area_emissions_Tg_Neq <- total_area_emissions_mg 
# values in mg--> need to convert to Tg : 1e15

calculate_present_day_n2o <- function(n2o_mean_stack,
                                      years = 1850:2099,
                                      ref_years = 1880:1900,
                                      pres_years = 2000:2020) {
  
  ref_idx  <- which(years %in% ref_years)
  pres_idx <- which(years %in% pres_years)
  
  # per-cell means for each period
  ref_mean  <- mean(n2o_mean_stack[[ref_idx]],  na.rm = TRUE)
  pres_mean <- mean(n2o_mean_stack[[pres_idx]], na.rm = TRUE)
  
  # convert difference to annual rate using midpoint-year difference
  ref_mid  <- mean(range(ref_years))
  pres_mid <- mean(range(pres_years))
  dt <- pres_mid - ref_mid
  
  annual_rate_anomaly <- (pres_mean - ref_mean) / dt   
  # land footprint
  land_mask <- !is.na(annual_rate_anomaly)
  
  # cell area, masked to land only (ocean becomes NA)
  area_m2 <- cellSize(annual_rate_anomaly, unit = "m")
  area_m2 <- mask(area_m2, land_mask, maskvalues = FALSE)
  
  # area-weighted mean over land
  area_total_emissions <- global(annual_rate_anomaly * area_m2, "sum", na.rm = TRUE)[[1]]
  area_total_emissions<-area_total_emissions* 0.636 / 1e15 # Tg n2o_n emissions
  
  return(area_total_emissions)
}

calculate_present_day_n2o <- function(n2o_mean_stack,
                                      years = 1850:2099,
                                      ref_years = 1880:1900,
                                      pres_years = 2000:2020) {
  
  ref_idx  <- which(years %in% ref_years)
  pres_idx <- which(years %in% pres_years)
  
  if (length(ref_idx) == 0 || length(pres_idx) == 0)
    stop("Reference or present-day years not found in 'years' vector.")
  
  # per-cell means for each period
  ref_mean  <- mean(n2o_mean_stack[[ref_idx]],  na.rm = TRUE)
  pres_mean <- mean(n2o_mean_stack[[pres_idx]], na.rm = TRUE)
  
  # convert difference to annual rate using midpoint-year difference
  ref_mid  <- mean(range(ref_years))
  pres_mid <- mean(range(pres_years))
  dt <- pres_mid - ref_mid
  
  annual_rate_anomaly <- (pres_mean - ref_mean) / dt   # kg N/m2/yr
  
  # land footprint
  land_mask <- !is.na(annual_rate_anomaly)
  
  # cell area, masked to land only (ocean becomes NA)
  area_m2 <- cellSize(annual_rate_anomaly, unit = "m")
  area_m2 <- mask(area_m2, land_mask, maskvalues = FALSE)
  
  # area-weighted mean over land
  num <- global(annual_rate_anomaly * area_m2, "sum", na.rm = TRUE)[[1]]
  den <- global(area_m2, "sum", na.rm = TRUE)[[1]]
  result<-num / den
  return(result)
}



# Calculate for each SSP
n2o_126_present_day_value_mean <- calculate_present_day_n2o(n2o_126_mean)
n2o_245_present_day_value_mean <- calculate_present_day_n2o(n2o_245_mean)
n2o_370_present_day_value_mean <- calculate_present_day_n2o(n2o_370_mean)
n2o_585_present_day_value_mean <- calculate_present_day_n2o(n2o_585_mean)

# Calculate for each SSP
n2o_126_present_day_value_std <- calculate_present_day_n2o(n2o_126_std)
n2o_245_present_day_value_std <- calculate_present_day_n2o(n2o_245_std)
n2o_370_present_day_value_std <- calculate_present_day_n2o(n2o_370_std)
n2o_585_present_day_value_std <- calculate_present_day_n2o(n2o_585_std)


# Create the dataframe
n2o_present_day_df <- data.frame(
  SSP = c("SSP126", "SSP245", "SSP370", "SSP585"),
  present_day_n2o_mean = c(
    n2o_126_present_day_value_mean,
    n2o_245_present_day_value_mean,
    n2o_370_present_day_value_mean,
    n2o_585_present_day_value_mean
  ),
  present_day_n2o_std = c(
    n2o_126_present_day_value_std,   # Assun2og you have std stacks too
    n2o_245_present_day_value_std,
    n2o_370_present_day_value_std,
    n2o_585_present_day_value_std
  )
)

# Convert to tibble for better printing
library(tibble)
n2o_present_day_df <- as_tibble(n2o_present_day_df)

# View the result
print(n2o_present_day_df)






# future n2o



calculate_future_n2o <- function(n2o_mean_stack,
                                      years = 1850:2099,
                                      ref_years = 2000:2020,
                                      pres_years = 2080:2099) {
  
  ref_idx  <- which(years %in% ref_years)
  pres_idx <- which(years %in% pres_years)
  
  if (length(ref_idx) == 0 || length(pres_idx) == 0)
    stop("Reference or present-day years not found in 'years' vector.")
  
  # per-cell means for each period
  ref_mean  <- mean(n2o_mean_stack[[ref_idx]],  na.rm = TRUE)
  pres_mean <- mean(n2o_mean_stack[[pres_idx]], na.rm = TRUE)
  
  # convert difference to annual rate using midpoint-year difference
  ref_mid  <- mean(range(ref_years))
  pres_mid <- mean(range(pres_years))
  dt <- pres_mid - ref_mid
  
  # linear change per year at each grid cell (same units as stack per year)
  rate_future <- (pres_mean - ref_mean) / dt
  
  # land footprint
  land_mask <- !is.na(rate_future)
  
  # cell area (m²), land only
  area_m2 <- cellSize(rate_future, unit = "m")
  area_m2 <- mask(area_m2, land_mask, maskvalues = FALSE)
  
  # area-weighted mean rate over whole area
  area_total_emissions <- global(rate_future * area_m2, "sum", na.rm = TRUE)[[1]] # total mg N2O emissions
  area_total_emissions<-area_total_emissions* 0.636 / 1e15 # Tg n2o_n emissions
  
  return(area_total_emissions)
}

n2o_126_future_value_mean <- calculate_future_n2o(n2o_126_mean)
n2o_245_future_value_mean <- calculate_future_n2o(n2o_245_mean)
n2o_370_future_value_mean <- calculate_future_n2o(n2o_370_mean)
n2o_585_future_value_mean <- calculate_future_n2o(n2o_585_mean)

# Calculate for each SSP
n2o_126_future_value_std <- calculate_future_n2o(n2o_126_std)
n2o_245_future_value_std <- calculate_future_n2o(n2o_245_std)
n2o_370_future_value_std <- calculate_future_n2o(n2o_370_std)
n2o_585_future_value_std <- calculate_future_n2o(n2o_585_std)


# Create the dataframe
n2o_future_df <- data.frame(
  SSP = c("SSP126", "SSP245", "SSP370", "SSP585"),
  future_n2o_mean = c(
    n2o_126_future_value_mean,
    n2o_245_future_value_mean,
    n2o_370_future_value_mean,
    n2o_585_future_value_mean
  ),
  future_n2o_std = c(
    n2o_126_future_value_std,   # Assun2og you have std stacks too
    n2o_245_future_value_std,
    n2o_370_future_value_std,
    n2o_585_future_value_std
  )
)

# Convert to tibble for better printing
library(tibble)
n2o_future_df <- as_tibble(n2o_future_df)

# View the result
print(n2o_future_df)

write.csv(n2o_future_df, "n2o_future_budget_3m_lim_2%.csv")
n2o<-read.csv("n2o_future_budget_3m_lim_2%.csv")
library(dplyr)
library(tidyr)
library(ggplot2)

df_plot<-n2o
df_plot<-df_plot%>%
  rename(mean = future_n2o_mean,
         sd   = future_n2o_std)



df_plot <- df_plot %>%
  mutate(
    SSP = recode(
      SSP,
      "SSP585" = "5-8.5",
      "SSP370" = "3-7.0",
      "SSP245" = "2-4.5",
      "SSP126" = "1-2.6"
    ),
    SSP = factor(
      SSP,
      levels = c("1-2.6", "2-4.5", "3-7.0", "5-8.5")
    )
  )
library(ggplot2)


ssp_colors <- c(
  "3-7.0" = "#D73027",
  "5-8.5" = "#7B3294",
  "2-4.5" = "orange",
  "1-2.6" = "blue"
)

n2o_increase<-ggplot(df_plot, aes(x = SSP, y = mean, fill = SSP)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.7) +
  geom_errorbar(
    aes(ymin = mean - sd, ymax = mean + sd),
    width = 0.2,
    position = position_dodge(width = 0.8)
  ) +
  scale_fill_manual(values = ssp_colors) +
  scale_y_continuous(n.breaks = 8) +   # <- more ticks
  labs(
    y = expression("N"[2]*"O - N emissions (Tg N yr"^{-1}*")"),
    x = NULL
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
    plot.subtitle = element_text(hjust = 0.5, size = 8),
    axis.text.x = element_text(angle = 0, vjust=-3, size= 8),
    legend.position = "bottom"
  )



taiga_NPP<-read.csv("total_biomass_taiga_mean_yearly.csv")
tundra_NPP<-read.csv("total_biomass_Tundra_mean_yearly.csv")

library(dplyr)
library(tidyr)

# Tundra
tundra_df <- tundra_NPP %>%
  pivot_longer(-X, names_to = "Scenario", values_to = "value") %>%
  pivot_wider(names_from = X, values_from = value) %>%
  mutate(Biome = "Tundra")

# Taiga
taiga_df <- taiga_NPP %>%
  pivot_longer(-X, names_to = "Scenario", values_to = "value") %>%
  pivot_wider(names_from = X, values_from = value) %>%
  mutate(Biome = "Taiga")

# Combine both
df <- bind_rows(tundra_df, taiga_df)

df <- df %>%
  mutate(
    mean = ifelse(Biome == "Tundra", agb_tundra_mean, agb_taiga_mean),
    upper = ifelse(Biome == "Tundra", agb_Tundra_upper, agb_taiga_upper),
    lower = ifelse(Biome == "Tundra", agb_Tundra_lower, agb_taiga_lower)
  )

df <- df %>%
  mutate(
    Scenario = recode(
      Scenario,
      "V1" = "5-8.5",
      "V2" = "3-7.0",
      "V3" = "2-4.5",
      "V4" = "1-2.6"
    ),
    Scenario = factor(
      Scenario,
      levels = c("1-2.6", "2-4.5", "3-7.0", "5-8.5")
    )
  )
library(ggplot2)



dodge <- position_dodge(width = 0.7)

NPP_increase<-ggplot(df, aes(x = Scenario, y = mean, fill = Biome)) +
  geom_col(position = dodge, width = 0.6) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    position = dodge,
    width = 0.2
  ) +
  scale_y_continuous(n.breaks = 10) +   # <- more ticks
  #geom_text(
  #  aes(label = round(mean, 2)),
  #  position = dodge,
  #  hjust = -0.3,   # move text to the left of bar
  #  vjust = -0.2,
  #  size = 3.5
  #) +
  theme_minimal() +
  scale_fill_manual(
    values = c(
      "Tundra" = "#C7D6C1", 
      "Taiga"  = "#2F6B4F" 
    ))+labs(
      x = "",
      y = expression("NPP increase (" * g * C * m^{-2} * " yr"^{-1} * ")"),
      fill = "Biome"
    )+
  theme(legend.position = "bottom")


library(patchwork)

 

n2o_increase <- n2o_increase +
  theme(
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.09, "cm"),
    legend.spacing.y = unit(0.09, "cm")
  )

NPP_increase <- NPP_increase +
  theme(
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.09, "cm"),
    legend.spacing.y = unit(0.09, "cm")
  )


combined_plot<-n2o_increase | NPP_increase + plot_annotation(tag_levels = "a", tag_suffix = ") ")
ggsave("Fig4.png", combined_plot,
       width = 15,
       height = 15,
       dpi=500,
       units = "cm")


################################################################################
install.packages("DiagrammeR")
library(DiagrammeR)

grViz("
digraph {
  graph [layout = dot, rankdir = TB]

  permafrost [label = 'Permafrost organic matter']
  thaw_pool  [label = 'N thaw']
  organic    [label = 'Organic N']
  bio        [label = 'Bioavailable N']
  plant      [label = 'Plant uptake']
  loss       [label = 'Losses (N2O, runoff)']

  permafrost -> thaw_pool [label = 'thaw']
  thaw_pool -> organic
  thaw_pool -> bio [label = 'mineralisation']
  bio -> plant
  bio -> loss
}
")

