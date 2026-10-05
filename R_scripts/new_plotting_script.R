library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(zoo)

# ============================================================
# Total thawed N anomaly plot
# Mean ± SD, relative to 2000–2020
# 20-year rolling mean
# ============================================================

# -----------------------------
# 1. Settings
# -----------------------------

ssps <- c("126", "245", "370", "585")

common_extent <- ext(-179.95, 179.95, 60, 90)

out_dir <- "total_thawed_plots"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ssp_colors <- c(
  "SSP3-7.0" = "#D73027",
  "SSP5-8.5" = "#7B3294",
  "SSP2-4.5" = "orange",
  "SSP1-2.6" = "blue"
)

# -----------------------------
# 2. Helper function
# -----------------------------

extract_total_pg <- function(mean_rast, sd_rast, ssp, years, area_rast) {

  w <- terra::mask(area_rast, mean_rast[[1]])

  total_mean_kg <- terra::global(
    mean_rast * w,
    "sum",
    na.rm = TRUE
  )[[1]]

  total_sd_kg <- terra::global(
    sd_rast * w,
    "sum",
    na.rm = TRUE
  )[[1]]

  tibble::tibble(
    Year = years,
    Mean = total_mean_kg / 1e12,
    SD = total_sd_kg / 1e12,
    SSP = ssp
  )
}

# extract_total_pg <- function(mean_rast, sd_rast, ssp, years, area_rast) {
#   
#   w <- terra::mask(area_rast, mean_rast[[1]])
#   
#   # Arctic total mean for every year
#   total_mean_kg <- terra::global(
#     mean_rast * w,
#     "sum",
#     na.rm = TRUE
#   )[[1]]
#   
#   # Arctic total SD, assuming independent grid-cell errors
#   total_variance_kg2 <- terra::global(
#     (sd_rast * w)^2,
#     "sum",
#     na.rm = TRUE
#   )[[1]]
#   
#   total_sd_kg <- sqrt(total_variance_kg2)
#   
#   tibble::tibble(
#     Year = years,
#     Mean = total_mean_kg / 1e12,
#     SD = total_sd_kg / 1e12,
#     SSP = ssp
#   )
# }
# -----------------------------
# 3. Read and extract all SSPs
# -----------------------------

df_list <- list()

for (ssp in ssps) {
  
  cat("Processing SSP", ssp, "\n")
  
  mean_file <- paste0(
    "total_thawed_extended/arctic_total_thawed_",
    ssp,
    "_60N_mean.nc"
  )
  
  sd_file <- paste0(
    "total_thawed_extended/arctic_total_thawed_",
    ssp,
    "_60N_std.nc"
  )
  
  mean_rast <- rast(mean_file)
  sd_rast <- rast(sd_file)
  
  mean_rast <- crop(mean_rast, common_extent)
  sd_rast <- crop(sd_rast, common_extent)
  
  years <- as.integer(format(time(mean_rast), "%Y"))
  
  area_rast <- cellSize(mean_rast[[1]], unit = "m")
  
  df_list[[ssp]] <- extract_total_pg(
    mean_rast = mean_rast,
    sd_rast = sd_rast,
    ssp = paste0("SSP", ssp),
    years = years,
    area_rast = area_rast
  )
  
  rm(mean_rast, sd_rast, area_rast)
  gc()
}

df_all <- bind_rows(df_list)

active_layer_N_historical <- df_all %>%
  filter(Year >= 2080, Year <= 2099) %>%
  group_by(SSP, Year) %>%
  summarise(
    Mean = mean(Mean, na.rm = TRUE),
    SD = mean(SD, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(SSP) %>%
  summarise(
    Mean_active_layer_N_Pg = mean(Mean, na.rm = TRUE),
    Mean_annual_SD_Pg = mean(SD, na.rm = TRUE),
    Interannual_SD_Pg = sd(Mean, na.rm = TRUE),
    n_years = n(),
    .groups = "drop"
  )

print(active_layer_N_historical)

df_all %>%
  filter(Year == 2099) %>%
  select(SSP, Year, Mean, SD)
# -----------------------------
# 4. Compute anomaly relative to 2000–2020
# -----------------------------

compute_anomaly <- function(df) {
  
  # 1850–1900 baseline
  baseline <- df %>%
    filter(Year >= 1850, Year <= 1900)
  
  baseline_mean <- mean(baseline$Mean, na.rm = TRUE)
  baseline_sd   <- mean(baseline$SD, na.rm = TRUE)
  
  # Baseline-corrected series
  df <- df %>%
    mutate(
      Mean_baseline = Mean - baseline_mean,
      SD_baseline   = SD - baseline_sd
    )
  
  # 2000–2020 reference on the baseline-corrected series
  ref <- df %>%
    filter(Year >= 2000, Year <= 2020)
  
  ref_mean <- mean(ref$Mean_baseline, na.rm = TRUE)
  ref_sd <- mean(ref$SD_baseline, na.rm = TRUE)
  df %>%
    mutate(
      Mean_anom = Mean_baseline - ref_mean,
      SD_anom   = SD_baseline - ref_sd
    )
}

df_all_anom <- df_all %>%
  group_by(SSP) %>%
  group_modify(~ compute_anomaly(.x)) %>%
  ungroup()

# compute_anomaly <- function(df) {
#   
#   ref <- df %>%
#     filter(Year >= 2000, Year <= 2020)
#   
#   ref_mean <- mean(ref$Mean, na.rm = TRUE)
#   ref_std <- mean(ref$SD, na.rm = TRUE)
#   
#   df %>%
#     mutate(
#       Mean_anom = Mean - ref_mean,
#       SD_anom = SD - ref_std
#     )
# }
# 
# df_all_anom <- df_all %>%
#   group_by(SSP) %>%
#   group_modify(~ compute_anomaly(.x)) %>%
#   ungroup()

# -----------------------------
# 5. Rolling 20-year mean
# -----------------------------

df_all_anom <- df_all_anom %>%
  arrange(SSP, Year) %>%
  group_by(SSP) %>%
  mutate(
    Rolling_Mean_anom = rollapply(
      Mean_anom,
      width = 20,
      FUN = mean,
      align = "center",
      fill = NA,
      na.rm = TRUE
    ),
    Rolling_SD_anom = rollapply(
      SD_anom,
      width = 20,
      FUN = mean,
      align = "center",
      fill = NA,
      na.rm = TRUE
    )
  ) %>%
  ungroup()

# -----------------------------
# 6. Recode SSP labels
# -----------------------------

df_all_anom <- df_all_anom %>%
  mutate(
    SSP = recode(
      SSP,
      "SSP126" = "SSP1-2.6",
      "SSP245" = "SSP2-4.5",
      "SSP370" = "SSP3-7.0",
      "SSP585" = "SSP5-8.5"
    ),
    Period = ifelse(Year < 2015, "Before 2015", "After 2015")
  )

# -----------------------------
# 7. Save data
# -----------------------------

write.csv(
  df_all_anom,
  file.path(out_dir, "total_thawed_N_anomaly_rolling20_all_ssps.csv"),
  row.names = FALSE
)

# -----------------------------
# 8. Plot
# -----------------------------

plot_total_thawed <- ggplot(
  df_all_anom,
  aes(x = Year, y = Mean_anom, group = SSP)
) +
  
  # Before 2015: black raw lines
  geom_line(
    data = filter(df_all_anom, Period == "Before 2015"),
    color = "black",
    linewidth = 0.3,
    alpha = 0.5
  ) +
  
  # After 2015: coloured raw lines
  geom_line(
    data = filter(df_all_anom, Period == "After 2015"),
    aes(color = SSP),
    linewidth = 0.3,
    alpha = 0.5
  ) +
  
  # Rolling SD before 2015
  geom_ribbon(
    data = filter(df_all_anom, Year < 2015),
    aes(
      ymin = Rolling_Mean_anom - Rolling_SD_anom,
      ymax = Rolling_Mean_anom + Rolling_SD_anom
    ),
    fill = "grey",
    alpha = 0.15
  ) +
  
  # Rolling SD after 2015
  geom_ribbon(
    data = filter(df_all_anom, Year >= 2015),
    aes(
      ymin = Rolling_Mean_anom - Rolling_SD_anom,
      ymax = Rolling_Mean_anom + Rolling_SD_anom,
      fill = SSP
    ),
    alpha = 0.15
  ) +
  
  # Rolling mean before 2015
  geom_line(
    data = filter(df_all_anom, Year < 2015),
    aes(y = Rolling_Mean_anom),
    color = "black",
    linewidth = 0.3
  ) +
  
  # Rolling mean after 2015
  geom_line(
    data = filter(df_all_anom, Year >= 2015),
    aes(y = Rolling_Mean_anom, color = SSP),
    linewidth = 0.3
  ) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  geom_vline(xintercept = 2015, linetype = "dashed", color = "black") +
  
  scale_color_manual(values = ssp_colors) +
  scale_fill_manual(values = ssp_colors) +
  scale_y_continuous(breaks = seq(-5, 35, by = 5)) +
  ylim(-7, 30)+
  # Theme settings
  theme_minimal(base_size = 7) +
  theme(
    axis.text = element_text(size = 6.5),
    axis.title = element_text(size = 7),
    legend.text = element_text(size = 6.5),
    legend.title = element_text(size = 7),
    plot.tag = element_text(size = 8), 
    legend.position = "bottom"
  )+
  labs(
    x = "Year",
    y = "Total N thaw [Pg N]",
    color = "SSP Scenario",
    fill = "SSP Scenario",
    title= ""
  )

plot_total_thawed




df<-read.csv("ALD_no_lim_bias_corr.csv")

df_period <- df %>%
  mutate(
    Period = case_when(
      Year >= 1880 & Year <= 1900 ~ "preindustrial",
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

df_change_preindustrial <- df_mean %>%
  pivot_wider(names_from = Period, values_from = starts_with("mean")) %>%
  mutate(
    delta_585 = mean_585_baseline - mean_585_preindustrial,
    delta_370 = mean_370_baseline - mean_370_preindustrial,
    delta_245 = mean_245_baseline - mean_245_preindustrial,
    delta_126 = mean_126_baseline - mean_126_preindustrial
  )

print(df_change %>% select(starts_with("delta")))

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

long_data <- long_data %>%
  mutate(
    Mean_ALD = -Mean_ALD,
    Rolling_Mean_ALD = -Rolling_Mean_ALD
  )
# Plot
ALD<-ggplot(long_data, aes(x = Year, y = Mean_ALD, group = SSP)) +
  
  # Before 2015: Black lines & grey ribbons
  geom_line(data = filter(long_data, Period == "Before 2015"), color = "black", linewidth = 0.3) +
  #geom_ribbon(data = filter(long_data, Period == "Before 2015"),
  #            aes(ymin = Mean_ALD - STD_ALD, ymax = Mean_ALD + STD_ALD),
  #            fill = "grey", alpha = 0.2) +
  
  # After 2015: Colored lines & ribbons
  geom_line(data = filter(long_data, Period == "After 2015"),
            aes(color = SSP), linewidth = 0.3) +
  #geom_ribbon(data = filter(long_data, Period == "After 2015"),
  #            aes(ymin = Mean_ALD - STD_ALD, ymax = Mean_ALD + STD_ALD, fill = SSP),
  #            alpha = 0.2) +
  
  # Rolling Mean before 2015 (Black line)
  geom_line(data = filter(long_data, Year < 2015),
            aes(y = Rolling_Mean_ALD), color = "black", linewidth = 0.3, linetype = "solid") +
  
  # Rolling Mean after 2015 (Colored lines)
  geom_line(data = filter(long_data, Year >= 2015),
            aes(y = Rolling_Mean_ALD, color = SSP), linewidth = 0.3, linetype = "solid") +
  
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
  theme_minimal(base_size = 7) +
  theme(
    axis.text = element_text(size = 6.5),
    axis.title = element_text(size = 7),
    legend.text = element_text(size = 6.5),
    legend.title = element_text(size = 7),
    plot.tag = element_text(size = 8), 
    legend.position = "bottom"
  )

ALD


library(patchwork)

ALD <- ALD + theme(legend.position = "none")

combined_plot <- ALD | plot_total_thawed +
  plot_annotation(tag_levels = "a", tag_suffix = ")")&
  
  theme(legend.position = "right") &
  
  plot_layout(guides = "collect")



combined_plot
ggsave("Fig1.png", combined_plot,
       width = 15,
       height = 6,
       units = "cm", 
       #scale = 1.4,
       dpi = 500,  #70
       #device = cairo_pdf
)



ggsave(
  file.path(out_dir, "total_thawed_N_anomaly_rolling20_all_ssps.png"),
  plot_total_thawed,
  width = 8,
  height = 5,
  dpi = 300
)

#########################
# remaining organic N plot

# -----------------------------
# 1. Settings
# -----------------------------

ssps <- c("126", "245", "370", "585")

common_extent <- ext(-179.95, 179.95, 60, 90)

out_dir <- "organic_N_plots"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ssp_colors <- c(
  "SSP3-7.0" = "#D73027",
  "SSP5-8.5" = "#7B3294",
  "SSP2-4.5" = "orange",
  "SSP1-2.6" = "blue"
)

# -----------------------------
# 2. Helper function
# -----------------------------

extract_total_pg <- function(mean_rast, sd_rast, ssp, years, area_rast) {
  
  w <- terra::mask(area_rast, mean_rast[[1]])
  
  total_mean_kg <- terra::global(
    mean_rast * w,
    "sum",
    na.rm = TRUE
  )[[1]]
  
  total_sd_kg <- terra::global(
    sd_rast * w,
    "sum",
    na.rm = TRUE
  )[[1]]
  
  tibble::tibble(
    Year = years,
    Mean = total_mean_kg / 1e12,
    SD = total_sd_kg / 1e12,
    SSP = ssp
  )
}

# -----------------------------
# 3. Read and extract all SSPs
# -----------------------------

df_list <- list()

for (ssp in ssps) {
  
  cat("Processing SSP", ssp, "\n")
  
  mean_file <- paste0(
    "monthly_mineralised/latest_2perc_baseline/arctic_organic_N_pool_remaining_yearly_",
    ssp,
    "_mean_1850_2099_w_temp.nc"
  )
  
  sd_file <- paste0(
    "monthly_mineralised/whole_arctic_std_8deg/arctic_organic_N_pool_remaining_yearly_",
    ssp,
    "_mean_1850_2099_w_temp.nc"
  )
  
  mean_rast <- rast(mean_file)
  sd_rast <- rast(sd_file)
  
  mean_rast <- crop(mean_rast, common_extent)
  sd_rast <- crop(sd_rast, common_extent)
  
  years <- as.integer(format(time(mean_rast), "%Y"))
  
  area_rast <- cellSize(mean_rast[[1]], unit = "m")
  
  df_list[[ssp]] <- extract_total_pg(
    mean_rast = mean_rast,
    sd_rast = sd_rast,
    ssp = paste0("SSP", ssp),
    years = years,
    area_rast = area_rast
  )
  
  rm(mean_rast, sd_rast, area_rast)
  gc()
}

df_all <- bind_rows(df_list)

# -----------------------------
# 4. Compute anomaly relative to 2000–2020
# -----------------------------
mean_file <- paste0(
  "monthly_mineralised/whole_arctic_mean_60degN/arctic_organic_N_pool_remaining_yearly_",
  ssp,
  "_mean_1850_2099_w_temp.nc"
)

sd_file <- paste0(
  "monthly_mineralised/whole_arctic_std_8deg/arctic_organic_N_pool_remaining_yearly_",
  ssp,
  "_mean_1850_2099_w_temp.nc"
)

mean_rast <- rast(mean_file)
sd_rast <- rast(sd_file)




compute_anomaly <- function(df) {
  
  # 1850–1900 baseline
  baseline <- df %>%
    filter(Year >= 1850, Year <= 1900)
  
  baseline_mean <- mean(baseline$Mean, na.rm = TRUE)
  baseline_sd   <- mean(baseline$SD, na.rm = TRUE)
  
  # Baseline-corrected series
  df <- df %>%
    mutate(
      Mean_baseline = Mean - baseline_mean,
      SD_baseline   = SD - baseline_sd
    )
  
  # 2000–2020 reference on the baseline-corrected series
  ref <- df %>%
    filter(Year >= 2000, Year <= 2020)
  
  ref_mean <- mean(ref$Mean_baseline, na.rm = TRUE)
  ref_sd <- mean(ref$SD_baseline, na.rm = TRUE)
  df %>%
    mutate(
      Mean_anom = Mean_baseline - ref_mean,
      SD_anom   = SD_baseline - ref_sd
    )
}

df_all_anom <- df_all %>%
  group_by(SSP) %>%
  group_modify(~ compute_anomaly(.x)) %>%
  ungroup()

# compute_anomaly <- function(df) {
#   
#   ref <- df %>%
#     filter(Year >= 2000, Year <= 2020)
#   
#   ref_mean <- mean(ref$Mean, na.rm = TRUE)
#   ref_std <- mean(ref$SD, na.rm = TRUE)
#   
#   df %>%
#     mutate(
#       Mean_anom = Mean - ref_mean,
#       SD_anom = SD - ref_std
#     )
# }
# 
# df_all_anom <- df_all %>%
#   group_by(SSP) %>%
#   group_modify(~ compute_anomaly(.x)) %>%
#   ungroup()

# -----------------------------
# 5. Rolling 20-year mean
# -----------------------------

df_all_anom <- df_all_anom %>%
  arrange(SSP, Year) %>%
  group_by(SSP) %>%
  mutate(
    Rolling_Mean_anom = rollapply(
      Mean_anom,
      width = 20,
      FUN = mean,
      align = "center",
      fill = NA,
      na.rm = TRUE
    ),
    Rolling_SD_anom = rollapply(
      SD_anom,
      width = 20,
      FUN = mean,
      align = "center",
      fill = NA,
      na.rm = TRUE
    )
  ) %>%
  ungroup()

# -----------------------------
# 6. Recode SSP labels
# -----------------------------

df_all_anom <- df_all_anom %>%
  mutate(
    SSP = recode(
      SSP,
      "SSP126" = "SSP1-2.6",
      "SSP245" = "SSP2-4.5",
      "SSP370" = "SSP3-7.0",
      "SSP585" = "SSP5-8.5"
    ),
    Period = ifelse(Year < 2015, "Before 2015", "After 2015")
  )

# -----------------------------
# 7. Save data
# -----------------------------

write.csv(
  df_all_anom,
  file.path(out_dir, "organic_N_anomaly_rolling20_all_ssps.csv"),
  row.names = FALSE
)

# -----------------------------
# 8. Plot
# -----------------------------

plot_organic <- ggplot(
  df_all_anom,
  aes(x = Year, y = Mean_anom, group = SSP)
) +
  
  # Before 2015: black raw lines
  geom_line(
    data = filter(df_all_anom, Period == "Before 2015"),
    color = "black",
    linewidth = 0.4,
    alpha = 0.5
  ) +
  
  # After 2015: coloured raw lines
  geom_line(
    data = filter(df_all_anom, Period == "After 2015"),
    aes(color = SSP),
    linewidth = 0.4,
    alpha = 0.5
  ) +
  
  # Rolling SD before 2015
  geom_ribbon(
    data = filter(df_all_anom, Year < 2015),
    aes(
      ymin = Rolling_Mean_anom - Rolling_SD_anom,
      ymax = Rolling_Mean_anom + Rolling_SD_anom
    ),
    fill = "grey",
    alpha = 0.15
  ) +
  
  # Rolling SD after 2015
  geom_ribbon(
    data = filter(df_all_anom, Year >= 2015),
    aes(
      ymin = Rolling_Mean_anom - Rolling_SD_anom,
      ymax = Rolling_Mean_anom + Rolling_SD_anom,
      fill = SSP
    ),
    alpha = 0.15
  ) +
  
  # Rolling mean before 2015
  geom_line(
    data = filter(df_all_anom, Year < 2015),
    aes(y = Rolling_Mean_anom),
    color = "black",
    linewidth = 0.7
  ) +
  
  # Rolling mean after 2015
  geom_line(
    data = filter(df_all_anom, Year >= 2015),
    aes(y = Rolling_Mean_anom, color = SSP),
    linewidth = 0.7
  ) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  geom_vline(xintercept = 2015, linetype = "dashed", color = "black") +
  
  scale_color_manual(values = ssp_colors) +
  scale_fill_manual(values = ssp_colors) +
  scale_y_continuous(breaks = seq(-45, 45, by = 5)) +
  
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 8),
    legend.spacing.y = unit(0.08, "cm"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 10)
  ) +
  labs(
    x = "Year",
    y = "Dynamic organic N pool [Pg N]",
    color = "SSP Scenario",
    fill = "SSP Scenario"
  )

plot_organic

ggsave(
  file.path(out_dir, "organic_N_anomaly_rolling20_all_ssps.png"),
  plot_organic,
  width = 8,
  height = 5,
  dpi = 300
)


######################################################################


library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

# ============================================================
# Cumulative mineralised and bioavailable N anomalies
# Mean and standard-deviation simulations
# Reference period: 2000–2020
# ============================================================

# ------------------------------------------------------------
# 1. Settings
# ------------------------------------------------------------

base_dir <- "monthly_mineralised"

mean_dir <- file.path(
  base_dir,
  "test_arctic_long"
)

sd_dir <- file.path(
  base_dir,
  "test_arctic_long_std"
)

out_dir <- file.path(
  "coding_first_chapter",
  "diagnostic_csv_analysis"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

ssps <- c("126", "245", "370", "585")

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

reference_start <- 2000
reference_end   <- 2020

# ------------------------------------------------------------
# 2. Find diagnostic CSV
# ------------------------------------------------------------

find_diagnostic_csv <- function(folder, ssp) {
  
  expected_names <- c(
    paste0(
      "arctic_monthly_diagnostics_w_temp",
      ssp,
      "_1850_2100.csv"
    ),
    paste0(
      "arctic_monthly_diagnostics_w_temp_",
      ssp,
      "_1850_2100.csv"
    )
  )
  
  expected_paths <- file.path(folder, expected_names)
  
  existing <- expected_paths[file.exists(expected_paths)]
  
  if (length(existing) > 0) {
    return(existing[1])
  }
  
  # Search recursively in case CSVs are inside SSP subfolders
  candidates <- list.files(
    folder,
    pattern = paste0(
      "arctic_monthly_diagnostics.*",
      ssp,
      ".*\\.csv$"
    ),
    recursive = TRUE,
    full.names = TRUE
  )
  
  if (length(candidates) == 0) {
    stop(
      "Could not find diagnostic CSV for SSP ",
      ssp,
      " in:\n",
      folder
    )
  }
  
  if (length(candidates) > 1) {
    warning(
      "Multiple files found for SSP ",
      ssp,
      ". Using:\n",
      candidates[1]
    )
  }
  
  candidates[1]
}

# ------------------------------------------------------------
# 3. Read one monthly diagnostic CSV
# ------------------------------------------------------------

read_monthly_diagnostics <- function(folder, ssp, simulation) {
  
  file <- find_diagnostic_csv(folder, ssp)
  
  cat(
    "Reading",
    simulation,
    "SSP",
    ssp,
    ":\n",
    file,
    "\n"
  )
  
  df <- read_csv(
    file,
    show_col_types = FALSE
  )
  
  required_columns <- c(
    "Year",
    "mineralised_pg_monthly",
    "bioavailable_pg_monthly"
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(df)
  )
  
  if (length(missing_columns) > 0) {
    stop(
      "The following columns are missing from ",
      basename(file),
      ":\n",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  df %>%
    mutate(
      SSP_code = ssp,
      SSP = unname(ssp_labels[ssp]),
      Simulation = simulation
    )
}

# ------------------------------------------------------------
# 4. Convert monthly mean simulations to annual totals
# ------------------------------------------------------------

monthly_mean <- bind_rows(
  lapply(
    ssps,
    function(ssp) {
      read_monthly_diagnostics(
        folder = mean_dir,
        ssp = ssp,
        simulation = "Mean"
      )
    }
  )
)

annual_mean <- monthly_mean %>%
  filter(
    Year >= 1850,
    Year <= 2099
  ) %>%
  group_by(
    SSP_code,
    SSP,
    Year
  ) %>%
  summarise(
    Mineralised_mean = sum(
      mineralised_pg_monthly,
      na.rm = TRUE
    ),
    Bioavailable_mean = sum(
      bioavailable_pg_monthly,
      na.rm = TRUE
    ),
    n_months_mineralised = sum(
      !is.na(mineralised_pg_monthly)
    ),
    n_months_bioavailable = sum(
      !is.na(bioavailable_pg_monthly)
    ),
    .groups = "drop"
  )

# Check for incomplete years
incomplete_mean <- annual_mean %>%
  filter(
    n_months_mineralised != 12 |
      n_months_bioavailable != 12
  )

if (nrow(incomplete_mean) > 0) {
  warning(
    "Some mean simulation years do not contain 12 valid months."
  )
  
  print(incomplete_mean)
}

# ------------------------------------------------------------
# 5. Convert monthly SD simulations to annual SD
# ------------------------------------------------------------

monthly_sd <- bind_rows(
  lapply(
    ssps,
    function(ssp) {
      read_monthly_diagnostics(
        folder = sd_dir,
        ssp = ssp,
        simulation = "SD"
      )
    }
  )
)

annual_sd <- monthly_sd %>%
  filter(
    Year >= 1850,
    Year <= 2099
  ) %>%
  group_by(
    SSP_code,
    SSP,
    Year
  ) %>%
  summarise(
    
    # Root-sum-of-squares propagation across months
    Mineralised_SD = sqrt(
      sum(
        mineralised_pg_monthly^2,
        na.rm = TRUE
      )
    ),
    
    Bioavailable_SD = sqrt(
      sum(
        bioavailable_pg_monthly^2,
        na.rm = TRUE
      )
    ),
    
    n_months_mineralised_sd = sum(
      !is.na(mineralised_pg_monthly)
    ),
    
    n_months_bioavailable_sd = sum(
      !is.na(bioavailable_pg_monthly)
    ),
    
    .groups = "drop"
  )

# Check for incomplete SD years
incomplete_sd <- annual_sd %>%
  filter(
    n_months_mineralised_sd != 12 |
      n_months_bioavailable_sd != 12
  )

if (nrow(incomplete_sd) > 0) {
  warning(
    "Some standard-deviation simulation years do not contain 12 valid months."
  )
  
  print(incomplete_sd)
}

# ------------------------------------------------------------
# 6. Combine annual mean and SD
# ------------------------------------------------------------

annual_fluxes <- annual_mean %>%
  left_join(
    annual_sd,
    by = c(
      "SSP_code",
      "SSP",
      "Year"
    )
  ) %>%
  select(
    SSP_code,
    SSP,
    Year,
    Mineralised_mean,
    Mineralised_SD,
    Bioavailable_mean,
    Bioavailable_SD
  ) %>%
  arrange(
    SSP,
    Year
  )

# ------------------------------------------------------------
# 7. Convert to long format
# ------------------------------------------------------------

annual_long <- bind_rows(
  
  annual_fluxes %>%
    transmute(
      SSP_code,
      SSP,
      Year,
      Variable = "Mineralised N",
      Annual_Mean = Mineralised_mean,
      Annual_SD = Mineralised_SD
    ),
  
  annual_fluxes %>%
    transmute(
      SSP_code,
      SSP,
      Year,
      Variable = "Bioavailable N",
      Annual_Mean = Bioavailable_mean,
      Annual_SD = Bioavailable_SD
    )
) %>%
  mutate(
    SSP = factor(
      SSP,
      levels = unname(ssp_labels[ssps])
    ),
    Variable = factor(
      Variable,
      levels = c(
        "Mineralised N",
        "Bioavailable N"
      )
    )
  ) %>%
  arrange(
    Variable,
    SSP,
    Year
  )

# ------------------------------------------------------------
# 8. Calculate cumulative totals and anomalies
# ------------------------------------------------------------

cumulative_anomaly <- annual_long %>%
  group_by(
    Variable,
    SSP
  ) %>%
  arrange(
    Year,
    .by_group = TRUE
  ) %>%
  mutate(
    
    # Cumulative mean flux
    Cumulative_Mean = cumsum(
      replace_na(Annual_Mean, 0)
    ),
    
    # Propagate annual SD through cumulative sum using RSS
    Cumulative_SD = sqrt(
      cumsum(
        replace_na(Annual_SD, 0)^2
      )
    )
  ) %>%
  mutate(
    
    # Mean cumulative value during 2000–2020
    Reference_Cumulative_Mean = mean(
      Cumulative_Mean[
        Year >= reference_start &
          Year <= reference_end
      ],
      na.rm = TRUE
    ),
    
    # Cumulative anomaly centred on the reference period
    Cumulative_Anomaly = Cumulative_Mean -
      Reference_Cumulative_Mean,
    
    Lower = Cumulative_Anomaly - Cumulative_SD,
    Upper = Cumulative_Anomaly + Cumulative_SD
  ) %>%
  ungroup()

# ------------------------------------------------------------
# 9. Save annual and cumulative data
# ------------------------------------------------------------

write_csv(
  annual_long,
  file.path(
    out_dir,
    "annual_mineralised_bioavailable_mean_sd_all_ssps.csv"
  )
)

write_csv(
  cumulative_anomaly,
  file.path(
    out_dir,
    paste0(
      "cumulative_mineralised_bioavailable_anomaly_",
      reference_start,
      "_",
      reference_end,
      "_all_ssps.csv"
    )
  )
)

# ------------------------------------------------------------
# 10. Plot both variables
# ------------------------------------------------------------

p_cumulative <- ggplot(
  cumulative_anomaly,
  aes(
    x = Year,
    y = Cumulative_Anomaly,
    color = SSP,
    fill = SSP
  )
) +
  geom_ribbon(
    aes(
      ymin = Lower,
      ymax = Upper
    ),
    alpha = 0.15,
    color = NA
  ) +
  geom_line(
    linewidth = 0.9
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_vline(
    xintercept = c(
      reference_start,
      reference_end
    ),
    linetype = "dotted",
    linewidth = 0.4
  ) +
  facet_wrap(
    ~Variable,
    ncol = 1,
    scales = "free_y"
  ) +
  scale_color_manual(
    values = ssp_colors
  ) +
  scale_fill_manual(
    values = ssp_colors
  ) +
  scale_x_continuous(
    breaks = seq(
      1850,
      2100,
      by = 25
    )
  ) +
  labs(
    title = "Cumulative Arctic nitrogen flux anomalies",
    subtitle = paste0(
      "Anomalies relative to the mean cumulative value during ",
      reference_start,
      "–",
      reference_end,
      "; ribbons show propagated SD"
    ),
    x = "Year",
    y = expression(
      "Cumulative N anomaly [Pg N]"
    ),
    color = NULL,
    fill = NULL
  ) +
  theme_bw(
    base_size = 12
  ) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.text = element_text(
      face = "bold"
    ),
    plot.title = element_text(
      face = "bold"
    )
  )

print(p_cumulative)

ggsave(
  filename = file.path(
    out_dir,
    paste0(
      "cumulative_N_anomaly_",
      reference_start,
      "_",
      reference_end,
      "_mean_sd_all_ssps.png"
    )
  ),
  plot = p_cumulative,
  width = 10,
  height = 8,
  dpi = 300
)





library(terra)
library(dplyr)
library(ggplot2)
library(zoo)

# -----------------------------
# Paths and settings
# -----------------------------

mean_base <- "monthly_mineralised/whole_arctic_mean_new_8deg"
out_dir <- "cumulative_mineralised_bioavailable_plots_no_sd"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ssps <- c("126", "245", "370", "585")
years <- 1850:2099

common_extent <- ext(-180, 180, 45, 90)

ssp_labels <- c(
  "126" = "SSP1-2.6",
  "245" = "SSP2-4.5",
  "370" = "SSP3-7.0",
  "585" = "SSP5-8.5"
)

ssp_colors <- c(
  "SSP3-7.0" = "#D73027",
  "SSP5-8.5" = "#7B3294",
  "SSP2-4.5" = "orange",
  "SSP1-2.6" = "blue"
)

# -----------------------------
# File helper
# -----------------------------

make_file <- function(base_dir, var, ssp, year) {
  file.path(
    base_dir,
    paste0("_w_temp_w_sm_yearly_nc_", ssp),
    paste0("region_", var, "_N_monthly_", ssp, "_", year, "_w_temp_sm.nc")
  )
}

# -----------------------------
# Extract annual Pg N
# -----------------------------

extract_annual_pg <- function(var, ssp, years) {
  
  out <- list()
  
  for (yr in years) {
    
    mean_file <- make_file(mean_base, var, ssp, yr)
    
    if (!file.exists(mean_file)) {
      warning("Missing file: ", mean_file)
      next
    }
    
    r <- rast(mean_file)
    r <- crop(r, common_extent)
    
    area_r <- cellSize(r[[1]], unit = "m")
    area_r <- mask(area_r, r[[1]])
    
    monthly_pg <- global(
      r * area_r,
      "sum",
      na.rm = TRUE
    )[, 1] / 1e12
    
    annual_pg <- sum(monthly_pg, na.rm = TRUE)
    
    out[[as.character(yr)]] <- tibble(
      Year = yr,
      Mean = annual_pg,
      SSP = ssp_labels[as.character(ssp)],
      Variable = var
    )
    
    rm(r, area_r, monthly_pg)
    gc()
  }
  
  bind_rows(out)
}

# -----------------------------
# Process and plot variable
# -----------------------------

process_variable <- function(var) {
  
  df <- bind_rows(lapply(ssps, function(ssp) {
    cat("Processing", var, "SSP", ssp, "\n")
    extract_annual_pg(var = var, ssp = ssp, years = years)
  }))
  
  df_anom <- df %>%
    arrange(SSP, Year) %>%
    group_by(SSP) %>%
    mutate(
      Cum_Mean = cumsum(Mean),
      ref_mean = mean(Cum_Mean[Year >= 2000 & Year <= 2020], na.rm = TRUE),
      Cum_Mean_anom = Cum_Mean - ref_mean,
      Rolling_Mean_anom = zoo::rollapply(
        Cum_Mean_anom,
        width = 20,
        FUN = mean,
        align = "center",
        fill = NA,
        na.rm = TRUE
      )
    ) %>%
    ungroup() %>%
    mutate(
      SSP = factor(
        SSP,
        levels = c("SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5")
      ),
      Period = ifelse(Year < 2015, "Before 2015", "After 2015")
    )
  
  y_lab <- if (var == "bioavailable") {
    "Cumulative bioavailable N anomaly [Pg N]"
  } else {
    "Cumulative mineralised N anomaly [Pg N]"
  }
  
  p <- ggplot(df_anom, aes(x = Year, y = Cum_Mean_anom, group = SSP)) +
    
    geom_line(
      data = filter(df_anom, Period == "Before 2015"),
      color = "black",
      linewidth = 0.4,
      alpha = 0.5
    ) +
    
    geom_line(
      data = filter(df_anom, Period == "After 2015"),
      aes(color = SSP),
      linewidth = 0.4,
      alpha = 0.5
    ) +
    
    geom_line(
      data = filter(df_anom, Year < 2015),
      aes(y = Rolling_Mean_anom),
      color = "black",
      linewidth = 0.8,
      na.rm = TRUE
    ) +
    
    geom_line(
      data = filter(df_anom, Year >= 2015),
      aes(y = Rolling_Mean_anom, color = SSP),
      linewidth = 0.8,
      na.rm = TRUE
    ) +
    
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "grey30"
    ) +
    
    geom_vline(
      xintercept = 2015,
      linetype = "dashed",
      color = "black"
    ) +
    
    scale_color_manual(values = ssp_colors, drop = FALSE) +
    
    theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 10)
    ) +
    
    labs(
      x = "Year",
      y = y_lab,
      color = "SSP Scenario"
    )
  
  write.csv(
    df_anom,
    file.path(out_dir, paste0("cumulative_", var, "_N_anomaly_no_sd.csv")),
    row.names = FALSE
  )
  
  ggsave(
    file.path(out_dir, paste0("cumulative_", var, "_N_anomaly_no_sd.png")),
    p,
    width = 8,
    height = 5,
    dpi = 300
  )
  
  return(p)
}

# -----------------------------
# Run
# -----------------------------

plot_bioavailable <- process_variable("bioavailable")
plot_mineralised <- process_variable("mineralised")

print(plot_bioavailable)
print(plot_mineralised)

plot_bioavailable
plot_mineralised


###################################
# mineralised files: 

library(terra)
library(dplyr)
library(ggplot2)
library(zoo)
library(tibble)

ssps <- c("126", "245", "370", "585")
years <- 1850:2099

common_extent <- ext(-179.95, 179.95, 60, 90)

out_dir <- "mineralised_N_plots"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

base_dir <- "monthly_mineralised/whole_arctic_mean"

ssp_colors <- c(
  "SSP3-7.0" = "#D73027",
  "SSP5-8.5" = "#7B3294",
  "SSP2-4.5" = "orange",
  "SSP1-2.6" = "blue"
)

extract_bioavailable_annual_pg <- function(ssp, years, base_dir, common_extent) {
  
  folder <- file.path(
    base_dir,
    paste0("_w_temp_w_sm_yearly_nc_", ssp)
  )
  
  out <- vector("list", length(years))
  
  for (ii in seq_along(years)) {
    
    yr <- years[ii]
    
    f <- file.path(
      folder,
      paste0(
        "arctic_mineralised_N_monthly_",
        ssp, "_", yr, "_w_temp_sm.nc"
      )
    )
    
    if (!file.exists(f)) {
      warning("Missing file: ", f)
      out[[ii]] <- tibble(
        Year = yr,
        Mean = NA_real_,
        SSP = paste0("SSP", ssp)
      )
      next
    }
    
    r <- rast(f)
    r <- crop(r, common_extent)
    
    area_rast <- cellSize(r[[1]], unit = "m")
    
    # yearly total = sum over 12 monthly bioavailable N layers
    r_year <- app(r, sum, na.rm = TRUE)
    
    total_pg <- global(
      r_year * area_rast,
      "sum",
      na.rm = TRUE
    )[1, 1] / 1e12
    
    out[[ii]] <- tibble(
      Year = yr,
      Mean = total_pg,
      SSP = paste0("SSP", ssp)
    )
    
    rm(r, r_year, area_rast)
    gc()
  }
  
  bind_rows(out)
}

df_list <- list()

for (ssp in ssps) {
  cat("Processing SSP", ssp, "\n")
  
  df_list[[ssp]] <- extract_bioavailable_annual_pg(
    ssp = ssp,
    years = years,
    base_dir = base_dir,
    common_extent = common_extent
  )
}

df_all <- bind_rows(df_list)

df_all <- df_all %>%
  mutate(
    SSP = recode(
      SSP,
      "SSP126" = "SSP1-2.6",
      "SSP245" = "SSP2-4.5",
      "SSP370" = "SSP3-7.0",
      "SSP585" = "SSP5-8.5"
    ),
    Period = ifelse(Year < 2015, "Before 2015", "After 2015")
  )

library(tidyr)

df_all <- df_all %>%
  arrange(SSP, Year) %>%
  group_by(SSP) %>%
  mutate(
    Cumulative_Pg = cumsum(replace_na(Mean, 0))
  ) %>%
  ungroup()

compute_cumulative_anomaly <- function(df) {
  
  ref_mean <- df %>%
    filter(Year >= 2000, Year <= 2020) %>%
    summarise(ref = mean(Cumulative_Pg, na.rm = TRUE)) %>%
    pull(ref)
  
  df %>%
    mutate(
      Cumulative_anom = Cumulative_Pg - ref_mean
    )
}

df_all <- df_all %>%
  group_by(SSP) %>%
  group_modify(~ compute_cumulative_anomaly(.x)) %>%
  ungroup()

df_all <- df_all %>%
  arrange(SSP, Year) %>%
  group_by(SSP) %>%
  mutate(
    Rolling_Cumulative_anom = rollapply(
      Cumulative_anom,
      width = 20,
      FUN = mean,
      align = "center",
      fill = NA,
      na.rm = TRUE
    )
  ) %>%
  ungroup()

write.csv(
  df_all,
  file.path(out_dir, "bioavailable_N_annual_rolling20_all_ssps.csv"),
  row.names = FALSE
)

plot_bioavailable_cumulative <- ggplot(
  df_all,
  aes(x = Year, y = Cumulative_anom, group = SSP)
) +
  geom_line(
    data = filter(df_all, Period == "Before 2015"),
    color = "black",
    linewidth = 0.4,
    alpha = 0.5
  ) +
  geom_line(
    data = filter(df_all, Period == "After 2015"),
    aes(color = SSP),
    linewidth = 0.4,
    alpha = 0.5
  ) +
  geom_line(
    data = filter(df_all, Year < 2015),
    aes(y = Rolling_Cumulative_anom),
    color = "black",
    linewidth = 0.8
  ) +
  geom_line(
    data = filter(df_all, Year >= 2015),
    aes(y = Rolling_Cumulative_anom, color = SSP),
    linewidth = 0.8
  ) +
  geom_vline(xintercept = 2015, linetype = "dashed", color = "black") +
  scale_color_manual(values = ssp_colors) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(
    x = "Year",
    y = "Cumulative mineralised N (Pg N)",
    color = "SSP Scenario"
  )

plot_bioavailable_cumulative

plot_bioavailable

ggsave(
  file.path(out_dir, "bioavailable_N_annual_rolling20_all_ssps.png"),
  plot_bioavailable,
  width = 8,
  height = 5,
  dpi = 300
)

# ============================================================
# Spatial plot of bioavailable N from monthly dynamic-depth files
# SSP370, 1850–2099
# ============================================================

library(terra)
library(dplyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(scales)

ssp <- "370"

in_dir <- "monthly_mineralised/mean_2perc_baseline/yearly_nc_370"


out_dir <- "monthly_mineralised/spatial_plots"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

years <- 1850:2099

# Choose period to plot
plot_years <- 2080:2099
# plot_years <- 2000:2020

# ============================================================
# 1. Read monthly files and calculate annual/period mean
# ============================================================

read_year_bio <- function(year) {
  
  f <- file.path(
    in_dir,
    paste0(
      "arctic_mineralised_N_monthly_w_temp_sm_",
      ssp, "_", year, ".nc"
    )
  )
  
  if (!file.exists(f)) {
    stop("Missing file: ", f)
  }
  
  r <- rast(f)
  
  # Sum 12 months to annual mineralised N
  # input: kg N m-2 month-1
  # output: kg N m-2 yr-1
  annual <- app(r, sum, na.rm = TRUE)
  
  names(annual) <- paste0("mineralised_", year)
  
  annual
}

bio_list <- lapply(plot_years, read_year_bio)

bio_stack <- rast(bio_list)

# Mean annual mineralised N over selected period
rate_per_year <- app(bio_stack, mean, na.rm = TRUE)

# Convert kg N m-2 yr-1 to g N m-2 yr-1
rate_per_year <- rate_per_year * 1000

# Remove negative values
#rate_per_year <- terra::ifel(rate_per_year < 0, 0, rate_per_year)

plot(rate_per_year)

# ============================================================
# 2. Convert raster to sf
# ============================================================

df_min_n_future <- as.data.frame(rate_per_year, xy = TRUE, na.rm = FALSE)
names(df_min_n_future)[3] <- "mineralised_N"

df_min_n_future <- df_min_n_future %>%
  filter(!is.na(mineralised_N))

arctic_sf <- st_as_sf(
  df_min_n_future,
  coords = c("x", "y"),
  crs = 4326
)


# North Pole LAEA projection
laea_crs <- "+proj=laea +lat_0=90 +lon_0=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

arctic_sf_laea <- st_transform(arctic_sf, crs = laea_crs)

arctic_sf_laea$mineralised_N <- pmax(arctic_sf_laea$mineralised_N, 0)

summary(arctic_sf_laea$mineralised_N)
bbox <- st_bbox(arctic_sf_laea)
print(bbox)

sum(arctic_sf_laea$mineralised_N > 0)

coastlines <- ne_coastline(scale = "medium", returnclass = "sf")
coastlines_laea <- st_transform(coastlines, crs = laea_crs)

# Plot limits
xlim <- c(-3294699, 3294699)
ylim <- c(-3076617, 3300207)

# ============================================================
# 3. Plot
# ============================================================

spatial_plot <- ggplot() +
  geom_sf(
    data = arctic_sf_laea,
    aes(color = mineralised_N),
    size = 0.01
  ) +scale_color_gradientn(
    colours = c("dodgerblue3", "#abd9e9", "#ffffbf", "orange", "darkred"),
    values = scales::rescale(c(0, 0.37, 1.27, 1.6357, 20)),   # ~min, Q1, median, Q3, then tail
    limits = c(0, 20),
    breaks = c(0, 0.37, 1.27, 1.6357, 20),
    labels = c("0", "0.37", "1.27", "1.6357", "20+"),
    oob = scales::squish,
    name = expression(N~(g~m^{-2}~yr^{-1})),
    guide = guide_colorbar(barheight = unit(2.5, "cm"), barwidth = unit(0.25, "cm"), ticks = TRUE)
  )+
  geom_sf(
    data = coastlines_laea,
    color = "black",
    linewidth = 0.05
  ) +
  coord_sf(
    crs = laea_crs,
    xlim = xlim,
    ylim = ylim
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(
    title = paste0("Gross mineralised N flux, SSP", ssp)
    #subtitle = paste0(min(plot_years), "–", max(plot_years), "")
  )
spatial_plot

spatial_plot <- ggplot() +
  geom_sf(data = arctic_sf_laea, aes(color = mineralised_N), size = 0.01) +
  scale_color_gradientn(
    colours = c("dodgerblue3", "#abd9e9", "#ffffbf", "orange", "darkred"),
    trans = "sqrt",             
    limits = c(0, 20),
    oob = scales::squish,
    name = expression(N~(g~m^{-2}~yr^{-1})),
    guide = guide_colorbar(barheight = unit(2.5, "cm"), barwidth = unit(0.25, "cm"), ticks = TRUE)
  ) +
  geom_sf(data = coastlines_laea, color = "black", linewidth = 0.05) +
  coord_sf(crs = laea_crs, xlim = xlim, ylim = ylim) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(title = paste0("Gross mineralised N flux, SSP", ssp))

spatial_plot

ggsave(
  file.path(
    out_dir,
    paste0(
      "spatial_bioavailable_N_SSP",
      ssp, "_",
      min(plot_years), "_", max(plot_years),
      "_mean.png"
    )
  ),
  spatial_plot,
  width = 6,
  height = 6,
  dpi = 300
)





# plot Palmtag dataset: 

palmtag<-rast("TN_30deg_corr.nc")
palmtag<-rast("total_thawed_extended/arctic_total_thawed_585_no_lim_mean.nc")
palmtag<-rast("monthly_mineralised/arctic_final_mean/yearly_nc_245/arctic_mineralised_N_monthly_w_temp_sm_245_2099.nc")
permafrost<-rast("ESACCI-PERMAFROST-L4-PFR-MODISLST_CRYOGRID-AREA4_PP-2023-fv05.0.nc")
df_plot <- function(r, period) {
  df <- as.data.frame(r, xy = TRUE)
  names(df)[3] <- "TN"
  df$Period <- period
  return(df)
}

df_thawedN<-df_plot(permafrost, "Total")

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
xlim <- c(-4539074, 4539080)
ylim <- c(-2989335, 4006331)
# Adjust these values based on your data
coastlines <- ne_coastline(scale = "medium", returnclass = "sf")
summary(arctic_sf_laea$TN)
# Step 5: Reproject the coastlines to match the LAEA projection
coastlines_laea <- st_transform(coastlines, crs = laea_crs)
arctic_sf_laea <- arctic_sf_laea[!is.na(arctic_sf_laea$TN), ]
# Step 6: Create the plot with the LAEA projection
ggplot() +
  geom_sf(data = arctic_sf_laea, aes(color = TN), size = 0.01) +
  scale_color_viridis_c(
    option = "viridis",   # or "magma", "plasma", "inferno", "cividis"
    name = expression(N~(kg~m^{-2}))
  ) +
  geom_sf(data = coastlines_laea, color = "black", size = 0.05) +
  coord_sf(crs = laea_crs, xlim = xlim, ylim = ylim) +
  theme_minimal() +
  theme(
    legend.position = "right"
  ) +
  labs(x = "", y = "")





library(terra)
library(sf)
library(ggplot2)
library(viridis)
library(rnaturalearth)
library(scico)
install.packages("scico")
# ------------------------------------------------------------
# 1. Read the raster
# ------------------------------------------------------------

permafrost <- rast(
  "ESACCI-PERMAFROST-L4-PFR-MODISLST_CRYOGRID-AREA4_PP-2023-fv05.0.nc"
)

# If the file contains several layers, choose the required one:
 permafrost <- permafrost[[1]]

names(permafrost) <- "PFR"

# Remove missing cells before processing
permafrost <- mask(permafrost, permafrost)

# ------------------------------------------------------------
# 2. Define Arctic LAEA projection
# ------------------------------------------------------------

laea_crs <- paste(
  "+proj=laea",
  "+lat_0=90",
  "+lon_0=30",
  "+x_0=0",
  "+y_0=0",
  "+datum=WGS84",
  "+units=m",
  "+no_defs"
)

# ------------------------------------------------------------
# 3. Reproject the raster directly
# ------------------------------------------------------------

permafrost_laea <- project(
  permafrost,
  laea_crs,
  method = "near"
)

permafrost_quick <- aggregate(
  permafrost,
  fact = 4,
  fun = mean,
  na.rm = TRUE
)

permafrost_laea <- project(
  permafrost_quick,
  laea_crs,
  method = "near"
)

# Convert the projected raster to an ordinary data frame,
# not an sf point object
permafrost_df <- as.data.frame(
  permafrost_laea,
  xy = TRUE,
  na.rm = TRUE
)

names(permafrost_df)[3] <- "PFR"


# ------------------------------------------------------------
# 4. Coastlines
# ------------------------------------------------------------

coastlines <- ne_coastline(
  scale = "medium",
  returnclass = "sf"
)

coastlines_laea <- st_transform(
  coastlines,
  crs = laea_crs
)

# ------------------------------------------------------------
# 5. Plot
# ------------------------------------------------------------

xlim <- c(-4539074, 4539080)
ylim <- c(-2989335, 4006331)

ggplot() +
  geom_tile(
    data = permafrost_df,
    aes(x = x, y = y, fill = PFR)
  ) +
  geom_sf(
    data = coastlines_laea,
    colour = "black",
    linewidth = 0.15,
    fill = NA,
    inherit.aes = FALSE
  ) +
  scale_fill_viridis_c(
    option = "mako",
    direction = -1,
    name = "Permafrost fraction",
    na.value = "transparent"
  ) +
  coord_sf(
    crs = laea_crs,
    xlim = xlim,
    ylim = ylim,
    expand = FALSE
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

ggplot() +
  geom_tile(
    data = permafrost_df,
    aes(x = x, y = y, fill = PFR)
  ) +
  geom_sf(
    data = coastlines_laea,
    colour = "black",
    linewidth = 0.15,
    fill = NA,
    inherit.aes = FALSE
  ) +
  scale_fill_scico(
    palette = "bilbao",
    limits = c(0, 100),
    direction = -1,
    name = "Permafrost fraction"
  ) +
  coord_sf(
    crs = laea_crs,
    xlim = xlim,
    ylim = ylim,
    expand = FALSE
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )


library(terra)
library(dplyr)
library(ggplot2)

ssps <- c("126", "245", "370", "585")

all_df <- list()

for (ssp in ssps) {
  
  thawed <- rast(
    paste0(
      "total_thawed_extended/arctic_total_thawed_",
      ssp,
      "_no_lim_mean.nc"
    )
  )
  
  years <- as.integer(format(time(thawed), "%Y"))
  
  baseline <- app(
    thawed[[years >= 1850 & years <= 1900]],
    mean,
    na.rm = TRUE
  )
  
  thawed_pf <- thawed - baseline
  #thawed_pf <- ifel(thawed_pf < 0, 0, thawed_pf)
  
  area <- cellSize(thawed_pf[[1]], unit = "m")
  
  total_pg <- global(
    thawed_pf * area,
    "sum",
    na.rm = TRUE
  )[,1] / 1e12
  
  all_df[[ssp]] <- data.frame(
    Year = years,
    SSP = paste0("SSP", substr(ssp,1,1), "-", substr(ssp,2,3)),
    thawed_pg = total_pg
  )
}

thawed_df <- bind_rows(all_df)

ggplot(thawed_df,
       aes(Year, thawed_pg, colour = SSP)) +
  geom_line(linewidth = 1.2) +
  geom_vline(xintercept = 2015,
             linetype = "dashed") +
  theme_bw(base_size = 16) +
  labs(
    x = "Year",
    y = "Baseline-corrected thawed permafrost N (Pg N)",
    colour = "SSP scenario"
  )




