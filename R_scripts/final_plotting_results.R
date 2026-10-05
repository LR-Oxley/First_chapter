# R script for plotting graphs

###### ------------------ figure 1 ---------------------- ######
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
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




changes <- df_all %>%
  group_by(SSP) %>%
  summarise(
    mean_pi = mean(Mean[Year >= 1850 & Year <= 1900], na.rm = TRUE),
    mean_pd = mean(Mean[Year >= 2000 & Year <= 2020], na.rm = TRUE),
    mean_fu = mean(Mean[Year >= 2080 & Year <= 2099], na.rm = TRUE),
    sd_pi   = mean(SD[Year >= 1850 & Year <= 1900], na.rm = TRUE),
    sd_pd   = mean(SD[Year >= 2000 & Year <= 2020], na.rm = TRUE),
    sd_fu   = mean(SD[Year >= 2080 & Year <= 2099], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # historical: present day (2000-2020) minus preindustrial (1850-1900)
    hist_change = mean_pd - mean_pi,
    hist_sd     = sd_pd - sd_pi,
    # future: 2080-2099 minus present day (2000-2020)
    fut_change  = mean_fu - mean_pd,
    fut_sd      = sd_fu - sd_pd,
    hist_text   = sprintf("%.1f (\u00b1 %.1f) Pg N", hist_change, hist_sd),
    fut_text    = sprintf("%.1f (\u00b1 %.1f) Pg N", fut_change, fut_sd)
  )

print(changes %>% select(SSP, hist_text, fut_text))

print(changes, width = Inf)

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
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 7),
    plot.tag = element_text(size = 8), 
    legend.position = "bottom"
  )+
  labs(
    x = "Year",
    y = "Total N thaw [Pg N]",
    color = "",
    fill = "",
    title= ""
  )

plot_total_thawed <- plot_total_thawed +
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
  "plot_total_thawed.pdf",
  plot_total_thawed,
  width = 8.5,
  height = 6,
  units = "cm", 
  dpi = 500
)




# ============================================================
# Historical increase: present-day (2000-2020) vs preindustrial (1850-1900)
# ============================================================

historical_increase <- df_all_anom %>%
  filter(Year >= 2000, Year <= 2020) %>%
  group_by(SSP) %>%
  summarise(
    Mean_increase_Pg = mean(Mean_baseline, na.rm = TRUE),
    SD_increase_Pg   = sqrt(mean(SD_baseline^2, na.rm = TRUE)),  # RMS pooling
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("%.2f \u00B1 %.2f Pg N", Mean_increase_Pg, SD_increase_Pg))

print(historical_increase)

# ============================================================
# Future increase: 2080-2099 relative to present day (2000-2020), by SSP
# ============================================================

future_increase <- df_all_anom %>%
  filter(Year >= 2080, Year <= 2099, Period == "After 2015") %>%
  group_by(SSP) %>%
  summarise(
    Mean_increase_Pg = mean(Mean_anom, na.rm = TRUE),
    SD_increase_Pg   = sqrt(mean(SD_anom^2, na.rm = TRUE)),  # RMS pooling
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("%.2f \u00B1 %.2f Pg N", Mean_increase_Pg, SD_increase_Pg))

print(future_increase)


# ============================================================
# Period-average version (e.g., mean over 2080-2099, matching
# your "future" period definition elsewhere in the analysis)
# ============================================================

summary_future_period <- df_all_anom %>%
  filter(Year >= 2080 & Year <= 2099, Period == "After 2015") %>%
  group_by(SSP) %>%
  summarise(
    period_mean_ALD = mean(Mean_anom, na.rm = TRUE),
    period_mean_SD  = mean(SD_anom, na.rm = TRUE),   # average of yearly inter-model SDs
    .groups = "drop"
  ) %>%
  mutate(
    label = sprintf("%.2f ± %.2f m", period_mean_ALD, period_mean_SD)
  )

print(summary_future_period)

########################################################################




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
# Plot
ALD<-ggplot(long_data, aes(x = Year, y = Mean_ALD, group = SSP)) +
  
  # Before 2015: Black lines & grey ribbons
  geom_line(data = filter(long_data, Period == "Before 2015"), color = "black", linewidth = 0.3) +
  #geom_ribbon(data = filter(long_data, Period == "Before 2015"),
  #            aes(ymin = Mean_ALD - STD_ALD, ymax = Mean_ALD + STD_ALD),
  #            fill = "grey", alpha = 0.2) +
  
  # After 2015: Colored lines & ribbons
  geom_line(data = filter(long_data, Period == "After 2015"),
            aes(color = SSP), linewidth = 0.15) +
  #geom_ribbon(data = filter(long_data, Period == "After 2015"),
  #            aes(ymin = Mean_ALD - STD_ALD, ymax = Mean_ALD + STD_ALD, fill = SSP),
  #            alpha = 0.2) +
  
  # Rolling Mean before 2015 (Black line)
  geom_line(data = filter(long_data, Year < 2015),
            aes(y = Rolling_Mean_ALD), color = "black", linewidth = 0.15, linetype = "solid") +
  
  # Rolling Mean after 2015 (Colored lines)
  geom_line(data = filter(long_data, Year >= 2015),
            aes(y = Rolling_Mean_ALD, color = SSP), linewidth = 0.15, linetype = "solid") +
  
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
              alpha = 0.15) +
  
  # Vertical line for 2015
  geom_vline(xintercept = 2015, linetype = "dashed", color = "black") +
  
  # Vertical line for 2015
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  
  # Labels and title
  labs(x = "Year", y = "Active layer depth [m]", title = "", 
       color = "", fill = "") +
  
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
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 7),
    plot.tag = element_text(size = 8), 
    legend.position = "bottom"
  )

ALD



plot_ALD <- ALD +
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
  "plot_ALD.pdf",
  plot_ALD,
  width = 8.5,
  height = 6,
  units = "cm", 
  dpi = 500
)









# =========================================================
# Paths
# =========================================================

mean_dir <- "monthly_mineralised/mean_2perc_baseline_new"
std_dir  <- "monthly_mineralised/std_2perc_baseline_redone"   
out_dir <- file.path(mean_dir, "diagnostic_csv_analysis_60N")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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

# Columns we expect to carry an SD counterpart (edit list if needed)
sd_value_cols <- c(
  "mineralised_pg_monthly",
  "total_inorg_pg_monthly",
  "inorg_rapid_available",
  "k_T_mean_monthly",
  "k_T_pool_weighted_mean",
  "k_env_pool_weighted_mean",
  "organic_pool_remaining_pg",
  "thawed_N_annual_change_pg"
)

# =========================================================
# Readers
# =========================================================

read_diag_csv <- function(dir, ssp) {
  
  f <- file.path(
    dir,
    paste0("arctic_monthly_diagnostics_w_temp", ssp, "_1850_2100.csv")
  )
  
  if (!file.exists(f)) {
    warning("Missing file: ", f)
    return(NULL)
  }
  
  read_csv(f, show_col_types = FALSE) %>%
    mutate(
      SSP_raw = ssp,
      SSP = ssp_labels[ssp]
    )
}

# --- Mean (existing) diagnostics ---
diag_monthly_mean <- bind_rows(lapply(ssps, function(s) read_diag_csv(mean_dir, s))) %>%
  mutate(
    SSP = factor(SSP, levels = c("SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5"))
  )

# --- Std diagnostics ---
diag_monthly_sd_raw <- bind_rows(lapply(ssps, function(s) read_diag_csv(std_dir, s))) %>%
  mutate(
    SSP = factor(SSP, levels = c("SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5"))
  )

# Rename the value columns in the SD table to *_sd so they don't collide on join
sd_cols_present <- intersect(sd_value_cols, names(diag_monthly_sd_raw))

diag_monthly_sd <- diag_monthly_sd_raw %>%
  select(SSP_raw, SSP, Year, any_of(c("Month", sd_cols_present))) %>%
  rename_with(~ paste0(.x, "_sd"), all_of(sd_cols_present))

# Join key: Year (+ Month if present in both)
join_keys <- intersect(c("SSP_raw", "SSP", "Year", "Month"), names(diag_monthly_mean))
join_keys <- intersect(join_keys, names(diag_monthly_sd))

diag_monthly <- diag_monthly_mean %>%
  left_join(diag_monthly_sd, by = join_keys)

# Check column names
print(names(diag_monthly))
print(head(diag_monthly))

# =========================================================
# Annual summaries (mean + SD)
# =========================================================
# Matching the reference approach exactly: the SD is aggregated with
# the SAME function used for the corresponding mean quantity (sum -> sum,
# mean -> mean). No quadrature/error-propagation anywhere.

diag_annual <- diag_monthly %>%
  group_by(SSP, Year) %>%
  summarise(
    mineralised_pg_yr = sum(mineralised_pg_monthly, na.rm = TRUE),
    total_inorg_pg_yr = sum(total_inorg_pg_monthly, na.rm = TRUE),
    inorg_rapid_pg_yr = sum(inorg_rapid_available, na.rm = TRUE),
    mean_k_t = mean(k_T_mean_monthly, na.rm = TRUE),
    mean_k_t_weighted = mean(k_T_pool_weighted_mean, na.rm = TRUE),
    organic_pool_remaining_pg = mean(organic_pool_remaining_pg, na.rm = TRUE),
    thawed_N_annual_change_pg = mean(thawed_N_annual_change_pg, na.rm = TRUE),
    
    # ---- SD: same aggregation function as its mean counterpart above ----
    mineralised_pg_yr_sd = sum(mineralised_pg_monthly_sd, na.rm = TRUE),
    total_inorg_pg_yr_sd = sum(total_inorg_pg_monthly_sd, na.rm = TRUE),
    inorg_rapid_pg_yr_sd = sum(inorg_rapid_available_sd, na.rm = TRUE),
    mean_k_t_sd = mean(k_T_mean_monthly_sd, na.rm = TRUE),
    mean_k_t_weighted_sd = mean(k_T_pool_weighted_mean_sd, na.rm = TRUE),
    organic_pool_remaining_pg_sd = mean(organic_pool_remaining_pg_sd, na.rm = TRUE),
    thawed_N_annual_change_pg_sd = mean(thawed_N_annual_change_pg_sd, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  group_by(SSP) %>%
  arrange(Year, .by_group = TRUE) %>%
  mutate(
    # Cumulative annual thawed-N change
    cumulative_thawed_N_pg = cumsum(replace_na(thawed_N_annual_change_pg, 0)),
    
    cumulative_mineralised_pg = cumsum(mineralised_pg_yr),
    cumulative_inorg_pg = cumsum(total_inorg_pg_yr),
    cumulative_inorg_rapid_pg = cumsum(inorg_rapid_pg_yr),
    # Cumulative SD: same function (cumsum) as used for the cumulative mean
    cumulative_mineralised_pg_sd = cumsum(mineralised_pg_yr_sd),
    cumulative_inorg_pg_sd = cumsum(total_inorg_pg_yr_sd),
    cumulative_inorg_rapid_pg_sd = cumsum(inorg_rapid_pg_yr_sd),
    mineralised_20yr = zoo::rollapply(mineralised_pg_yr, 20, mean, align = "center", fill = NA),
    inorg_20yr = zoo::rollapply(total_inorg_pg_yr, 20, mean, align = "center", fill = NA),
    
    mineralised_20yr_sd = zoo::rollapply(mineralised_pg_yr_sd, 20, mean, align = "center", fill = NA),
    inorg_20yr_sd = zoo::rollapply(total_inorg_pg_yr_sd, 20, mean, align = "center", fill = NA)
  ) %>%
  ungroup()

# =========================================================
# Anomaly relative to present-day (2000-2020) reference period
# =========================================================
# Single-stage: subtract the 2000-2020 reference-period mean directly
# from the raw series. Applied identically to Mean and SD (SD is just
# walked through the same subtraction -- no quadrature, no baseline
# re-centering step).

compute_anomaly_pair <- function(df, value_col, sd_col,
                                 ref_range = c(2000, 2020)) {
  
  ref <- df %>% filter(Year >= ref_range[1], Year <= ref_range[2])
  
  ref_mean <- mean(ref[[value_col]], na.rm = TRUE)
  ref_sd   <- mean(ref[[sd_col]], na.rm = TRUE)
  
  df %>%
    mutate(
      "{value_col}_anom" := .data[[value_col]] - ref_mean,
      "{sd_col}_anom"     := .data[[sd_col]] - ref_sd
    )
}

diag_annual <- diag_annual %>%
  group_by(SSP) %>%
  group_modify(~ compute_anomaly_pair(.x, "mineralised_pg_yr", "mineralised_pg_yr_sd")) %>%
  group_modify(~ compute_anomaly_pair(.x, "total_inorg_pg_yr", "total_inorg_pg_yr_sd")) %>%
  group_modify(~ compute_anomaly_pair(.x, "inorg_rapid_pg_yr", "inorg_rapid_pg_yr_sd")) %>%
  group_modify(~ compute_anomaly_pair(.x, "organic_pool_remaining_pg", "organic_pool_remaining_pg_sd")) %>%
  group_modify(~ compute_anomaly_pair(.x, "cumulative_mineralised_pg", "cumulative_mineralised_pg_sd")) %>%
  group_modify(~ compute_anomaly_pair(.x, "cumulative_inorg_pg", "cumulative_inorg_pg_sd")) %>%
  group_modify(~ compute_anomaly_pair(.x, "cumulative_inorg_rapid_pg", "cumulative_inorg_rapid_pg_sd")) %>%
  ungroup() %>%
  rename(
    annual_mineralised_anom_pg       = mineralised_pg_yr_anom,
    annual_mineralised_anom_pg_sd    = mineralised_pg_yr_sd_anom,
    
    annual_inorg_anom_pg             = total_inorg_pg_yr_anom,
    annual_inorg_anom_pg_sd          = total_inorg_pg_yr_sd_anom,
    
    annual_rapid_inorg_anom_pg       = inorg_rapid_pg_yr_anom,
    annual_rapid_inorg_anom_pg_sd    = inorg_rapid_pg_yr_sd_anom,
    
    organic_pool_anom_pg             = organic_pool_remaining_pg_anom,
    organic_pool_anom_pg_sd          = organic_pool_remaining_pg_sd_anom,
    
    cumulative_mineralised_anom_pg    = cumulative_mineralised_pg_anom,
    cumulative_mineralised_anom_pg_sd = cumulative_mineralised_pg_sd_anom,
    
    cumulative_inorg_anom_pg          = cumulative_inorg_pg_anom,
    cumulative_inorg_anom_pg_sd       = cumulative_inorg_pg_sd_anom,
    
    cumulative_rapid_inorg_anom_pg    = cumulative_inorg_rapid_pg_anom,
    cumulative_rapid_inorg_anom_pg_sd = cumulative_inorg_rapid_pg_sd_anom
  )

# Rolling 20-year mean of the anomalies (plain mean, same as the reference)
# =========================================================
# 20-year running means of anomalies
# =========================================================

diag_annual <- diag_annual %>%
  group_by(SSP) %>%
  arrange(Year, .by_group = TRUE) %>%
  mutate(
    
    #  mineralised N
    mineralised_pg_yr_anom_20yr =
      zoo::rollapply(
        annual_mineralised_anom_pg,
        20,
        mean,
        align = "center",
        fill = NA
      ),
    
    mineralised_pg_yr_anom_sd_20yr =
      zoo::rollapply(
        annual_mineralised_anom_pg_sd,
        20,
        mean,
        align = "center",
        fill = NA
      ),
    
    # Cumulative mineralised N
    cumulative_mineralised_anom_20yr =
      zoo::rollapply(
        cumulative_mineralised_anom_pg,
        20,
        mean,
        align = "center",
        fill = NA
      ),
    
    cumulative_mineralised_anom_20yr_sd =
      zoo::rollapply(
        cumulative_mineralised_anom_pg_sd,
        20,
        mean,
        align = "center",
        fill = NA
      ),
    
    # Cumulative inorganic N
    cumulative_inorg_anom_20yr =
      zoo::rollapply(
        cumulative_inorg_anom_pg,
        20,
        mean,
        align = "center",
        fill = NA
      ),
    
    cumulative_inorg_anom_20yr_sd =
      zoo::rollapply(
        cumulative_inorg_anom_pg_sd,
        20,
        mean,
        align = "center",
        fill = NA
      ),
    
    # Organic N pool
    organic_pool_anom_20yr =
      zoo::rollapply(
        organic_pool_anom_pg,
        20,
        mean,
        align = "center",
        fill = NA
      ),
    
    organic_pool_anom_20yr_sd =
      zoo::rollapply(
        organic_pool_anom_pg_sd,
        20,
        mean,
        align = "center",
        fill = NA
      ),
    # Annual rapid inorganic N
    rapid_inorg_pg_yr_anom_20yr =
      zoo::rollapply(
        annual_rapid_inorg_anom_pg,
        20,
        mean,
        align = "center",
        fill = NA
      ),
    
    rapid_inorg_pg_yr_anom_sd_20yr =
      zoo::rollapply(
        annual_rapid_inorg_anom_pg_sd,
        20,
        mean,
        align = "center",
        fill = NA
      ),
    
    # Cumulative rapid inorganic N
    cumulative_rapid_inorg_anom_20yr =
      zoo::rollapply(
        cumulative_rapid_inorg_anom_pg,
        20,
        mean,
        align = "center",
        fill = NA
      ),
    
    cumulative_rapid_inorg_anom_20yr_sd =
      zoo::rollapply(
        cumulative_rapid_inorg_anom_pg_sd,
        20,
        mean,
        align = "center",
        fill = NA
      )
  ) %>%
  ungroup()
# =========================================================
# Plot function
# Historical (1850-2014) = black + grey ribbon
# Future (2015-2099)      = SSP colours + coloured ribbons
# =========================================================

plot_diag <- function(df,
                      yvar,
                      ylab,
                      title,
                      filename,
                      sd_var = NULL,
                      running_var = NULL,
                      running_sd_var = NULL) {
  
  df <- df %>%
    arrange(SSP, Year)
  
  # Historical period
  df_hist <- df %>%
    filter(Year <= 2014)
  
  # Future period
  df_future <- df %>%
    filter(Year >= 2015)
  
  # SSP used to represent the shared historical period
  historical_ssp <- unique(df$SSP)[1]
  
  df_hist <- df_hist %>%
    filter(SSP == historical_ssp)
  
  # =======================================================
  # Base plot
  # =======================================================
  
  p <- ggplot(df, aes(x = Year, group = SSP))
  
  # =======================================================
  # RAW VALUES — HISTORICAL
  # =======================================================
  
  p <- p +
    geom_line(
      data = df_hist,
      aes(y = .data[[yvar]]),
      color = "black",
      linewidth = 0.2,
      alpha = 0.45,
      na.rm = TRUE
    )
  
  # =======================================================
  # RAW VALUES — FUTURE
  # =======================================================
  
  p <- p +
    geom_line(
      data = df_future,
      aes(
        y = .data[[yvar]],
        color = SSP
      ),
      linewidth = 0.2,
      alpha = 0.45,
      na.rm = TRUE
    )
  
  # =======================================================
  # RAW SD RIBBON — OPTIONAL
  # =======================================================
  
  if (!is.null(sd_var) && sd_var %in% names(df)) {

    p <- p +
      geom_ribbon(
        data = df_hist,
        aes(
          ymin = .data[[yvar]] - .data[[sd_var]],
          ymax = .data[[yvar]] + .data[[sd_var]]
        ),
        fill = "grey",
        color = NA,
        alpha = 0.03
      ) +

      geom_ribbon(
        data = df_future,
        aes(
          ymin = .data[[yvar]] - .data[[sd_var]],
          ymax = .data[[yvar]] + .data[[sd_var]],
          fill = SSP
        ),
        color = NA,
        alpha = 0.03
      )
  }
  
  # =======================================================
  # RUNNING MEAN — HISTORICAL
  # =======================================================
  
  if (!is.null(running_var) && running_var %in% names(df)) {
    
    p <- p +
      geom_line(
        data = df_hist,
        aes(y = .data[[running_var]]),
        color = "black",
        linewidth = 0.2,
        na.rm = TRUE
      )
  }
  
  # =======================================================
  # RUNNING MEAN — FUTURE
  # =======================================================
  
  if (!is.null(running_var) && running_var %in% names(df)) {
    
    p <- p +
      geom_line(
        data = df_future,
        aes(
          y = .data[[running_var]],
          color = SSP
        ),
        linewidth = 0.2,
        na.rm = TRUE
      )
  }
  
  # =======================================================
  # RUNNING MEAN SD RIBBON
  # =======================================================
  
  if (!is.null(running_var) &&
      !is.null(running_sd_var) &&
      running_var %in% names(df) &&
      running_sd_var %in% names(df)) {
    
    p <- p +
      geom_ribbon(
        data = df_hist,
        aes(
          ymin = .data[[running_var]] - .data[[running_sd_var]],
          ymax = .data[[running_var]] + .data[[running_sd_var]]
        ),
        fill = "grey",
        color = NA,
        alpha = 0.5
      ) +
      
      geom_ribbon(
        data = df_future,
        aes(
          ymin = .data[[running_var]] - .data[[running_sd_var]],
          ymax = .data[[running_var]] + .data[[running_sd_var]],
          fill = SSP
        ),
        color = NA,
        alpha = 0.2
      )
  }
  
  # =======================================================
  # 2015 boundary
  # =======================================================
  
  p <- p +
    geom_vline(
      xintercept = 2015,
      linetype = "dashed",
      color = "black"
    )
  
  # =======================================================
  # Scales
  # =======================================================
  
  p <- p +
    scale_color_manual(
      values = ssp_colors,
      drop = FALSE
    ) +
    
    scale_fill_manual(
      values = ssp_colors,
      drop = FALSE,
      guide = "none"
    ) +
    
    # =====================================================
  # Theme
  # =====================================================
  
  theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 8),
      legend.text = element_text(size = 5),
      legend.title = element_text(size = 4)
    ) +
    
    labs(
      x = "Year",
      y = ylab,
      color = "",
      title = title
    )
    #ylim(-0.25, 1.5)

  
  return(p)
}
# =========================================================
# Make plots (with SD ribbons where available)
# =========================================================
p_mineralised <- plot_diag(
  diag_annual,
  
  # Data column
  
  yvar = "annual_mineralised_anom_pg",
  
  # Axis label
  
  ylab = expression("Mineralised N [Pg N yr"^-1*"]"),
  title = "",
  filename = "annual_mineralised_anomaly.png",
  
  # Raw SD
  sd_var = "annual_mineralised_anom_pg_sd",
  
  # 20-year running mean
  running_var = "mineralised_pg_yr_anom_20yr",
  
  # Running-mean SD
  running_sd_var = "mineralised_pg_yr_anom_sd_20yr"
)

p_mineralised <- p_mineralised +
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
  "p_mineralised.pdf",
  p_mineralised,
  width = 8.5,
  height = 6,
  units = "cm", 
  dpi = 500
)

p_cum_mineralised <- plot_diag(
  diag_annual,
  
  # Data column
  
  yvar = "cumulative_mineralised_anom_pg",
  
  # Axis label
  
  ylab = expression("Cumulative Mineralised N [Pg N yr"^-1*"]"),
  title = "",
  filename = "cumulative_mineralised_pg.png",
  
  # Raw SD
  sd_var = "cumulative_mineralised_anom_pg_sd",
  
  # 20-year running mean
  running_var = "cumulative_mineralised_anom_20yr",
  
  # Running-mean SD
  running_sd_var = "cumulative_mineralised_anom_20yr_sd"
)

cum_mineral_N_release <- p_cum_mineralised +
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
  "cum_mineral_N_release.pdf",
  cum_mineral_N_release,
  width = 8.5,
  height = 6,
  units = "cm", 
  dpi = 500
)



p3 <- plot_diag(
  diag_annual,
  "annual_inorg_anom_pg",
  "Annual total inorganic N [Pg N yr\u207b\u00b9]",
  "Total inorganic N, annual",
  "annual_total_inorg_N.png",
  sd_var = "annual_inorg_anom_pg_sd"
)

p_cum_inorg <- plot_diag(
  diag_annual,
  
  # Data column
  
  yvar = "cumulative_inorg_anom_pg",
  
  # Axis label
  
  ylab = expression("Total inorganic N [Pg N]"),
  title = "",
  filename = "cumulative_inorg_N_anomaly.png",
  
  # Raw SD
  sd_var = "cumulative_inorg_anom_pg_sd",
  
  # 20-year running mean
  running_var = "cumulative_inorg_anom_20yr",
  
  # Running-mean SD
  running_sd_var = "cumulative_inorg_anom_20yr_sd"
)

p_cum_inorg <- p_cum_inorg +
  ylim(-0.15, 1.5)+
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
  "cum_inorg_N_release.pdf",
  p_cum_inorg,
  width = 8.5,
  height = 6,
  units = "cm", 
  dpi = 500
)




p_organic <- plot_diag(
  diag_annual,
  
  # Raw annual anomaly
  yvar = "organic_pool_anom_pg",
  
  ylab = "Organic N [Pg N]",
  title = "",
  filename = "organic_N_pool_anomaly.png",
  
  # Raw SD
  sd_var = "organic_pool_anom_pg_sd",
  
  # 20-year running mean
  running_var = "organic_pool_anom_20yr",
  
  # Running-mean SD
  running_sd_var = "organic_pool_anom_20yr_sd"
)

p_organic <- p_organic+
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
  "cum_organic_N.pdf",
  p_organic,
  width = 8.5,
  height = 6,
  units = "cm", 
  dpi = 500
)


p_rapid_inorg <- plot_diag(
  diag_annual,
  yvar = "annual_rapid_inorg_anom_pg",
  ylab = expression("Pre-thaw inorganic N [Pg N yr"^-1*"]"),
  title = "",
  filename = "annual_rapid_inorg_anomaly.png",
  sd_var = "annual_rapid_inorg_anom_pg_sd",
  running_var = "rapid_inorg_pg_yr_anom_20yr",
  running_sd_var = ""
)

p_rapid_abs <- plot_diag(
  diag_annual,
  yvar           = "inorg_rapid_pg_yr",
  ylab           = expression("Pre-thaw inorganic N [Pg N yr"^-1*"]"),
  title          = "",
  filename       = "annual_rapid_inorg_absolute.png",
  sd_var         = "inorg_rapid_pg_yr_sd",
  running_var    = "inorg_rapid_20yr",
  running_sd_var = "inorg_rapid_20yr_sd"
)
print(p_rapid_abs)

p_rapid_inorg <- p_rapid_inorg +
  theme(
    legend.margin = margin(t = 0),
    legend.box.margin = margin(1, 1, 1, 1),
    plot.margin = margin(3, 3, 2, 3),
    legend.position = "bottom",
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 5),
    axis.title.x = element_text(margin = margin(t = 2)),
    axis.title.y = element_text(margin = margin(r = 2)),
    axis.text.x = element_text(margin = margin(t = 1)),
    axis.text.y = element_text(margin = margin(r = 1)),
    axis.title = element_text(size = 6),
    axis.text  = element_text(size = 6)
  )

ggsave(
  "p_rapid_inorg.pdf",
  p_rapid_inorg,
  width = 8.5,
  height = 6,
  units = "cm",
  dpi = 500
)

p_cum_rapid_inorg <- plot_diag(
  diag_annual,
  yvar = "cumulative_rapid_inorg_anom_pg",
  ylab = expression("Cumulative rapid inorganic N [Pg N]"),
  title = "",
  filename = "cumulative_rapid_inorg_pg.png",
  sd_var = "cumulative_rapid_inorg_anom_pg_sd",
  running_var = "cumulative_rapid_inorg_anom_20yr",
  running_sd_var = "cumulative_rapid_inorg_anom_20yr_sd"
)

p_cum_rapid_inorg <- p_cum_rapid_inorg +
  theme(
    legend.margin = margin(t = 0),
    legend.box.margin = margin(1, 1, 1, 1),
    plot.margin = margin(3, 3, 2, 3),
    legend.position = "bottom",
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 5),
    legend.title = element_text(size = 5),
    axis.title.x = element_text(margin = margin(t = 2)),
    axis.title.y = element_text(margin = margin(r = 2)),
    axis.text.x = element_text(margin = margin(t = 1)),
    axis.text.y = element_text(margin = margin(r = 1)),
    axis.title = element_text(size = 6),
    axis.text  = element_text(size = 6)
  )

ggsave(
  "p_cum_rapid_inorg.pdf",
  p_cum_rapid_inorg,
  width = 8.5,
  height = 6,
  units = "cm",
  dpi = 500
)







## spatial plotting

# ============================================================
# Spatial plot of mineralised N from monthly dynamic-depth files
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
      "arctic_total_inorg_N_monthly_w_temp_sm_",
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
  geom_sf(data = arctic_sf_laea, aes(color = mineralised_N), size = 0.005) +
  scale_color_gradientn(
    colours = c("dodgerblue3", "#abd9e9", "#ffffbf", "orange", "darkred"),
    trans = "log1p",
    limits = c(0, 22),
    oob = scales::squish,
    name = expression(N~(g~m^{-2}~yr^{-1})),
    guide = guide_colorbar(barheight = unit(2.5, "cm"), barwidth = unit(0.3, "cm"), ticks = TRUE)
  ) +
  geom_sf(data = coastlines_laea, color = "black", linewidth = 0.09) +
  coord_sf(
    crs = laea_crs,
    xlim = xlim, ylim = ylim,
    datum = sf::st_crs(4326)        # <-- graticule drawn in lat/long, not projected meters
    #label_graticule = "SW"           # optional: adds lat/long tick labels on south & west edges
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    axis.text = element_text(size = 6),     # needed if you want the lat/long labels visible
    axis.title = element_blank(),
    panel.grid = element_line(color = "grey60", linewidth = 0.15, linetype = "dashed")
  )
spatial_plot


#-----------#
# Figure 2: spatial plot of total mineral N release (right column) + cumulative total mineral N release (left column)

spatial_plot <- spatial_plot +
  theme(legend.position = "right")


spatial_plot <- spatial_plot +
  theme(
    legend.position = "right",
    legend.margin = margin(t = 0),
    legend.box.margin = margin(1, 1, 1, 1),
    plot.margin = margin(3, 3, 2, 3),
    
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
  "spatial_mineralised_N.pdf",
  spatial_plot,
  width = 8.5,
  height = 6,
  units = "cm", 
  dpi = 500
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




library(terra)
library(dplyr)
library(ggplot2)

out_dir <- "monthly_mineralised/mean_2perc_baseline_new"
ssp_list <- c("126", "245", "370", "585")   # adjust to match your actual SSP codes
years_to_plot <- 1850:2099

# -----------------------------
# 1. Read k_env yearly files and compute area-weighted spatial mean
#    both monthly and annually, for each SSP
# -----------------------------

read_k_env_ssp <- function(ssp, years) {
  
  out_dir_yearly <- file.path(out_dir, paste0("yearly_nc_", ssp))
  
  monthly_rows <- vector("list", length(years))
  
  for (i in seq_along(years)) {
    yr <- years[i]
    
    f <- file.path(out_dir_yearly, paste0("arctic_k_env_monthly_w_temp_sm_", ssp, "_", yr, ".nc"))
    
    if (!file.exists(f)) {
      warning("Missing file: ", f)
      next
    }
    
    r <- rast(f)   # 12 layers, Jan-Dec
    area_rast <- cellSize(r[[1]], unit = "m")
    
    monthly_means <- terra::global(r, "mean", weights = area_rast, na.rm = TRUE)[, 1]
    
    monthly_rows[[i]] <- data.frame(
      SSP = ssp,
      Year = yr,
      Month = 1:12,
      Date = seq(as.Date(paste0(yr, "-01-16")), by = "month", length.out = 12),
      k_env_monthly = monthly_means
    )
  }
  
  bind_rows(monthly_rows)
}

k_env_monthly_all <- bind_rows(lapply(ssp_list, read_k_env_ssp, years = years_to_plot))

# -----------------------------
# 2. Annual mean (average across the 12 months within each year)
# -----------------------------

k_env_annual_all <- k_env_monthly_all %>%
  group_by(SSP, Year) %>%
  summarise(
    k_env_annual = mean(k_env_monthly, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    SSP = recode(
      SSP,
      "126" = "SSP1-2.6",
      "245" = "SSP2-4.5",
      "370" = "SSP3-7.0",
      "585" = "SSP5-8.5"
    )
  )
# -----------------------------
# 3. Define the three time windows and label the monthly data
# -----------------------------

time_windows <- list(
  "1880-1900" = 1880:1900,
  "2000-2020" = 2000:2020,
  "2080-2099" = 2080:2099
)

assign_window <- function(yr) {
  for (w in names(time_windows)) {
    if (yr %in% time_windows[[w]]) return(w)
  }
  NA_character_
}

k_env_monthly_windows <- k_env_monthly_all %>%
  mutate(Window = sapply(Year, assign_window)) %>%
  filter(!is.na(Window)) %>%
  mutate(Window = factor(Window, levels = names(time_windows)))

# -----------------------------
# 3b. Average across years WITHIN each window, per month, per SSP
#     -> one 12-month seasonal cycle per SSP per window
# -----------------------------

k_env_seasonal <- k_env_monthly_windows %>%
  group_by(SSP, Window, Month) %>%
  summarise(k_env_mean = mean(k_env_monthly, na.rm = TRUE), .groups = "drop")
# -----------------------------
# 4. Plot styling
# -----------------------------

ssp_labels <- c(
  "126" = "SSP1-2.6",
  "245" = "SSP2-4.5",
  "370" = "SSP3-7.0",
  "585" = "SSP5-8.5"
)

window_colors <- c(
  "1880-1900" = "steelblue",
  "2000-2020" = "darkorange",
  "2080-2099" = "firebrick"
)


ssp_colors <- c(
  "SSP1-2.6" = "blue",
  "SSP2-4.5" = "orange",
  "SSP3-7.0" = "#D73027",
  "SSP5-8.5" = "#7B3294"
)



# -----------------------------
# 5. Plot: one line per time window, faceted by SSP
# -----------------------------

p_seasonal <- ggplot(k_env_seasonal, aes(x = Month, y = k_env_mean, colour = Window)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1) +
  facet_wrap(~ SSP, nrow = 2, labeller = as_labeller(ssp_labels)) +
  scale_colour_manual(values = window_colors, name = "Time period") +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        strip.text = element_text(size = 5),
        title = element_text(size = 5),
        plot.subtitle = element_text(size = 5)) +
  labs(
    title = "",
    subtitle = "",
    x = "Month",
    y = expression(k[env]~"[-]")
  )

print(p_seasonal)

p_seasonal <- p_seasonal +
  theme(
    legend.margin = margin(t = 0),
    legend.box.margin = margin(1, 1, 1, 1),
    plot.margin = margin(2, 2, 1, 2),
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
  "p_seasonal_kenv.pdf",
  p_seasonal,
  width = 8.5,
  height = 6,
  units = "cm", 
  dpi = 500
)
# -----------------------------
# 6. Plot: annual mean time series, full 1850-2099, coloured by SSP
# -----------------------------

p_annual <- ggplot(k_env_annual_all, aes(x = Year, y = k_env_annual, colour = SSP)) +
  geom_line(linewidth = 0.5) +
  scale_colour_manual(values = ssp_colors, labels = ssp_labels, name = "SSP scenario") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        strip.text = element_text(size = 5),
        title = element_text(size = 5),
        plot.subtitle = element_text(size = 5)) +
  labs(
    title = "",
    x = "Year",
    y = expression(k[env]~"[-]")
  )

print(p_annual)


p_annual <- p_annual +
  theme(
    legend.margin = margin(t = 0),
    legend.box.margin = margin(1, 1, 1, 1),
    plot.margin = margin(2, 2, 1, 2),
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
  "p_annual.pdf",
  p_annual,
  width = 8.5,
  height = 6,
  units = "cm", 
  dpi = 500
)







library(dplyr)
library(readr)
library(zoo)
library(ggplot2)

mean_dir <- "monthly_mineralised/mean_2perc_baseline_new"
std_dir  <- "monthly_mineralised/std_2perc_baseline_new"

ssps <- c("126", "245", "370", "585")
ssp_labels <- c("126" = "SSP1-2.6", "245" = "SSP2-4.5", "370" = "SSP3-7.0", "585" = "SSP5-8.5")
ssp_colors <- c("SSP1-2.6" = "blue", "SSP2-4.5" = "orange", "SSP3-7.0" = "#D73027", "SSP5-8.5" = "#7B3294")

read_diag_csv <- function(dir, ssp) {
  f <- file.path(dir, paste0("arctic_monthly_diagnostics_w_temp", ssp, "_1850_2100.csv"))
  if (!file.exists(f)) { warning("Missing file: ", f); return(NULL) }
  read_csv(f, show_col_types = FALSE) %>% mutate(SSP_raw = ssp, SSP = ssp_labels[ssp])
}

diag_monthly_mean <- bind_rows(lapply(ssps, function(s) read_diag_csv(mean_dir, s))) %>%
  select(SSP_raw, SSP, Year, Month, inorg_rapid_available)

diag_monthly_sd <- bind_rows(lapply(ssps, function(s) read_diag_csv(std_dir, s))) %>%
  select(SSP_raw, SSP, Year, Month, inorg_rapid_available) %>%
  rename(inorg_rapid_available_sd = inorg_rapid_available)

diag_monthly <- diag_monthly_mean %>%
  left_join(diag_monthly_sd, by = c("SSP_raw", "SSP", "Year", "Month"))

# ---- Annual totals + cumulative sum ----
diag_annual <- diag_monthly %>%
  group_by(SSP, Year) %>%
  summarise(
    inorg_rapid_pg_yr    = sum(inorg_rapid_available, na.rm = TRUE),
    inorg_rapid_pg_yr_sd = sum(inorg_rapid_available_sd, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(SSP = factor(SSP, levels = c("SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5"))) %>%
  group_by(SSP) %>%
  arrange(Year, .by_group = TRUE) %>%
  mutate(
    cumulative_inorg_rapid_pg    = cumsum(inorg_rapid_pg_yr),
    cumulative_inorg_rapid_pg_sd = cumsum(inorg_rapid_pg_yr_sd)
  ) %>%
  ungroup()

# ---- Anomaly relative to 2000-2020, for both annual and cumulative ----
compute_anomaly_pair <- function(df, value_col, sd_col, ref_range = c(2000, 2020)) {
  ref <- df %>% filter(Year >= ref_range[1], Year <= ref_range[2])
  ref_mean <- mean(ref[[value_col]], na.rm = TRUE)
  ref_sd   <- mean(ref[[sd_col]], na.rm = TRUE)
  df %>%
    mutate(
      "{value_col}_anom" := .data[[value_col]] - ref_mean,
      "{sd_col}_anom"     := .data[[sd_col]] - ref_sd
    )
}

diag_annual <- diag_annual %>%
  group_by(SSP) %>%
  group_modify(~ compute_anomaly_pair(.x, "inorg_rapid_pg_yr", "inorg_rapid_pg_yr_sd")) %>%
  group_modify(~ compute_anomaly_pair(.x, "cumulative_inorg_rapid_pg", "cumulative_inorg_rapid_pg_sd")) %>%
  ungroup() %>%
  rename(
    annual_rapid_inorg_anom_pg        = inorg_rapid_pg_yr_anom,
    annual_rapid_inorg_anom_pg_sd     = inorg_rapid_pg_yr_sd_anom,
    cumulative_rapid_inorg_anom_pg    = cumulative_inorg_rapid_pg_anom,
    cumulative_rapid_inorg_anom_pg_sd = cumulative_inorg_rapid_pg_sd_anom
  )

# ---- 20-year rolling means ----
diag_annual <- diag_annual %>%
  group_by(SSP) %>%
  arrange(Year, .by_group = TRUE) %>%
  mutate(
    rapid_inorg_pg_yr_anom_20yr    = zoo::rollapply(annual_rapid_inorg_anom_pg, 20, mean, align = "center", fill = NA),
    rapid_inorg_pg_yr_anom_sd_20yr = zoo::rollapply(annual_rapid_inorg_anom_pg_sd, 20, mean, align = "center", fill = NA),
    cumulative_rapid_inorg_anom_20yr    = zoo::rollapply(cumulative_rapid_inorg_anom_pg, 20, mean, align = "center", fill = NA),
    cumulative_rapid_inorg_anom_20yr_sd = zoo::rollapply(cumulative_rapid_inorg_anom_pg_sd, 20, mean, align = "center", fill = NA)
  ) %>%
  ungroup()

# ---- Generic plotting function ----
plot_anom <- function(df, raw_var, roll_var, roll_sd_var, ylab, filename) {
  
  df_hist   <- df %>% filter(Year <= 2014, SSP == levels(SSP)[1])
  df_future <- df %>% filter(Year >= 2015)
  
  p <- ggplot(df, aes(x = Year, group = SSP)) +
    
    # Raw values - historical
    geom_line(
      data = df_hist,
      aes(y = .data[[raw_var]]),
      color = "black", linewidth = 0.2, alpha = 0.45, na.rm = TRUE
    ) +
    
    # Raw values - future
    geom_line(
      data = df_future,
      aes(y = .data[[raw_var]], color = SSP),
      linewidth = 0.2, alpha = 0.45, na.rm = TRUE
    ) +
    
    # Rolling mean - historical
    geom_line(
      data = df_hist,
      aes(y = .data[[roll_var]]),
      color = "black", linewidth = 0.2, na.rm = TRUE
    ) +
    
    # Rolling mean - future
    geom_line(
      data = df_future,
      aes(y = .data[[roll_var]], color = SSP),
      linewidth = 0.2, na.rm = TRUE
    ) +
    
    # Rolling SD ribbon - historical
    geom_ribbon(
      data = df_hist,
      aes(
        ymin = .data[[roll_var]] - .data[[roll_sd_var]],
        ymax = .data[[roll_var]] + .data[[roll_sd_var]]
      ),
      fill = "grey", color = NA, alpha = 0.2
    ) +
    
    # Rolling SD ribbon - future
    geom_ribbon(
      data = df_future,
      aes(
        ymin = .data[[roll_var]] - .data[[roll_sd_var]],
        ymax = .data[[roll_var]] + .data[[roll_sd_var]],
        fill = SSP
      ),
      color = NA, alpha = 0.15
    ) +
    
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.2) +
    geom_vline(xintercept = 2015, linetype = "dashed", color = "black", linewidth = 0.2) +
    
    scale_color_manual(values = ssp_colors, drop = FALSE) +
    scale_fill_manual(values = ssp_colors, drop = FALSE, guide = "none") +
    
    theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 8),
      legend.text = element_text(size = 5),
      legend.title = element_text(size = 4)
    ) +
    ylim(-0.075, 0.3)+
    
    labs(x = "Year", y = ylab, color = "", title = "")
  
  ggsave(filename, p, width = 8.5, height = 6, units = "cm", dpi = 500)
  return(p)
}

p_rapid_inorg <- plot_anom(
  diag_annual,
  raw_var = "annual_rapid_inorg_anom_pg",
  roll_var = "rapid_inorg_pg_yr_anom_20yr",
  roll_sd_var = "rapid_inorg_pg_yr_anom_sd_20yr",
  ylab = expression("Rapid inorganic N anomaly [Pg N yr"^-1*"]"),
  filename = "p_rapid_inorg.pdf"
)

p_cum_rapid_inorg <- plot_anom(
  diag_annual,
  raw_var = "cumulative_rapid_inorg_anom_pg",
  roll_var = "cumulative_rapid_inorg_anom_20yr",
  roll_sd_var = "cumulative_rapid_inorg_anom_20yr_sd",
  ylab = expression("Pre-thaw inorganic N [Pg N]"),
  filename = "p_cum_rapid_inorg.pdf"
)



p_cum_rapid_inorg <- p_cum_rapid_inorg +
  theme(
    legend.margin = margin(t = 0),
    legend.box.margin = margin(1, 1, 1, 1),
    plot.margin = margin(2, 2, 1, 2),
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
  "p_cum_pre-thaw_inorg.pdf",
  p_cum_rapid_inorg,
  width = 8.5,
  height = 6,
  units = "cm", 
  dpi = 500
)



