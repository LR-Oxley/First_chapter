library(terra)

# -----------------------------
# SETTINGS
# -----------------------------
terraOptions(progress = 1, memfrac = 0.8)
sc <- "585"
message("Running SSP", sc)
temp_offsets <- c(2,2.5)
years <- 1850:2099

ref_start <- 2000
ref_end   <- 2020

k_base_mineralisation <- 0.019
Ea <- 80000
Ed <- 200000
t_opt_C <- 28
f_inorg_rapid <- 0.1235

common_extent <- ext(-179.95, 179.95, 30, 90)

# -----------------------------
# FUNCTION
# -----------------------------
peaked_arrhenius <- function(temp_C, Ea, Ed, t_opt_C = 28) {
  R <- 8.314462618
  
  temp_K <- temp_C + 273.15
  t_opt_K <- t_opt_C + 273.15
  
  hlp1 <- temp_K - t_opt_K
  hlp2 <- temp_K * t_opt_K * R
  
  num <- Ed * exp(Ea * hlp1 / hlp2)
  den <- Ed - Ea * (1 - exp(Ed * hlp1 / hlp2))
  
  num / den
}

# -----------------------------
# LOAD DATA (SSP585)
# -----------------------------
temp_data <- rast("mean_ST_585_clean_final.nc")
thawed_total <- rast("1.9%_min/arctic_total_thawed_585_no_lim.tif")

ext(temp_data) <- common_extent
temp_data <- temp_data[[1:250]] - 273.15


# -----------------------------
# BASELINE MINERALISATION
# -----------------------------
message("Computing baseline mineralisation...")

k_T <- app(temp_data,
           fun = function(x) peaked_arrhenius(x, Ea, Ed, t_opt_C))

ref_idx <- which(years >= ref_start & years <= ref_end)
k_ref <- mean(k_T[[ref_idx]], na.rm = TRUE)

k_factor <- k_T / k_ref
k_t <- k_base_mineralisation * k_factor

f_org <- 1 - f_inorg_rapid
org_pool <- thawed_total * f_org

mineralised_baseline <- org_pool * k_t

# -----------------------------
# OUTPUT DIR
# -----------------------------
out_dir <- "ssp585_temp_sensitivity_mineralisation"
dir.create(out_dir, showWarnings = FALSE)

writeRaster(mineralised_baseline,
            file.path(out_dir, "mineralised_baseline.tif"),
            overwrite = TRUE)

# -----------------------------
# TEMPERATURE LOOP
# -----------------------------
for (offset in temp_offsets) {
  
  message("Offset: ", offset, " °C")
  
  temp_mod <- temp_data + offset
  
  k_T_mod <- app(temp_mod,
                 fun = function(x) peaked_arrhenius(x, Ea, Ed, t_opt_C))
  
  k_factor_mod <- k_T_mod / k_ref
  k_t_mod <- k_base_mineralisation * k_factor_mod
  
  mineralised_mod <- org_pool * k_t_mod
  
  sensitivity <- mineralised_mod / mineralised_baseline
  
  # -----------------------------
  # SAVE OUTPUTS
  # -----------------------------
  writeRaster(
    mineralised_mod,
    file.path(out_dir, paste0("mineralised_offset_", offset, "C.tif")),
    overwrite = TRUE
  )
  
  writeRaster(
    sensitivity,
    file.path(out_dir, paste0("sensitivity_offset_", offset, "C.tif")),
    overwrite = TRUE
  )
  
  # Late century (2080–2100)
  future_idx <- which(years >= 2080 & years <= 2100)
  late_mean <- mean(sensitivity[[future_idx]], na.rm = TRUE)
  
  writeRaster(
    late_mean,
    file.path(out_dir, paste0("late_century_sensitivity_", offset, "C.tif")),
    overwrite = TRUE
  )
  
  rm(temp_mod, k_T_mod, k_factor_mod, k_t_mod,
     mineralised_mod, sensitivity, late_mean)
  gc()
}

message("Done.")


# check outputs

library(terra)
library(dplyr)
library(ggplot2)

# -----------------------------
# SETTINGS
# -----------------------------
years <- 1850:2099
future_years<- years
#future_years <- 2080:2099

offsets <- c(-2.5,-2,-1.5,-1,-0.5,0.5,1,1.5,2,2.5)

data_dir <- "ssp585_temp_sensitivity_mineralisation"

# -----------------------------
# HELPER FUNCTIONS
# -----------------------------
weighted_mean_timeseries <- function(r, weights) {
  sapply(1:nlyr(r), function(i) {
    num <- global(r[[i]] * weights, "sum", na.rm = TRUE)[[1]]
    den <- global(weights, "sum", na.rm = TRUE)[[1]]
    num / den
  })
}

# -----------------------------
# LOAD BASELINE
# -----------------------------
baseline <- rast(file.path(data_dir, "mineralised_baseline.tif"))
plot(baseline[[249]])
# area weights (cell area in m²)
weights <- cellSize(baseline[[1]], unit = "m")
compareGeom(baseline, weights, stopOnError = FALSE)
nlyr(baseline)
length(years)
# baseline time series
baseline_ts <- global(baseline * weights, "sum", na.rm = TRUE)[,1] 
baseline_ts<- baseline_ts / 1e12

# -----------------------------
# LOOP OVER OFFSETS
# -----------------------------
results_list <- list()

# -----------------------------
# PRECOMPUTE
# -----------------------------
future_idx <- which(years %in% future_years)
den <- global(weights, "sum", na.rm = TRUE)[[1]]

# ---- baseline late century ----
baseline_late <- mean(baseline[[future_idx]], na.rm = TRUE)

baseline_late_mean <- global(baseline_late * weights, "sum", na.rm = TRUE)[[1]] / den

# -----------------------------
# LOOP
# -----------------------------
results_list <- list()

for (offset in offsets) {
  
  cat("Processing offset:", offset, "°C\n")
  
  # --- load raster ---
  mineral_file <- file.path(data_dir, paste0("mineralised_offset_", offset, "C.tif"))
  mineralised <- rast(mineral_file)
  
  # -----------------------------
  # LATE CENTURY MEAN
  # -----------------------------
  miner_late <- mean(mineralised[[future_idx]], na.rm = TRUE)
  
  late_miner <- global(miner_late * weights, "sum", na.rm = TRUE)[[1]] / den
  
  # -----------------------------
  # ✅ CORRECT SENSITIVITY
  # -----------------------------
  relative_change <- late_miner / baseline_late_mean
  pct_change <- (relative_change - 1) * 100
  
  # optional: Q10 only for positive ΔT
  Q10 <- if (offset > 0) relative_change^(10 / offset) else NA
  
  results_list[[as.character(offset)]] <- data.frame(
    offset = offset,
    relative_change = relative_change,
    percent_change = pct_change,
    Q10 = Q10,
    late_mineralisation = late_miner
  )
  
  rm(mineralised, miner_late)
  gc()
}

results <- bind_rows(results_list)
print(results)

# -----------------------------
# PLOT 1: Sensitivity curve
# -----------------------------
ggplot(results, aes(x = offset, y = percent_change)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Temperature change (°C)",
    y = "Change in mineralisation (%)",
    title = "Temperature sensitivity of mineralisation (SSP585)"
  ) +
  theme_minimal()

# -----------------------------
# PLOT 2: Q10 behaviour
# -----------------------------
ggplot(results, aes(x = offset, y = Q10)) +
  geom_line() +
  geom_point() +
  labs(
    x = "Temperature change (°C)",
    y = "Effective Q10",
    title = "Emergent Q10 from peaked Arrhenius"
  ) +
  theme_minimal()

# -----------------------------
# OPTIONAL: MAP VISUALISATION
# -----------------------------
# Example: plot +2°C late-century sensitivity

example_map <- rast(file.path(data_dir, "late_century_sensitivity_2C.tif"))

plot(example_map, main = "+2°C sensitivity (2080–2100)")

# Convert to % change for interpretation
example_pct <- (example_map - 1) * 100
plot(example_pct, main = "+2°C % change in mineralisation")

############################################################################

future_idx <- which(years %in% future_years)

# mean late-century raster
baseline_late <- mean(baseline[[future_idx]], na.rm = TRUE)

# total kg N
baseline_total_kg <- global(baseline_late * weights, "sum", na.rm = TRUE)[[1]]

# convert to Pg N
baseline_total_pg <- baseline_total_kg / 1e12

global(baseline_late * weights, "sum", na.rm = TRUE)[[1]] / 1e12
results_pg$total_pg_N
results_list <- list()

for (offset in offsets) {
  
  cat("Processing offset:", offset, "°C\n")
  
  mineral_file <- file.path(data_dir, paste0("mineralised_offset_", offset, "C.tif"))
  mineralised <- rast(mineral_file)
  
  # --- late century mean ---
  miner_late <- mean(mineralised[[future_idx]], na.rm = TRUE)
  
  # --- total kg N ---
  total_kg <- global(miner_late * weights, "sum", na.rm = TRUE)[[1]]
  
  # --- convert to Pg N ---
  total_pg <- total_kg / 1e12
  
  # --- absolute change ---
  delta_pg <- total_pg - baseline_total_pg
  
  results_list[[as.character(offset)]] <- data.frame(
    offset = offset,
    total_pg_N = total_pg,
    delta_pg_N = delta_pg
  )
  
  rm(mineralised, miner_late)
  gc()
}

results_pg <- do.call(rbind, results_list)
print(results_pg)





### with ALD temp

mineralised_585<-rast("1.9%_min/arctic_bioavailable_N_pool_585_no_lim_mean.tif")[[1:250]]
mineralised_370<-rast("1.9%_min/arctic_bioavailable_N_pool_370_no_lim_mean.tif")[[1:250]]
mineralised_245<-rast("1.9%_min/arctic_bioavailable_N_pool_245_no_lim_mean.tif")[[1:250]]
mineralised_126<-rast("1.9%_min/arctic_bioavailable_N_pool_126_no_lim_mean.tif")[[1:250]]
mineralised_585_std<-rast("1.9%_min/arctic_bioavailable_N_pool_585_no_lim_std.tif")[[1:250]]
mineralised_370_std<-rast("1.9%_min/arctic_bioavailable_N_pool_370_no_lim_std.tif")[[1:250]]
mineralised_245_std<-rast("1.9%_min/arctic_bioavailable_N_pool_245_no_lim_std.tif")[[1:250]]
mineralised_126_std<-rast("1.9%_min/arctic_bioavailable_N_pool_126_no_lim_std.tif")[[1:250]]

mineralised_585_ALD<-rast("1.9%_min/ALD_temp/arctic_bioavailable_N_pool_585_no_lim_mean.tif")[[1:250]]
mineralised_370_ALD<-rast("1.9%_min/ALD_temp/arctic_bioavailable_N_pool_370_no_lim_mean.tif")[[1:250]]
mineralised_245_ALD<-rast("1.9%_min/ALD_temp/arctic_bioavailable_N_pool_245_no_lim_mean.tif")[[1:250]]
mineralised_126_ALD<-rast("1.9%_min/ALD_temp/arctic_bioavailable_N_pool_126_no_lim_mean.tif")[[1:250]]
mineralised_585_ALD_std<-rast("1.9%_min/ALD_temp/arctic_bioavailable_N_pool_585_no_lim_std.tif")[[1:250]]
mineralised_370_ALD_std<-rast("1.9%_min/ALD_temp/arctic_bioavailable_N_pool_370_no_lim_std.tif")[[1:250]]
mineralised_245_ALD_std<-rast("1.9%_min/ALD_temp/arctic_bioavailable_N_pool_245_no_lim_std.tif")[[1:250]]
mineralised_126_ALD_std<-rast("1.9%_min/ALD_temp/arctic_bioavailable_N_pool_126_no_lim_std.tif")[[1:250]]


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


# Define years (adjust if needed)
years <- 1850:(1850 + nlyr(mineralised_585) - 1)

# Calculate cell area
cell_area <- cellSize(mineralised_585, unit = "m")



# Extract time series for each scenario
df_585 <- extract_ts(mineralised_585, mineralised_585_std, "SSP585", years, cell_area)
df_370 <- extract_ts(mineralised_370, mineralised_370_std, "SSP370", years, cell_area)
df_245 <- extract_ts(mineralised_245, mineralised_245_std, "SSP245", years, cell_area)
df_126 <- extract_ts(mineralised_126, mineralised_126_std, "SSP126", years, cell_area)

df_585_ALD <- extract_ts(mineralised_585_ALD, mineralised_585_ALD_std, "SSP585", years, cell_area)
df_370_ALD <- extract_ts(mineralised_370_ALD, mineralised_370_ALD_std, "SSP370", years, cell_area)
df_245_ALD <- extract_ts(mineralised_245_ALD, mineralised_245_ALD_std, "SSP245", years, cell_area)
df_126_ALD <- extract_ts(mineralised_126_ALD, mineralised_126_ALD_std, "SSP126", years, cell_area)


df_all <- bind_rows(df_585, df_370, df_245, df_126)
df_all_ALD <- bind_rows(df_585_ALD , df_370_ALD , df_245_ALD , df_126_ALD )


### compare mineralisation flux of surface and ald temperature

library(dplyr)

df_compare <- df_all_ALD %>%
  rename(Mean_ALD_temp = Mean) %>%
  left_join(
    df_all %>% rename(Mean_surface_temp = Mean),
    by = c("Year", "SSP")
  )

df_compare <- df_compare %>%
  mutate(
    pct_diff = 100 * (Mean_ALD_temp- Mean_surface_temp) / Mean_surface_temp
  )

# how much smaller is the mineralisation rate under ALD temperatures?

df_summary <- df_compare %>%
  group_by(SSP) %>%
  summarise(
    pct_diff_baseline = mean(pct_diff[Year >= 2000 & Year <= 2020], na.rm = TRUE),
    pct_diff_future   = mean(pct_diff[Year >= 2080 & Year <= 2099], na.rm = TRUE),
    change_in_effect  = pct_diff_future - pct_diff_baseline,
    .groups = "drop"
  )



print(df_summary)

library(ggplot2)

ggplot(df_compare, aes(x = Year, y = pct_diff, color = SSP)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    y = "% difference",
    title = "Decrease in mineralised N due to ALD temperature"
  ) +
  theme_minimal()

df_compare <- df_compare %>%
  mutate(
    abs_diff = Mean_surface_temp - Mean_ALD_temp
  )



ssp_colors <- c(
  "SSP585" = "#7B3294",
  "SSP370" = "#D73027",
  "SSP245" = "orange",
  "SSP126" = "blue"
)

ggplot(df_compare, aes(x = Year, y = pct_diff, color = SSP)) +
  
  geom_line(linewidth = 1) +
  
  # smooth signal (optional but recommended)
  #geom_smooth(se = FALSE, linewidth = 1.2, linetype = "solid") +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  scale_color_manual(values = ssp_colors) +
  
  labs(
    y = "Change in mineralisation due to ALD (%)",
    x = "Year",
    title = "Effect of ALD temperature on mineralisation",
    subtitle = "Relative to surface temperature"
  ) +
  
  theme_minimal() +
  theme(
    legend.position = "right"
  )







library(dplyr)

df_all$type <- "Surface temperature"
df_all_ALD$type <- "ALD temperature"
df_plot <- bind_rows(df_all, df_all_ALD)
df_plot <- df_plot %>%
  mutate(
    SSP = recode(
      SSP,
      "SSP126" = "SSP1-2.6",
      "SSP245" = "SSP2-4.5",
      "SSP370" = "SSP3-7.0",
      "SSP585" = "SSP5-8.5"
    )
  )
library(ggplot2)

ssp_colors <- c(
  "SSP3-7.0" = "#D73027",
  "SSP5-8.5" = "#7B3294",
  "SSP2-4.5" = "orange",
  "SSP1-2.6" = "blue"
)
ggplot(df_plot, aes(x = Year, y = Mean, color = SSP, linetype = type)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = ssp_colors) +
  labs(
    y = "[Pg N]",
    title = "Cumulative bio-available N"
  ) +
  theme_minimal()


library(ggplot2)
library(dplyr)
library(zoo)

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

# ---- STANDARD ----
df_std <- df_all %>%
  mutate(type = "Surface layer temperature") %>%
  group_by(SSP) %>%
  mutate(
    baseline_mean = mean(Mean[Year >= 2000 & Year <= 2020], na.rm = TRUE),
    baseline_sd = mean(SD[Year >= 2000 & Year <= 2020], na.rm = TRUE),
    Mean_anom = Mean - baseline_mean,
    SD_anom = SD-baseline_sd,
    Rolling_Mean_anom = rollapply(Mean_anom, width = 20, fill = NA,FUN = mean, align = "center"),
    Rolling_SD_anom   = rollapply(SD_anom, width = 20, FUN = mean, fill = NA, align = "center"),
    Period = ifelse(Year < 2015, "Before 2015", "After 2015")
  ) %>%
  ungroup()

# ---- ALD ----
df_ald <- df_all_ALD %>%
  mutate(type = "Active layer temperature") %>%
  group_by(SSP) %>%
  mutate(
    baseline_mean = mean(Mean[Year >= 2000 & Year <= 2020], na.rm = TRUE),
    baseline_sd = mean(SD[Year >= 2000 & Year <= 2020], na.rm = TRUE),
    Mean_anom = Mean - baseline_mean,
    SD_anom = SD-baseline_sd,
    Rolling_Mean_anom = rollapply(Mean_anom, width = 20, fill = NA,FUN = mean, align = "center"),
    Rolling_SD_anom   = rollapply(SD_anom, width = 20, FUN = mean, fill = NA, align = "center"),
    Period = ifelse(Year < 2015, "Before 2015", "After 2015")
  ) %>%
  ungroup()

# ---- COMBINE ----
df_plot <- bind_rows(df_std, df_ald)



ssp_colors <- c(
  "SSP585" = "#7B3294",
  "SSP370" = "#D73027",
  "SSP245" = "orange",
  "SSP126" = "blue"
)

library(ggplot2)

ggplot(df_plot, aes(x = Year, y = Mean_anom, group = interaction(SSP, type))) +
  
  # ---- BEFORE 2015 (black, both types) ----
geom_line(
  data = filter(df_plot, Period == "Before 2015"),
  aes(linetype = type),
  color = "black",
  linewidth = 0.8
) +
  
  # ---- AFTER 2015 (colored by SSP, linetype = type) ----
geom_line(
  data = filter(df_plot, Period == "After 2015"),
  aes(color = SSP, linetype = type),
  linewidth = 0.8
) +
  
  # ---- Rolling mean BEFORE 2015 ----
geom_line(
  data = filter(df_plot, Year < 2015),
  aes(y = Rolling_Mean_anom, linetype = type),
  color = "black",
  linewidth = 1.2
) +
  
  # ---- Rolling mean AFTER 2015 ----
geom_line(
  data = filter(df_plot, Year >= 2015),
  aes(y = Rolling_Mean_anom, color = SSP, linetype = type),
  linewidth = 1.2
) +
  
  # ---- Rolling SD BEFORE 2015 ----
geom_ribbon(
  data = filter(df_plot, Year < 2015),
  aes(
    ymin = Rolling_Mean_anom - Rolling_SD_anom,
    ymax = Rolling_Mean_anom + Rolling_SD_anom
  ),
  fill = "grey",
  alpha = 0.15
) +
  
  # ---- Rolling SD AFTER 2015 ----
geom_ribbon(
  data = filter(df_plot, Year >= 2015),
  aes(
    ymin = Rolling_Mean_anom - Rolling_SD_anom,
    ymax = Rolling_Mean_anom + Rolling_SD_anom,
    fill = SSP
  ),
  alpha = 0.15
) +
  
  # ---- Reference lines ----
geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 2015, linetype = "dashed") +
  
  # ---- Scales ----
scale_color_manual(values = ssp_colors) +
  scale_fill_manual(values = ssp_colors) +
  
  # ---- Linetype legend (IMPORTANT) ----
scale_linetype_manual(values = c("Surface layer temperature" = "solid", "Active layer temperature" = "dashed")) +
  
  # ---- Labels ----
labs(
  x = "Year",
  y = "Cumulative mineralised N (Pg N)",
  title = "Cumulative mineralised N",
  color = "SSP",
  fill = "SSP",
  linetype = "Model"
) +
  
  # ---- Theme ----
theme_minimal() +
  theme(
    legend.position = "below",
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9)
  )



# ALD

ALD_plot<-ggplot(df_ald, aes(x = Year, y = Mean_anom, group = interaction(SSP, type))) +
  
  # ---- BEFORE 2015 (black, both types) ----
geom_line(
  data = filter(df_ald, Period == "Before 2015"),
  aes(linetype = type),
  color = "black",
  linewidth = 0.8
) +
  
  # ---- AFTER 2015 (colored by SSP, linetype = type) ----
geom_line(
  data = filter(df_ald, Period == "After 2015"),
  aes(color = SSP, linetype = type),
  linewidth = 0.8
) +
  
  # ---- Rolling mean BEFORE 2015 ----
geom_line(
  data = filter(df_ald, Year < 2015),
  aes(y = Rolling_Mean_anom, linetype = type),
  color = "black",
  linewidth = 1.2
) +
  
  # ---- Rolling mean AFTER 2015 ----
geom_line(
  data = filter(df_ald, Year >= 2015),
  aes(y = Rolling_Mean_anom, color = SSP, linetype = type),
  linewidth = 1.2
) +
  
  # ---- Rolling SD BEFORE 2015 ----
geom_ribbon(
  data = filter(df_ald, Year < 2015),
  aes(
    ymin = Rolling_Mean_anom - Rolling_SD_anom,
    ymax = Rolling_Mean_anom + Rolling_SD_anom
  ),
  fill = "grey",
  alpha = 0.15
) +
  
  # ---- Rolling SD AFTER 2015 ----
geom_ribbon(
  data = filter(df_ald, Year >= 2015),
  aes(
    ymin = Rolling_Mean_anom - Rolling_SD_anom,
    ymax = Rolling_Mean_anom + Rolling_SD_anom,
    fill = SSP
  ),
  alpha = 0.15
) +
  
  # ---- Reference lines ----
geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 2015, linetype = "dashed") +
  
  # ---- Scales ----
scale_color_manual(values = ssp_colors) +
  scale_fill_manual(values = ssp_colors) +
  

  # ---- Labels ----
labs(
  x = "Year",
  y = "Pg N",
  title = "Bio-available N, ALD temp",
  color = "SSP",
  fill = "SSP",
  linetype = "Model"
) +
  
  # ---- Theme ----
theme_minimal() +
  theme(
    legend.position = "below",
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9)
  )


# surface temp

plot_surface<-ggplot(df_std, aes(x = Year, y = Mean_anom, group = interaction(SSP, type))) +
  
  # ---- BEFORE 2015 (black, both types) ----
geom_line(
  data = filter(df_std, Period == "Before 2015"),
  aes(linetype = type),
  color = "black",
  linewidth = 0.8
) +
  
  # ---- AFTER 2015 (colored by SSP, linetype = type) ----
geom_line(
  data = filter(df_std, Period == "After 2015"),
  aes(color = SSP, linetype = type),
  linewidth = 0.8
) +
  
  # ---- Rolling mean BEFORE 2015 ----
geom_line(
  data = filter(df_std, Year < 2015),
  aes(y = Rolling_Mean_anom, linetype = type),
  color = "black",
  linewidth = 1.2
) +
  
  # ---- Rolling mean AFTER 2015 ----
geom_line(
  data = filter(df_std, Year >= 2015),
  aes(y = Rolling_Mean_anom, color = SSP, linetype = type),
  linewidth = 1.2
) +
  
  # ---- Rolling SD BEFORE 2015 ----
geom_ribbon(
  data = filter(df_std, Year < 2015),
  aes(
    ymin = Rolling_Mean_anom - Rolling_SD_anom,
    ymax = Rolling_Mean_anom + Rolling_SD_anom
  ),
  fill = "grey",
  alpha = 0.15
) +
  
  # ---- Rolling SD AFTER 2015 ----
geom_ribbon(
  data = filter(df_std, Year >= 2015),
  aes(
    ymin = Rolling_Mean_anom - Rolling_SD_anom,
    ymax = Rolling_Mean_anom + Rolling_SD_anom,
    fill = SSP
  ),
  alpha = 0.15
) +
  
  # ---- Reference lines ----
geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 2015, linetype = "dashed") +
  
  # ---- Scales ----
scale_color_manual(values = ssp_colors) +
  scale_fill_manual(values = ssp_colors) +
  
  
  # ---- Labels ----
labs(
  x = "Year",
  y = "Pg N",
  title = "Bio-available N, surface temp",
  color = "SSP",
  fill = "SSP",
  linetype = "Model"
) +
  
  # ---- Theme ----
theme_minimal() +
  theme(
    legend.position = "below",
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9)
  )


combined_plot
combined_plot<-plot_surface | ALD_plot+ plot_annotation(tag_levels = "a", tag_suffix = ") ")
ggsave("Fig2.png", combined_plot,
       width = 15,
       height = 15,
       units = "cm", 
       #scale = 1.4,
       dpi = 500,  #70
       #device = cairo_pdf
)


ylim_range <- range(-3,8)

# apply same y-limits
plot_surface <- plot_surface + coord_cartesian(ylim = ylim_range)
ALD_plot     <- ALD_plot     + coord_cartesian(ylim = ylim_range)

combined_plot <- (plot_surface | ALD_plot) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", tag_suffix = ") ") +
  theme(legend.position = "below")
combined_plot
