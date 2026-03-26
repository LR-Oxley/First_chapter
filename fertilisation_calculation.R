# fertilisation calculation
#mineralised N:
library(terra)
library(here)
min_126_mean<-rast("1.9%_min/3m_lim/mineralised_N_pool_126_3m_lim_mean.tif")
min_245_mean<-rast("1.9%_min/3m_lim/mineralised_N_pool_245_3m_lim_mean.tif")
min_370_mean<-rast("1.9%_min/3m_lim/mineralised_N_pool_370_3m_lim_mean.tif")
min_585_mean<-rast("1.9%_min/3m_lim/mineralised_N_pool_585_3m_lim_mean.tif")

min_126_std<-rast("1.9%_min/3m_lim/mineralised_N_pool_126_3m_lim_std.tif")
min_245_std<-rast("1.9%_min/3m_lim/mineralised_N_pool_245_3m_lim_std.tif")
min_370_std<-rast("1.9%_min/3m_lim/mineralised_N_pool_370_3m_lim_std.tif")
min_585_std<-rast("1.9%_min/3m_lim/mineralised_N_pool_585_3m_lim_std.tif")

plot(min_585_std[[2]])

# differentiate between tundra and taiga biomes
LC <- rast("LC_remapnn_corr.nc")
# Reproject LC if needed
LC <- project(LC, min_126_std, method = "near")
plot(min_126_mean[[2]])

# 2) Force EXACT same grid as template (resolution + origin + nrow/ncol)
LC_aligned <- resample(LC, min_126_std, method = "near")

# 3) Crop/mask just in case (usually resample already matches)
LC_aligned <- crop(LC_aligned, min_126_std)
# Crop to identical extent
LC <- crop(LC, min_126_std)
plot(LC_aligned)
# Final geometry check
compareGeom(min_126_std, LC, stopOnError = TRUE)
# Ensure masks are numeric (not logical)
taiga_mask <- as.numeric(LC %in% c(1,2,3,4,5,8,9))
tundra_mask <- as.numeric(LC %in% c(6,7,10, 11, 12, 13, 14, 15, 16))

tundra_min_126 <- mask(min_126_std, tundra_mask, maskvalue = 0)
tundra_min_245 <- mask(min_245_std, tundra_mask, maskvalue = 0)
tundra_min_370 <- mask(min_370_std, tundra_mask, maskvalue = 0)
tundra_min_585 <- mask(min_585_std, tundra_mask, maskvalue = 0)
plot(tundra_min_585[[2]])
r1 <- tundra_min_585[[1]]

# area per cell (m²) in lon/lat
area_m2 <- cellSize(r1, unit = "m")

# keep only cells that are land (non-NA in layer 1)
area_land_m2 <- mask(area_m2, r1)
plot(area_land_m2)
# sum area over land cells
total_land_area_m2 <- global(area_land_m2, "sum", na.rm = TRUE)[1,1]

total_land_area_m2
tundra<-total_land_area_m2/(10^12)

#rm(tundra_min_126, tundra_min_370, tundra_min_370, tundra_min_585)
# first tundra: 
area_rast <- cellSize(area_land_m2[[1]], unit = "m")
tundra_mean_kg_m2_yr <- sapply(1:nlyr(tundra_min_126), function(i) {
  
  min_n <- values(tundra_min_126[[i]], mat = FALSE)
  area <- values(area_rast, mat = FALSE)
  
  remove_missing_pixels <- !is.na(min_n) & !is.na(area)
  
  sum(min_n[remove_missing_pixels] * area[remove_missing_pixels]) / sum(area[remove_missing_pixels])
})

tundra_126<-tundra_mean_kg_m2_yr

rm(tundra_min_126)
tundra_min<-cbind(tundra_585, tundra_370, tundra_245, tundra_126)
tundra_min<-data.frame(tundra_min)
tundra_min$Year <- 1850:2099
tail(tundra_min)
colnames(tundra_min)<-c("Min_370", "Min_585","Min_245", "Min_126", "Year")




results <- list(
  period_1900 = colMeans(tundra_min[tundra_min$Year %in% 1900, 1:4]),
  period_2000 = colMeans(tundra_min[tundra_min$Year %in% 2000, 1:4]),
  period_2099 = colMeans(tundra_min[tundra_min$Year %in% 2099, 1:4])
)

# Print the results
results

# Calculate the difference (additional Nin soil due to permafrost thawing (difference 2080-2100 to 1990-2010)):
# this is in kg / m2
#additional_N <- mean_2080_2100 - mean_2000_2020
additional_N <-results$period_2099-results$period_2000
additional_N<-additional_N / 100
add_N_gm2<-additional_N*1000 # in g / m2* yr

write.csv(add_N_gm2, "min_n_bioavailable_tundra_std.csv")

permafrost_inorg_N_mean<-read.csv("min_n_bioavailable_tundra_mean.csv") # g_N_yr
permafrost_inorg_N_std<-read.csv("min_n_bioavailable_tundra_std.csv") # g_N_yr

n_sensitivity_tundra<-1.11 # %_biomass_increase_per_g_N

potential_increase<-n_sensitivity_tundra*permafrost_inorg_N$x #percent_biomass_increase_per_yr


#



# aboveground biomass productivity increase through N fertilisation:
# 2.66 g biomass per gram N increased productivity 
add_N_gm2_mean<-read.csv("min_n_bioavailable_tundra_mean.csv")
add_N_gm2_std<-read.csv("min_n_bioavailable_tundra_std.csv")


agb_tundra_mean<-add_N_gm2_mean[,2]*2.59
agb_tundra_std<-add_N_gm2_std[,2]*2.59
agb_Tundra_upper<-agb_tundra_mean+agb_tundra_std
agb_Tundra_lower<-agb_tundra_mean-agb_tundra_std
agb_Tundra<-rbind(agb_Tundra_upper, agb_tundra_mean, agb_Tundra_lower)

bgb_Tundra_mean<-add_N_gm2_mean[,2]*0.5
bgb_Tundra_std<-add_N_gm2_std[,2]*0.5

bgb_Tundra_upper<-bgb_Tundra_mean+bgb_Tundra_std
bgb_Tundra_lower<-bgb_Tundra_mean-bgb_Tundra_std
bgb_Tundra<-rbind(bgb_Tundra_upper, bgb_Tundra_mean, bgb_Tundra_lower)
rownames(bgb_Tundra)<-c("upper","mean","lower")


#mean biomass in g C / m2
increase_Tundra<-agb_Tundra+bgb_Tundra

total_biomass_Tundra_mean_yearly<-increase_Tundra*0.5
write.csv(total_biomass_Tundra_mean_yearly, "total_biomass_Tundra_mean_yearly.csv")

total_biomass_Tundra_mean_yearly<-read.csv("total_biomass_Tundra_mean_yearly.csv")

### taiga: 

taiga_min_126 <- mask(min_126_mean, taiga_mask, maskvalue = 0)
taiga_min_245 <- mask(min_245_mean, taiga_mask, maskvalue = 0)
taiga_min_370 <- mask(min_370_mean, taiga_mask, maskvalue = 0)
taiga_min_585 <- mask(min_585_mean, taiga_mask, maskvalue = 0)
r1 <- taiga_min_585[[1]]
# area per cell (m²) in lon/lat
area_m2 <- cellSize(r1, unit = "m")

# keep only cells that are land (non-NA in layer 1)
area_land_m2 <- mask(area_m2, r1)

# sum area over land cells
total_land_area_m2 <- global(area_land_m2, "sum", na.rm = TRUE)[1,1]
total_land_area_m2
taiga<-total_land_area_m2/(10^12)
plot(taiga_min_585[[2]])
area_rast <- cellSize(taiga_min_126[[1]], unit = "m")
taiga_mean_kg_m2_yr <- sapply(1:nlyr(taiga_min_126), function(i) {
  
  min_n <- values(taiga_min_126[[i]], mat = FALSE)
  area <- values(area_rast, mat = FALSE)
  
  remove_missing_pixels <- !is.na(min_n) & !is.na(area)
  
  sum(min_n[remove_missing_pixels] * area[remove_missing_pixels]) / sum(area[remove_missing_pixels])
})

taiga_126<-taiga_mean_kg_m2_yr
#

taiga_min<-cbind(taiga_585, taiga_370, taiga_245, taiga_126)
taiga_min<-data.frame(taiga_min)
taiga_min$Year <- 1850:2099

colnames(taiga_min)<-c("Min_585","Min_370", "Min_245", "Min_126", "Year")


results <- list(
  period_1900 = colMeans(taiga_min[taiga_min$Year %in% 1900, 1:4]),
  period_2000 = colMeans(taiga_min[taiga_min$Year %in% 2000, 1:4]),
  period_2099 = colMeans(taiga_min[taiga_min$Year %in% 2099, 1:4])
)

# Print the results
results

# Calculate the difference (additional Nin soil due to permafrost thawing (difference 2080-2100 to 1990-2010)):
# this is in kg / m2
#additional_N <- mean_2080_2100 - mean_2000_2020
additional_N <-results$period_2099-results$period_2000
additional_N<-additional_N / 100
add_N_gm2<-additional_N*1000 # in g / m2


write.csv(add_N_gm2, "min_n_bioavailable_taiga_mean.csv")

# aboveground biomass productivity increase through N fertilisation:
add_N_gm2_mean<-read.csv("min_n_bioavailable_taiga_mean.csv")
add_N_gm2_std<-read.csv("min_n_bioavailable_taiga_std.csv")

agb_taiga_mean<-add_N_gm2_mean[,2]*14.1
add_N_gm2_std <- rev(add_N_gm2_std[,2])
agb_taiga_std<-add_N_gm2_std*14.1
agb_taiga_upper<-agb_taiga_mean+agb_taiga_std
agb_taiga_lower<-agb_taiga_mean-agb_taiga_std
agb_taiga<-rbind(agb_taiga_upper, agb_taiga_mean, agb_taiga_lower)

bgb_taiga_mean<-agb_taiga_mean*0.2
bgb_taiga_std<-agb_taiga_std*0.2

bgb_taiga_upper<-bgb_taiga_mean+bgb_taiga_std
bgb_taiga_lower<-bgb_taiga_mean-bgb_taiga_std
bgb_taiga<-rbind(bgb_taiga_upper, bgb_taiga_mean, bgb_taiga_lower)
rownames(bgb_taiga)<-c("upper","mean","lower")

# no multiplication by 0.5 as 14.1 value is 14.1 g C per g N;

#mean biomass in g C / m2
total_biomass_taiga_mean_yearly<-agb_taiga+bgb_taiga

write.csv(total_biomass_taiga_mean_yearly, "total_biomass_taiga_mean_yearly.csv")



# permafrost-n induced c sink

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
      "V1" = "SSP5-8.5",
      "V2" = "SSP3-7.0",
      "V3" = "SSP2-4.5",
      "V4" = "SSP1-2.6"
    ),
    Scenario = factor(
      Scenario,
      levels = c("SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5")
    )
  )
library(ggplot2)


ggplot(df, aes(x = Scenario, y = mean, fill = Biome)) +
  geom_col(position = "dodge") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  facet_wrap(~Biome)+
  theme_minimal()+
  labs(
    x = "SSP",
    y = "Permafrost- inorganic N induced NPP",
    fill = "Biome"
  )


dodge <- position_dodge(width = 0.7)


ggplot(df, aes(x = Scenario, y = mean, fill = Biome)) +
  geom_col(position = dodge, width = 0.6) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    position = dodge,
    width = 0.2
  ) +
  geom_text(
    aes(label = round(mean, 2)),
    position = dodge,
    hjust = -0.3,   # move text to the left of bar
    vjust = -0.2,
    size = 3.5
  ) +
  theme_minimal() +
  scale_fill_manual(
    values = c(
      "Tundra" = "#C7D6C1", 
      "Taiga"  = "#2F6B4F" 
    ))+
  labs(
    x = "SSP",
    y = expression("g C m"^{-2}*" yr"^{-1}),
    fill = "Biome"
  )




### to get Tg C / yr

min_585_mean<-rast("1.9%_min/3m_lim/mineralised_N_pool_585_3m_lim_mean.tif")
plot(min_585_mean[[2]])

# differentiate between tundra and taiga biomes
LC <- rast("LC_remapnn_corr.nc")
# Reproject LC if needed
LC <- project(LC, min_585_mean, method = "near")
plot(min_126_mean[[2]])

# 2) Force EXACT same grid as template (resolution + origin + nrow/ncol)
LC_aligned <- resample(LC, min_585_mean, method = "near")

# 3) Crop/mask just in case (usually resample already matches)
LC_aligned <- crop(LC_aligned, min_585_mean)
# Crop to identical extent
LC <- crop(LC, min_585_mean)
plot(LC_aligned)
# Final geometry check
# Ensure masks are numeric (not logical)
taiga_mask <- as.numeric(LC %in% c(1,2,3,4,5,8,9))
tundra_mask <- as.numeric(LC %in% c(6,7,10, 11, 12, 13, 14, 15, 16))

tundra_min_585 <- mask(min_585_mean, tundra_mask, maskvalue = 0)
plot(tundra_min_585[[2]])
area_rast <- cellSize(tundra_min_585[[1]], unit="m")
tundra_area_m2 <- global(area_rast * !is.na(tundra_min_585[[1]]), 
                         "sum", na.rm = TRUE)[1]


# conversion factor
g_to_Tg <- 1e-12

tundra_gCyr <- as.matrix(tundra_NPP[, 2:5]) * tundra_area_m2$sum 
tundra_TgCyr<-tundra_gCyr* g_to_Tg


taiga_min_585 <- mask(min_585_mean, taiga_mask, maskvalue = 0)

area_rast <- cellSize(taiga_min_585[[1]], unit="m")
# taiga mask

taiga_area_m2 <- global(area_rast * !is.na(taiga_min_585[[1]]), 
                        "sum", na.rm = TRUE)[1]
# conversion factor gram to Teragram
g_to_Tg <- 1e-12

taiga_gCyr <- as.matrix(taiga_NPP[, 2:5]) * taiga_area_m2$sum 
taiga_TgCyr<-taiga_gCyr * g_to_Tg

NPP_arctic_TgC_yr<-taiga_TgCyr+ tundra_TgCyr

colnames(NPP_arctic_TgC_yr)<-c("SSP5-8.5", "SSP3-7.0", "SSP2-4.5", "SSP1-2.6")
colnames(taiga_TgCyr)<-c("SSP5-8.5", "SSP3-7.0", "SSP2-4.5", "SSP1-2.6")
colnames(tundra_TgCyr)<-c("SSP5-8.5", "SSP3-7.0", "SSP2-4.5", "SSP1-2.6")

rownames(tundra_TgCyr) <- c("Upper", "Mean", "Lower")
rownames(taiga_TgCyr)  <- c("Upper", "Mean", "Lower")
rownames(NPP_arctic_TgC_yr) <- c("Upper", "Mean", "Lower")

tundra_df <- as.data.frame(tundra_TgCyr) |>
  tibble::rownames_to_column("Quantile") |>
  tidyr::pivot_longer(
    cols = -Quantile,
    names_to = "SSP",
    values_to = "Tundra_TgC_yr"
  )

taiga_df <- as.data.frame(taiga_TgCyr) |>
  tibble::rownames_to_column("Quantile") |>
  tidyr::pivot_longer(
    cols = -Quantile,
    names_to = "SSP",
    values_to = "Taiga_TgC_yr"
  )

arctic_df <- as.data.frame(NPP_arctic_TgC_yr) |>
  tibble::rownames_to_column("Quantile") |>
  tidyr::pivot_longer(
    cols = -Quantile,
    names_to = "SSP",
    values_to = "Arctic_TgC_yr"
  )

NPP_arctic_df <- tundra_df |>
  left_join(taiga_df,  by = c("SSP", "Quantile")) |>
  left_join(arctic_df, by = c("SSP", "Quantile"))

NPP_arctic_df |>
  arrange(SSP, Quantile)
                                                                                                       
write.csv(NPP_arctic_df, "NPP_arctic.csv")

tundra<-tundra_area_m2/10^12
taiga<-taiga_area_m2 / 10^12


library(dplyr)
library(tidyr)
library(ggplot2)

NPP_arctic_df<-read.csv("NPP_arctic.csv")

# data from Ramage et al 2024
present_day <- tibble(
  Quantile = c("Upper", "Mean", "Lower"),
  SSP = "Present day",
  Tundra_TgC_yr = c(295.7, 69.27,14.2),
  Taiga_TgC_yr  = c(539.8, 270.32, 0.9),
  Arctic_TgC_yr = c(835.5,339.59,  156.3)
)

NPP_arctic_df2 <- dplyr::bind_rows(
  dplyr::select(NPP_arctic_df, -X),
  present_day
)

NPP_arctic_df2 <- NPP_arctic_df2 %>%
  mutate(
    Quantile = factor(Quantile, levels = c("Lower", "Mean", "Upper")),
    SSP = factor(SSP, levels = c("Present day" , "SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5"))
  ) %>%
  arrange(Quantile, SSP)

# 1. Reshape data from wide to long format for plotting
df_long <- NPP_arctic_df2 %>%
  pivot_longer(
    cols = c(Tundra_TgC_yr, Taiga_TgC_yr, Arctic_TgC_yr),
    names_to = "Biome",
    values_to = "Value"
  ) %>%
  mutate(Biome = gsub("_TgC_yr", "", Biome)) %>%  # Clean biome names
  mutate(Biome = factor(Biome, levels = c("Tundra", "Taiga", "Arctic")))

# 2. Separate mean values from upper/lower bounds
df_mean <- df_long %>% 
  filter(Quantile == "Mean") %>%
  dplyr::select(SSP, Biome, Mean_Value = Value)

df_bounds <- df_long %>%
  filter(Quantile %in% c("Upper", "Lower")) %>%
  pivot_wider(
    id_cols = c(SSP, Biome),
    names_from = Quantile,
    values_from = Value
  )

# 3. Combine mean with bounds
plot_data <- df_mean %>%
  left_join(df_bounds, by = c("SSP", "Biome")) %>%
  mutate(SSP = factor(SSP, levels = c("Present day", "SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5")))


plot_data_future <- plot_data %>%
  dplyr::filter(SSP != "Present day")

# Alternative: Grouped bars (all SSPs together)
ggplot(plot_data_future, aes(x = SSP, y = Mean_Value, fill = Biome)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.7) +
  geom_errorbar(
    aes(ymin = Lower, ymax = Upper),
    position = position_dodge(width = 0.9),
    width = 0.25,
    linewidth = 0.5
  ) +
  labs(
    title = "",
    subtitle = "",
    x = "",
    y = expression("Tg C yr"^{-1}),
    fill = "Biome"
  ) +scale_fill_manual(
    values = c(
      "Tundra" = "#C7D6C1", 
      "Taiga"  = "#2F6B4F", 
      "Arctic" = "#556270"
    ),
    labels = c(
      "Tundra" = "Tundra",
      "Taiga"  = "Taiga",
      "Arctic" = "Pan-Arctic Total"
    )
  )+
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    axis.text.x = element_text(angle = 0, vjust=-3, size= 14),
    legend.position = "bottom"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))


######## 
# ============================================================
# Produces two plots:
#   1) ΔNPP as % of present-day MEAN baseline NPP
#   2) ΔNPP as % of present-day baseline SD (variability)

# ============================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

# ------------------------------------------------------------
# 1) Read Rodal et al. data & compute baseline mean/sd for tundra/taiga
# ------------------------------------------------------------
excel_sheets("Global_NetPrimaryProduction_rodal_et_al.xlsx")   # inspect

npp  <- read_xlsx("Global_NetPrimaryProduction_rodal_et_al.xlsx", sheet = "NPP_estimates")
sites <- read_xlsx("Global_NetPrimaryProduction_rodal_et_al.xlsx", sheet = "Site_Info")

npp_geo <- npp %>%
  left_join(sites, by = "site_ID")

npp_arctic <- npp_geo %>%
  filter(latitude > 60)

# Tundra-like biomes (grassland excluded)
tundra_df <- npp_arctic %>%
  filter(grepl("tundra|shurbland|marsh|northern|peatland|shrubland|savanna",
               biome, ignore.case = TRUE))

# Taiga-like biomes
taiga_df <- npp_arctic %>%
  filter(grepl("forest", biome, ignore.case = TRUE))

# Use NPP_tot2 if present, else NPP_tot1
npp_tundra_values <- tundra_df %>%
  mutate(NPP = coalesce(NPP_tot2, NPP_tot1)) %>%
  dplyr::select(site_ID, latitude, NPP)

npp_taiga_values <- taiga_df %>%
  mutate(NPP = coalesce(NPP_tot2, NPP_tot1)) %>%
  dplyr::select(site_ID, latitude, NPP)

stats_tundra <- npp_tundra_values %>%
  summarise(
    mean_NPP = mean(NPP, na.rm = TRUE),
    sd_NPP   = sd(NPP, na.rm = TRUE),
    n        = sum(!is.na(NPP))
  )

stats_taiga <- npp_taiga_values %>%
  summarise(
    mean_NPP = mean(NPP, na.rm = TRUE),
    sd_NPP   = sd(NPP, na.rm = TRUE),
    n        = sum(!is.na(NPP))
  )

# Convert to C units (your factor 0.5)
mean_NPP_C_tundra <- stats_tundra[, 1:2] * 0.5
mean_NPP_C_taiga  <- stats_taiga[, 1:2] * 0.5

baseline_tundra_mean <- mean_NPP_C_tundra$mean_NPP
baseline_taiga_mean  <- mean_NPP_C_taiga$mean_NPP

baseline_tundra_sd <- mean_NPP_C_tundra$sd_NPP
baseline_taiga_sd  <- mean_NPP_C_taiga$sd_NPP


# ------------------------------------------------------------
# 2) Read your CMIP6-derived ΔNPP tables (absolute units!)
#    These should NOT be pre-normalized to %.
# ------------------------------------------------------------
taiga_NPP  <- read.csv("total_biomass_taiga_mean_yearly.csv")
tundra_NPP <- read.csv("total_biomass_Tundra_mean_yearly.csv")

# Ensure column names match expected format
colnames(tundra_NPP) <- c("x", "SSP585", "SSP370", "SSP245", "SSP126")
colnames(taiga_NPP)  <- c("x", "SSP585", "SSP370", "SSP245", "SSP126")


# ------------------------------------------------------------
# 3) Reshape wide -> long with columns: biome, scenario, mean, lower, upper
# ------------------------------------------------------------
to_rel_long <- function(df_wide, biome_name) {
  df_wide %>%
    pivot_longer(
      cols = starts_with("SSP"),
      names_to = "scenario",
      values_to = "value"
    ) %>%
    mutate(
      stat = case_when(
        grepl("lower", x, ignore.case = TRUE) ~ "lower",
        grepl("mean",  x, ignore.case = TRUE) ~ "mean",
        grepl("upper", x, ignore.case = TRUE) ~ "upper",
        TRUE ~ NA_character_
      ),
      biome = biome_name
    ) %>%
    dplyr::select(biome, scenario, stat, value) %>%
    pivot_wider(names_from = stat, values_from = value) %>%
    filter(!is.na(mean), !is.na(lower), !is.na(upper))
}

rel_df <- bind_rows(
  to_rel_long(tundra_NPP, "Tundra"),
  to_rel_long(taiga_NPP,  "Taiga")
)

# Keep consistent ordering
rel_df$biome <- factor(rel_df$biome, levels = c("Tundra", "Taiga"))
rel_df$scenario <- factor(rel_df$scenario, levels = c("SSP126", "SSP245", "SSP370", "SSP585"))


# ------------------------------------------------------------
# 4) Denominators per biome
# ------------------------------------------------------------
denoms <- tibble(
  biome = c("Tundra", "Taiga"),
  baseline_mean = c(baseline_tundra_mean, baseline_taiga_mean),
  baseline_sd   = c(baseline_tundra_sd,   baseline_taiga_sd)
)

# ------------------------------------------------------------
# 5) Make TWO datasets (simple scaling; no uncertainty propagation)
# ------------------------------------------------------------

# A) % of present-day MEAN baseline NPP
rel_to_mean_df <- rel_df %>%
  left_join(denoms, by = "biome") %>%
  mutate(
    mean  = (mean  / baseline_mean) * 100,
    lower = (lower / baseline_mean) * 100,
    upper = (upper / baseline_mean) * 100
  ) %>%
  dplyr::select(biome, scenario, mean, lower, upper)

# B) % of present-day baseline SD
rel_to_sd_df <- rel_df %>%
  left_join(denoms, by = "biome") %>%
  mutate(
    mean  = (mean  / baseline_sd) * 100,
    lower = (lower / baseline_sd) * 100,
    upper = (upper / baseline_sd) * 100
  ) %>%
  dplyr::select(biome, scenario, mean, lower, upper)


# ------------------------------------------------------------
# 6) Plot function (your grouped-bar styling)
# ------------------------------------------------------------
plot_grouped <- function(df, ylab_text, title_text = "") {
  dodge <- position_dodge(width = 0.9)
  
  ggplot(df, aes(x = scenario, y = mean, fill = biome)) +
    geom_col(position = dodge, width = 0.7) +
    geom_errorbar(
      aes(ymin = lower, ymax = upper),
      position = dodge,
      width = 0.25,
      linewidth = 0.5
    ) +
    labs(
      title = title_text,
      subtitle = "",
      x = "",
      y = ylab_text,
      fill = "Biome"
    ) +
    scale_fill_manual(
      values = c(
        "Tundra" = "#C7D6C1",
        "Taiga"  = "#2F6B4F"
      ),
      labels = c(
        "Tundra" = "Tundra",
        "Taiga"  = "Taiga"
      )
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.text.x = element_text(angle = 0, vjust = -0.5, size = 14),
      legend.position = "bottom"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
}

# ------------------------------------------------------------
# 7) Make the TWO plots
# ------------------------------------------------------------
p_mean <- plot_grouped(
  rel_to_mean_df,
  ylab_text = "Potential NPP increase (% of present-day mean baseline)"
)

p_sd <- plot_grouped(
  rel_to_sd_df,
  ylab_text = "Potential NPP increase (% of present-day baseline SD)"
)

p_mean
p_sd

### make one plot

# find global y range across BOTH datasets
ymax <- max(
  rel_to_mean_df$upper,
  rel_to_sd_df$upper,
  na.rm = TRUE
)

ymin <- min(
  rel_to_mean_df$lower,
  rel_to_sd_df$lower,
  na.rm = TRUE
)

shared_ylim <- c(ymin, ymax)

plot_grouped <- function(df, ylab_text, title_text = "", ylim = NULL) {
  dodge <- position_dodge(width = 0.9)
  
  ggplot(df, aes(x = scenario, y = mean, fill = biome)) +
    geom_col(position = dodge, width = 0.7) +
    geom_errorbar(
      aes(ymin = lower, ymax = upper),
      position = dodge,
      width = 0.25,
      linewidth = 0.5
    ) +
    labs(
      title = title_text,
      x = "",
      y = ylab_text,
      fill = "Biome"
    ) +
    scale_fill_manual(
      values = c("Tundra"="#C7D6C1","Taiga"="#2F6B4F")
    ) +
    coord_cartesian(ylim = ylim) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.text.x = element_text(angle = 0, vjust = -0.5, size = 14),
      axis.text.y = element_text(angle = 0, vjust = -0.5, size = 14),
      axis.title.y = element_text(size = 14),
      legend.position = "bottom"
    ) 
}

p_mean <- plot_grouped(
  rel_to_mean_df,
  ylab_text = "% of present-day mean NPP",
  title_text = ""
  #ylim = shared_ylim
)

plot(p_mean)

p_sd <- plot_grouped(
  rel_to_sd_df,
  ylab_text = "Standardized response (N-induced potential NPP increase / present day NPP SD)",
  title_text = "",
  ylim = shared_ylim
)

#install.packages("patchwork")
library(patchwork)

(p_mean | p_sd) + plot_layout(guides = "collect")
