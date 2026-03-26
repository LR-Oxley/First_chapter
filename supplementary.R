## Supplementary Information

### spatial mapping of land cover data
library(terra)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(here)
library(viridis)

# Load the landcover raster
landcover_map <- rast("LC_remapnn_corr.nc")

# Define Lambert Azimuthal Equal-Area (LAEA) projection centered on North Pole
laea_proj <- "+proj=laea +lat_0=90 +lon_0=0 +datum=WGS84"

coastlines <- ne_coastline(scale = "medium", returnclass = "sf")

# Convert raster to a data frame
df <- as.data.frame(landcover_map, xy = TRUE)
colnames(df) <- c("x", "y", "landcover")

# Step 1: Create a new column for land cover categories
df$landcover_category <- NA

# Step 2: Map the original land cover values to the new categories
df$landcover_category[df$landcover %in% c(1, 2, 3, 4, 5, 8, 9)] <- "Taiga"
df$landcover_category[df$landcover %in% c(6, 7, 10)] <- "Tundra"
df$landcover_category[df$landcover == 11] <- "Wetlands"
df$landcover_category[df$landcover %in% c(15, 16)] <- "Barren"


# Step 3: Filter the data to include only the Arctic region
arctic_df <- subset(df, y >= 30)

# Step 4: Convert the filtered data to an sf object
arctic_sf <- st_as_sf(arctic_df, coords = c("x", "y"), crs = 4326)  # WGS84

# Step 5: Define the LAEA projection centered on the North Pole
laea_crs <- "+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# Step 6: Reproject the data into the LAEA projection
arctic_sf_laea <- st_transform(arctic_sf, crs = laea_crs)
bbox <- st_bbox(arctic_sf_laea)
print(bbox)

# Define custom limits for the plot
xlim <- c(-4500000, 4500000)  # Adjust these values based on your data
ylim <- c(-4500000, 4500000)  # Adjust these values based on your data

# Step 7: Reproject the coastlines to match the LAEA projection
coastlines_laea <- st_transform(coastlines, crs = laea_crs)

# Step 8: Remove rows with NA landcover values
arctic_sf_laea <- arctic_sf_laea[!is.na(arctic_sf_laea$landcover_category), ]

landcover_colors <- c(
  "Tundra" = "#C7D6C1",  # tundra
  "Taiga" = "#2F6B4F",  # taiga
  "Barren" = "#9E9E9E",  # barren
  "Wetlands" = "#4FA3A5"  # wetlands
)

# Step 10: Create the plot with the LAEA projection and custom colors
lc <-ggplot() +
  geom_sf(data = arctic_sf_laea, aes(color = landcover_category), size = 0.05) +
  scale_color_manual(
    name = "Arctic Biome",
    values = landcover_colors  # Use the defined colors
  ) +
  geom_sf(data = coastlines_laea, color = "black", size = 0.05) +  # Add reprojected coastlines
  coord_sf(crs = laea_crs, xlim = xlim, ylim = ylim) +  # LAEA projection with custom limits
  theme_minimal() +
  theme(
    legend.text = element_text(size = 8),  # Increase legend text size
    legend.title = element_text(size = 8, face = "bold"),  # Increase legend title size and make it bold
    legend.position = "bottom"
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 3)))  # Increase the size of legend color keys







# with ggplot
install.packages("systemfonts")
library(ggplot2)
library(shadowtext)

# Equation
eq <- function(depth, a, b, k) {
  a + b * exp(-k * abs(depth))
}

# Parameters
params <- data.frame(
  biome = c("Taiga", "Tundra", "Barren"),
  a = c(0.007, 0.01, 0),
  b = c(0.097, 0.017, 0.0161),
  k = c(0.027, 0.019, 0.016)
)

# Colors and line widths
colors <- c(
  Taiga = "#2F6B4F",
  Tundra = "#C7D6C1",
  Barren = "#9E9E9E",
  Wetlands = "#4FA3A5"
)

lwds <- c(
  Taiga = 1.1,
  Tundra = 0.9,
  Barren = 0.7,
  Wetlands = 0.5
)

# Depth grid
depth_seq <- seq(0, 300, length.out = 300)

# Build curves for fitted biomes
curve_df <- do.call(rbind, lapply(seq_len(nrow(params)), function(i) {
  data.frame(
    biome = params$biome[i],
    depth = -depth_seq,
    TN = eq(depth_seq, params$a[i], params$b[i], params$k[i])
  )
}))

# Add wetlands as constant line
wetlands_df <- data.frame(
  biome = "Wetlands",
  depth = -depth_seq,
  TN = 0.022
)

# Combine all curves
curve_df <- rbind(curve_df, wetlands_df)

# Split by biome
curve_list <- split(curve_df, curve_df$biome)

# Compute label positions for Taiga, Tundra, Barren
label_df <- do.call(rbind, lapply(setdiff(names(curve_list), "Wetlands"), function(biome) {
  
  df_i <- curve_list[[biome]]
  others <- curve_list[setdiff(names(curve_list), c(biome, "Wetlands"))]
  
  min_dist <- sapply(seq_len(nrow(df_i)), function(i) {
    min(abs(df_i$TN[i] - sapply(others, function(o) o$TN[i])))
  })
  
  idx <- which.max(min_dist)
  
  data.frame(
    biome = biome,
    TN = df_i$TN[idx],
    depth = df_i$depth[idx]
  )
}))

# Wetlands label
wetlands_label <- data.frame(
  biome = "Wetlands",
  TN = 0.022,
  depth = -150
)

# Final label dataframe
label_df_all <- rbind(label_df, wetlands_label)
# ---- Plot ----
TN <- ggplot() +
  # Curves
  geom_path(data = curve_df,
            aes(TN, depth, color = biome, linewidth = biome)) +
  
  # Wetlands curve
  geom_path(data = wetlands_df,
            aes(TN, depth),
            color = colors["Wetlands"],
            linewidth = lwds["Wetlands"]) +
  
  # Straight labels with subtle halo
  # shadowtext::geom_shadowtext(data = label_df_all,
  #                             aes(TN, depth, label = biome, color = biome),
  #                             bg.colour = "white",
  #                             bg.r = 0.12,
  #                             size = 2,
  #                             hjust = -0.1,
  #                             show.legend = FALSE) +
  
  scale_color_manual(values = colors) +
  scale_linewidth_manual(values = lwds) +
  
  coord_cartesian(xlim = c(0, 0.12), ylim = c(-300, 0)) +
  labs(x = "Total Nitrogen [kg/m²]", y = "Depth [cm]") +
  theme_minimal(base_size = 8) +
  theme(legend.position = "none")


combined_plot<-lc /  TN + plot_annotation(tag_levels = "a", tag_suffix = ") ")



ggsave("Fig_S1.png", combined_plot,
       width = 15,
       height = 20,
       units = "cm", 
       #scale = 1.4,
       dpi = 400,  #70
       #device = cairo_pdf
)


combined_plot <- (lc / TN) +
  plot_layout(
    heights = c(1.3, 1),
    guides = "collect"
  ) +
  plot_annotation(tag_levels = "a", tag_suffix = ") ")

lc <- lc + theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5))
TN <- TN + theme(plot.margin = margin(5.5, 5.5, 5.5, 5.5))







# # to see % of total N that has mineralised
# N_total<-rast("no_lim/arctic_total_thawed_585_no_lim.tif")
# N_fast<-rast("arctic_mineralised_N_585_no_lim_mean.tif")
# 
# library(terra)
# 
# Ntot_time <- app(N_total, sum, na.rm = TRUE)
# Nfast_time <- app(N_fast, sum, na.rm = TRUE)
# 
# percent_mineralised_time <- ifel(
#   Ntot_time > 0,
#   (Nfast_time / Ntot_time) * 100,
#   NA
# )
# 
# cell_area <- cellSize(Ntot_time, unit = "km")
# 
# total_N_aw <- global(Ntot_time * cell_area, "sum", na.rm = TRUE)
# fast_N_aw  <- global(Nfast_time * cell_area, "sum", na.rm = TRUE)
# 
# 
# percent_mineralised_arctic_aw <- (fast_N_aw / total_N_aw) * 100
# 
# # Use the LAST layer only
# Ntot_end  <- N_total[[nlyr(N_total)]]
# Nfast_end <- N_fast[[nlyr(N_fast)]]
# 
# cell_area <- cellSize(Ntot_end, unit = "m")
# 
# total_N_aw <- global(Ntot_end * cell_area, "sum", na.rm = TRUE)
# fast_N_aw  <- global(Nfast_end * cell_area, "sum", na.rm = TRUE)
# total_N_pg<- total_N_aw * (10^(-12))
# fast_N_pg<- fast_N_aw * (10^(-12))
# percent_mineralised_arctic_aw <- (fast_N_aw / total_N_aw) * 100
# 
# #8.29 (+/- 7)% for 5.85
# #
# #5% for 3.70
# # 4.4% for 2.45
# #4.39%


# calculate what amount of N is in present-day active layer: 
#thawed N (state) = f(ALD_present, TN_profile)
# Load NetCDF files
ALD <- rast("mean_ssp585_corr_new.nc")
common_extent <- ext(-179.95, 179.95, 30, 90)
ALD_arctic<-crop(ALD,common_extent)
test_ALD <- ALD_arctic[[1:250]] # 1980 - 2099
years <- 1850:2099

present_idx <- which(years >= 2000 & years <= 2020)
#present_idx<-years
ALD_present <- mean(test_ALD[[present_idx]], na.rm = TRUE)
plot(ALD_present)
ALD_mean <- global(ALD_present, mean, na.rm = TRUE)
ALD_mean

# full thaw down to 5m depth: 
#values(ALD_present) <- 5

N_data <- rast("TN_30deg_corr.nc", lyr = 1)
LC <- rast("LC_remapnn_corr.nc")
# Set common extent
ext(ALD_present) <- ext(N_data) <- ext(LC) <- common_extent
LC <- resample(LC, N_data, method = "near")


#plot(test_temp[[2]])
# Create masks
taiga_mask    <- LC %in% c(1,2,3,4,5,8,9)
tundra_mask   <- LC %in% c(6,7,10)
wetlands_mask <- LC == 11
barren_mask   <- LC %in% c(15,16)

# Land cover parameters
params <- list(
  taiga = c(a = 0.007, b = 0.097, k = 0.027),
  tundra = c(a = 0.01, b = 0.017, k = 0.019),
  barren = c(a = 0, b = 0.0161, k = 0.016)
)


# Normalize Nitrogen
normalize_N <- function(N, a, b, k) {
  A_3m <- 3 * a + (b / k) * (1 - exp(-3 * k))
  N / A_3m
}

taiga_N <- normalize_N(N_data * taiga_mask, params$taiga['a'], params$taiga['b'], params$taiga['k'])
tundra_N <- normalize_N(N_data * tundra_mask, params$tundra['a'], params$tundra['b'], params$tundra['k'])
barren_N <- normalize_N(N_data * barren_mask, params$barren['a'], params$barren['b'], params$barren['k'])
wetlands_N <- (N_data * wetlands_mask) / 3

# --- 1️⃣ Compute thawed N ---

compute_thawed_N <- function(ALD, N, a, b, k) {
  A_ALD <- a * ALD + (b / k) * (1 - exp(-k * ALD))
  total_thawed <-ifel(N == 0, NA, N * A_ALD)
}
# --- Wetlands variant ---

compute_thawed_N_wetlands <- function(ALD, N) {
  
  total_thawed <- ifel(N == 0, NA, N * ALD)
  
  
  list(total_thawed = total_thawed)
}



# Taiga
thawed_taiga_present <- compute_thawed_N(
  ALD_present,
  taiga_N,
  a = params$taiga['a'],
  b = params$taiga['b'],
  k = params$taiga['k']
)
# Tundra
thawed_tundra_present <- compute_thawed_N(
  ALD_present,
  tundra_N,
  a = params$tundra['a'],
  b = params$tundra['b'],
  k = params$tundra['k']
)

# Barren
thawed_barren_present <- compute_thawed_N(
  ALD_present,
  barren_N,
  a = params$barren['a'],
  b = params$barren['b'],
  k = params$barren['k']
)

# Wetlands (linear)
thawed_wetlands_present <- compute_thawed_N_wetlands(
  ALD_present,
  wetlands_N
)$total_thawed
thawed_taiga_present[is.na(thawed_taiga_present)] <- 0
thawed_tundra_present[is.na(thawed_tundra_present)] <- 0
thawed_barren_present[is.na(thawed_barren_present)] <- 0
thawed_wetlands_present[is.na(thawed_wetlands_present)] <- 0
# Combine the rasters for each year
combined_thawed <- thawed_tundra_present + thawed_barren_present + thawed_wetlands_present+thawed_taiga_present
#combined_thawed <- tundra_thawed_total + barren_thawed_total + wetlands_thawed_total
TN_active_layer_present <- ifel(combined_thawed == 0, NA, combined_thawed)

weights <- cellSize(TN_active_layer_present, unit = "m")

total_N_present <- global(
  TN_active_layer_present * weights,
  "sum",
  na.rm = TRUE
)[[1]]


# convert kg → Pg (if your units are kg/m²)
total_N_present_Pg <- total_N_present / 1e12
total_N_present_Pg

# present day 2000-2020, 585: 28.3 + / - 11.5 Pg; 
# present day 2000-2020, 370: 28.0 + / - 11.3 Pg; 
# 2-4.5: 25.4+ / - 9.1 Pg; 
# 1-2.6: 25.2+ / - 9.1 Pg; 

# full thaw down to 5m depth: 93.4 Pg; 
x<-(100 / 93.4) * 28.3
y<-93.4*0.3


ALD_df <- as.data.frame(ALD_present, xy = TRUE, na.rm = TRUE)

## for a flat circular projection
# Step 1: Filter the data to include only the Arctic region
arctic_df <- subset(ALD_df)

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
summary(arctic_sf_laea$mean)
# Step 5: Reproject the coastlines to match the LAEA projection
coastlines_laea <- st_transform(coastlines, crs = laea_crs)
arctic_sf_laea <- arctic_sf_laea[!is.na(arctic_sf_laea$mean), ]
# Step 6: Create the plot with the LAEA projection

spatial_plot <- ggplot() +
  geom_sf(data = arctic_sf_laea, aes(color = TN), size = 0.01) +
  
  scale_color_gradientn(
    colours = c("dodgerblue3", "#abd9e9", "orange"),
    #values = scales::rescale(c(0,0.5,1,2,2.5,3,3.5, 4,4.5,5,5.5, 6, 10,12)),
    limits = c(0, 8),
    breaks = seq(0, 8, 1),
    #labels = c(0,0.5,1,2,2.5,3,3.5, 4,4.5,5,5.5, 6, 10,12),
    name = expression("kg N m-2"),
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

# ALD
spatial_plot <- ggplot() +
  geom_sf(data = arctic_sf_laea, aes(color = mean), size = 0.01) +
  
  scale_color_viridis_c(
    name = "ALD [m]",
    limits = c(0, 5),
    breaks = seq(0, 5, 1),
    option = "viridis",
    oob = scales::squish,
    guide = guide_colorbar(
      barheight = unit(2.5, "cm"),
      barwidth  = unit(0.25, "cm"),
      ticks = FALSE
    )
  )  +
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

ggsave("ALD_present_day_mean.png", spatial_plot,
       width = 15,
       height = 15,
       dpi=500,
       units = "cm")




library(terra)
library(here)
library(ggplot2)
library(dplyr)
library(tidyr)
## finding out when ALD reaches 5m depth
rast_mean_585 <- rast("mean_ssp585_corr_new.nc")
rast_mean_370 <- rast("mean_ssp370_corr_new.nc")
rast_mean_245 <- rast("mean_ssp245_corr_new.nc")
rast_mean_126 <- rast("mean_ssp126_corr_new.nc")
rast_std_126<- rast("std_ssp126_corr_new.nc")
rast_std_245<- rast("std_ssp245_corr_new.nc")
rast_std_370<- rast("std_ssp370_corr_new.nc")
rast_std_585<- rast("std_ssp585_corr_new.nc")

ald_5m <- rast_mean_585 >= 5

years <- 1850:2099

first_5m <- app(ald_5m, function(x) {
  if (any(x, na.rm = TRUE)) {
    years[which(x)[1]]
  } else {
    NA
  }
})

year_all_5m <- global(first_5m, fun = "max", na.rm = TRUE)
year_all_5m

sum(is.na(values(first_5m)))

plot(first_5m,
     main = "First year ALD reaches 5 m (SSP5-8.5)",
     col = hcl.colors(20, "YlOrRd"))


n_cells <- global(!is.na(rast_mean_585[[1]]), "sum", na.rm = TRUE)[1,1]
n_cells
plot(rast_mean_585[[2]])

n_5m <- global(ald_5m, "sum", na.rm = TRUE)[,1]

df_frac <- data.frame(
  year = years,
  frac_5m = n_5m / n_cells
)

ggplot(df_frac, aes(x = year, y = frac_5m)) +
  geom_line(linewidth = 1) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent
  ) +
  labs(
    x = "Year",
    y = "Fraction of Arctic area (ALD ≥ 5 m)",
    title = "Fraction of Arctic area with active layer depth ≥ 5 m",
    subtitle = "SSP5–8.5"
  ) +
  theme_bw()


df_frac |>
  subset(frac_5m >= 0.5) |>
  head(1)

df_frac |>
  subset(frac_5m >= 0.75) |>
  head(1)

df_frac |>
  subset(frac_5m >= 0.9) |>
  head(1)


cell_area_585 <- cellSize(rast_mean_585, unit = "m")
cell_area_370 <- cellSize(rast_mean_370, unit = "m")
cell_area_245<- cellSize(rast_mean_245, unit = "m")
cell_area_126 <- cellSize(rast_mean_126, unit = "m")




compute_fraction_area_mean <- function(ald, cell_area, threshold = 5) {
  
  # Threshold
  ald_5m <- ald >= threshold
  
  # Baseline permafrost mask (1850)
  pm_mask <- ald[[1]] < threshold
  
  ald_5m <- mask(ald_5m, pm_mask)
  cell_area_pm <- mask(cell_area, pm_mask)
  
  total_area <- global(cell_area_pm, "sum", na.rm = TRUE)[1,1]
  
  area_5m <- global(ald_5m * cell_area_pm, "sum", na.rm = TRUE)[,1]
  
  data.frame(
    year = 1850:2099,
    frac = area_5m / total_area,
    area_km2 = area_5m
  )
}

df_126 <- compute_fraction_area_mean(rast_mean_126, cell_area_126)
df_126$ssp <- "SSP1-2.6"
rm(rast_mean_126, cell_area_126)
df_245 <- compute_fraction_area_mean(rast_mean_245, cell_area_245)
df_245$ssp <- "SSP2-4.5"
rm(rast_mean_245, cell_area_245)
df_370 <- compute_fraction_area_mean(rast_mean_370, cell_area_370)
df_370$ssp <- "SSP3-7.0"
rm(rast_mean_370, cell_area_370)
df_585 <- compute_fraction_area_mean(rast_mean_585, cell_area_585)
df_585$ssp <- "SSP5-8.5"
rm(rast_mean_585, cell_area_585)
df_all <- rbind(df_126, df_245, df_370, df_585)

library(ggplot2)

ggplot(df_all, aes(year, frac, color = ssp)) +
  geom_line(linewidth = 1) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1)
  ) +
  labs(
    x = "Year",
    y = "Fraction of baseline permafrost area (ALD ≥ 5 m)",
    title = "Permafrost loss across SSP scenarios",
    subtitle = "Multi-model mean active layer depth"
  ) +
  theme_bw()

ggplot(df_all, aes(year, area_km2 / 1e6, color = ssp)) +
  geom_line(linewidth = 1) +
  labs(
    x = "Year",
    y = "Area with ALD ≥ 5 m (million km²)",
    title = "Permafrost thaw area across SSP scenarios"
  ) +
  theme_bw()

library(dplyr)

df_rates <- df_all |>
  filter(year >= 2000) |>
  group_by(ssp) |>
  summarise(
    rate_pct_per_decade =
      coef(lm(frac ~ year))[2] * 10 * 100,
    .groups = "drop"
  )

df_rates






###### with uncertainty
library(terra)


prob_ge_5m <- function(mean_ald, sd_ald, threshold = 5) {
  
  lapp(c(mean_ald, sd_ald), fun = function(mu, sd) {
    
    if (is.na(mu) || is.na(sd)) return(NA_real_)
    
    if (sd == 0) {
      return(as.numeric(mu >= threshold))
    }
    
    1 - pnorm((threshold - mu) / sd)
  })
}

compute_fraction_area_with_sd <- function(mean_ald, sd_ald, cell_area, threshold = 5) {
  
  n_years <- nlyr(mean_ald)
  prob_list <- vector("list", n_years)
  
  # Loop over each year (layer)
  for (i in 1:n_years) {
    mu <- mean_ald[[i]]
    s  <- sd_ald[[i]]
    
    # Compute probability per cell
    prob_list[[i]] <- app(c(mu, s), fun = function(x) {
      m <- x[1]
      sd <- x[2]
      if (is.na(m) || is.na(sd)) return(NA_real_)
      if (sd == 0) return(as.numeric(m >= threshold))
      1 - pnorm((threshold - m)/sd)
    })
  }
  
  # Stack all probability layers
  p <- rast(prob_list)
  
  # Compute expected area and variance per year
  mean_area <- numeric(n_years)
  sd_area   <- numeric(n_years)
  
  for (i in 1:n_years) {
    mean_area[i] <- global(p[[i]] * cell_area, "sum", na.rm = TRUE)[,1]
    var_area     <- global(p[[i]] * (1 - p[[i]]) * (cell_area^2), "sum", na.rm = TRUE)[,1]
    sd_area[i]   <- sqrt(var_area)
  }
  
  # Total permafrost area
  total_area <- global(cell_area, "sum", na.rm = TRUE)[1,1]
  
  data.frame(
    year = 1850:(1850 + n_years - 1),
    mean_area_km2 = mean_area,
    sd_area_km2   = sd_area,
    mean_frac     = mean_area / total_area,
    sd_frac       = sd_area / total_area
  )
}


# Compute thawed fraction over time
df_126 <- compute_fraction_area_with_sd(
  mean_ald  = rast_mean_126,
  sd_ald    = rast_std_126,
  cell_area = cell_area_126,
  threshold = 5
)

df_126$ssp <- "SSP1-2.6"
global(rast_mean_126[[1]] - rast_mean_126[[2]], range, na.rm = TRUE)
global(prob_ge_5m(rast_mean_126, rast_std_126), range, na.rm = TRUE)
plot(df_126$year, df_126$mean_frac, type = "l")
plot(values(rast_mean_126[1, 1, ]), type = "l")


df_245 <- compute_fraction_area_with_sd(
  rast_mean_245,
  rast_std_245,
  cell_area_245
)
df_245$ssp <- "SSP2-4.5"

df_370 <- compute_fraction_area_with_sd(
  rast_mean_370,
  rast_std_370,
  cell_area_370
)
df_370$ssp <- "SSP3-7.0"

df_585 <- compute_fraction_area_with_sd(
  rast_mean_585,
  rast_std_585,
  cell_area_585
)
df_585$ssp <- "SSP5-8.5"

df_all <- rbind(df_126, df_245, df_370, df_585)


library(ggplot2)
# fraction
ggplot(df_all, aes(year, mean_frac, color = ssp, fill = ssp)) +
  geom_ribbon(
    aes(
      ymin = mean_frac - sd_frac,
      ymax = mean_frac + sd_frac
    ),
    alpha = 0.25,
    color = NA
  ) +
  geom_line(linewidth = 1) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(
    x = "Year",
    y = "Fraction of baseline permafrost area (ALD ≥ 5 m)",
    title = "Permafrost thaw with inter-model uncertainty"
  ) +
  theme_bw()


# area
ggplot(df_all, aes(year, mean_area_km2 / 1e6, color = ssp, fill = ssp)) +
  geom_ribbon(
    aes(
      ymin = (mean_area_km2 - sd_area_km2) / 1e6,
      ymax = (mean_area_km2 + sd_area_km2) / 1e6
    ),
    alpha = 0.25,
    color = NA
  ) +
  geom_line(linewidth = 1) +
  labs(
    x = "Year",
    y = "Area with ALD ≥ 5 m (million km²)"
  ) +
  theme_bw()

library(dplyr)

df_rates <- df_all |>
  filter(year >= 2000) |>
  group_by(ssp) |>
  summarise(
    rate_pct_per_decade =
      coef(lm(mean_frac ~ year))[2] * 10 * 100,
    .groups = "drop"
  )

df_rates