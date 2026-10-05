## Supplementary Information

### spatial mapping of land cover data
library(terra)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(here)
library(viridis)
setwd("/storage/homefs/lo21q255/coding_first_chapter")
thawed<-rast("total_thawed_extended/arctic_total_thawed_245_no_lim_mean.nc")
plot(thawed[[250]])
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
# Equation



# Equation (exponential profiles)

eq <- function(depth, a, b, k) {
  
  a + b * exp(-k * abs(depth))
  
}

# Parameters

params <- list(
  
  Taiga  = c(a = 0.007, b = 0.097,  k = 2.7),
  
  Tundra = c(a = 0.01,  b = 0.017,  k = 1.9),
  
  Barren = c(a = 0,     b = 0.0161, k = 1.6)
  
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

# Depth (0 at surface, -3 downward)

depths <- seq(0, -3, length.out = 200)

# Exponential profiles

df_exp <- lapply(names(params), function(type) {
  
  p <- params[[type]]
  
  data.frame(
    
    depth = depths,
    
    value = eq(depths, p["a"], p["b"], p["k"]),
    
    type = type
    
  )
  
}) %>% bind_rows()

# Wetlands: uniform distribution

df_wet <- data.frame(
  
  depth = depths,
  
  value = rep(0.022, length(depths)),
  
  type = "Wetlands"
  
)

# Combine all

df <- bind_rows(df_exp, df_wet)

# Plot

ggplot(df, aes(x = value, y = depth, color = type, size = type)) +
  
  geom_line() +
  
  scale_color_manual(values = colors) +
  
  scale_size_manual(values = lwds) +
  
  scale_y_continuous(limits = c(-3, 0), breaks = seq(0, -3, -0.5)) +
  
  labs(x = "Normalized N", y = "Depth (m)") +
  
  theme_minimal() +
  
  theme(legend.title = element_blank())






TN

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







# to see % of total N that has mineralised
N_total<-rast("total_thawed_extended/arctic_total_thawed_585_no_lim_mean.nc")
N_fast<-rast("mineralised_results/no_ALD_limit/1perc/ALD_temperature/na_filled/arctic_bioavailable_N_pool_245_no_lim_mean.nc")
plot(N_fast[[250]])
library(terra)


# Use the LAST layer only
Ntot_end  <- N_total[[nlyr(N_total)]]
Nfast_end <- N_fast[[nlyr(N_fast)]]

cell_area <- cellSize(Ntot_end, unit = "m")

total_N_aw <- global(Ntot_end * cell_area, "sum", na.rm = TRUE)
fast_N_aw  <- global(Nfast_end * cell_area, "sum", na.rm = TRUE)
total_N_pg<- total_N_aw * (10^(-12))
fast_N_pg<- fast_N_aw * (10^(-12))
percent_mineralised_arctic_aw <- (fast_N_aw / total_N_aw) * 100
percent_mineralised_arctic_aw
# 4.35 %


# calculate what amount of N is in present-day active layer: 
#thawed N (state) = f(ALD_present, TN_profile)
# Load NetCDF files
library(terra)
ALD <- rast("60N/mean_ALD_ssp370_60deg.nc")
getwd()
common_extent <- ext(-179.95, 179.95, 60, 90)
ALD_arctic<-crop(ALD,common_extent)
test_ALD <- ALD_arctic[[1:250]] # 1980 - 2099
years <- 1850:2099

present_idx <- which(years >= 2000 & years <= 2020)
#present_idx<-years
ALD_present <- mean(test_ALD[[present_idx]], na.rm = TRUE)

plot(ALD_present)
ALD_mean <- global(ALD_present, mean, na.rm = TRUE)
ALD_mean

# full thaw down to x m depth: 
ALD_present <- ifel(is.na(ALD_present), NA, 5)
plot(ALD_present)
N_data <- rast("TN_30deg_corr.nc", lyr = 1)
N_data <-crop(N_data,common_extent)
plot(N_data)
LC <- rast("LC_remapnn_corr.nc")
LC <-crop(LC,common_extent)
plot(LC)
# Set common extent
ext(ALD_present) <- ext(N_data) <- ext(LC) <- common_extent
LC <- resample(LC, N_data, method = "near")
plot(LC)



#plot(test_temp[[2]])
# Create masks
taiga_mask    <- LC %in% c(1,2,3,4,5,8,9)
tundra_mask   <- LC %in% c(6,7,10, 12, 14)
wetlands_mask <- LC == 11
barren_mask   <- LC %in% c(13,15,16)

# Land cover parameters
params <- list(
  taiga = c(a = 0.007, b = 0.097, k = 2.7),
  tundra = c(a = 0.01, b = 0.017, k = 1.9),
  barren = c(a = 0, b = 0.0161, k = 1.6)
)



f_inorg <- 0.00789


# Normalize Nitrogen
normalize_N <- function(N, a, b, k) {
  A_3m <- 3 * a + (b / k) * (1 - exp(-3 * k))
  N / A_3m
}

taiga_N <- normalize_N(N_data * taiga_mask, params$taiga['a'], params$taiga['b'], params$taiga['k'])
tundra_N <- normalize_N(N_data * tundra_mask, params$tundra['a'], params$tundra['b'], params$tundra['k'])
barren_N <- normalize_N(N_data * barren_mask, params$barren['a'], params$barren['b'], params$barren['k'])
wetlands_N <- (N_data * wetlands_mask) / 3


compute_thawed_N <- function(ALD, N, a, b, k, f_inorg) {
  A_ALD <- a * ALD + (b / k) * (1 - exp(-k * ALD))
  total_thawed <- ifel(N == 0, NA, N * A_ALD)
  total_thaw_corr <- total_thawed / (1 - f_inorg)
  return(total_thaw_corr)
}


compute_thawed_N_wetlands <- function(ALD, N, f_inorg) {
  total_thawed <- ifel(N == 0, NA, N * ALD)
  total_thaw_corr <- total_thawed / (1 - f_inorg)
  return(total_thaw_corr)
}

# Taiga
thawed_taiga_present <- compute_thawed_N(
  ALD_present,
  taiga_N,
  a = params$taiga['a'],
  b = params$taiga['b'],
  k = params$taiga['k'], 
  f_inorg = 0.00789
)
#plot(thawed_taiga_present)
# Tundra
thawed_tundra_present <- compute_thawed_N(
  ALD_present,
  tundra_N,
  a = params$tundra['a'],
  b = params$tundra['b'],
  k = params$tundra['k'], 
  f_inorg = 0.00789
)

plot(thawed_tundra_present)
# Barren
thawed_barren_present <- compute_thawed_N(
  ALD_present,
  barren_N,
  a = params$barren['a'],
  b = params$barren['b'],
  k = params$barren['k'], 
  f_inorg = 0.00789
  
)

# Wetlands (linear)
thawed_wetlands_present <- compute_thawed_N_wetlands(
  ALD_present,
  wetlands_N, 
  f_inorg = 0.00789
)

thawed_taiga_present[is.na(thawed_taiga_present)] <- 0
thawed_tundra_present[is.na(thawed_tundra_present)] <- 0
thawed_barren_present[is.na(thawed_barren_present)] <- 0
thawed_wetlands_present[is.na(thawed_wetlands_present)] <- 0
# Combine the rasters for each year
combined_thawed <- thawed_tundra_present + thawed_barren_present + thawed_wetlands_present+thawed_taiga_present
TN_active_layer_present <- ifel(combined_thawed == 0, NA, combined_thawed)
plot(TN_active_layer_present)
weights <- cellSize(TN_active_layer_present, unit = "m")
total_N_present <- global(
  TN_active_layer_present * weights,
  "sum",
  na.rm = TRUE
)[[1]]
total_N_present_Pg <- total_N_present / 1e12
total_N_present_Pg

# mean 5-8.5: 34.95
# 3-7.0: 34.75
# 2-4.5: 32.56
# 1-2.6: 32.5 
mean(32.5,32.56, 34.75, 34.95)
x<- sum(32.5, 32.56, 34.75, 34.95)
y<-x/4

# down to 5m: 69.67 Pg; 
# down to 3m: 49.58

# fill nans: down to 5m: 80.76
# down to 3m: 58 Pg




palmtag<- rast("TN_30deg_corr.nc",lyr=1)
weights <- cellSize(palmtag, unit = "m", transform=TRUE)
?cellSize
plot(palmtag)

palmtag_masked<-mask(palmtag, TN_active_layer_present)
plot(palmtag_masked)

total_N_present <- global(
  palmtag * weights,
  "sum",
  na.rm = TRUE
)[[1]]
total_N_present_Pg <- total_N_present / 1e12
total_N_present_Pg
# 49.19
# 61.81238
plot(TN_active_layer_present)
plot(thawed_wetlands_present)

weights <- cellSize(thawed_wetlands_present, unit = "m")
total_N_present <- global(
  thawed_wetlands_present * weights,
  "sum",
  na.rm = TRUE
)[[1]]
# convert kg → Pg 
total_N_wet <- total_N_present / 1e12
total_N_wet

weights <- cellSize(thawed_taiga_present, unit = "m")
total_N_present <- global(
  thawed_taiga_present * weights,
  "sum",
  na.rm = TRUE
)[[1]]
# convert kg → Pg 
total_N_tai <- total_N_present / 1e12
total_N_tai

weights <- cellSize(thawed_tundra_present, unit = "m")
total_N_present <- global(
  thawed_tundra_present * weights,
  "sum",
  na.rm = TRUE
)[[1]]
# convert kg → Pg 
total_N_tun <- total_N_present / 1e12
total_N_tun

weights <- cellSize(thawed_barren_present, unit = "m")
total_N_present <- global(
  thawed_barren_present * weights,
  "sum",
  na.rm = TRUE
)[[1]]

# convert kg → Pg 
total_N_bar <- total_N_present / 1e12
total_N_bar
plot(thawed_barren_present)

weights_palmtag <- cellSize(palmtag_masked, unit = "m", transform = TRUE)
weights_wet <- cellSize(thawed_wetlands_present, unit = "m", transform = TRUE)
weights_tai <- cellSize(thawed_taiga_present, unit = "m", transform = TRUE)
weights_tun <- cellSize(thawed_tundra_present, unit = "m", transform = TRUE)
weights_bar <- cellSize(thawed_barren_present, unit = "m", transform = TRUE)

area_palmtag<-global(weights_palmtag * !is.na(palmtag_masked), "sum", na.rm = TRUE)[[1]]
area_wet <- global(weights_wet * !is.na(thawed_wetlands_present), "sum", na.rm = TRUE)[[1]]
area_tai <- global(weights_tai * !is.na(thawed_taiga_present), "sum", na.rm = TRUE)[[1]]
area_tun <- global(weights_tun * !is.na(thawed_tundra_present), "sum", na.rm = TRUE)[[1]]
area_bar <- global(weights_bar * !is.na(thawed_barren_present), "sum", na.rm = TRUE)[[1]]

density_wet <- total_N_wet / area_wet
density_tai <- total_N_tai / area_tai
density_tun <- total_N_tun / area_tun
density_bar <- total_N_bar / area_bar

df <- tibble(
  biome = c("wetlands", "taiga", "tundra", "barren"),
  total_N = c(total_N_wet, total_N_tai, total_N_tun, total_N_bar),
  area = c(area_wet, area_tai, area_tun, area_bar),
  density = c(density_wet, density_tai, density_tun, density_bar)
) %>%
  mutate(
    frac_N = total_N / sum(total_N),
    frac_area = area / sum(area),
    disproportionality = frac_N / frac_area
  )

df
total_area<-sum(df$area)/10^12
# convert kg → Pg 
total_N_present_Pg <- total_N_present / 1e12
total_N_present_Pg

# present day 2000-2020, 585: 28.3 + / - 11.5 Pg; 
# present day 2000-2020, 370: 28.0 + / - 11.3 Pg; 
# 2-4.5: 25.4+ / - 9.1 Pg; 
# 1-2.6: 25.2+ / - 9.1 Pg; 

# full thaw down to 5m depth: 80.5 Pg; 



p <- params$taiga

a <- p["a"]
b <- p["b"]
k <- p["k"]

A_3m <- 3 * a + (b / k) * (1 - exp(-3 * k))
A_5m <- 5 * a + (b / k) * (1 - exp(-5 * k))

fraction_3_to_5 <- (A_5m - A_3m) / A_5m
fraction_3_to_5









ALD_df <- as.data.frame(TN_active_layer_present, xy = TRUE, na.rm = TRUE)

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
summary(arctic_sf_laea$TN)

# Step 5: Reproject the coastlines to match the LAEA projection
coastlines_laea <- st_transform(coastlines, crs = laea_crs)
arctic_sf_laea <- arctic_sf_laea[!is.na(arctic_sf_laea$TN), ]
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

ggsave("TN_present_day_mean.png", spatial_plot,
       width = 20,
       height = 20,
       dpi=500,
       units = "cm")




library(terra)
library(here)
library(ggplot2)
library(dplyr)
library(tidyr)

ESA<-rast("ESACCI-PERMAFROST-L4-ALT-MODISLST_CRYOGRID-AREA4_PP-fv04.0_regrid_fillednans.nc", lyr=1)
palmtag<-rast("TN_30deg_corr.nc", lyr=1)
total_thawed<-rast("total_thawed_extended/arctic_total_thawed_585_no_lim_mean.nc")
ALD<-rast("mean_ssp585_fillnans.nc", lyr=1)
ESA_ALD<-rast("ESA_ALD_present_day.nc", lyr=1)
plot(ESA)
plot(ALD)
plot(total_thawed[[250]])
plot(palmtag)

compareGeom(palmtag, total_thawed, stopOnError = FALSE)
compareGeom(palmtag, ALD, stopOnError = FALSE)
compareGeom(palmtag, ESA_ALD, stopOnError = FALSE)

missing_tt <- is.na(palmtag) & !is.na(total_thawed)
plot(missing_tt, main = "Cells present in Palmtag but missing in Total_thawed")

missing_tt <- is.na(total_thawed) & !is.na(palmtag)
plot(missing_tt, main = "Cells present in Palmtag but missing in Active layer depth file")

total_thawed_original<-rast("total_thawed_extended/arctic_total_thawed_585_no_lim_mean.nc")
total_thawed<-rast("total_thawed_extended/arctic_total_thawed_585_no_lim_mean_try.nc")

area_raster_original <- cellSize(total_thawed_original, unit = "m")
area_raster <- cellSize(total_thawed, unit = "m")
plot(area_raster)
ts_original <- global(total_thawed_original * area_raster_original, 
                      "sum", na.rm = TRUE)

ts_try <- global(total_thawed * area_raster, 
                 "sum", na.rm = TRUE)
ts_try_pg<-ts_try/10^12
df <- data.frame(
  year = 1850:(1850 + nlyr(total_thawed_original) - 1),
  original = ts_original,
  try = ts_try
)





n_missing <- global(missing_tt, "sum", na.rm = TRUE)
n_total <- global(!is.na(palmtag), "sum", na.rm = TRUE)
percent_missing <- 100 * n_missing / n_total
percent_missing

n_total <- ncell(palmtag)
diff_tt <- global(palmtag != total_thawed, "sum", na.rm = TRUE)
diff_ald <- global(palmtag != ALD, "sum", na.rm = TRUE)
pct_diff_tt <- 100 * diff_tt / n_total
pct_diff_ald <- 100 * diff_ald / n_total
agreement_tt <- 100 - pct_diff_tt
agreement_ald <- 100 - pct_diff_ald

c(
  Palmtag = 100,
  Total_Thawed_agreement = agreement_tt,
  ALD_agreement = agreement_ald
)
global(diff_pt_tt, "mean", na.rm = TRUE)
global(diff_pt_ald, "mean", na.rm = TRUE)
