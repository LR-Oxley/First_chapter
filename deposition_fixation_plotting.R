# deposition and fixation
library(here)
library(dplyr)
# Define an extent covering longitudes -180 to 180, latitudes 60 to 90
library(terra)
library(ggplot2)
ext_sub <- ext(-179.95, 179.95, 60, 90)
dep_126<-rast("dep_1850_2100_ssp126_totalN.nc")
plot(dep_126[[2]])
dep_370<-rast("dep_1850_2100_ssp370_totalN.nc")
dep_585<-rast("dep_1850_2100_ssp585_totalN.nc")
# Crop the raster
dep_126 <- crop(dep_126, ext_sub)
dep_370 <- crop(dep_370, ext_sub)
dep_585 <- crop(dep_585, ext_sub)
plot(dep_126[[170]])

# N fixation
fix_126<-rast("fix_1850_2100_ssp126.nc")
fix_370<-rast("fix_1850_2100_ssp370.nc")
fix_585<-rast("fix_1850_2100_ssp585.nc")
plot(fix_585[[251]])
# Crop the raster
fix_126 <- crop(fix_126, ext_sub)
fix_370 <- crop(fix_370, ext_sub)
fix_585 <- crop(fix_585, ext_sub)
fix_585<-(fix_585)


# time series , deposition
weights <- terra::cellSize(dep_126, unit = "m")
wmean <- function(r, w) {
  num <- global(r * w, "sum", na.rm = TRUE)[[1]]
  den <- global(w,      "sum", na.rm = TRUE)[[1]]
  num / den
}

yearly_weighted <- sapply(1:nlyr(dep_126), function(i) {
  wmean(dep_126[[i]], weights)
})

dep_126<-yearly_weighted

yearly_weighted <- sapply(1:nlyr(dep_370), function(i) {
  wmean(dep_370[[i]], weights)
})

dep_370<-yearly_weighted

yearly_weighted <- sapply(1:nlyr(dep_585), function(i) {
  wmean(dep_585[[i]], weights)
})

dep_585<-yearly_weighted

dep_yr<-cbind(dep_126,dep_370,dep_585)
dep_yr<-data.frame(dep_yr)

dep_yr<- dep_yr* 3600 * 24* 365
dep_yr$Year<-1850:2100

colnames(dep_yr)<-c("SSP 1-2.6", "SSP 3-7.0", "SSP 5-8.5", "Year")
ggplot(dep_yr, aes(x = Year)) +
  geom_line(aes(y =`SSP 1-2.6`, colour = "SSP 1-2.6"), linewidth = 1.1) +
  geom_line(aes(y =`SSP 3-7.0`, colour= "SSP 3-7.0"), linewidth = 0.9) +
  geom_line(aes(y =`SSP 5-8.5`, colour="SSP 5-8.5"), linewidth = 1.1) +
  labs(
    y = "N deposition [kg N / m2 * yr]",
    x = "Year",
    color = "",
    title = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")
# r_stack = fix_585 (or any SpatRaster)
# dates = seq of dates for each layer
years_vec <- 2015:2100


# time series , fixosition
weights <- terra::cellSize(fix_126, unit = "m")
wmean <- function(r, w) {
  num <- global(r * w, "sum", na.rm = TRUE)[[1]]
  den <- global(w,      "sum", na.rm = TRUE)[[1]]
  num / den
}

yearly_weighted <- sapply(1:nlyr(fix_126), function(i) {
  wmean(fix_126[[i]], weights)
})

fix_126<-yearly_weighted

fix_yr<-cbind(fix_126,fix_370,fix_585)
fix_yr<-data.frame(fix_yr)
fix_yr<- fix_yr* 3600 * 24 * 365
fix_yr$Year<-1850:2100

colnames(fix_yr)<-c("SSP 1-2.6", "SSP 3-7.0", "SSP 5-8.5", "Year")

library( dplyr)
ggplot(fix_yr, aes(x = Year)) +
  geom_line(aes(y =`SSP 1-2.6`, colour = "SSP 1-2.6"), linewidth = 1.1) +
  geom_line(aes(y =`SSP 3-7.0`, colour= "SSP 3-7.0"), linewidth = 0.9) +
  geom_line(aes(y =`SSP 5-8.5`, colour="SSP 5-8.5"), linewidth = 1.1) +
  labs(
    y = "N fixation [kg N / m2 * yr]",
    x = "Year",
    color = "",
    title = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")






# mean fixation for present day: 
fix_present<-fix_yr[140:160]

years_vec <- 2015:2100


library(terra, tidyr, dplyr, tibble)
library(ggplot2)
slow_pool_r       <- rast("thawed_organic_N_pool_370_test.tif")
fast_pool_r       <- rast("mineralised_N_pool_370_test.tif")
denit_loss_r <- rast("denitrification_loss_cumulative_370_mean.tif")
runoff_loss_r      <- rast("runoff_loss_cumulative_370_mean.tif")
years<-1850:2099
# Compute cell areas (m²)
cell_area <- cellSize(slow_pool_r, unit = "m")

get_totals_Pg <- function(r_stack) {
  # r_stack = SpatRaster with one layer per year (mask already applied)
  
  # cell area in m²
  cell_area <- terra::cellSize(r_stack, unit = "m")
  
  # multiply N flux (kg N/m²/year) by area → kg N/year
  total_kg <- global(r_stack * cell_area, "sum", na.rm = TRUE)[,1]
  
  # convert kg → Tg
  total_Pg <- total_kg / 1e12
  
  return(as.numeric(total_Pg))  # MUST be numeric vector
}

df_pools<- tibble(
  year = years,
  slow_pool = get_totals_Pg(slow_pool_r),
  mineralised_pool = get_totals_Pg(fast_pool_r),
  denit_loss = get_totals_Pg(denit_loss_r), 
  runoff_loss = get_totals_Pg(runoff_loss_r)
)
names(df_pools)<-c("Year", "Thawed_organic_N_pool", "Mineralised_N_Pool","denit_loss", "runoff_loss")
#names(df_pools)<-c("Year","Total_loss" )

ggplot(df_pools, aes(x = Year)) +
  #geom_line(aes(y = Thawed_organic_N_pool, colour = "Thawed_organic_N_pool"), linewidth = 1.1) +
  geom_line(aes(y = Mineralised_N_Pool, colour= "Mineralised_N_Pool"), linewidth = 0.9) +
  geom_line(aes(y = denit_loss, colour="denit_loss"), linewidth = 1.1) +
  geom_line(aes(y = runoff_loss, colour="runoff_loss"), linewidth = 1.1) +
  labs(
    y = "Nitrogen (Pg N)",
    x = "Year",
    color = "",
    title = "Cumulative N dynamics for SSP585 (1950-2099)"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")
tail(df_pools)

write.csv(df_pools, "final_results/n_pools_losses_126.csv")





### spatial plotting

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
library(viridisLite)



# Fig 1a and b)

# Load the NetCDF file
fBNF <- rast("fix_1850_2100_ssp585.nc")  # replace with your actual file name
#ext_sub <- ext(-179.95, 179.95, 60, 90)
#fBNF <- crop(fBNF, ext_sub)

fBNF_future <- mean(fBNF[[231:251]], na.rm = TRUE)

seconds_per_year <- 60 * 60 * 24 * 365
fBNF_future_yr <- fBNF_future * seconds_per_year

df_fBNF <- as.data.frame(
  fBNF_future_yr,
  xy = TRUE,
  na.rm = TRUE
)

colnames(df_fBNF) <- c("x", "y", "fBNF")

arctic_sf <- st_as_sf(
  df_fBNF,
  coords = c("x", "y"),
  crs = 4326
)
arctic_sf$fBNF[arctic_sf$fBNF == 0] <- NA
laea_crs <- "+proj=laea +lat_0=90 +lon_0=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

arctic_sf_laea <- st_transform(arctic_sf, crs = laea_crs)

coastlines <- ne_coastline(scale = "medium", returnclass = "sf")
coastlines_laea <- st_transform(coastlines, crs = laea_crs)
bbox <- st_bbox(arctic_sf_laea)
print(bbox)
# Define custom limits for the plot
xlim <- c(-3304438, 3304438)  # Adjust these values based on your data
ylim <- c(-3304438, 3304438)  # Adjust these values based on your data
#xlim <- c(-4.5e6, 4.5e6)
#ylim <- c(-4.5e6, 4.5e6)

ggplot() +
  geom_sf(
    data = arctic_sf_laea,
    aes(color = fBNF),
    size = 0.2
  ) +
  scale_color_viridis_c(
    name = "BNF (kg N m⁻² yr⁻¹)",
    option = "viridis", 
    na.value="white"
  ) +
  geom_sf(
    data = coastlines_laea,
    color = "black",
    size = 0.4
  ) +
  coord_sf(
    crs = laea_crs,
    xlim = xlim,
    ylim = ylim,
    expand = FALSE
  ) +
  theme_minimal() +
  labs(
    title = "Future Arctic biological nitrogen fixation (2080–2100, SSP5-8.5)",
    x = NULL, y = NULL
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

summary(arctic_sf$fBNF)
# fig 1 a)


## for a flat circular projection
# Step 1: Filter the data to include only the Arctic region
arctic_df <- subset(df_2080)

# Step 2: Convert the filtered data to an sf object
arctic_sf <- st_as_sf(arctic_df, coords = c("x", "y"), crs = 4326)  # WGS84

# Step 3: Define the LAEA projection centered on the North Pole
laea_crs <- "+proj=laea +lat_0=90 +lon_0=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# Step 4: Reproject the data into the LAEA projection
arctic_sf_laea <- st_transform(arctic_sf, crs = laea_crs)
bbox <- st_bbox(arctic_sf_laea)
print(bbox)
# Define custom limits for the plot
xlim <- c(-4500000, 4500000)  # Adjust these values based on your data
ylim <- c(-4500000, 4500000)  # Adjust these values based on your data
coastlines <- ne_coastline(scale = "medium", returnclass = "sf")

# Step 5: Reproject the coastlines to match the LAEA projection
coastlines_laea <- st_transform(coastlines, crs = laea_crs)
arctic_sf_laea <- arctic_sf_laea[!is.na(arctic_sf_laea$fBNF), ]
# Step 6: Create the plot with the LAEA projection
ggplot() +
  geom_sf(data = arctic_sf_laea, aes(color = fBNF), size = 0.2) +
  scale_color_viridis_c(
    name = "Soil N (kg N / m2)",
  )+
  geom_sf(data = coastlines_laea, color = "black", size = 0.5) +  
  coord_sf(crs = laea_crs, xlim = xlim, ylim = ylim) +  
  theme_minimal() +
  theme(legend.position = "right") +  # Remove legend
  labs(
    x = "Longitude", y = "Latitude",
    title="Estimated present day pan-Arctic nitrogen storage (kg / m2)"
  )
