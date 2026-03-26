# thawing and mineralisation
# clean version of computing thawed N, the N losses through denitrification and runoff, 
# and the increased mineralisation due to higher soil temp;
library(terra)
options(repos = c(CRAN = "https://cloud.r-project.org"))
#install.packages("tidyverse")
library(ggplot2)
#.rs.restartR()
# Load NetCDF files
ALD <- rast("mean_ssp585_corr_new.nc")

# For tundra, barren, and wetlands: set to NA where ALD < 1
#capped_ALD <- ifel(ALD > 1, 1, ALD)
#capped_ALD <- ifel(ALD > 3, 3, ALD)

#ALD<-capped_ALD
#rm(capped_ALD)
plot(ALD[[250]])
N_data <- rast("TN_30deg_corr.nc", lyr = 1)
coords <- crds(N_data, na.rm = TRUE)
min(coords[,2])
max(coords[,1])
plot(N_data)
LC <- rast("LC_remapnn_corr.nc")
# Set common extent
common_extent <- ext(-179.95, 179.95, 30, 90)
ext(ALD) <- ext(N_data) <- ext(LC) <- common_extent
LC <- resample(LC, N_data, method = "near")

test_ALD <- ALD[[1:250]] # 1980 - 2099


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


#######
thawed_taiga <- compute_thawed_N(test_ALD, taiga_N,
                                 a=params$taiga['a'],
                                 b=params$taiga['b'],
                                 k=params$taiga['k'])

#taiga_thawed_total <- writeRaster(thawed_taiga, "taiga_total_thawed_N_585_3m_lim_mean.tif")
#rm(taiga_thawed_total,thawed_taiga )
thawed_tundra <- compute_thawed_N(test_ALD, tundra_N,
                                  a=params$tundra['a'],
                                  b=params$tundra['b'],
                                  k=params$tundra['k'])

#tundra_thawed_total <- writeRaster(thawed_tundra, "tundra_total_thawed_N_585_3m_lim_mean.tif")
#plot(tundra_thawed_total)
#rm(thawed_tundra,tundra_thawed_total )
thawed_barren <- compute_thawed_N(test_ALD, barren_N,
                                  a=params$barren['a'],
                                  b=params$barren['b'],
                                  k=params$barren['k'])

#barren_thawed_total <- writeRaster(thawed_barren, "barren_total_thawed_N_585_3m_lim_mean.tif")
#rm(barren_thawed_total, thawed_barren, tundra_N_mat_2d, barren_N_mat_2d, taiga_N_mat_2d)
thawed_wetlands <- compute_thawed_N_wetlands(test_ALD, wetlands_N)

#wetlands_thawed_total <- writeRaster(thawed_wetlands$total_thawed, "wetlands_total_thawed_N_585_3m_lim_mean.tif")
#rm(wetlands_thawed_total, thawed_wetlands)

library(terra)
tundra_thawed_total<-thawed_tundra
taiga_thawed_total<-thawed_taiga
barren_thawed_total<-thawed_barren
wetlands_thawed_total<-thawed_wetlands$total_thawed

#tundra_thawed_total<-rast("tundra_total_thawed_N_585_3m_lim_mean.tif")
#taiga_thawed_total<-rast("taiga_total_thawed_N_585_3m_lim_mean.tif")
#barren_thawed_total<-rast("barren_total_thawed_N_585_3m_lim_mean.tif")
#wetlands_thawed_total<-rast("wetlands_total_thawed_N_585_3m_lim_mean.tif")

taiga_thawed_total[is.na(taiga_thawed_total)] <- 0
tundra_thawed_total[is.na(tundra_thawed_total)] <- 0
barren_thawed_total[is.na(barren_thawed_total)] <- 0
wetlands_thawed_total[is.na(wetlands_thawed_total)] <- 0
# Combine the rasters for each year
combined_thawed <- tundra_thawed_total + barren_thawed_total + wetlands_thawed_total+taiga_thawed_total
rm(tundra_thawed_total,barren_thawed_total,wetlands_thawed_total, taiga_thawed_total)
#combined_thawed <- tundra_thawed_total + barren_thawed_total + wetlands_thawed_total
combined_thawed <- ifel(combined_thawed == 0, NA, combined_thawed)

writeRaster(combined_thawed, "3m_lim/arctic_total_thawed_585_3m_lim_mean.tif")
#rm(combined_thawed)


