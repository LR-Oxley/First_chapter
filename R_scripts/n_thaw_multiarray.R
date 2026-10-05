library(terra)
terra::gdal(drivers=TRUE)
# -----------------------------
# Configuration
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)
ssp <- args[1]
#ssp<-585 # if just one ssp
ALD_path <- sprintf("60N/std_ALD_ssp%s_60deg.nc", ssp)
output_dir <- "total_thawed_extended"
dir.create(output_dir, showWarnings = FALSE)

N_data <- rast("TN_30deg_corr.nc", lyr = 1)
LC <- rast("60N/LC_60deg.nc")
common_extent <- ext(-179.95, 179.95, 60, 90)
N_data <- crop(N_data, common_extent)
LC <- crop(LC, common_extent)
plot(LC)
# -----------------------------
# Land cover masks
# -----------------------------
taiga_classes <- c(1,2,3,4,5,8,9)
tundra_classes <- c(6,7,10,12,14)
wetlands_classes <- 11
barren_classes <- c(13,15,16)

# parameters in m
# Land cover parameters
params <- list(
     taiga = c(a = 0.007, b = 0.097, k = 2.7),
     tundra = c(a = 0.01, b = 0.017, k = 1.9),
     barren = c(a = 0, b = 0.0161, k = 1.6)
  )


f_inorg <- 0.01

# -----------------------------
# Functions
# -----------------------------
normalize_N <- function(N, a, b, k) {
  A_3m <- 3 * a + (b / k) * (1 - exp(-3 * k))
  N / A_3m
}

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

# N (kg/m² over 3 m) × A_ALD (unitless) = kg/m²

# -----------------------------
# Process SSP
# -----------------------------
ALD_CMIP6 <- rast(ALD_path)
ALD_ESA_REF <- rast("ESA_ALD_present_day.nc", lyr=1) # lyr 1 for average, lyr 2 for std
ALD_ESA_REF <- crop(ALD_ESA_REF, common_extent)
plot(ALD_ESA_REF)
present_idx <- which(time(ALD_CMIP6) >= 2000 & time(ALD_CMIP6) <= 2020)
ALD_present <- ALD_CMIP6[[present_idx]]
ALD_present_mean <- mean(ALD_present, na.rm = TRUE)
bias <- ALD_present_mean - ALD_ESA_REF
plot(bias)

# bias correction of present-day ALD
ALD_corr <- ALD_CMIP6 - bias
ALD_corr[ALD_corr < 0] <- 0


plot(ALD_corr[[2]])
# limit ALD_corr to 3 m depth:
#ALD_corr <- ifel(ALD_corr > 3, 3, ALD_corr)
ALD_corr <- crop(ALD_corr, common_extent)
LC_res <- resample(LC, N_data, method = "near")

taiga_mask <- LC_res %in% taiga_classes
tundra_mask <- LC_res %in% tundra_classes
wetlands_mask <- LC_res == wetlands_classes
barren_mask <- LC_res %in% barren_classes

taiga_N <- normalize_N(N_data * taiga_mask, params$taiga["a"], params$taiga["b"], params$taiga["k"])
tundra_N <- normalize_N(N_data * tundra_mask, params$tundra["a"], params$tundra["b"], params$tundra["k"])
barren_N <- normalize_N(N_data * barren_mask, params$barren["a"], params$barren["b"], params$barren["k"])
wetlands_N <- (N_data * wetlands_mask) / 3

thawed_taiga <- compute_thawed_N(ALD_corr, taiga_N, params$taiga["a"], params$taiga["b"], params$taiga["k"], f_inorg)
thawed_tundra <- compute_thawed_N(ALD_corr, tundra_N, params$tundra["a"], params$tundra["b"], params$tundra["k"], f_inorg)
thawed_barren <- compute_thawed_N(ALD_corr, barren_N, params$barren["a"], params$barren["b"], params$barren["k"], f_inorg)
thawed_wetlands <- compute_thawed_N_wetlands(ALD_corr, wetlands_N, f_inorg)

thawed_taiga[is.na(thawed_taiga)] <- 0
thawed_tundra[is.na(thawed_tundra)] <- 0
thawed_barren[is.na(thawed_barren)] <- 0
thawed_wetlands[is.na(thawed_wetlands)] <- 0

combined_thawed <- thawed_taiga + thawed_tundra + thawed_barren + thawed_wetlands

# Keep true land/N-data cells, but allow thawed N to be 0
domain_mask <- ifel(!is.na(N_data) & N_data > 0, 1, NA)
plot(domain_mask)
combined_thawed <- terra::mask(
  combined_thawed,
  domain_mask
)

plot(combined_thawed[[250]])
time_dates <- seq(
  as.Date("1850-07-01"),
  by = "1 year",
  length.out = nlyr(combined_thawed)
)

time(combined_thawed) <- time_dates

writeCDF(
  combined_thawed,
  file.path(output_dir, paste0("arctic_total_thawed_", ssp, "_60N_std.nc")),
  overwrite = TRUE,
  varname="total_thawed_N",
  longname="Total thawed soil nitrogen", unit="kgN_m2"
  )



cat("SSP", ssp, "finished.\n")
