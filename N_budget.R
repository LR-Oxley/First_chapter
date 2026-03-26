# N budget
# 
library(terra)

# min N 


library(terra)
min_126_mean<-rast("0.01%_min/arctic_min_N_pool_126_no_lim_mean_new.tif")
min_245_mean<-rast("0.01%_min/arctic_min_N_pool_245_no_lim_mean_new.tif")
min_370_mean<-rast("0.01%_min/arctic_min_N_pool_370_no_lim_mean_new.tif")
min_585_mean<-rast("0.01%_min/arctic_min_N_pool_585_no_lim_mean_new.tif")

min_126_std<-rast("0.01%_min/arctic_min_N_pool_126_no_lim_std_new.tif")
min_245_std<-rast("0.01%_min/arctic_min_N_pool_245_no_lim_std_new.tif")
min_370_std<-rast("0.01%_min/arctic_min_N_pool_370_no_lim_std_new.tif")
min_585_std<-rast("0.01%_min/arctic_min_N_pool_585_no_lim_std_new.tif")


calculate_present_day_min <- function(min_mean_stack,
                                      years = 1850:2099,
                                      ref_years = 1880:1900,
                                      pres_years = 2000:2020) {
  
  ref_idx  <- which(years %in% ref_years)
  pres_idx <- which(years %in% pres_years)
  
  if (length(ref_idx) == 0 || length(pres_idx) == 0)
    stop("Reference or present-day years not found in 'years' vector.")
  
  # per-cell means for each period
  ref_mean  <- mean(min_mean_stack[[ref_idx]],  na.rm = TRUE)
  pres_mean <- mean(min_mean_stack[[pres_idx]], na.rm = TRUE)
  
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
  
  num / den
}
# Calculate for each SSP
min_126_present_day_value_mean <- calculate_present_day_min(min_126_mean)
min_245_present_day_value_mean <- calculate_present_day_min(min_245_mean)
min_370_present_day_value_mean <- calculate_present_day_min(min_370_mean)
min_585_present_day_value_mean <- calculate_present_day_min(min_585_mean)

# Calculate for each SSP
min_126_present_day_value_std <- calculate_present_day_min(min_126_std)
min_245_present_day_value_std <- calculate_present_day_min(min_245_std)
min_370_present_day_value_std <- calculate_present_day_min(min_370_std)
min_585_present_day_value_std <- calculate_present_day_min(min_585_std)


# Create the dataframe
min_present_day_df <- data.frame(
  SSP = c("SSP126", "SSP245", "SSP370", "SSP585"),
  present_day_min_mean = c(
    min_126_present_day_value_mean,
    min_245_present_day_value_mean,
    min_370_present_day_value_mean,
    min_585_present_day_value_mean
  ),
  present_day_min_std = c(
    min_126_present_day_value_std,   # Assuming you have std stacks too
    min_245_present_day_value_std,
    min_370_present_day_value_std,
    min_585_present_day_value_std
  )
)

# Convert to tibble for better printing
library(tibble)
min_present_day_df <- as_tibble(min_present_day_df)

# View the result
print(min_present_day_df)

write.csv(min_present_day_df, "min_present_day_budget_no_lim_0001.csv")





# future min

calculate_future_min <- function(min_mean_stack, years = 1850:2099,
                                 y0 = 2000, y1 = 2099) {
  
  i0 <- which(years == y0)
  i1 <- which(years == y1)
  if (length(i0) == 0 || length(i1) == 0) stop("y0/y1 not found in years")
  
  rate_future <- (min_mean_stack[[i1]] - min_mean_stack[[i0]]) / (y1 - y0)
  
  # land mask: TRUE where rate has data
  land_mask <- !is.na(rate_future)
  
  # cell area, then mask to land footprint
  area_m2 <- cellSize(rate_future, unit = "m")
  area_m2 <- mask(area_m2, land_mask, maskvalues = FALSE)  # keep only land cells
  
  # area-weighted mean (equivalent to weighted mean over land)
  num <- global(rate_future * area_m2, "sum", na.rm = TRUE)[[1]]
  den <- global(area_m2, "sum", na.rm = TRUE)[[1]]
  
  num / den
}


area_m2<- cellSize(min_126_mean, unit = "m")
land_mask <- !is.na(min_126_mean)
area_m2 <- mask(area_m2, land_mask, maskvalues = FALSE)
total_area <- global(area_m2, "sum", na.rm = TRUE)
# Calculate for each SSP
min_126_future_value_mean <- calculate_future_min(min_126_mean)
min_245_future_value_mean <- calculate_future_min(min_245_mean)
min_370_future_value_mean <- calculate_future_min(min_370_mean)
min_585_future_value_mean <- calculate_future_min(min_585_mean)

# Calculate for each SSP
min_126_future_value_std <- calculate_future_min(min_126_std)
min_245_future_value_std <- calculate_future_min(min_245_std)
min_370_future_value_std <- calculate_future_min(min_370_std)
min_585_future_value_std <- calculate_future_min(min_585_std)


# Create the dataframe
min_future_df <- data.frame(
  SSP = c("SSP126", "SSP245", "SSP370", "SSP585"),
  future_min_mean = c(
    min_126_future_value_mean,
    min_245_future_value_mean,
    min_370_future_value_mean,
    min_585_future_value_mean
  ),
  future_min_std = c(
    min_126_future_value_std,   # Assuming you have std stacks too
    min_245_future_value_std,
    min_370_future_value_std,
    min_585_future_value_std
  )
)

# Convert to tibble for better printing
library(tibble)
min_future_df <- as_tibble(min_future_df)

# View the result
print(min_future_df)

write.csv(min_future_df, "min_future_budget_no_lim_new_0001.csv")


####### ----------- N deposition and fixation -----------
# Define an extent covering longitudes -180 to 180, latitudes 60 to 90
library(terra)
ext_sub <- ext(-179.95, 179.95, 60, 90)
dep_126<-rast("dep_1850_2100_ssp126_totalN.nc")
dep_370<-rast("dep_1850_2100_ssp370_totalN.nc")
dep_585<-rast("dep_1850_2100_ssp585_totalN.nc")
plot(dep_585[[2]])
# Crop the raster
dep_126 <- crop(dep_126, ext_sub)
dep_370 <- crop(dep_370, ext_sub)
dep_585 <- crop(dep_585, ext_sub)
# convert kg N / m2 * s to kg N / m2 * yr
dep_126<- dep_126*31556926
dep_370<- dep_370*31556926
dep_585<- dep_585*31556926

# N fixation
fix_126<-rast("fix_1850_2100_ssp126.nc")
fix_370<-rast("fix_1850_2100_ssp370.nc")
fix_585<-rast("fix_1850_2100_ssp585.nc")
# Crop the raster
fix_126 <- crop(fix_126, ext_sub)
fix_370 <- crop(fix_370, ext_sub)
fix_585 <- crop(fix_585, ext_sub)
# convert to kg N / m2 * yr
fix_126<- fix_126*31556926
fix_370<- fix_370*31556926
fix_585<- fix_585*31556926


# r_stack = fix_585 (or any SpatRaster)
# dates = seq of dates for each layer
years_vec <- 1850:2100

fix_present<-0.0004 # kg N / m2 * yr

## for present day values
dep_370_present<-dep_370[[150:170]]
fix_370_present<-fix_370[[150:170]]

weights <- terra::cellSize(dep_370_present, unit = "m")
wmean <- function(r, w) {
  num <- global(r * w, "sum", na.rm = TRUE)[[1]]
  den <- global(w,      "sum", na.rm = TRUE)[[1]]
  num / den
}

yearly_weighted_dep <- sapply(1:nlyr(dep_370_present), function(i) {
  wmean(dep_370_present[[i]], weights)
})

weights <- terra::cellSize(fix_370_present, unit = "m")
wmean <- function(r, w) {
  num <- global(r * w, "sum", na.rm = TRUE)[[1]]
  den <- global(w,      "sum", na.rm = TRUE)[[1]]
  num / den
}

yearly_weighted_fix <- sapply(1:nlyr(fix_370_present), function(i) {
  wmean(fix_370_present[[i]], weights)
})

dep_370_present_value_mean <- mean(yearly_weighted_dep, na.rm = TRUE)
fix_370_present_value_mean<- mean(yearly_weighted_fix, na.rm = TRUE)
present_day_fix<-cbind(fix_370_present_value_mean)
present_day_dep<-cbind(dep_370_present_value_mean)
write.csv(present_day_fix, "fix_budget_present_day.csv")
write.csv(present_day_dep, "dep_budget_present_day.csv")



fix_126_present<-fix_126[[150:170]]
fix_370_present<-fix_370[[150:170]]
fix_585_present<-fix_585[[150:170]]

fix_126_future<-fix_126[[231:251]]
fix_370_future<-fix_370[[231:251]]
fix_585_future<-fix_585[[231:251]]


weights <- terra::cellSize(fix_126_present, unit = "m")
wmean <- function(r, w) {
  num <- global(r * w, "sum", na.rm = TRUE)[[1]]
  den <- global(w,      "sum", na.rm = TRUE)[[1]]
  num / den
}

yearly_weighted <- sapply(1:nlyr(fix_126_present), function(i) {
  wmean(fix_126_present[[i]], weights)
})

fix_126_present_value <- mean(yearly_weighted, na.rm = TRUE)


yearly_weighted_370 <- sapply(1:nlyr(fix_370_present), function(i) {
  wmean(fix_370_present[[i]], weights)
})

fix_370_present_value <- mean(yearly_weighted_370, na.rm = TRUE)

yearly_weighted_585 <- sapply(1:nlyr(fix_585_present), function(i) {
  wmean(fix_585_present[[i]], weights)
})

fix_585_present_value <- mean(yearly_weighted_585, na.rm = TRUE)

# future: 

weights <- terra::cellSize(fix_126_future, unit = "m")
wmean <- function(r, w) {
  num <- global(r * w, "sum", na.rm = TRUE)[[1]]
  den <- global(w,      "sum", na.rm = TRUE)[[1]]
  num / den
}

yearly_weighted <- sapply(1:nlyr(fix_126_future), function(i) {
  wmean(fix_126_future[[i]], weights)
})

fix_126_future_value <- mean(yearly_weighted, na.rm = TRUE)


yearly_weighted_370 <- sapply(1:nlyr(fix_370_future), function(i) {
  wmean(fix_370_future[[i]], weights)
})

fix_370_future_value <- mean(yearly_weighted_370, na.rm = TRUE)

yearly_weighted_585 <- sapply(1:nlyr(fix_585_future), function(i) {
  wmean(fix_585_future[[i]], weights)
})

fix_585_future_value <- mean(yearly_weighted_585, na.rm = TRUE)

fix_N <- data.frame(
  fix_126_present = fix_126_present_value,
  fix_126_future  = fix_126_future_value,
  fix_370_present = fix_370_present_value,
  fix_370_future  = fix_370_future_value,
  fix_585_present = fix_585_present_value,
  fix_585_future  = fix_585_future_value
)


library(tibble)
library(dplyr)

fix_added_N <- tibble(
  SSP  = c("SSP126", "SSP370", "SSP585"),
  present = c(fix_126_present_value,
              fix_370_present_value,
              fix_585_present_value),
  future = c(fix_126_future_value,
             fix_370_future_value,
             fix_585_future_value)
) %>%
  mutate(
    additional_N_future  = future - present)

fix_added_Ng <- fix_added_N %>%
  mutate(
    present = present * 1000,      # kg/m² → g/m²
    future  = future  * 1000,
    additional_N = additional_N * 1000
  )

write.csv(fix_added_N, "fix_N_budget.csv")


# total fix


# deposition

dep_126_present<-dep_126[[150:170]]
dep_370_present<-dep_370[[150:170]]
dep_585_present<-dep_585[[150:170]]

dep_126_future<-dep_126[[231:251]]
dep_370_future<-dep_370[[231:251]]
dep_585_future<-dep_585[[231:251]]


weights <- terra::cellSize(dep_126_present, unit = "m")
wmean <- function(r, w) {
  num <- global(r * w, "sum", na.rm = TRUE)[[1]]
  den <- global(w,      "sum", na.rm = TRUE)[[1]]
  num / den
}

yearly_weighted <- sapply(1:nlyr(dep_126_present), function(i) {
  wmean(dep_126_present[[i]], weights)
})

dep_126_present_value <- mean(yearly_weighted, na.rm = TRUE)


yearly_weighted_370 <- sapply(1:nlyr(dep_370_present), function(i) {
  wmean(dep_370_present[[i]], weights)
})

dep_370_present_value <- mean(yearly_weighted_370, na.rm = TRUE)

yearly_weighted_585 <- sapply(1:nlyr(dep_585_present), function(i) {
  wmean(dep_585_present[[i]], weights)
})

dep_585_present_value <- mean(yearly_weighted_585, na.rm = TRUE)

# future: 

weights <- terra::cellSize(dep_126_future, unit = "m")
wmean <- function(r, w) {
  num <- global(r * w, "sum", na.rm = TRUE)[[1]]
  den <- global(w,      "sum", na.rm = TRUE)[[1]]
  num / den
}

yearly_weighted <- sapply(1:nlyr(dep_126_future), function(i) {
  wmean(dep_126_future[[i]], weights)
})

dep_126_future_value <- mean(yearly_weighted, na.rm = TRUE)


yearly_weighted_370 <- sapply(1:nlyr(dep_370_future), function(i) {
  wmean(dep_370_future[[i]], weights)
})

dep_370_future_value <- mean(yearly_weighted_370, na.rm = TRUE)

yearly_weighted_585 <- sapply(1:nlyr(dep_585_future), function(i) {
  wmean(dep_585_future[[i]], weights)
})

dep_585_future_value <- mean(yearly_weighted_585, na.rm = TRUE)

dep_N <- data.frame(
  dep_126_present = dep_126_present_value,
  dep_126_future  = dep_126_future_value,
  dep_370_present = dep_370_present_value,
  dep_370_future  = dep_370_future_value,
  dep_585_present = dep_585_present_value,
  dep_585_future  = dep_585_future_value
)

library(tibble)
library(dplyr)

dep_added_N <- tibble(
  SSP  = c("SSP126", "SSP370", "SSP585"),
  present = c(dep_126_present_value,
              dep_370_present_value,
              dep_585_present_value),
  future = c(dep_126_future_value,
             dep_370_future_value,
             dep_585_future_value)
) %>%
  mutate(
    additional_N_future  = future - present)

dep_added_Ng <- dep_added_N %>%
  mutate(
    present = present * 1000,      # kg/m² → g/m²
    future  = future  * 1000,
    additional_N = additional_N * 1000
  )

write.csv(dep_added_N, "dep_N_budget.csv")






####### ---------- denitrification loss -------

denitrification_loss_126_mean<-rast("3%_min/denitrification_loss_cumulative_126_no_lim_mean_new.tif")
denitrification_loss_245_mean<-rast("3%_min/denitrification_loss_cumulative_245_no_lim_mean_new.tif")
denitrification_loss_370_mean<-rast("3%_min/denitrification_loss_cumulative_370_no_lim_mean_new.tif")
denitrification_loss_585_mean<-rast("3%_min/denitrification_loss_cumulative_585_no_lim_mean_new.tif")

denitrification_loss_126_std<-rast("3%_min/denitrification_loss_cumulative_126_std.tif")
denitrification_loss_245_std<-rast("3%_min/denitrification_loss_cumulative_245_std.tif")
denitrification_loss_370_std<-rast("3%_min/denitrification_loss_cumulative_370_std.tif")
denitrification_loss_585_std<-rast("3%_min/denitrification_loss_cumulative_585_std.tif")

calculate_present_day_runoff <- function(runoff_mean_stack,
                                      years = 1850:2099,
                                      ref_years = 1880:1900,
                                      pres_years = 2000:2020) {
  
  ref_idx  <- which(years %in% ref_years)
  pres_idx <- which(years %in% pres_years)
  
  if (length(ref_idx) == 0 || length(pres_idx) == 0)
    stop("Reference or present-day years not found in 'years' vector.")
  
  # per-cell means for each period
  ref_mean  <- mean(runoff_mean_stack[[ref_idx]],  na.rm = TRUE)
  pres_mean <- mean(denit_mean_stack[[pres_idx]], na.rm = TRUE)
  
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
  
  num / den
}


# Calculate for each SSP
denit_loss_126_present_day_value_mean <- calculate_present_day_denit(denitrification_loss_126_mean)
denit_loss_245_present_day_value_mean <- calculate_present_day_denit(denitrification_loss_245_mean)
denit_loss_370_present_day_value_mean <- calculate_present_day_denit(denitrification_loss_370_mean)
denit_loss_585_present_day_value_mean <- calculate_present_day_denit(denitrification_loss_585_mean)

# Calculate for each SSP
denit_loss_126_present_day_value_std <- calculate_present_day_denit(denitrification_loss_126_std)
denit_loss_245_present_day_value_std <- calculate_present_day_denit(denitrification_loss_245_std)
denit_loss_370_present_day_value_std <- calculate_present_day_denit(denitrification_loss_370_std)
denit_loss_585_present_day_value_std <- calculate_present_day_denit(denitrification_loss_585_std)

# Create the dataframe
denit_loss_present_day_df <- data.frame(
  SSP = c("SSP126", "SSP245", "SSP370", "SSP585"),
  present_day_denit_loss_mean = c(
    denit_loss_126_present_day_value_mean,
    denit_loss_245_present_day_value_mean,
    denit_loss_370_present_day_value_mean,
    denit_loss_585_present_day_value_mean
  ))
  #present_day_denit_loss_std = c(
   # denit_loss_126_present_day_value_std,   
    #denit_loss_245_present_day_value_std,
    #denit_loss_370_present_day_value_std,
    #denit_loss_585_present_day_value_std
  #)
#)

# Convert to tibble for better printing
library(tibble)
denit_loss_present_day_df <- as_tibble(denit_loss_present_day_df)

# View the result
print(denit_loss_present_day_df)

write.csv(denit_loss_present_day_df, "denit_loss_present_day_budget.csv")





# future denit
calculate_future_denit <- function(denit_mean_stack, years = 1850:2099,
                                 y0 = 2000, y1 = 2099) {
  
  i0 <- which(years == y0)
  i1 <- which(years == y1)
  if (length(i0) == 0 || length(i1) == 0) stop("y0/y1 not found in years")
  
  rate_future <- (denit_mean_stack[[i1]] - denit_mean_stack[[i0]]) / (y1 - y0)
  
  # land mask: TRUE where rate has data
  land_mask <- !is.na(rate_future)
  
  # cell area, then mask to land footprint
  area_m2 <- cellSize(rate_future, unit = "m")
  area_m2 <- mask(area_m2, land_mask, maskvalues = FALSE)  # keep only land cells
  
  # area-weighted mean (equivalent to weighted mean over land)
  num <- global(rate_future * area_m2, "sum", na.rm = TRUE)[[1]]
  den <- global(area_m2, "sum", na.rm = TRUE)[[1]]
  
  num / den
}



# Calculate for each SSP
denit_loss_126_future_value_mean <- calculate_future_denit(denitrification_loss_126_mean)
denit_loss_245_future_value_mean <- calculate_future_denit(denitrification_loss_245_mean)
denit_loss_370_future_value_mean <- calculate_future_denit(denitrification_loss_370_mean)
denit_loss_585_future_value_mean <- calculate_future_denit(denitrification_loss_585_mean)

# Calculate for each SSP
denit_loss_126_future_value_std <- calculate_future_denit(denitrification_loss_126_std)
denit_loss_245_future_value_std <- calculate_future_denit(denitrification_loss_245_std)
denit_loss_370_future_value_std <- calculate_future_denit(denitrification_loss_370_std)
denit_loss_585_future_value_std <- calculate_future_denit(denitrification_loss_585_std)


denit_loss_future_df <- data.frame(
  SSP = c("SSP126", "SSP245", "SSP370", "SSP585"),
  future_denit_loss_mean = c(
    denit_loss_126_future_value_mean,
    denit_loss_245_future_value_mean,
    denit_loss_370_future_value_mean,
    denit_loss_585_future_value_mean
  ))
#future_denit_loss_std = c(
# denit_loss_126_future_value_std,   
#denit_loss_245_future_value_std,
#denit_loss_370_future_value_std,
#denit_loss_585_future_value_std
#)
#)


# Convert to tibble for better printing
library(tibble)
denit_loss_future_df <- as_tibble(denit_loss_future_df)

# View the result
print(denit_loss_future_df)

write.csv(denit_loss_future_df, "denit_loss_future_budget.csv")




####### ----------- runoff loss -----------


runoff_loss_126_mean<-rast("3%_min/runoff_loss_cumulative_126_no_lim_mean_new.tif")
runoff_loss_245_mean<-rast("3%_min/runoff_loss_cumulative_245_no_lim_mean_new.tif")
runoff_loss_370_mean<-rast("3%_min/runoff_loss_cumulative_370_no_lim_mean_new.tif")
runoff_loss_585_mean<-rast("3%_min/runoff_loss_cumulative_585_no_lim_mean_new.tif")

runoff_loss_126_std<-rast("3%_min/runoff_loss_cumulative_126_std.tif")
runoff_loss_245_std<-rast("3%_min/runoff_loss_cumulative_245_std.tif")
runoff_loss_370_std<-rast("3%_min/runoff_loss_cumulative_370_std.tif")
runoff_loss_585_std<-rast("3%_min/runoff_loss_cumulative_585_std.tif")

calculate_present_day_runoff <- function(runoff_mean_stack,
                                        years = 1850:2099,
                                        ref_years = 1880:1900,
                                        pres_years = 2000:2020) {
  
  ref_idx  <- which(years %in% ref_years)
  pres_idx <- which(years %in% pres_years)
  
  if (length(ref_idx) == 0 || length(pres_idx) == 0)
    stop("Reference or present-day years not found in 'years' vector.")
  
  # per-cell means for each period
  ref_mean  <- mean(runoff_mean_stack[[ref_idx]],  na.rm = TRUE)
  pres_mean <- mean(runoff_mean_stack[[pres_idx]], na.rm = TRUE)
  
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
  
  num / den
}


# Calculate for each SSP
runoff_loss_126_present_day_value_mean <- calculate_present_day_runoff(runoff_loss_126_mean)
runoff_loss_245_present_day_value_mean <- calculate_present_day_runoff(runoff_loss_245_mean)
runoff_loss_370_present_day_value_mean <- calculate_present_day_runoff(runoff_loss_370_mean)
runoff_loss_585_present_day_value_mean <- calculate_present_day_runoff(runoff_loss_585_mean)

# Calculate for each SSP
runoff_loss_126_present_day_value_std <- calculate_present_day_runoff(runoff_loss_126_std)
runoff_loss_245_present_day_value_std <- calculate_present_day_runoff(runoff_loss_245_std)
runoff_loss_370_present_day_value_std <- calculate_present_day_runoff(runoff_loss_370_std)
runoff_loss_585_present_day_value_std <- calculate_present_day_runoff(runoff_loss_585_std)

# Create the dataframe
runoff_loss_present_day_df <- data.frame(
  SSP = c("SSP126", "SSP245", "SSP370", "SSP585"),
  present_day_runoff_loss_mean = c(
    runoff_loss_126_present_day_value_mean,
    runoff_loss_245_present_day_value_mean,
    runoff_loss_370_present_day_value_mean,
    runoff_loss_585_present_day_value_mean
  ))
#present_day_runoff_loss_std = c(
# runoff_loss_126_present_day_value_std,   
#runoff_loss_245_present_day_value_std,
#runoff_loss_370_present_day_value_std,
#runoff_loss_585_present_day_value_std
#)
#)

# Convert to tibble for better printing
library(tibble)
runoff_loss_present_day_df <- as_tibble(runoff_loss_present_day_df)

# View the result
print(runoff_loss_present_day_df)

write.csv(runoff_loss_present_day_df, "runoff_loss_present_day_budget.csv")





# future runoff
calculate_future_runoff <- function(runoff_mean_stack, years = 1850:2099,
                                   y0 = 2000, y1 = 2099) {
  
  i0 <- which(years == y0)
  i1 <- which(years == y1)
  if (length(i0) == 0 || length(i1) == 0) stop("y0/y1 not found in years")
  
  rate_future <- (runoff_mean_stack[[i1]] - runoff_mean_stack[[i0]]) / (y1 - y0)
  
  # land mask: TRUE where rate has data
  land_mask <- !is.na(rate_future)
  
  # cell area, then mask to land footprint
  area_m2 <- cellSize(rate_future, unit = "m")
  area_m2 <- mask(area_m2, land_mask, maskvalues = FALSE)  # keep only land cells
  
  # area-weighted mean (equivalent to weighted mean over land)
  num <- global(rate_future * area_m2, "sum", na.rm = TRUE)[[1]]
  den <- global(area_m2, "sum", na.rm = TRUE)[[1]]
  
  num / den
}



# Calculate for each SSP
runoff_loss_126_future_value_mean <- calculate_future_runoff(runoff_loss_126_mean)
runoff_loss_245_future_value_mean <- calculate_future_runoff(runoff_loss_245_mean)
runoff_loss_370_future_value_mean <- calculate_future_runoff(runoff_loss_370_mean)
runoff_loss_585_future_value_mean <- calculate_future_runoff(runoff_loss_585_mean)

# Calculate for each SSP
runoff_loss_126_future_value_std <- calculate_future_runoff(runoff_loss_126_std)
runoff_loss_245_future_value_std <- calculate_future_runoff(runoff_loss_245_std)
runoff_loss_370_future_value_std <- calculate_future_runoff(runoff_loss_370_std)
runoff_loss_585_future_value_std <- calculate_future_runoff(runoff_loss_585_std)


runoff_loss_future_df <- data.frame(
  SSP = c("SSP126", "SSP245", "SSP370", "SSP585"),
  future_runoff_loss_mean = c(
    runoff_loss_126_future_value_mean,
    runoff_loss_245_future_value_mean,
    runoff_loss_370_future_value_mean,
    runoff_loss_585_future_value_mean
  ))
#future_runoff_loss_std = c(
# runoff_loss_126_future_value_std,   
#runoff_loss_245_future_value_std,
#runoff_loss_370_future_value_std,
#runoff_loss_585_future_value_std
#)
#)


# Convert to tibble for better printing
library(tibble)
runoff_loss_future_df <- as_tibble(runoff_loss_future_df)

# View the result
print(runoff_loss_future_df)

write.csv(runoff_loss_future_df, "runoff_loss_future_budget.csv")







####### ----------- creating n budget graph -----------

# present day fixation estimate: 4 kg N / ha * yr == 0.0004 kg N / m2⋅yr
fix <-0.0004 * 1000
# present day deposition estimate: 2 kg N / ha *yr == 0.0002 kg N / m2 * yr

# present day: 2000 - 2020: 

min_n_present_day<-read.csv("min_present_day_budget.csv")[,2:3]
min_n_present_day <- min_n_present_day %>%
  summarise(
    mean = mean(present_day_min_mean, na.rm = TRUE),
    sd   = sd(present_day_min_mean, na.rm = TRUE)
  )

fix_present_day<-read.csv("fix_budget_present_day.csv")[,2]
fix_present_day$obs<-c(0.0004)
fix_present_day <- data.frame(
  mean = mean(unlist(fix_present_day), na.rm = TRUE),
  sd   = sd(unlist(fix_present_day), na.rm = TRUE)
)
dep_present_day<-read.csv("dep_budget_present_day.csv")[,2]
dep_present_day$obs<-c(0.0002)
dep_present_day <- data.frame(
  mean = mean(unlist(dep_present_day), na.rm = TRUE),
  sd   = sd(unlist(dep_present_day), na.rm = TRUE)
)
# calculate present day runoff based on Zhao et al 2025: 
# 1.04 Tg N from permafrost region runoff into Arctic ocean; 
runoff_Tg<-1.04 # Tg N / yr 
area_km2 <- 20.4e6 # area from Zhao et al paper
# Convert units
area_m2 <- area_km2 * 1e6
runoff_kg<- runoff_Tg*10^9
runoff_kgN_m2_yr<-runoff_kg/area_m2
mean_min_N_present_day<-min_n_present_day[,2]
f_runoff<-runoff_kgN_m2_yr/mean_min_N_present_day
n_loss_runoff_present_day<- runoff_kgN_m2_yr

# calculate present day denitrification rate: 
denitrification_data<-read.csv("Guo_mengje_denitrification.csv")
mean_rate_nmol_n_g_h<-mean(denitrification_data$nmol_N_g_h)
sd_rate_nmol_n_g_h<-sd(denitrification_data$nmol_N_g_h)

# mean denitrification rate = 3.05 nmol N / g soil * h 
#  die Bodenmasseneinheit umrechnen und wegbekommen 
mean_bulk_density_g_cm3<-0.1 #estimated for the high arctic tundra soils in Guo study

mean_rate_nmol_n_cm3_h<-mean_rate_nmol_n_g_h*mean_bulk_density_g_cm3
sd_rate_nmol_n_cm3_h<-sd_rate_nmol_n_g_h*mean_bulk_density_g_cm3
# get it in m3
mean_rate_nmol_n_m3_h<-mean_rate_nmol_n_cm3_h*(1*10^6)
sd_rate_nmol_n_m3_h<-sd_rate_nmol_n_cm3_h*(1*10^6)
# 1 mol = 1*10^9 nanomol
mean_rate_mol_n_m3_h<-mean_rate_nmol_n_m3_h*(1*10^-9)
sd_rate_mol_n_m3_h<-sd_rate_nmol_n_m3_h*(1*10^-9)
# molecular weight N = 14 g / mol

mean_rate_g_N_m3_h<-mean_rate_mol_n_m3_h*14
sd_rate_g_N_m3_h<-sd_rate_mol_n_m3_h*14
# in kg
mean_rate_kg_N_m3_h<-mean_rate_g_N_m3_h/1000
sd_rate_kg_N_m3_h<-sd_rate_g_N_m3_h/1000
# in square m: multiply by average active layer depth
mean_rate_kg_N_m2_h<-mean_rate_kg_N_m3_h*1.3
sd_rate_kg_N_m2_h<-sd_rate_kg_N_m3_h*1.3
#  h-1 in yr-1 umrechnen
mean_rate_kg_N_m2_day<-mean_rate_kg_N_m2_h*12
sd_rate_kg_N_m2_day<-sd_rate_kg_N_m2_h*12
mean_rate_kg_N_m2_yr<-mean_rate_kg_N_m2_day*120
sd_rate_kg_N_m2_yr<-sd_rate_kg_N_m2_day*120

# denitrification rate acts on the nitrate part (10%) of mineralised N (which is 10% of total N) --> 1% of total N; 
# mean_rate_kg_N_m2_yr is for total N; need to divide by 100 to get the actual rate 
n_loss_denit_present_day<-mean_rate_kg_N_m2_yr*0.01
n_loss_denit_present_day_sd<-sd_rate_kg_N_m2_yr*0.01
n_loss_denit_present_day_g_nitrate<-n_loss_denit_present_day*1000

# mean present day denitrification: 0.8 g N / m2 * yr from Guo Menje et al 2025


n_loss_runoff_present_day <- data.frame(
  mean = n_loss_runoff_present_day,
  sd   = NA
)

n_loss_denit_present_day <- data.frame(
  mean = n_loss_denit_present_day,
  sd   = n_loss_denit_present_day_sd
)
min_n_present_day$Pool  <- "Mineralised N due to permafrost thaw"
dep_present_day$Pool    <- "N Deposition"
fix_present_day$Pool    <- "N Fixation"
n_loss_runoff_present_day$Pool   <- "N loss due to runoff"
n_loss_denit_present_day$Pool   <- "N loss due to denitrification"

present_day <- bind_rows(min_n_present_day, dep_present_day, fix_present_day, n_loss_denit_present_day, n_loss_runoff_present_day)
rownames(present_day)<- c("min_N", "dep", "fix", "denit", "runoff")
colnames(present_day)<-c("present_day", "std", "Pool")
present_day<-present_day[,1:3]
present_day_long <- present_day %>%
  dplyr::transmute(
    Pool  = Pool,
    SSP   = "Present day",
    Type  = "Present day reference period",
    value = present_day, 
    value_std = std
  )

write.csv(present_day_long,"present_day_long_budget.csv")
#### for all ssps
# units denit and runoff already converted

# future anomaly: divide values by 90 
library(ggplot2)
library(dplyr)
rm(min_N, dep, runoff, denit, fix)
min_N<- read.csv("min_future_budget.csv") # in kg / m2 
dep<- read.csv("dep_N_budget.csv") # in kg / m2 
runoff<- read.csv("runoff_loss_future_budget.csv") # in kg / m2 
denit<-read.csv("denit_loss_future_budget.csv") # in kg / m2 
fix<-read.csv("fix_N_budget.csv") # in kg / m2 

library(dplyr)

# Fixation: only SSP + additional_N_future
fix_df <- fix %>%
  dplyr::select(SSP, fix_additional = additional_N_future)

# Deposition: only SSP + additional_N_future
dep_df <- dep %>%
  dplyr::select(SSP, dep_additional = additional_N_future)

# Denitrification: mean + std
denit_df <- denit %>%
  dplyr::select(SSP,
         denit_mean = future_denit_mean,
         denit_std  = future_denit_std)

# Runoff: mean + std
runoff_df <- runoff %>%
  dplyr::select(SSP,
         runoff_mean = future_runoff_mean,
         runoff_std  = future_runoff_std)

# Mineralisation: mean + std
min_df <- min_N %>%
  dplyr::select(SSP,
         min_mean = future_min_mean,
         min_std  = future_min_std)

# Combine everything
N_budget_df <- fix_df %>%
  full_join(dep_df,    by = "SSP") %>%
  full_join(denit_df,  by = "SSP") %>%
  full_join(runoff_df, by = "SSP") %>%
  full_join(min_df,    by = "SSP")

pool_lookup <- c(
  fix    = "N Fixation",
  dep    = "N Deposition",
  denit  = "N loss due to denitrification",
  runoff = "N loss due to runoff",
  min    = "Mineralised N due to permafrost thaw"
)

N_budget_long <- N_budget_df %>%
  pivot_longer(
    cols = -SSP,
    names_to = c("process", "stat"),
    names_sep = "_",
    values_to = "value"
  ) %>%
  mutate(
    Pool = pool_lookup[process],
    Type = case_when(
      stat == "additional" ~ "Future change relative to present",
      stat == "mean"       ~ "Future mean",
      stat == "std"        ~ "Future standard deviation"
    )
  ) %>%
  filter(!is.na(value)) %>%
  dplyr::select(Pool, SSP, Type, value, ) %>%
  arrange(Pool, SSP, Type) %>%
  mutate(X = row_number()) %>%
  dplyr::select(X, Pool, SSP, Type, value)

write.csv(N_budget_long, "N_budget_future.csv")
# one big graph: 

library(tidyr)

present_day_long<-read.csv("present_day_long_budget.csv")
present_day_long<-present_day_long[,2:6]
N_budget_long<-N_budget_long[,2:5]
N_budget <- bind_rows(
  present_day_long,
  N_budget_long
)

future_long <- N_budget_long %>%
  filter(SSP != "Present day") %>%
  filter(Type %in% c(
    "Future mean",
    "Future standard deviation",
    "Future change relative to present"
  )) %>%
  mutate(
    Type = "Future anomaly"
  ) %>%
  group_by(Pool, SSP, Type) %>%
  reframe(
    value = case_when(
      any(Type == "Future anomaly" & cur_group()$Pool %in%
            N_budget$Pool[N_budget$Type == "Future change relative to present"]) ~
        value[match("Future change relative to present", N_budget$Type)],
      any(N_budget$Type == "Future mean") ~
        value[N_budget$Type == "Future mean"],
      TRUE ~ NA_real_
    ),
    value_std = value[N_budget$Type == "Future standard deviation"][1],
    .groups = "drop"
  )

N_budget_long <- bind_rows(
  present_day_long,
  future_long
) %>%
  arrange(Pool, SSP)

N_budget <- N_budget %>%  # <-- your full long table
  pivot_wider(
    names_from  = Type,
    values_from = value
  ) %>%
  mutate(
    value = case_when(
      SSP == "Present day" ~ `Present day reference period`,
      !is.na(`Future mean`) ~ `Future mean`,
      !is.na(`Future change relative to present`) ~ `Future change relative to present`
    ),
    value_std = `Future standard deviation`,
    Type = ifelse(
      SSP == "Present day",
      "Present day reference period",
      "Future anomaly"
    )
  ) %>%
  dplyr::select(Pool, SSP, Type, value, value_std) %>%
  arrange(Pool, SSP)


# Factor order
N_budget$Pool <- factor(
  N_budget$Pool,
  levels = c("N loss due to denitrification", "N loss due to runoff","N Fixation","N Deposition","Mineralised N due to permafrost thaw")
)



N_budget$Type <- factor(
  N_budget$Type,
  levels = c("Present day reference period", "Future anomaly")
)

N_budget <- N_budget %>%
  mutate(value_g = value * 1000, value_g_std = value_std *1000)




library(dplyr)
library(tidyr)

N_budget<-read.csv("n_budget_final_no_lim.csv")

N_present <- N_budget |> 
  dplyr::filter(Type == "Present day")

N_future <- N_budget |> 
  dplyr::filter(Type == "Future")

dodge_future <- position_dodge(width = 0.8)

# Reorder Type factor so Present day comes first
N_budget <- N_budget %>%
  mutate(Type = factor(Type, levels = c("Present day", "Future")))

N_budget <- N_budget %>%
  mutate(Type = factor(Type, labels = c("Present day", "Future anomaly")))
N_present <- N_budget |> 
  dplyr::filter(Type == "Present day")

N_future <- N_budget |> 
  dplyr::filter(Type == "Future anomaly")

# Create SSP_label for consistent coloring
N_budget <- N_budget %>%
  mutate(SSP_label = ifelse(Type == "Present day", "Present", SSP))

N_budget <- N_budget %>%
  mutate(mean = mean*1000, std = std*1000)

N_budget <- N_budget %>%
  mutate(
    SSP = recode(
      SSP,
      "Present day" = "Present day", 
      "SSP126" = "SSP1-2.6",
      "SSP245" = "SSP2-4.5",
      "SSP370" = "SSP3-7.0",
      "SSP585" = "SSP5-8.5"
    )
  )

# ssp order:
N_budget$SSP <- factor(
  N_budget$SSP,
  levels = c("SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5", "Present day")
)


pool_levels <- c(
  "Ecosystem imbalance",
  "N loss due to denitrification",
  "N loss due to runoff",
  "N Deposition",
  "N Fixation", 
  "Mineralised N due to permafrost thaw"
)


ssp_cols <- c(
  "SSP5-8.5" = "#7B3294",
  "SSP3-7.0" = "#D73027",
  "SSP2-4.5" = "orange",
  "SSP1-2.6" = "blue", 
  "Present day" = "grey"
)

loss_pools <- c(
  "N loss due to denitrification",
  "N loss due to runoff"
)


N_budget <- N_budget %>%
  mutate(
    mean_plot = ifelse(Pool %in% loss_pools, -abs(mean), abs(mean))
  )

N_present <- N_budget |> 
  dplyr::filter(Type == "Present day")

N_future <- N_budget |> 
  dplyr::filter(Type == "Future anomaly")

N_present <- N_present %>%
  dplyr::mutate(Pool = factor(Pool, levels = pool_levels))

N_future <- N_future %>%
  dplyr::mutate(Pool = factor(Pool, levels = pool_levels))

ggplot() +
  
  ## ---- Present-day bars (no dodge) ----
geom_col(
  data = N_present,
  aes(x = mean_plot, y = Pool, fill = SSP_label),
  width = 0.35
) +
  
  ## ---- Error bars for present day ----
geom_errorbar(
  data = N_present %>% dplyr::filter(!is.na(std)),
  aes(
    x = mean_plot,
    y = Pool,
    xmin = mean_plot - std,
    xmax = mean_plot + std
  ),
  width = 0.2,
  linewidth = 0.6,
  color = "black"
) +
  
  ## ---- Future bars ----
geom_col(
  data = N_future,
  aes(x = mean_plot, y = Pool, fill = SSP),
  width = 0.75,
  position = dodge_future
) +
  
  ## ---- Error bars for future ----
geom_errorbar(
  data = N_future %>% dplyr::filter(!is.na(std)),
  aes(
    x = mean_plot,
    y = Pool,
    xmin = mean_plot - std,
    xmax = mean_plot + std,
    group = SSP
  ),
  position = dodge_future,
  width = 0.2,
  linewidth = 0.6
) +
  
  ## ---- Present labels ----
geom_text(
  data = N_present,
  aes(
    x = mean_plot,
    y = Pool,
    label = sprintf("%.3f", mean_plot),
    hjust = ifelse(mean_plot >= 0, -0.2, 1.2),
    vjust = 1.6,        # <-- moves text below the bar
    group = SSP
  ),
  position = position_dodge2(width = 0.8),
  size = 3
) +
  
  ## ---- Future labels (MATCH dodge exactly) ----
geom_text(
  data = N_future,
  aes(
    x = mean_plot,
    y = Pool,
    label = sprintf("%.3f", mean_plot),
    hjust = ifelse(mean_plot >= 0, -0.2, 1.2),
    vjust = 1.6,
    group = SSP
  ),
  position = position_dodge2(width = 0.8),
  size = 3
) +
  
  ## ---- Asterisks for selected pools (Future anomaly ONLY) ----
geom_text(
  data = N_future %>%
    dplyr::filter(
      Pool %in% c(
        "N loss due to runoff",
        "N loss due to denitrification"
      )
    ),
  aes(
    x = -0.07,
    y = Pool,
    label = " * "
  ),
  inherit.aes = FALSE,
  hjust = -0.6,
  size = 6,
  fontface = "bold"
) +
  
  geom_vline(xintercept = 0, linewidth = 1) +
  
  facet_wrap(~ Type, nrow = 2) +
  
  # Define colors for scenarios
  scale_fill_manual(values = ssp_cols, name="Scenario") +
  labs(
    x = expression("N [g m"^{-2}*" yr"^{-1}*"]"),
    y = "",
    title = "Nitrogen budget: Present baseline and future anomalies"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )+ labs(
    caption = "* future losses only account for linearly-derived changes"
  )
####

