
###############################################################################################################################
# calculate potential n2o emissions for bioavailable mineralised N
# either mean or std
library(terra)
arctic<-rast("1.9%_min/arctic_min_N_pool_585_no_lim_std.tif")

LC <- rast("LC_remapnn_corr.nc")
common_extent <- ext(-179.95, 179.95, 30, 90)
ext(LC) <- common_extent
LC <- resample(LC, arctic, method = "near")
# Create masks
taiga_mask    <- LC %in% c(1,2,3,4,5,8,9)
tundra_mask   <- LC %in% c(6,7,10)
wetlands_mask <- LC == 11
barren_mask   <- LC %in% c(15,16)

# mask() keeps values where mask is TRUE/1, sets others to NA
taiga_min <- mask(arctic, taiga_mask, maskvalue = 0)
tundra_min<- mask(arctic, tundra_mask, maskvalue = 0)
wetlands_min<- mask(arctic, wetlands_mask, maskvalue = 0)
barren_min<- mask(arctic, barren_mask, maskvalue = 0)


# temperature (K -> °C)
temp_data <- rast("mean_ST_585_clean_final.nc")      # adapt filename to different SSPs
ext(temp_data) <- ext(-179.95, 179.95, 30, 90)
crs(temp_data) <- crs(arctic)

temp_C <- (temp_data - 273.15)

years <- 1850:2099
ref_idx <- which(years >= 2000 & years <= 2020)


# fixed land footprint (use your mineralised raster or LC mask)
land_mask <- !is.na(arctic[[1]])

# Parameters

Ea_denitrification <- 47000   # J/mol
Ed_denitrification <- 200000  # J/mol
k_base_nitrification <- 0.1235 # from QUINCY output
k_base_denitrification<- 0.091
t_opt_C <- 28
# Define peaked_arrhenius function with your specific formulation
peaked_arrhenius <- function(temp_C, Ea, Ed, t_opt_C = 28) {
  # Gas constant
  R_gas <- 8.314462618
  
  # Convert to Kelvin
  temp_K <- temp_C + 273.15
  t_opt_K <- t_opt_C + 273.15
  
  # Helper terms
  hlp1 <- temp_K - t_opt_K
  hlp2 <- temp_K * t_opt_K * R_gas
  
  # Core Arrhenius terms
  numerator <- Ed * exp(Ea * hlp1 / hlp2)
  denominator <- Ed - Ea * (1 - exp(Ed * hlp1 / hlp2))
  
  # Calculate result
  result <- numerator / denominator
  
  
  return(result)
}

# Arrhenius response (k_T) for every year, every pixel
k_T_year <- peaked_arrhenius(temp_C, Ea_denitrification, Ed_denitrification, t_opt_C=28)
k_T_year <- mask(k_T_year, land_mask, maskvalues = FALSE)

# per-pixel baseline (mean response during 2000–2020)
k_T_ref <- mean(k_T_year[[ref_idx]], na.rm = TRUE)

# normalized multiplier (~1 during baseline)
k_factor_stack <- k_T_year / k_T_ref

# make NA temps on land become 1, keep ocean NA
k_factor_stack <- ifel(land_mask, ifel(is.na(k_factor_stack), 1, k_factor_stack), NA)

#cell_area <- cellSize(k_factor_stack, unit = "m")

#arctic_mean <- function(r, cell_area) {
 # global(r, fun = "mean", weights = cell_area, na.rm = TRUE)[,1]
#}

#kt_n2o <- arctic_mean(k_factor_stack, cell_area)


#now convert thawed N kg --> mg 
thawed_tundra_bioavailable_mg<-tundra_min* 1e6
thawed_taiga_bioavailable_mg<-taiga_min* 1e6
thawed_wetlands_bioavailable_mg<-wetlands_min* 1e6
thawed_barren_bioavailable_mg<-barren_min* 1e6
rm(barren_min, wetlands_min, taiga_min, tundra_min)
# Compute potential N₂O emissions
# Voigt factor: mg N₂O per m depth thaw per m² fläche for dry vegetated
dry_vegetated_rate_1m_mean<-0.0001113896
dry_vegetated_rate_1m_std<-4.785834e-05
wet_vegetated_rate_1m_mean<-0.0001314326
wet_vegetated_rate_1m_std<-2.02204e-05
dry_bare_rate_1m_mean<-0.00369822
dry_bare_rate_1m_std<-0.002184002

emission_factor_dry_vegetated_mean <- dry_vegetated_rate_1m_mean
emission_factor_dry_vegetated_std <- dry_vegetated_rate_1m_std
emission_factor_wet_vegetated_mean<-wet_vegetated_rate_1m_mean
emission_factor_wet_vegetated_std<-wet_vegetated_rate_1m_std
emission_factor_dry_bare_mean<-dry_bare_rate_1m_mean
emission_factor_dry_bare_std<-dry_bare_rate_1m_std

# Voigt factor: mg N₂O per m depth thaw per m² fläche for dry vegetated per year (100 season days. Voigt et al 2020)
emission_factor_year_dry_vegetated <- emission_factor_dry_vegetated_mean*100
emission_factor_year_wet_vegetated<-emission_factor_wet_vegetated_mean*100
emission_factor_year_dry_bare<-emission_factor_dry_bare_mean*100

# in N2O-N: * 0.636
n2o_n<-emission_factor_year_dry_bare*0.636

n2o_emissions_per_yr_dry_vegetated_tundra <- thawed_tundra_bioavailable_mg * emission_factor_year_dry_vegetated * k_factor_stack # mg N₂O / m²
n2o_emissions_per_yr_dry_vegetated_taiga <- thawed_taiga_bioavailable_mg * emission_factor_year_dry_vegetated * k_factor_stack # mg N₂O / m²
n2o_emissions_per_yr_wet_vegetated <- thawed_wetlands_bioavailable_mg * emission_factor_year_wet_vegetated  * k_factor_stack# mg N₂O / m²
n2o_emissions_per_yr_dry_bare <- thawed_barren_bioavailable_mg * emission_factor_year_dry_bare  * k_factor_stack# mg N₂O / m²
rm(thawed_tundra_bioavailable_mg, thawed_taiga_bioavailable_mg, thawed_wetlands_bioavailable_mg, thawed_barren_bioavailable_mg)

#plot(n2o_emissions_per_yr_dry_bare[[250]])

n2o_emissions_per_yr_dry_vegetated<-merge(n2o_emissions_per_yr_dry_vegetated_tundra,n2o_emissions_per_yr_dry_vegetated_taiga)
rm(n2o_emissions_per_yr_dry_vegetated_tundra, n2o_emissions_per_yr_dry_vegetated_taiga)
#plot(n2o_emissions_per_yr_dry_vegetated[[250]])

### weighted calculation
# Step 1: Multiply raster by cell area to get mg emissions per cell
cell_area_m2_dry_veg <- cellSize(n2o_emissions_per_yr_dry_vegetated,mask = TRUE, unit = "m")
cell_area_m2_wet_veg <- cellSize(n2o_emissions_per_yr_wet_vegetated,mask = TRUE, unit = "m")
cell_area_m2_dry_bare <- cellSize(n2o_emissions_per_yr_dry_bare,mask = TRUE, unit = "m")

n2o_emissions_per_yr_dry_vegetated <- n2o_emissions_per_yr_dry_vegetated * cell_area_m2_dry_veg
n2o_emissions_per_yr_wet_vegetated <- n2o_emissions_per_yr_wet_vegetated * cell_area_m2_wet_veg
n2o_emissions_per_yr_dry_bare <- n2o_emissions_per_yr_dry_bare * cell_area_m2_dry_bare

# --- Per-cell totals (mg N2O / cell / yr) ---
cell_area_m2_dry_veg  <- cellSize(n2o_emissions_per_yr_dry_vegetated, mask=TRUE, unit="m")
cell_area_m2_wet_veg  <- cellSize(n2o_emissions_per_yr_wet_vegetated, mask=TRUE, unit="m")
cell_area_m2_dry_bare <- cellSize(n2o_emissions_per_yr_dry_bare,       mask=TRUE, unit="m")

n2o_cell_dry_veg  <- n2o_emissions_per_yr_dry_vegetated * cell_area_m2_dry_veg
n2o_cell_wet_veg  <- n2o_emissions_per_yr_wet_vegetated * cell_area_m2_wet_veg
n2o_cell_dry_bare <- n2o_emissions_per_yr_dry_bare       * cell_area_m2_dry_bare

# Total Arctic N2O emissions per grid cell per year (mg/cell/yr)
n2o_total <- merge(n2o_emissions_per_yr_dry_vegetated, n2o_emissions_per_yr_wet_vegetated, n2o_emissions_per_yr_dry_bare)
plot(n2o_total)
# Optional: if you want "annual incl winter" like your *2 step:
n2o_cell_total_annual <- n2o_total * 2

writeRaster(n2o_cell_total_annual, "n2o/arctic_total_n2o_585_std.tif")
#writeRaster(k_factor_stack, "n2o/k_factor_stack_585_n2o.tif")

# check k_t value

# Step 2: Sum all emissions (mg)
total_emissions_dry_vegetated <- global(n2o_emissions_per_yr_dry_vegetated, fun = "sum", na.rm = TRUE)
total_emissions_wet_vegetated <- global(n2o_emissions_per_yr_wet_vegetated, fun = "sum", na.rm = TRUE)

total_emissions_dry_bare <- global(n2o_emissions_per_yr_dry_bare, fun = "sum", na.rm = TRUE)

rm(n2o_emissions_per_yr_dry_bare, n2o_emissions_per_yr_wet_vegetated,n2o_emissions_per_yr_dry_vegetated )
# Step 3: Sum and double for annual total (winter included)
total_area_emissions_mg_annual <- (total_emissions_dry_vegetated + total_emissions_wet_vegetated + total_emissions_dry_bare) * 2
total_area_emissions_mg<- (total_emissions_dry_vegetated + total_emissions_wet_vegetated + total_emissions_dry_bare)
# Step 4: Convert to Tg
total_area_emissions_Tg <- total_area_emissions_mg / 1e15

#convert to Tg N equivalent
# conversion factor: 
# M (N2O) = 44.013 g/mol
# M (N) = 14.007 g / mol
# Atomic weight of N (2 N atoms in N₂O) = 2 × 14.007 = 28.014 g/mol
# So, the fraction of N in N₂O is: 28.014 / 44.013 = 0.636

total_area_emissions_Tg_Neq <- total_area_emissions_mg * 0.636 / 1e15
total_area_emissions_Tg_Neq$years<-(1850:2099)
baseline_mean <- mean(
  total_area_emissions_Tg_Neq$sum[total_area_emissions_Tg_Neq$years %in% 2000:2020],
  na.rm = TRUE
)

# Subtract the baseline mean from all years
total_area_emissions_Tg_Neq$anomaly <- total_area_emissions_Tg_Neq$sum - baseline_mean

write.csv(total_area_emissions_Tg_Neq, "n2o-n_emissions_Tg_std_585.csv")

