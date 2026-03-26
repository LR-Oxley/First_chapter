# calculating denitrifiaction rate

#1.die Mean Denitrification Rate direkt aus Guo et al. nehmen:
denitrification_data<-read.csv("Guo_mengje_denitrification.csv")
mean_rate_nmol_n_g_h<-mean(denitrification_data$nmol_N_g_h)
sd_rate<-sd(denitrification_data$nmol_N_g_h)

# mean denitrification rate = 3.05 nmol N / g soil * h 
#  die Bodenmasseneinheit umrechnen und wegbekommen 
mean_bulk_density_g_cm3<-0.1 #estimated for the high arctic tundra soils in Guo study

mean_rate_nmol_n_cm3_h<-mean_rate_nmol_n_g_h*mean_bulk_density_g_cm3
sd_rate_nmol_n_cm3_h<-sd_rate*mean_bulk_density_g_cm3
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

####
library(terra)
mineralised_N<-rast("1.9%_min/arctic_min_N_pool_585_no_lim_mean.tif")
present_day_subset <- mineralised_N[[150:170]]

# Calculate cell areas (in square meters)
cell_areas <- cellSize(present_day_subset, unit = "m")
cell_areas <- mask(cell_areas, present_day_subset)
plot(cell_areas)
# area-weighted mean (equivalent to weighted mean over land)
num <- global(present_day_subset * cell_areas, "sum", na.rm = TRUE)[[1]]
den <- global(cell_areas, "sum", na.rm = TRUE)[[1]]

result<-num / den
result<-mean(result)

n_pool<-result
# mean denitrification per year (Guo et al., 2025)
# mean_rate_kg_N_m2_yr: 0.008 (+/- ) kg N / m2 * yr 
k_base_denitrification <- mean_rate_kg_N_m2_yr / n_pool
print(k_base_denitrification)



#### Run off estimation: 

library(terra)
stock_raster<-rast("1.9%_min/arctic_total_thawed_after_mineralisation_585_no_lim_mean.tif")
runoff_Tg_yr <- 1.04  # Zhao et al. 2025 (make sure this is for your same domain)
runoff_kg_yr <- runoff_Tg_yr * 1e9

# cell area for every grid cell
area_m2 <- cellSize(stock_raster, unit = "m")

# mask ocean / non-permafrost cells using your stock raster
area_m2 <- mask(area_m2, stock_raster)
plot(area_m2)
stock_kg <- stock_raster                         # kg N m-2
stock_total_kg <- global(stock_kg[[150]] * area_m2[[150]], "sum", na.rm=TRUE)[[1]]

f_runoff <- runoff_kg_yr / stock_total_kg        # yr^-1
f_runoff

# f_runoff = 0.000051 / yr


# the below code is incorrect, as Zhao et al and palmtag don't cover same areal extent
# (Zhao whole pan-Arctic, incl. non-permafrost areas)
runoff_Tg_yr <- 1.04 # from Zhao et al 2025

stock_Pg <- 55 # palmtag
stock_Pg_sd <- 15

f_runoff <- runoff_Tg_yr / (stock_Pg * 1000)  # yr^-1

f_min <- runoff_Tg_yr / ((stock_Pg + stock_Pg_sd) * 1000)  # using 70 Pg
f_max <- runoff_Tg_yr / ((stock_Pg - stock_Pg_sd) * 1000)  # using 40 Pg

f_runoff; f_min; f_max
