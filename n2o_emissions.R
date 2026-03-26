terra::tmpFiles(remove = TRUE, current = TRUE, orphan = TRUE, old = TRUE)

tempdir()
### calculating N2O emissions
# calculate ratio N20 emission per total dissolved N (sum of organic and inorganic N) from permafrost thaw (data from Voigt et al 2017)

#Total dissolved N: mg N / kg DW
# dry - bare: 307 +/- 111 
# dry - vegetated: 570 +/- 199

# wet - bare: 219 +/- 40
# wet - vegetated: 314 +/- 75

# mass = density * volume
# convert total dissolved N: mg N / kg DW into g / cm'3 
# need bulk density: bare = 0.11 g / cm3 = 0.00011kg/cm3
# vegetated = 0.14 g / cm3 = 0.00014 kg/cm3

# dry bare: 307 mg N / kg DW --> how much is that in mg N / cm'3
# 307 mg N / kg * 0.00011 kg/cm'3 = 0.03377 mg N / cm'3
# dry vegetated: 570 mg N / kg DW * 0.00014 kg/ cm3 = 0.0798 mg N / cm3

# wet bare: 219 mg N / kg DW * 0.00011 kg / cm'3 = 0.02409 mg N / cm3
# wet vegetated: 314 mg N / kg DW * 0.00014 kg / cm3 = 0.04396 mg N / cm3

# to integrate this concentration of total dissolved N over 3cm deep soil slice that's 1m2
# mg N/ cm3 --> mg N/m2

# 1 m2 * 3cm ==> 10'000cm2 * 3 cm = 30'000cm3

# total dissolved N in /m2 for 3 cm deep slice
# dry bare = 0.03377 mg N / cm3 * 30'000 cm3 = 1'013.1 mg N / m2
# dry vegetated = 0.0798 mg N / cm3 * 30'000cm3 = 2'394 mg N / m2
# wet bare = 0.02409 mg N / cm3 * 30'000 cm3 = 722.7 mg N / m2 
# wet vegetated = 0.04396 mg N / cm3 * 30'000 cm3 = 1'318.8 mg N / m2

# to compare with 15 cm deep permafrost core: multiply by 5
# dry bare = 1'013.1 mg N / m2 * 5 = 5065.5 mg N / m2 
# dry vegetated = 2'394 mg N / m2 * 5= 11970 mg N / m2
# wet bare = 722.7 mg N / m2 * 5= 3613.5 mg N / m2
# wet vegetated = 1'318.8 mg N / m2 * 5 = 6594 mg N / m2

# palmtag et al: 0.5 kg N / m2 == 500'000 mg N / m2 (but this includes particulate organic N)
# for 50 cm depth increment

# to compare this with Voigt et al 15 cm deep permafrost core: multiply by 3
# dry bare = 5065.5 mg N / m2 * 3 = 15'195 mg N / m2
# dry vegetated = 11970 mg N / m2 * 3 = 35'910 mg N / m2
# wet bare = 3613.5 mg N / m2 *3 = 10'840.5 mg N / m2
# wet vegetated = 6590 mg N / m2 *3 = 19'770 mg N / m2

# as Voigt et al only included Total dissolved N, which is roughly 10% of the total N (particulate and dissolved) in the soil, the rest would be: 
# dry bare = 15195 / 0.1 = 151'950 mg N / m2
# dry vegetated = 35'910 /0.1 = 359'100 mg N / m2
# wet bare = 10'840 / 0.1 = 108'400 mg N / m2
# wet vegetated = 19'770 = 197'700 mg N / m2


# mean N2O emissions: mg N2O / m2*day for a 15cm slice of permafrost
# dry - bare: 2.81 +/- 0.6
# dry - vegetated: 0.20 +/- 0.03

# wet - bare: 0.21 +/- 0.03
# wet - vegetated: 0.13 +/- 0.02

# n20 emission rate per thawed total dissolved N for 15 cm deep permafrost slice

# dry bare = 2.81 mg N20 / m2 * day / 5065.5 mg N / m2 = 0.0005547 mg N20 / m2 *day per mg N / m2 
# dry vegetated = 0.2 mg N20 / m2 * day / 11970 mg N / m2 = 0.0000167 mg N20 / m2 *day per mg N / m2
# wet bare = 0.21 mg N20 / m2 * day / 3'613.5 mg N / m2 = 0.00005811 mg N20 / m2 *day per mg N / m2
# wet vegetated = 0.13 mg N20 / m2 / 6'594 mg N / m2 = 0.00001971 mg N20 / m2 *day per mg N / m2

# scale n20 emission rate per thawed total dissolved N for 100 cm deep permafrost slice
# dry bare = 0.0005547 mg N20 / m2 *day per mg N / m2 * 6.6 =  0.00366124 mg N20 / m2 * day per mg N / m2
# dry vegetated = 0.0000167 mg N20 / m2 *day per mg N / m2 * 6.6 = 0.00011022 
# wet bare = 0.00005811 mg N20 / m2 *day per mg N / m2 * 6.6 = 0.00038353 mg N20 / m2 * day per mg N / m2
# wet vegetated = 0.00001971 mg N20 / m2 *day per mg N / m2 * 6.6 = 0.00013009 mg N20 / m2 * day per mg N / m2


########
### calculating N2O emissions
# calculate ratio N20 emission per total dissolved N from permafrost thaw (data from Voigt et al 2017)

#Total dissolved N: mg N / kg DW
# dry - bare: 307 +/- 111 
# dry - vegetated: 570 +/- 199

# wet - bare: 219 +/- 40
# wet - vegetated: 314 +/- 75

# mass = density * volume
# convert total dissolved N: mg N / kg DW into g / cm'3 
# need bulk density: bare = 0.11 g / cm3 = 0.00011kg/cm3
# vegetated = 0.14 g / cm3 = 0.00014 kg/cm3

# dry bare: 307 mg N / kg DW --> how much is that in mg N / cm'3
# volume = mass / bulk density
dry_bare_mean<-307 *  0.00011 
dry_bare_std<-111 *  0.00011 
dry_vegetated_mean<- 570 * 0.00014 
dry_vegetated_std<- 199 * 0.00014 
wet_bare_mean<-219 * 0.00011 
wet_bare_std<-40 * 0.00011 
wet_vegetated_mean<-314 * 0.00014 
wet_vegetated_std<-75 * 0.00014 


# to integrate this concentration of total dissolved N over 3cm deep soil slice that's 1m2
# mg N/ cm3 --> mg N/m2

# 1 m2 * 3cm ==> 10'000cm2 * 3 cm = 30'000cm3

# total dissolved N in /m2 for 3 cm deep slice
dry_bare_m2_mean <- dry_bare_mean * 30000
dry_bare_m2_std <- dry_bare_std * 30000
dry_vegetated_m2_mean <-dry_vegetated_mean *30000
dry_vegetated_m2_std <-dry_vegetated_std *30000
wet_bare_m2_mean <- wet_bare_mean* 30000
wet_bare_m2_std <- wet_bare_std* 30000
wet_vegetated_m2_mean <- wet_vegetated_mean * 30000 
wet_vegetated_m2_std <- wet_vegetated_std * 30000 


# to compare with 15 cm deep permafrost core: multiply by 5
diss_N_dry_bare_15cm_mean<- dry_bare_m2_mean * 5 
diss_N_dry_bare_15cm_std<- dry_bare_m2_std * 5 
diss_N_dry_vegetated_15cm_mean<- dry_vegetated_m2_mean* 5
diss_N_dry_vegetated_15cm_std<- dry_vegetated_m2_std* 5
diss_N_wet_bare_15cm_mean<- wet_bare_m2_mean * 5
diss_N_wet_bare_15cm_std<- wet_bare_m2_std * 5
diss_N_wet_vegetated_15cm_mean <-wet_vegetated_m2_mean * 5 
diss_N_wet_vegetated_15cm_std <-wet_vegetated_m2_std * 5 

# palmtag et al: 0.5 kg N / m2 == 500'000 mg N / m2 (but this includes particulate organic N)
# for 50 cm depth increment

# to compare this with Voigt et al 15 cm deep permafrost core: multiply by 3 to scale it to 45 cm which is comparable to the 50 cm of palmtag
dry_bare_voigt_mean <- diss_N_dry_bare_15cm_mean * 3 
dry_bare_voigt_std <- diss_N_dry_bare_15cm_std * 3 
dry_vegetated_voigt_mean <- diss_N_dry_vegetated_15cm_mean * 3 
dry_vegetated_voigt_std <- diss_N_dry_vegetated_15cm_std * 3 
wet_bare_voigt_mean <- diss_N_wet_bare_15cm_mean*3 
wet_bare_voigt_std <- diss_N_wet_bare_15cm_std*3 
wet_vegetated_voigt_mean <- diss_N_wet_vegetated_15cm_mean *3 
wet_vegetated_voigt_std <- diss_N_wet_vegetated_15cm_std *3 

# as Voigt et al only included Total dissolved N, which is roughly 10% of the total N (particulate and dissolved) in the soil, the rest would be: 
dry_bare_comparison_mean<- dry_bare_voigt_mean / 0.1 
dry_bare_comparison_std<- dry_bare_voigt_std / 0.1 
dry_vegetated_comparison_mean <- dry_vegetated_voigt_mean/0.1 
dry_vegetated_comparison_std <- dry_vegetated_voigt_std/0.1 
wet_bare_comparison_mean <- wet_bare_voigt_mean / 0.1 
wet_bare_comparison_std <- wet_bare_voigt_std / 0.1 
wet_vegetated_mean <- wet_vegetated_voigt_mean / 0.1 
wet_vegetated_std <- wet_vegetated_voigt_std / 0.1 
comparison_with_palmtag<-dry_bare_comparison_mean+dry_vegetated_comparison_mean+wet_bare_comparison_mean+wet_vegetated_mean
print(comparison_with_palmtag)

# 817'290 mg N /m2 which is a bit more than the palmtag estimate


#### mean N2O emissions: mg N2O / m2*day for a 15cm slice of permafrost
dry_bare_n2o_mean<- 2.81 
dry_bare_n2o_std<- 0.6 
dry_vegetated_n2o_mean<-0.20 
dry_vegetated_n2o_std<-0.03 

wet_bare_n2o_mean<- 0.21 
wet_bare_n2o_std<- 0.03 
wet_vegetated_n2o_mean<- 0.13 
wet_vegetated_n2o_std<- 0.02 


# compare with Voigt, Carolina; van Delden, Lona; Marushchak, Maija E; Biasi, Christina; Abbott, Benjamin W; Elberling, Bo;
# Siciliano, Steven D; Sonnentag, Oliver; Stewart, Katherine J; Yang, Yuanhe; Martikainen, Pertti J (2020): Nitrous oxide fluxes from permafrost regions: worldwide synthesis dataset
# mean N2O emissions in mg N2O/m2*day
# dry bare: 3.88
# dry vegetated: 0.87
# wet vegetated: 0.602

# divide n20 emission rate from a 15cm core slice by thawed total dissolved N in a 15 cm deep permafrost slice == mg N20 / day per mg N 

dry_bare_rate_15_mean <-dry_bare_n2o_mean / diss_N_dry_bare_15cm_mean
dry_bare_rate_15_std <-dry_bare_n2o_std / diss_N_dry_bare_15cm_std
dry_vegetated_rate_15_mean <-dry_vegetated_n2o_mean / diss_N_dry_vegetated_15cm_mean
dry_vegetated_rate_15_std <-dry_vegetated_n2o_std / diss_N_dry_vegetated_15cm_std
wet_bare_rate_15_mean<-wet_bare_n2o_mean/ diss_N_wet_bare_15cm_mean
wet_bare_rate_15_std<-wet_bare_n2o_std/ diss_N_wet_bare_15cm_std
wet_vegetated_rate_15_mean<- wet_vegetated_n2o_mean / diss_N_wet_vegetated_15cm_mean 
wet_vegetated_rate_15_std<- wet_vegetated_n2o_std / diss_N_wet_vegetated_15cm_std

# n20 emission rate per thawed total dissolved N for 1 cm deep permafrost slice == mg N20 / day per mg N 

dry_bare_rate_1cm_mean <-dry_bare_rate_15_mean / 15 
dry_bare_rate_1cm_std <-dry_bare_rate_15_std / 15 
dry_vegetated_rate_1cm_mean <-dry_vegetated_rate_15_mean / 15
dry_vegetated_rate_1cm_std <-dry_vegetated_rate_15_std / 15
wet_bare_rate_1cm_mean<-wet_bare_rate_15_mean/ 15
wet_bare_rate_1cm_std<-wet_bare_rate_15_std/ 15
wet_vegetated_rate_1cm_mean<- wet_vegetated_rate_15_mean / 15 
wet_vegetated_rate_1cm_std<- wet_vegetated_rate_15_std / 15 

# scale n20 emission rate per day per thawed total dissolved N for 100 cm deep permafrost slice= *100 (as ALD is in m)
dry_bare_rate_1m_mean<- dry_bare_rate_1cm_mean* 100 
dry_bare_rate_1m_std<- dry_bare_rate_1cm_std* 100 
dry_vegetated_rate_1m_mean<- dry_vegetated_rate_1cm_mean* 100 
dry_vegetated_rate_1m_std<- dry_vegetated_rate_1cm_std* 100
wet_bare_rate_1m_mean<- wet_bare_rate_1cm_mean * 100 
wet_bare_rate_1m_std<- wet_bare_rate_1cm_std * 100
wet_vegetated_rate_1m_mean<- wet_vegetated_rate_1cm_mean * 100 
wet_vegetated_rate_1m_std<- wet_vegetated_rate_1cm_std * 100

###############################################################################################################################
# calculate potential n2o emissions for bioavailable mineralised N
# either mean or std

arctic<-rast("mineralised_N_pool_126_std.tif")

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
plot(taiga_min[[2]])
tundra_min<- mask(arctic, tundra_mask, maskvalue = 0)
wetlands_min<- mask(arctic, wetlands_mask, maskvalue = 0)
barren_min<- mask(arctic, barren_mask, maskvalue = 0)
plot(barren_min[[2]])
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

n2o_emissions_per_yr_dry_vegetated_tundra <- thawed_tundra_bioavailable_mg * emission_factor_year_dry_vegetated  # mg N₂O / m²
n2o_emissions_per_yr_dry_vegetated_taiga <- thawed_taiga_bioavailable_mg * emission_factor_year_dry_vegetated  # mg N₂O / m²
n2o_emissions_per_yr_wet_vegetated <- thawed_wetlands_bioavailable_mg * emission_factor_year_wet_vegetated  # mg N₂O / m²
n2o_emissions_per_yr_dry_bare <- thawed_barren_bioavailable_mg * emission_factor_year_dry_bare  # mg N₂O / m²
rm(thawed_tundra_bioavailable_mg, thawed_taiga_bioavailable_mg, thawed_wetlands_bioavailable_mg, thawed_barren_bioavailable_mg)

plot(n2o_emissions_per_yr_dry_vegetated_tundra[[250]])
#plot(n2o_emissions_per_yr_dry_bare[[250]])

n2o_emissions_per_yr_dry_vegetated<-merge(n2o_emissions_per_yr_dry_vegetated_tundra,n2o_emissions_per_yr_dry_vegetated_taiga)
rm(n2o_emissions_per_yr_dry_vegetated_tundra, n2o_emissions_per_yr_dry_vegetated_taiga)
#plot(n2o_emissions_per_yr_dry_vegetated[[250]])
### weighted calculation
# Step 1: Multiply raster by cell area to get mg emissions per cell
cell_area_m2_dry_veg <- cellSize(n2o_emissions_per_yr_dry_vegetated,mask = TRUE, unit = "m")
plot(cell_area_m2_dry_bare)
cell_area_m2_wet_veg <- cellSize(n2o_emissions_per_yr_wet_vegetated,mask = TRUE, unit = "m")
cell_area_m2_dry_bare <- cellSize(n2o_emissions_per_yr_dry_bare,mask = TRUE, unit = "m")

n2o_emissions_per_yr_dry_vegetated <- n2o_emissions_per_yr_dry_vegetated * cell_area_m2_dry_veg
n2o_emissions_per_yr_wet_vegetated <- n2o_emissions_per_yr_wet_vegetated * cell_area_m2_wet_veg
n2o_emissions_per_yr_dry_bare <- n2o_emissions_per_yr_dry_bare * cell_area_m2_dry_bare

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

write.csv(total_area_emissions_Tg_Neq, "n2o-n_emissions_Tg_std_126.csv")

plot(total_area_emissions_Tg_Neq$anomaly)
# estimates (Repo et al 2009): 0.1 Tg per year = 100'000'000 kg

total_area_emissions_Tg_Neq<-read.csv("n2o-n_emissions_total.csv")

# annual: double the emissions (as wintertime contributes to 30-60% of annual emissions, Voigt et al 2020)
total_area_emissions_Tg_Neq<-total_area_emissions_Tg_Neq[,2:9]*2
total_area_emissions_Tg_Neq$years<-c(1850:2099)
present_day<-total_area_emissions_Tg_Neq[total_area_emissions_Tg_Neq$years %in% 2000:2020, ]
mean_present_day <- colMeans(present_day, na.rm = TRUE)

future<-total_area_emissions_Tg_Neq[total_area_emissions_Tg_Neq$years %in% 2080:2099, ]
mean_future<- colMeans(future, na.rm = TRUE)

increase<-mean_future-mean_present_day
increase<-increase/90
library(ggplot2)

df <- data.frame(
  scenario = c("SSP 1-2.6", "SSP 2-4.5", "SSP 3-7.0", "SSP 5-8.5"),
  mean = increase[c("mean_126", "mean_245", "mean_370", "mean_585")],
  std  = increase[c("std_126",  "std_245",  "std_370",  "std_585")]
)


ggplot(df, aes(x = scenario, y = mean, fill = scenario)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean - std, ymax = mean + std), width = 0.2) +
  xlab("") +
  ylab(expression(N[2]*O~"["*Tg~N[2]*O-N~yr^{-1}*"]")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, vjust=-3, size= 14))+
  ggtitle("") +
  scale_fill_manual(
    name = "SSP scenario",
    values = c(
      "SSP 1-2.6" = "blue", 
      "SSP 2-4.5" = "orange",
      "SSP 3-7.0" = "#D73027",
      "SSP 5-8.5" = "#7B3294"
    )
  )

