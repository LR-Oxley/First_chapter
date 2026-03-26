# meta-analysis of tundra N- fertilisation experiments

# meta analysis for absolute response ratio

library(dplyr)
library(metafor)

tundra<-read.csv("n_fertilisation_experiments/AG_tundra_new.csv")

# absolute effect: mean_treatment - mean_control
df_meta <- tundra %>%
  mutate(
    # Convert SE to SD
    sd_control   = agb_control_SE * sqrt(rep_control),
    sd_treatment = agb_treatment_SE * sqrt(rep_treatment),
    
    # Absolute effect (mean difference)
    absolute_effect = agb_treatment - agb_control,
    
    # Variance of the absolute effect
    var_absolute_effect =
      (sd_treatment^2 / rep_treatment) +
      (sd_control^2 / rep_control),
    
    se_absolute_effect = sqrt(var_absolute_effect)
  )

summary(df_meta$absolute_effect)
summary(df_meta$se_absolute_effect)
any(df_meta$var_absolute_effect <= 0)

res_abs <- rma(
  yi  = absolute_effect,
  vi  = var_absolute_effect,
  data = df_meta,
  method = "REML"
)

summary(res_abs)

forest(
  res_abs,
  slab = paste(df_meta$citation, df_meta$Biome, sep = " – "),
  xlab = expression("Absolute effect (Treatment − Control)")
)




### absolute sensitivity to N 

# effect size
df_meta <- df_meta %>%
  mutate(
    AS = absolute_effect / gN_total,
    var_AS = var_absolute_effect / gN_total^2
  )

res_AS <- rma(
  yi = AS,
  vi = var_AS,
  #mods = ~ fertiliser,
  data = df_meta,
  method = "REML"
)

summary(res_AS)

forest(
  res_AS,
  slab = paste(df_meta$citation, sep = ". – "),
  xlab = expression("g C per g N")
)

summary(res_AS)


# mutliple factors: 
res_mod <- rma(
  yi = AS,
  vi = var_AS,
  mods = ~ gN_yr,
  data = df_meta,
  method = "REML"
)

summary(res_mod)



funnel(res_abs)
regtest(res_abs)




## relative ratio:
tundra<-read.csv("n_fertilisation_experiments/AG_tundra_new.csv")

df_meta <- tundra %>%
  mutate(
    # Convert SE to SD
    sd_control   = agb_control_SE * sqrt(rep_control),
    sd_treatment = agb_treatment_SE * sqrt(rep_treatment),
    
    # Log response ratio
    lnRR = log(agb_treatment / agb_control),
    
    # Variance of lnRR
    var_lnRR =
      (sd_treatment^2 / (rep_treatment * agb_treatment^2)) +
      (sd_control^2   / (rep_control   * agb_control^2)),
    
    se_lnRR = sqrt(var_lnRR)
  )

summary(df_meta$lnRR)
summary(df_meta$se_lnRR)

any(!is.finite(df_meta$lnRR))
any(df_meta$var_lnRR <= 0)

res_lnRR <- rma(
  yi = lnRR,
  vi = var_lnRR,
  data = df_meta,
  method = "REML"
)

summary(res_lnRR)

exp(res_lnRR$b)          # multiplicative response
(exp(res_lnRR$b) - 1)*100  # % change

forest(
  res_lnRR,
  slab = paste(df_meta$citation, df_meta$Biome, sep = " – "),
  xlab = "Log response ratio (ln[Treatment / Control])"
)

forest(
  res_lnRR,
  slab = paste(df_meta$citation, df_meta$dom_species, sep = " – "),
  transf = exp,
  refline = 1,
  xlab = "Relative response ratio (Treatment / Control)")

forest(
  res_lnRR,
  slab = paste(df_meta$citation, df_meta$Biome, sep = " – "),
  transf = function(x) (exp(x) - 1) * 100,
  refline = 0,
  xlab = "Biomass change (%)")





##### N-normalised

tundra<-read.csv("n_fertilisation_experiments/AG_tundra_new.csv")

# Prepare your data with the calculated effect sizes
tundra <- tundra %>%
  mutate(
    # Calculate SD from SE
    sd_control = agb_control_SE * sqrt(rep_control),
    sd_treatment = agb_treatment_SE * sqrt(rep_treatment),
    
    # Calculate log response ratio
    yi = log(agb_treatment / agb_control),  # effect size (lnRR)
    
    # Calculate variance of lnRR
    vi = (sd_treatment^2 / (rep_treatment * agb_treatment^2)) +
      (sd_control^2 / (rep_control * agb_control^2)),
    
    # N-normalised effect size (per g N)
    yi_per_gN = yi / gN_total,
    vi_per_gN = vi / (gN_total^2)
  )


# Run meta-analysis on N-normalised effect sizes
res <- rma(yi = yi_per_gN, vi = vi_per_gN, data = tundra, 
           method = "REML", slab = citation)

# Print results
summary(res)

# Back-transform to percent change
percent_change <- (exp(coef(res)) - 1) * 100
percent_ci_lower <- (exp(res$ci.lb) - 1) * 100
percent_ci_upper <- (exp(res$ci.ub) - 1) * 100

cat(sprintf("\nOverall relative response: %.2f%% per g N [%.2f%%, %.2f%%]",
            percent_change, percent_ci_lower, percent_ci_upper))


# Create forest plot
forest(res, 
       slab = tundra$citation,
       xlab = "Response ratio (Treatment / Control)",
       header = "Study",
       transf = exp,  # Transform to ratio scale
       refline = 1,
       main = "N-normalised Biomass Response to Fertilization")
par(cex = 1)

forest(
  res,
  slab = paste(tundra$citation, tundra$dom_species, sep = " – "),
  transf = function(x) (exp(x) - 1) * 100,
  refline = 0,
  xlab = "Biomass change (% per g N m^-2)"
)

par(
  cex = 1.3,        # overall text
  cex.axis = 1.2,   # axis numbers
  cex.lab = 1.3     # axis labels
)

forest(
  res,
  slab = paste(tundra$citation, tundra$dom_species, sep = " – "),
  transf = function(x) (exp(x) - 1) * 100,
  refline = 0,
  xlab = "Biomass change (% per g N m^-2)"
)


### belowground: 
tundra<-read.csv("n_fertilisation_experiments/BG_tundra_new.csv")

# N-normalised log response ratio: 
df_meta <- tundra %>%
  mutate(
    # Convert SE to SD
    sd_control   = bgb_control_SE * sqrt(rep_control),
    sd_treatment = bgb_treatment_SE * sqrt(rep_treatment),
    
    # Log response ratio
    lnRR = log(bgb_treatment / bgb_control),
    
    # Variance of lnRR
    var_lnRR =
      (sd_treatment^2 / (rep_treatment * bgb_treatment^2)) +
      (sd_control^2   / (rep_control   * bgb_control^2)),
    
    se_lnRR = sqrt(var_lnRR),
    
    # ----- N-normalised effect size -----
    lnRR_per_gN = lnRR / gN_total,
    
    var_lnRR_per_gN = var_lnRR / (gN_total^2),
    
    se_lnRR_per_gN = sqrt(var_lnRR_per_gN)
  )

# random-effects model using N-normalised effect size
res_Nsens <- rma(
  yi = lnRR_per_gN,
  vi = var_lnRR_per_gN,
  data = df_meta,
  method = "REML"
)

forest(res_Nsens, 
       slab = tundra$citation,
       xlab = "ln(Response Ratio) per g N/yr",
       header = "Study",
       transf = exp,  # Transform to ratio scale
       refline = 1,
       main = "N-normalised Biomass Response to Fertilization")



forest(
  res_Nsens,
  slab = paste(df_meta$citation, df_meta$dom_species, sep = " – "),
  refline = 0,
  xlab = expression("lnRR per g N added")
)

forest(
  res_Nsens,
  slab = paste(df_meta$citation, df_meta$dom_species, sep = " – "),
  transf = function(x) (exp(x) - 1) * 100,
  refline = 0,
  xlab = "% biomass increase per g N added"
)



# absolute effect
df_meta <- tundra %>%
  mutate(
    # Convert SE to SD
    sd_control   = bgb_control_SE * sqrt(rep_control),
    sd_treatment = bgb_treatment_SE * sqrt(rep_treatment),
    
    # Absolute effect (mean difference)
    absolute_effect = bgb_treatment - bgb_control,
    
    # Variance of the absolute effect
    var_absolute_effect =
      (sd_treatment^2 / rep_treatment) +
      (sd_control^2 / rep_control),
    
    se_absolute_effect = sqrt(var_absolute_effect)
  )

df_meta <- df_meta %>%
  mutate(
    AS = absolute_effect / gN_total,
    var_AS = var_absolute_effect / gN_total^2
  )

res_AS <- rma(
  yi = AS,
  vi = var_AS,
  data = df_meta,
  method = "REML"
)


forest(
  res_AS,
  slab = paste(df_meta$citation, df_meta$dom_type, sep = ". – "),
  xlab = expression("g C per g N")
)

summary(res_AS)

## relative ratio:
tundra<-read.csv("n_fertilisation_experiments/BG_tundra_new.csv")

df_meta <- tundra %>%
  mutate(
    # Convert SE to SD
    sd_control   = bgb_control_SE * sqrt(rep_control),
    sd_treatment = bgb_treatment_SE * sqrt(rep_treatment),
    
    # Log response ratio
    lnRR = log(bgb_treatment / bgb_control),
    
    # Variance of lnRR
    var_lnRR =
      (sd_treatment^2 / (rep_treatment * bgb_treatment^2)) +
      (sd_control^2   / (rep_control   * bgb_control^2)),
    
    se_lnRR = sqrt(var_lnRR)
  )

summary(df_meta$lnRR)
summary(df_meta$se_lnRR)

any(!is.finite(df_meta$lnRR))
any(df_meta$var_lnRR <= 0)

res_lnRR <- rma(
  yi = lnRR,
  vi = var_lnRR,
  data = df_meta,
  method = "REML"
)

summary(res_lnRR)

exp(res_lnRR$b)          # multiplicative response
(exp(res_lnRR$b) - 1)*100  # % change

forest(
  res_lnRR,
  slab = paste(df_meta$citation, df_meta$Biome, sep = ". – "),
  xlab = "Log response ratio (ln[Treatment / Control])"
)

forest(
  res_lnRR,
  slab = paste(df_meta$citation, df_meta$dom_species, sep = " – "),
  transf = exp,
  refline = 1,
  xlab = "Relative response ratio (Treatment / Control)")

forest(
  res_lnRR,
  slab = paste(df_meta$citation, df_meta$dom_species, sep = " – "),
  transf = function(x) (exp(x) - 1) * 100,
  refline = 0,
  xlab = "Biomass change (%)")



## for taiga: 
# absolute C-N response: 14.1 kg C per kg N; 
# relative response: 1.20;
# average N addition rate: 2 g N / m2* yr
# average length of experiment: 17.41 years
# total N added: 35 g N / m2
library(dplyr)
library(metafor)

taiga <- data.frame(
  yi = 1.2,
  vi = 0.1,
  gN_total = 35
)

taiga$yi_per_gN <- taiga$yi / taiga$gN_total
taiga$vi_per_gN <- taiga$vi / (taiga$gN_total^2)

res <- rma(yi = yi_per_gN, vi = vi_per_gN, data = taiga)
summary(res)
# Back-transform to percent change
percent_change <- (exp(coef(res)) - 1) * 100
percent_ci_lower <- (exp(res$ci.lb) - 1) * 100
percent_ci_upper <- (exp(res$ci.ub) - 1) * 100

cat(sprintf("\nOverall relative response: %.2f%% per g N [%.2f%%, %.2f%%]",
            percent_change, percent_ci_lower, percent_ci_upper))






# differentiate between fertiliser types
tundra<-read.csv("n_fertilisation_experiments/AG_tundra_long_new.csv")

# absolute effect: mean_treatment - mean_control
df_meta <- tundra %>%
  mutate(
    # Convert SE to SD
    sd_control   = agb_control_SE * sqrt(rep_control),
    sd_treatment = agb_treatment_SE * sqrt(rep_treatment),
    
    # Absolute effect (mean difference)
    absolute_effect = agb_treatment - agb_control,
    
    # Variance of the absolute effect
    var_absolute_effect =
      (sd_treatment^2 / rep_treatment) +
      (sd_control^2 / rep_control),
    
    se_absolute_effect = sqrt(var_absolute_effect)
  )

# effect size
df_meta <- df_meta %>%
  mutate(
    AS = absolute_effect / gN_total,
    var_AS = var_absolute_effect / gN_total^2
  )

# keep only the fertiliser groups you want
df_forest <- df_meta %>%
  filter(fertiliser %in% c("N", "P", "NP", "NPK")) %>%
  mutate(study_label = paste(citation, Biome, sep = " | "))

# split data
df_N  <- df_forest %>% filter(fertiliser == "N")
df_P  <- df_forest %>% filter(fertiliser == "P")
df_NP <- df_forest %>%
  filter(fertiliser %in% c("NP", "NPK"))

# fit separate models
res_N  <- rma(yi = AS, vi = var_AS, data = df_N, method = "REML")
res_P  <- rma(yi = AS, vi = var_AS, data = df_P, method = "REML")
res_NP <- rma(yi = AS, vi = var_AS, data = df_NP, method = "REML")

# plot all three in one figure
par(mfrow = c(3, 1), mar = c(4, 4, 2, 2))
par(mfrow = c(1, 1), mar = c(4, 4, 2, 2))
forest(
  res_N,
  slab = df_N$citation,
  xlab = "Absolute sensitivity (g biomass increase per g N added)",
  main = "Fertiliser: N"
)
addpoly(res_N, row = -1, mlab = "RE Model")

forest(
  res_P,
  slab = df_P$study_label,
  xlab = "Absolute sensitivity (g biomass increase per g N added)",
  main = "Fertiliser: P"
)
addpoly(res_P, row = -1, mlab = "RE Model")

forest(
  res_NP,
  slab = df_NP$study_label,
  xlab = "Absolute sensitivity (g biomass increase per g N added)",
  main = "Fertiliser: NP"
)
addpoly(res_NP, row = -1, mlab = "RE Model")


# n-normalized relative ratio

# Prepare your data with the calculated effect sizes
tundra <- tundra %>%
  mutate(
    # Calculate SD from SE
    sd_control = agb_control_SE * sqrt(rep_control),
    sd_treatment = agb_treatment_SE * sqrt(rep_treatment),
    
    # Calculate log response ratio
    yi = log(agb_treatment / agb_control),  # effect size (lnRR)
    
    # Calculate variance of lnRR
    vi = (sd_treatment^2 / (rep_treatment * agb_treatment^2)) +
      (sd_control^2 / (rep_control * agb_control^2)),
    
    # N-normalised effect size (per g N)
    yi_per_gN = yi / gN_total,
    vi_per_gN = vi / (gN_total^2)
  )


# Run meta-analysis on N-normalised effect sizes
res <- rma(yi = yi_per_gN, vi = vi_per_gN, data = tundra, 
           method = "REML", slab = citation)

# Print results
summary(res)

# Back-transform to percent change
percent_change <- (exp(coef(res)) - 1) * 100
percent_ci_lower <- (exp(res$ci.lb) - 1) * 100
percent_ci_upper <- (exp(res$ci.ub) - 1) * 100

cat(sprintf("\nOverall relative response: %.2f%% per g N [%.2f%%, %.2f%%]",
            percent_change, percent_ci_lower, percent_ci_upper))


# Create forest plot
forest(res, 
       slab = tundra$citation,
       xlab = "Response ratio (Treatment / Control)",
       header = "Study",
       transf = exp,  # Transform to ratio scale
       refline = 1,
       main = "N-normalised Biomass Response to Fertilization")
par(cex = 1)

forest(
  res,
  slab = paste(tundra$citation, tundra$dom_species, sep = " – "),
  transf = function(x) (exp(x) - 1) * 100,
  refline = 0,
  xlab = "Biomass change (% per g N m^-2)"
)



## distinct N, P and NPK
# Keep relevant groups
tundra <- tundra %>%
  filter(fertiliser %in% c("N", "P", "NP", "NPK")) %>%
  mutate(
    fert_group = case_when(
      fertiliser == "N" ~ "N",
      fertiliser == "P" ~ "P",
      fertiliser %in% c("NP", "NPK") ~ "NP"
    ),
    study_label = paste(citation, Biome, sep = " | ")
  )

# Split
df_N  <- tundra %>% filter(fert_group == "N")
df_P  <- tundra %>% filter(fert_group == "P")
df_NP <- tundra %>% filter(fert_group == "NP")

# Models
res_N  <- rma(yi = yi_per_gN, vi = vi_per_gN, data = df_N,  method = "REML")
res_P  <- rma(yi = yi_per_gN, vi = vi_per_gN, data = df_P,  method = "REML")
res_NP <- rma(yi = yi_per_gN, vi = vi_per_gN, data = df_NP, method = "REML")

forest(
  res_N,
  slab = df_N$study_label,
  xlab = "ln(response ratio) per g N added",
  main = "Fertiliser: N"
)
addpoly(res_N, row = -1, mlab = "RE Model")

forest(
  res_P,
  slab = df_P$study_label,
  xlab = "ln(response ratio) per g N added",
  main = "Fertiliser: P"
)
addpoly(res_P, row = -1, mlab = "RE Model")

forest(
  res_NP,
  slab = df_NP$study_label,
  xlab = "ln(response ratio) per g N added",
  main = "Fertiliser: NP + NPK"
)
addpoly(res_NP, row = -1, mlab = "RE Model")




#### visualize the locations on a map

library(dplyr)
library(readr)
library(stringr)

tundra<-read.csv("n_fertilisation_experiments/AG_tundra_long_new.csv")

df_clean <- tundra %>%
  mutate(
    lat = str_replace(latitude, "_", "."),
    lon = str_replace(longitude, "_", "."),
    
    lat = as.numeric(str_remove(lat, "_[NS]$")),
    lon = as.numeric(str_remove(lon, "_[EW]$")),
    
    lat = ifelse(str_detect(latitude, "S"), -lat, lat),
    lon = ifelse(str_detect(longitude, "W"), -lon, lon)
  )

library(ggplot2)
library(maps)

world <- map_data("world")

ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),
               fill = "grey90", color = "white") +
  geom_point(data = df_clean,
             aes(x = lon, y = lat, color = Biome),
             size = 3) +
  coord_fixed(1.3, xlim = c(-180, 180), ylim = c(50, 90)) +
  theme_minimal() +
  labs(title = "Tundra Sites",
       x = "Longitude", y = "Latitude")

ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group),
               fill = "grey95", color = "white") +
  geom_point(data = df_clean,
             aes(x = lon, y = lat, color = Biome),
             size = 3) +
  coord_map("ortho", orientation = c(90, 0, 0)) +  # Arctic view
  theme_void() +
  labs(title = "Arctic Tundra Sampling Sites")
