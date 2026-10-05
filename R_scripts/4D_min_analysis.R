
library(terra)
library(dplyr)
library(ggplot2)

# ============================================================
# Settings
# ============================================================

base_dir <- "monthly_mineralised/mean_2perc_baseline_new"
# Replace with the folder containing the NetCDF files

ssps <- c("126", "245", "370", "585")

ssp_labels <- c(
  "126" = "SSP1-2.6",
  "245" = "SSP2-4.5",
  "370" = "SSP3-7.0",
  "585" = "SSP5-8.5"
)

ssp_colours <- c(
  "SSP1-2.6" = "blue",
  "SSP2-4.5" = "orange",
  "SSP3-7.0" = "#D73027",
  "SSP5-8.5" = "#7B3294"
)

# ============================================================
# Read and summarise one SSP
# ============================================================

read_net_pool_change <- function(ssp) {
  
  nc_file <- file.path(
    base_dir,
    paste0(
      "arctic_net_organic_N_pool_change_yearly_",
      ssp,
      "_mean_1850_2099_w_temp.nc"
    )
  )
  
  if (!file.exists(nc_file)) {
    stop("File not found: ", nc_file)
  }
  
  cat("Reading:", nc_file, "\n")
  
  net_change <- rast(nc_file)
  
  # Extract years from the NetCDF time dimension
  net_time <- time(net_change)
  
  if (
    is.null(net_time) ||
    length(net_time) != nlyr(net_change)
  ) {
    years <- 1850:(1850 + nlyr(net_change) - 1)
  } else {
    years <- as.integer(format(as.Date(net_time), "%Y"))
  }
  
  # Grid-cell areas in m²
  area_rast <- cellSize(
    net_change[[1]],
    unit = "m"
  )
  
  # Spatial sum:
  # kg N m-2 × m2 = kg N
  # kg N / 1e12 = Pg N
  net_change_pg <- global(
    net_change * area_rast,
    "sum",
    na.rm = TRUE
  )[, 1] / 1e12
  
  data.frame(
    Year = years,
    SSP_raw = ssp,
    SSP = unname(ssp_labels[ssp]),
    net_organic_pool_change_pg = net_change_pg
  )
}

# ============================================================
# Combine all SSPs
# ============================================================

net_pool_change_df <- bind_rows(
  lapply(ssps, read_net_pool_change)
) %>%
  mutate(
    SSP = factor(
      SSP,
      levels = unname(ssp_labels[ssps])
    )
  ) %>%
  arrange(SSP, Year)

print(head(net_pool_change_df))
print(summary(net_pool_change_df))

# ============================================================
# Plot
# ============================================================

p_net_pool_change <- ggplot(
  net_pool_change_df,
  aes(
    x = Year,
    y = net_organic_pool_change_pg,
    colour = SSP
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_line(
    linewidth = 0.9,
    na.rm = TRUE
  ) +
  scale_colour_manual(
    values = ssp_colours,
    drop = FALSE
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Annual net change in the organic N pool",
    subtitle = "Annual organic N input minus annual mineralisation",
    x = "Year",
    y = expression("Net organic N pool change (Pg N yr"^{-1}*")"),
    colour = "SSP scenario"
  )

print(p_net_pool_change)

ggsave(
  filename = file.path(
    base_dir,
    "annual_net_organic_N_pool_change_all_SSPs.png"
  ),
  plot = p_net_pool_change,
  width = 10,
  height = 5.5,
  dpi = 300
)

write.csv(
  net_pool_change_df,
  file.path(
    base_dir,
    "annual_net_organic_N_pool_change_all_SSPs.csv"
  ),
  row.names = FALSE
)

########## compare org pool remaining and cumulative newly thawed org 

library(terra)
library(ggplot2)

terraOptions(memfrac = 0.4)

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

ssp <- "126"   # change to 126, 245, 370 or 585
getwd()
base_dir <- "monthly_mineralised/mean_2perc_baseline"

yearly_dir <- file.path(
  base_dir,
  paste0("_w_temp_w_sm_yearly_nc_", ssp)
)



# ------------------------------------------------------------
# Load yearly rasters
# ------------------------------------------------------------

new_thawed_organic_N <- rast(
  file.path(
    base_dir,
    paste0(
      "arctic_new_thawed_organic_N_yearly_",
      ssp,
      "_mean_1850_2099_w_temp.nc"
    )
  )
)

organic_pool_remaining <- rast(
  file.path(
    base_dir,
    paste0(
      "arctic_organic_N_pool_remaining_yearly_",
      ssp,
      "_mean_1850_2099_w_temp.nc"
    )
  )
)

# ------------------------------------------------------------
# Convert to Pg N
# ------------------------------------------------------------

area_rast <- cellSize(new_thawed_organic_N[[1]], unit = "m")

annual_input_pg <-
  global(new_thawed_organic_N * area_rast,
         "sum",
         na.rm = TRUE)[,1] / 1e12


cum_input_pg <- cumsum(annual_input_pg)

pool_pg <-
  global(organic_pool_remaining * area_rast,
         "sum",
         na.rm = TRUE)[,1] / 1e12

# ------------------------------------------------------------
# Data frame
# ------------------------------------------------------------

compare_df <- data.frame(
  Year = as.integer(format(time(new_thawed_organic_N), "%Y")),
  cumulative_input_pg = cum_input_pg,
  pool_remaining_pg = pool_pg
)

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

ggplot(compare_df, aes(x = Year)) +
  geom_line(
    aes(y = cumulative_input_pg,
        colour = "Cumulative newly thawed organic N"),
    linewidth = 1.2
  ) +
  geom_line(
    aes(y = pool_remaining_pg,
        colour = "Organic pool remaining"),
    linewidth = 1.2
  ) +
  scale_colour_manual(
    values = c(
      "Cumulative newly thawed organic N" = "black",
      "Organic pool remaining" = "#D73027"
    ),
    name = NULL
  ) +
  labs(
    x = "Year",
    y = "Pg N",
    title = paste("Organic N mass balance (SSP", ssp, ")", sep = "")
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  )



library(dplyr)
library(ggplot2)
library(readr)
library(zoo)

# -----------------------------
# Settings
# -----------------------------

base_dir <- "monthly_mineralised/mean_2perc_baseline_new"

out_dir <- file.path(base_dir, "diagnostic_csv_analysis_60N")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ssps <- c("126", "245", "370", "585")

ssp_labels <- c(
  "126" = "SSP1-2.6",
  "245" = "SSP2-4.5",
  "370" = "SSP3-7.0",
  "585" = "SSP5-8.5"
)

ssp_colors <- c(
  "SSP1-2.6" = "blue",
  "SSP2-4.5" = "orange",
  "SSP3-7.0" = "#D73027",
  "SSP5-8.5" = "#7B3294"
)

# -----------------------------
# Read diagnostic CSVs
# -----------------------------

read_diag_csv <- function(ssp) {
  
  f <- file.path(
    base_dir,
    paste0("arctic_monthly_diagnostics_w_temp", ssp, "_1850_2100.csv")
  )
  
  if (!file.exists(f)) {
    warning("Missing file: ", f)
    return(NULL)
  }
  
  read_csv(f, show_col_types = FALSE) %>%
    mutate(
      SSP_raw = ssp,
      SSP = ssp_labels[ssp]
    )
}


diag_monthly <- bind_rows(lapply(ssps, read_diag_csv)) %>%
  mutate(
    SSP = factor(
      SSP,
      levels = c("SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5")
    )
  )

# Check column names
print(names(diag_monthly))
print(head(diag_monthly))

# -----------------------------
# Annual summaries
# -----------------------------
# from sum or mean over all month values to produce annual

diag_annual <- diag_monthly %>%
  group_by(SSP, Year) %>%
  summarise(
    mineralised_pg_yr = sum(mineralised_pg_monthly, na.rm = TRUE),
    total_inorg_pg_yr = sum(total_inorg_pg_monthly, na.rm = TRUE),
    mean_k_t = mean(k_T_mean_monthly, na.rm = TRUE),
    mean_k_t_weighted = mean(k_T_pool_weighted_mean, na.rm = TRUE),
    k_env_weighted_mean =
      mean(k_env_pool_weighted_mean, na.rm = TRUE),
    
    organic_pool_remaining_pg =
      mean(organic_pool_remaining_pg, na.rm = TRUE),
    
    thawed_N_annual_change_pg =
      mean(thawed_N_annual_change_pg, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  group_by(SSP) %>%
  arrange(Year, .by_group = TRUE) %>%
  mutate(
    # Cumulative annual thawed-N change
    cumulative_thawed_N_pg =
      cumsum(replace_na(thawed_N_annual_change_pg, 0)),
    
    cumulative_mineralised_pg =
      cumsum(mineralised_pg_yr),
    
    cumulative_inorg_pg =
      cumsum(total_inorg_pg_yr),
    
    ref_mineralised = mean(
      mineralised_pg_yr[Year >= 2000 & Year <= 2020],
      na.rm = TRUE
    ),
    
    ref_inorg = mean(
      total_inorg_pg_yr[Year >= 2000 & Year <= 2020],
      na.rm = TRUE
    ),
    
    ref_cum_mineralised = mean(
      cumulative_mineralised_pg[Year >= 2000 & Year <= 2020],
      na.rm = TRUE
    ),
    
    ref_cum_inorg = mean(
      cumulative_inorg_pg[Year >= 2000 & Year <= 2020],
      na.rm = TRUE
    ),
    
    annual_mineralised_anom_pg =
      mineralised_pg_yr - ref_mineralised,
    
    annual_inorg_anom_pg =
      total_inorg_pg_yr - ref_inorg,
    
    cumulative_mineralised_anom_pg =
      cumulative_mineralised_pg - ref_cum_mineralised,
    
    cumulative_inorg_anom_pg =
      cumulative_inorg_pg - ref_cum_inorg,
    
    mineralised_20yr =
      zoo::rollapply(
        mineralised_pg_yr,
        20,
        mean,
        align = "center",
        fill = NA
      ),
    
    inorg_20yr =
      zoo::rollapply(
        total_inorg_pg_yr,
        20,
        mean,
        align = "center",
        fill = NA
      ),
    
    cumulative_mineralised_anom_20yr =
      zoo::rollapply(
        cumulative_mineralised_anom_pg,
        20,
        mean,
        align = "center",
        fill = NA
      ),
    
    cumulative_inorg_anom_20yr =
      zoo::rollapply(
        cumulative_inorg_anom_pg,
        20,
        mean,
        align = "center",
        fill = NA
      )
  ) %>%
  ungroup()


# -----------------------------
# Plot function
# -----------------------------

plot_diag <- function(df, yvar, ylab, title, filename) {
  
  p <- ggplot(df, aes(x = Year, y = .data[[yvar]], color = SSP)) +
    geom_line(linewidth = 0.8, na.rm = TRUE) +
    geom_vline(xintercept = 2015, linetype = "dashed", color = "black") +
    scale_color_manual(values = ssp_colors, drop = FALSE) +
    theme_minimal() +
    #ylim(-0.2, 1.5)+
    theme(
      legend.position = "bottom",
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12)
    ) +
    labs(
      x = "Year",
      y = ylab,
      color = "SSP scenario"
      #title = title
    )
  
  ggsave(
    file.path(out_dir, filename),
    p,
    width = 8,
    height = 5,
    dpi = 300
  )
  
  return(p)
}


# -----------------------------
# Make plots
# -----------------------------

p1 <- plot_diag(
  diag_annual,
  "mineralised_pg_yr",
  "Annual mineralised N [Pg N yr⁻¹]",
  "Mineralised N, annual",
  "annual_mineralised_N.png"
)
p1
p1 <- plot_diag(
  diag_annual,
  "annual_mineralised_anom_pg",
  "Annual mineralised N [Pg N yr⁻¹]",
  "Mineralised N, annual",
  "annual_mineralised_N.png"
)
p1
p1 <- plot_diag_facet(
  diag_annual,
  "mineralised_pg_yr",
  "Annual mineralised N [Pg N yr⁻¹]",
  "annual_mineralised_N_faceted.png"
)

p1
ggsave(
  file.path(out_dir, "annual_mineralised_N_with_20yr_mean.png"),
  p1,
  width = 8,
  height = 5,
  dpi = 300
)
p1

p2 <- plot_diag(
  diag_annual,
  "annual_inorg_anom_pg",
  "Annual inorganic N [Pg N yr⁻¹]",
  "Total inorganic N, annual",
  "annual_bioavailable_N.png"
)



p2 <- ggplot(
  diag_annual,
  aes(x = Year, color = SSP, fill = SSP)
) +
  # Annual bioavailable N as thin bars, one per year
  geom_col(
    aes(y = bioavailable_pg_yr),
    width = 0.8,
    alpha = 0.35,
    position = "identity",
    na.rm = TRUE
  ) +
  
  # 20-year centred moving average
  geom_line(
    aes(y = bioavailable_20yr),
    linewidth = 1.2,
    na.rm = TRUE
  ) +
  
  geom_vline(
    xintercept = 2015,
    linetype = "dashed",
    color = "black"
  ) +
  
  scale_color_manual(values = ssp_colors, drop = FALSE) +
  scale_fill_manual(values = ssp_colors, drop = FALSE) +
  
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 10)
  ) +
  
  labs(
    x = "Year",
    y = "Annual bioavailable N [Pg N yr⁻¹]",
    color = "SSP scenario",
    fill = "SSP scenario"
  )

ggsave(
  file.path(out_dir, "annual_bioavailable_N_with_20yr_mean.png"),
  p2,
  width = 8,
  height = 5,
  dpi = 300
)

p2

p3 <- plot_diag(
  diag_annual,
  "cumulative_mineralised_anom_pg",
  "Cumulative mineralised N [Pg N]",
  "Mineralised N, cumulative",
  "cumulative_mineralised_N_anomaly.png"
)
p3
p4 <- plot_diag(
  diag_annual,
  "cumulative_inorg_anom_pg",
  "Cumulative total inorganic N [Pg N]",
  "Total inorganic N, cumulative",
  "cumulative_bioavailable_N_anomaly.png"
)

p5 <- plot_diag(
  diag_annual,
  "k_env_weighted_mean",
  "Mean k env. - weighted",
  "Mean k environment",
  "mean_kT.png"
)


p6 <- plot_diag(
  diag_annual,
  "organic_pool_remaining_pg",
  "Organic N pool [Pg N]",
  "Organic N, cumulative",
  "organic_pool_remaining.png"
)

p7 <- plot_diag(
  diag_annual,
  "thawed_N_annual_change_pg",
  "Newly thawed N annual change [Pg N yr⁻¹]",
  "thawed_N_annual_change_pg.png"
)

p7 <- ggplot(
  diag_annual,
  aes(x = Year, color = SSP, fill = SSP)
) +
  # Annual bioavailable N as thin bars, one per year
  geom_col(
    aes(y = thawed_N_annual_change_pg),
    width = 0.8,
    alpha = 0.35,
    position = "identity",
    na.rm = TRUE
  ) +
  
  # 20-year centred moving average
  geom_line(
    aes(y = thawed_N_annual_change_pg),
    linewidth = 1.2,
    na.rm = TRUE
  ) +
  
  geom_vline(
    xintercept = 2015,
    linetype = "dashed",
    color = "black"
  ) +
  
  scale_color_manual(values = ssp_colors, drop = FALSE) +
  scale_fill_manual(values = ssp_colors, drop = FALSE) +
  
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 10)
  ) +
  
  labs(
    x = "Year",
    y = "Annual bioavailable N [Pg N yr⁻¹]",
    color = "SSP scenario",
    fill = "SSP scenario"
  )

p8 <- plot_diag(
  diag_annual,
  "k_env_weighted_mean",
  "k_env_pool_weighted [-]",
  "Mean environmental scaling factor, temperature + moisture",
  "mean_moisture.png"
)

print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
print(p7)


diag_annual_roll <- diag_annual %>%
  arrange(SSP, Year) %>%
  group_by(SSP) %>%
  mutate(
    mineralised_pg_yr_roll20     = rollmean(mineralised_pg_yr, k = 20, fill = NA, align = "center"),
    thawed_N_annual_change_roll20 = rollmean(thawed_N_annual_change_pg, k = 20, fill = NA, align = "center")
  ) %>%
  ungroup()

ggplot(diag_annual_roll, aes(x = Year, colour = SSP)) +
  geom_line(aes(y = mineralised_pg_yr), alpha = 0.25, linewidth = 0.4) +
  geom_line(aes(y = mineralised_pg_yr_roll20), linewidth = 1.1) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey40") +
  scale_colour_manual(values = c(
    "SSP1-2.6" = "blue",
    "SSP2-4.5" = "orange",
    "SSP3-7.0" = "red",
    "SSP5-8.5" = "purple"
  )) +
  labs(
    x = "Year",
    y = expression("Mineralised N [Pg " * yr^-1 * "]"),
    colour = "SSP scenario",
    title = "Annual mineralised N: raw vs. 20-yr rolling mean"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")


ggplot(diag_annual_roll, aes(x = Year, colour = SSP)) +
  geom_line(aes(y = thawed_N_annual_change_pg), alpha = 0.25, linewidth = 0.4) +
  geom_line(aes(y = thawed_N_annual_change_roll20), linewidth = 1.1) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey40") +
  scale_colour_manual(values = c(
    "SSP1-2.6" = "blue",
    "SSP2-4.5" = "orange",
    "SSP3-7.0" = "red",
    "SSP5-8.5" = "purple"
  )) +
  labs(
    x = "Year",
    y = expression("Thawed N annual change [Pg " * yr^-1 * "]"),
    colour = "SSP scenario",
    title = "Newly thawed N: raw vs. 20-yr rolling mean"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")


diag_long <- diag_annual_roll %>%
  select(SSP, Year, mineralised_pg_yr_roll20, thawed_N_annual_change_roll20) %>%
  pivot_longer(
    cols = c(mineralised_pg_yr_roll20, thawed_N_annual_change_roll20),
    names_to = "variable",
    values_to = "value"
  )

ggplot(diag_long, aes(x = Year, y = value, colour = SSP)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey40") +
  facet_wrap(~ variable, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = c(
    "SSP1-2.6" = "blue",
    "SSP2-4.5" = "orange",
    "SSP3-7.0" = "red",
    "SSP5-8.5" = "purple"
  )) +
  labs(x = "Year", y = NULL, colour = "SSP scenario") +
  theme_bw() +
  theme(legend.position = "bottom")


### effective mineralisation rate per month temperature and moisture modifiers"
k_base <- 0.02

annual_effective_rate <- diag_monthly %>%
  arrange(SSP_raw, Year, Month) %>%
  group_by(SSP_raw, Year) %>%
  summarise(
    
    active_months = sum(
      k_env_pool_weighted_mean > 0,
      na.rm = TRUE
    ),
    
    sum_k_env = sum(
      k_env_pool_weighted_mean,
      na.rm = TRUE
    ),
    
    # calculating the effective fraction of the organic N pool that is mineralised over an entire year, 
    # taking into account that mineralisation happens month after month and the pool gets smaller after each month.
    # (1 - k_base * k_env..) = how much organic N pool remains at the end of the year after monthly mineralisation;
    # 1- prod (1-k_base * k_env) )= fraction that is mineralised from organic pool = effective annual mineralisation rate
   
     effective_annual_rate = 1 - prod(
      1 - k_base * k_env_pool_weighted_mean,
      na.rm = TRUE
    ),
    
    approximate_annual_rate =
      k_base * sum_k_env,
    
    .groups = "drop"
  ) %>%
  mutate(
    effective_annual_percent = effective_annual_rate * 100,
    approximate_annual_percent = approximate_annual_rate * 100,
    
    SSP = factor(
      SSP_raw,
      levels = c("126", "245", "370", "585"),
      labels = c("1-2.6", "2-4.5", "3-7.0", "5-8.5")
    )
  )

ssp_colors <- c(
  "1-2.6" = "blue",
  "2-4.5" = "orange",
  "3-7.0" = "#D73027",
  "5-8.5" = "#7B3294"
)

annual_effective_rate %>%
  filter(Year == 2099) %>%
  select(SSP, Year, effective_annual_percent)

ggplot(
  annual_effective_rate,
  aes(
    x = Year,
    y = effective_annual_percent,
    colour = SSP
  )
) +
  geom_line(linewidth = 0.8) +
  geom_hline(
    
    yintercept = 2,
    
    linetype = "dashed",
    
    colour = "black"
    
  ) +
  scale_colour_manual(values = ssp_colors) +
  theme_bw() +
  labs(
    x = "Year",
    y = "Effective annual mineralisation rate (%)",
    colour = "SSP",
    title = "Effective annual mineralisation rate",
    subtitle = "Calculated from monthly temperature and moisture modifiers"
  )+
  
  annotate(
    
    "text",
    
    x = 1860,
    
    y = 2.03,
    
    label = "Base mineralisation parameter (2% yr⁻¹)",
    
    hjust = 0,
    
    size = 3
    
  )



range(
  k_base * diag_monthly$k_env_pool_weighted_mean,
  na.rm = TRUE
)





###monthly turnover rate: 

# ============================================================
# Mean monthly organic N turnover rate, all SSPs
# turnover_rate = k_base * k_env_pool_weighted (pool-weighted,
# as actually applied in mineralise_depth_pools())
# Reads the per-SSP diagnostic CSVs already written by the
# monthly mineralisation script.
# ============================================================

library(dplyr)
library(readr)
library(ggplot2)

ssps <- c("126", "245", "370", "585")

ssp_labels <- c(
  "126" = "SSP1-2.6",
  "245" = "SSP2-4.5",
  "370" = "SSP3-7.0",
  "585" = "SSP5-8.5"
)

ssp_colors <- c(
  "SSP1-2.6" = "blue",
  "SSP2-4.5" = "orange",
  "SSP3-7.0" = "#D73027",
  "SSP5-8.5" = "#7B3294"
)

k_base <- 0.02  # match whatever value was actually used in the runs

out_dir <- "monthly_mineralised/mean_2perc_baseline"
plot_out_dir <- file.path(out_dir, "turnover_rate_all_ssps")
dir.create(plot_out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Read and combine per-SSP diagnostic CSVs
# ------------------------------------------------------------

diagnostic_list <- lapply(ssps, function(ssp) {
  
  file <- file.path(
    out_dir,
    paste0("arctic_monthly_diagnostics_w_temp", ssp, "_1850_2100.csv")
  )
  
  if (!file.exists(file)) {
    stop("Missing file: ", file)
  }
  
  df <- read_csv(file, show_col_types = FALSE)
  
  df$SSP_raw <- ssp
  df$SSP <- unname(ssp_labels[ssp])
  
  df
})

diagnostic_all <- bind_rows(diagnostic_list) %>%
  mutate(
    turnover_rate_monthly = k_base * k_env_pool_weighted_mean,
    SSP = factor(SSP, levels = unname(ssp_labels))
  )

cat("Combined diagnostic rows:", nrow(diagnostic_all), "\n")



annual_rate_df <- diagnostic_all %>%
  group_by(Year) %>%
  summarise(
    effective_annual_rate =
      1 - prod(
        1 - k_base * k_env_pool_weighted_mean,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

annual_rate_df %>%
  filter(Year >= 2080, Year <= 2099) %>%
  summarise(
    mean_rate = mean(effective_annual_rate, na.rm = TRUE),
    median_rate = median(effective_annual_rate, na.rm = TRUE),
    min_rate = min(effective_annual_rate, na.rm = TRUE),
    max_rate = max(effective_annual_rate, na.rm = TRUE)
  )

# ------------------------------------------------------------
# 2. Annual mean monthly turnover rate, per SSP
# ------------------------------------------------------------

turnover_annual_all <- diagnostic_all %>%
  group_by(SSP_raw, SSP, Year) %>%
  summarise(
    mean_monthly_turnover_rate = mean(turnover_rate_monthly, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  turnover_annual_all,
  file.path(plot_out_dir, "all_ssps_annual_mean_monthly_turnover_rate.csv")
)

p_turnover_annual <- ggplot(
  turnover_annual_all,
  aes(x = Year, y = mean_monthly_turnover_rate * 100, color = SSP)
) +
  geom_line(linewidth = 1.0) +
  geom_hline(
    yintercept = k_base * 100,
    linetype = "dashed",
    color = "black"
  ) +
  annotate(
    "text",
    x = 1860,
    y = k_base * 100 + 0.02,
    label = paste0("k_base (", k_base * 100, "% mo\u207b\u00b9)"),
    hjust = 0,
    size = 4
  ) +
  scale_color_manual(values = ssp_colors) +
  theme_bw(base_size = 15) +
  labs(
    title = "Mean annual organic N turnover rate",
    subtitle = "Pool-weighted, k_base \u00d7 k_env, all SSPs",
    x = "Year",
    y = "Turnover rate (% of pool per month)",
    color = "Scenario"
  )

print(p_turnover_annual)

ggsave(
  file.path(plot_out_dir, "all_ssps_annual_mean_monthly_turnover_rate.png"),
  p_turnover_annual,
  width = 9,
  height = 5.5,
  dpi = 300
)

# ------------------------------------------------------------
# 3. Monthly seasonal cycle of turnover rate, per SSP
#    (mean seasonal cycle over a recent reference period)
# ------------------------------------------------------------

turnover_seasonal_all <- diagnostic_all %>%
  filter(Year >= 2080, Year <= 2099) %>%
  group_by(SSP_raw, SSP, Month) %>%
  summarise(
    mean_monthly_turnover_rate = mean(turnover_rate_monthly, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  turnover_seasonal_all,
  file.path(plot_out_dir, "all_ssps_seasonal_turnover_rate_2000_2020.csv")
)

p_turnover_seasonal <- ggplot(
  turnover_seasonal_all,
  aes(x = Month, y = mean_monthly_turnover_rate * 100, color = SSP)
) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  scale_color_manual(values = ssp_colors) +
  theme_bw(base_size = 15) +
  labs(
    title = "Seasonal cycle of organic N turnover rate",
    subtitle = "2080-2099 mean, all SSPs",
    x = "Month",
    y = "Turnover rate (% of pool per month)",
    color = "Scenario"
  )

print(p_turnover_seasonal)

ggsave(
  file.path(plot_out_dir, "all_ssps_seasonal_turnover_rate_2000_2020.png"),
  p_turnover_seasonal,
  width = 9,
  height = 5.5,
  dpi = 300
)

cat("Done. Outputs saved in:\n", plot_out_dir, "\n")



library(dplyr)
library(readr)
library(ggplot2)

# ============================================================
# Effective annual mineralisation rate for all SSPs
# ============================================================

base_dir <- "monthly_mineralised/mean_2perc_baseline"

ssps <- c("126", "245", "370", "585")

ssp_labels <- c(
  "126" = "1-2.6",
  "245" = "2-4.5",
  "370" = "3-7.0",
  "585" = "5-8.5"
)

ssp_colors <- c(
  "1-2.6" = "blue",
  "2-4.5" = "orange",
  "3-7.0" = "#D73027",
  "5-8.5" = "#7B3294"
)

k_base <- 0.02


# ------------------------------------------------------------
# 1. Read diagnostics and calculate annual effective rate
# ------------------------------------------------------------

annual_rate_df <- lapply(ssps, function(ssp) {
  
  file <- file.path(
    base_dir,
    paste0(
      "arctic_monthly_diagnostics_w_temp",
      ssp,
      "_1850_2100.csv"
    )
  )
  
  df <- read_csv(file, show_col_types = FALSE)
  
  df %>%
    arrange(Year, Month) %>%
    group_by(Year) %>%
    summarise(
      
      # Environmental activity summed across the year
      sum_k_env = sum(
        k_env_pool_weighted_mean,
        na.rm = TRUE
      ),
      
      # Approximation without monthly compounding
      approximate_annual_rate =
        k_base * sum_k_env,
      
      # Exact effective annual fraction of the pool mineralised
      effective_annual_rate =
        1 - prod(
          1 - k_base * k_env_pool_weighted_mean,
          na.rm = TRUE
        ),
      
      .groups = "drop"
    ) %>%
    mutate(
      SSP = ssp_labels[ssp]
    )
  
}) %>%
  bind_rows()

head(annual_rate_df)

annual_rate_df %>%
  group_by(SSP) %>%
  summarise(
    mean_rate = mean(effective_annual_rate, na.rm = TRUE),
    min_rate  = min(effective_annual_rate, na.rm = TRUE),
    max_rate  = max(effective_annual_rate, na.rm = TRUE)
  )


p_annual_rate <- ggplot(
  annual_rate_df,
  aes(
    x = Year,
    y = effective_annual_rate * 100,
    colour = SSP
  )
) +
  geom_hline(
    yintercept = 2,
    linetype = "dotted",
    linewidth = 0.8,
    colour = "black"
  ) +
  
  geom_line(linewidth = 0.9) +
  
  geom_vline(
    xintercept = 2015,
    linetype = "dashed",
    colour = "grey50"
  ) +
  
  scale_colour_manual(
    values = ssp_colors
  ) +
  
  labs(
    x = "Year",
    y = expression("Effective mineralisation rate (% yr"^{-1}*")"),
    colour = "SSP"
  ) +
  
  theme_bw() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

print(p_annual_rate)

##################################################################

library(dplyr)
library(tidyr)
library(ggplot2)

pool_comparison_df <- diag_annual %>%
  select(
    SSP,
    Year,
    cumulative_thawed_N_pg,
    organic_pool_remaining_pg
  ) %>%
  pivot_longer(
    cols = c(
      cumulative_thawed_N_pg,
      organic_pool_remaining_pg
    ),
    names_to = "Variable",
    values_to = "Pg_N"
  ) %>%
  mutate(
    SSP = recode(
      as.character(SSP),
      "SSP1-2.6" = "1-2.6",
      "SSP2-4.5" = "2-4.5",
      "SSP3-7.0" = "3-7.0",
      "SSP5-8.5" = "5-8.5"
    ),
    Variable = recode(
      Variable,
      cumulative_thawed_N_pg = "Cumulative thawed N",
      organic_pool_remaining_pg = "Organic N pool remaining"
    )
  )

ggplot(
  pool_comparison_df,
  aes(
    x = Year,
    y = Pg_N,
    colour = SSP,
    linetype = Variable
  )
) +
  geom_line(linewidth = 0.9, na.rm = TRUE) +
  scale_colour_manual(values = ssp_colors) +
  scale_linetype_manual(
    values = c(
      "Cumulative thawed N" = "dashed",
      "Organic N pool remaining" = "solid"
    )
  ) +
  labs(
    title = "Cumulative thawed N and remaining organic N pool",
    x = "Year",
    y = "Pg N",
    colour = "SSP",
    linetype = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )


# just 370: 

unique(pool_comparison_df$SSP)
names(ssp_colors)
pool_comparison_df <- diag_annual %>%
  filter(SSP == "SSP3-7.0") %>%   # or "370" depending on your column
  select(
    Year,
    cumulative_thawed_N_pg,
    organic_pool_remaining_pg
  ) %>%
  pivot_longer(
    cols = c(
      cumulative_thawed_N_pg,
      organic_pool_remaining_pg
    ),
    names_to = "Variable",
    values_to = "Pg_N"
  ) %>%
  mutate(
    Variable = recode(
      Variable,
      cumulative_thawed_N_pg = "Cumulative thawed N",
      organic_pool_remaining_pg = "Organic N pool remaining"
    )
  )

ggplot(pool_comparison_df,
       aes(x = Year, y = Pg_N, linetype = Variable)) +
  geom_line(linewidth = 1) +
  scale_linetype_manual(
    values = c(
      "Cumulative thawed N" = "solid",
      "Organic N pool remaining" = "dashed"
    )
  ) +
  labs(
    title = "SSP3-7.0",
    x = "Year",
    y = "Pg N",
    linetype = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )
# ==============================================================================
# Script: Diagnostic Analysis of k-Factors & Mineralisation Dynamics
# Description: Compares temperature, moisture, and env k-factors (weighted 
#              vs. unweighted) and their relationships to mineralisation flux.
# ==============================================================================

library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork) # For combining plots

# ------------------------------------------------------------------------------
# 1. Summary Statistics Table
# ------------------------------------------------------------------------------
cat("--- Summary Statistics of k-Factors ---\n")

k_summary <- diag_monthly %>%
  select(
    starts_with("k_")
  ) %>%
  summarise(across(everything(), list(
    mean = ~mean(.x, na.rm = TRUE),
    sd   = ~sd(.x, na.rm = TRUE),
    min  = ~min(.x, na.rm = TRUE),
    max  = ~max(.x, na.rm = TRUE)
  ))) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("k_variable", "metric"),
    names_pattern = "(.*)_(mean|sd|min|max)$"
  ) %>%
  pivot_wider(names_from = metric, values_from = value)

print(k_summary)

# ------------------------------------------------------------------------------
# 2. Time-Series Comparison: Temperature vs. env k-Factors
# ------------------------------------------------------------------------------
# Group by year to see multi-decadal trends across k-factors
k_annual_trends <- diag_monthly %>%
  group_by(Year, SSP) %>%
  filter(abs(k_env_unweighted_monthly) < 10) %>%
  summarise(
    k_T_spatial          = mean(k_T_mean_monthly, na.rm = TRUE),
    k_T_pool_weighted    = mean(k_T_pool_weighted_mean, na.rm = TRUE),
    k_comb_pool_weighted       = mean(k_env_pool_weighted_mean, na.rm = TRUE),
    k_comb_spatial = mean(k_env_unweighted_monthly, na.rm = TRUE),
    k_moisture_spatial   = mean(k_moisture_mean_monthly, na.rm = TRUE),
    .groups = "drop"
  )

# Plot 1: Pool-Weighted vs. Spatial Unweighted k-Factors
p1 <- ggplot(k_annual_trends, aes(x = Year, color = SSP)) +
  geom_line(aes(y = k_comb_pool_weighted, linetype = "Pool-Weighted (Actual Driver)"), linewidth = 1) +
  geom_line(aes(y = k_comb_spatial, linetype = "Spatial Average (Unweighted)"), linewidth = 0.8, alpha = 0.7) +
  scale_linetype_manual(values = c("Pool-Weighted (Actual Driver)" = "solid", 
                                   "Spatial Average (Unweighted)" = "dashed")) +
  scale_colour_manual(
    values = c(
      "SSP1-2.6" = "blue",
      "SSP2-4.5" = "orange",
      "SSP3-7.0" = "red",
      "SSP5-8.5" = "purple"
    ),
    labels = c(
      "SSP1-2.6" = "1-2.6",
      "SSP2-4.5" = "2-4.5",
      "SSP3-7.0" = "3-7.0",
      "SSP5-8.5" = "5-8.5"
    )
  ) +
  labs(
    title = "env k-Factor: Pool-Weighted vs. Spatial Mean",
    subtitle = "Divergence shows where organic N spatial distribution differs from max decay rates",
    y = "k_env. [-]",
    linetype = "Weighting Scheme",
    color = "SSP"
  ) +
  theme_minimal(base_size = 9) +
  theme(legend.position = "bottom")

p1
# ------------------------------------------------------------------------------
# 3. Moisture Limitation Analysis
# ------------------------------------------------------------------------------
# Plot 2: Moisture factor suppression effect
p2 <- ggplot(k_annual_trends, aes(x = Year, color = SSP)) +
  geom_line(
    aes(y = k_T_pool_weighted, linetype = "k_T (Temperature Only)"),
    linewidth = 0.9
  ) +
  geom_line(
    aes(y = k_comb_pool_weighted, linetype = "k_env (Temp * Moisture)"),
    linewidth = 0.9
  ) +
  geom_line(
    aes(y = k_moisture_spatial, linetype = "k_moisture (Moisture)"),
    linewidth = 0.9
  ) +
  scale_linetype_manual(
    values = c(
      "k_T (Temperature Only)" = "dashed",
      "k_env (Temp * Moisture)" = "solid",
      "k_moisture (Moisture)" = "dotdash"
    )
  ) +
  labs(
    title = "Impact of Soil Moisture Limitation on k-Factor",
    subtitle = "Gap between lines represents moisture suppression",
    y = "Rate Constant Multiplier [-]",
    linetype = "Factor Type",
    color = "SSP"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

p2

print(p1 / p2)

# ------------------------------------------------------------------------------
# 4. Correlation Matrix & Driver Verification
# ------------------------------------------------------------------------------
cat("\n--- Correlation with Monthly Mineralisation Flux ---\n")

k_correlations <- diag_monthly %>%
  select(
    mineralised_pg_monthly,
    organic_pool_remaining_pg,
    k_T_mean_monthly,
    k_T_pool_weighted_mean,
    k_env_mean_monthly,
    k_env_pool_weighted_mean,
    k_moisture_mean_monthly
  ) %>%
  cor(use = "complete.obs")

print(round(k_correlations["mineralised_pg_monthly", ], 3))

# ------------------------------------------------------------------------------
# 5. Diagnostic Sanity Checks
# ------------------------------------------------------------------------------
cat("\n--- Diagnostic Sanity Checks ---\n")

# Check A: Is moisture factor staying strictly within expected bounds [0, 1]?
invalid_moisture <- sum(diag_monthly$k_moisture_mean_monthly < 0 | diag_monthly$k_moisture_mean_monthly > 1, na.rm = TRUE)
cat("1. Out-of-bounds k_moisture records:", invalid_moisture, "\n")

# Check B: Is pool-weighted k_env higher or lower than unweighted spatial mean?
diag_monthly <- diag_monthly %>%
  mutate(weighting_bias = k_env_pool_weighted_mean - k_env_mean_monthly)

mean_bias <- mean(diag_monthly$weighting_bias, na.rm = TRUE)
cat("2. Average Weighting Bias (Pool-Weighted minus Spatial Mean):", round(mean_bias, 6), "\n")
if (mean_bias < 0) {
  cat("   -> Interpretation: Organic N pools are concentrated in colder/drier depths or regions than average.\n")
} else {
  cat("   -> Interpretation: Organic N pools are concentrated in warmer/wetter depths or regions than average.\n")
}


k_ref <- mean(
  k_annual_trends$k_comb_pool_weighted[
    k_annual_trends$Year >= 1850 &
      k_annual_trends$Year <= 1900
  ],
  na.rm = TRUE
)

pool_ref <- mean(
  k_annual_trends$organic_pool_remaining_pg[
    k_annual_trends$Year >= 1850 &
      k_annual_trends$Year <= 1900
  ],
  na.rm = TRUE
)

diagnostic <- k_annual_trends %>%
  mutate(
    # Actual estimated mineralisation driver
    M_actual_driver =
      organic_pool_remaining_pg * k_comb_pool_weighted,
    
    # Pool effect only: k held constant
    M_pool_only =
      organic_pool_remaining_pg * k_ref,
    
    # k-factor effect only: pool held constant
    M_k_only =
      pool_ref * k_comb_pool_weighted
  )




################################################################

library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(zoo)
library(lubridate)

# -----------------------------
# Settings
# -----------------------------
base_dir <- "monthly_mineralised/mean_2perc_baseline"
out_dir <- file.path(base_dir, "yearly_nc_and_seasonal_analysis")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ssps <- c("126", "245", "370", "585")

years <- 1850:2099

ssp_labels <- c(
  "126" = "SSP1-2.6",
  "245" = "SSP2-4.5",
  "370" = "SSP3-7.0",
  "585" = "SSP5-8.5"
)

ssp_colors <- c(
  "SSP1-2.6" = "blue",
  "SSP2-4.5" = "orange",
  "SSP3-7.0" = "#D73027",
  "SSP5-8.5" = "#7B3294"
)

common_extent <- ext(-180, 180, 60, 90)

# -----------------------------
# 1. Seasonal mineralisation from diagnostics CSVs
# -----------------------------

read_diag_csv <- function(ssp) {
  
  f <- file.path(
    base_dir,
    paste0("arctic_monthly_diagnostics_w_temp", ssp, "_1850_2100.csv")
  )
  
  if (!file.exists(f)) {
    warning("Missing file: ", f)
    return(NULL)
  }
  
  read_csv(f, show_col_types = FALSE) %>%
    mutate(
      SSP_raw = ssp,
      SSP = ssp_labels[ssp]
    )
}

diag_monthly <- bind_rows(lapply(ssps, read_diag_csv)) %>%
  mutate(
    SSP = factor(SSP, levels = names(ssp_colors)),
    Month = as.integer(Month),
    Season = case_when(
      Month %in% c(12, 1, 2) ~ "DJF",
      Month %in% c(3, 4, 5) ~ "MAM",
      Month %in% c(6, 7, 8) ~ "JJA",
      Month %in% c(9, 10, 11) ~ "SON"
    ),
    Season = factor(Season, levels = c("DJF", "MAM", "JJA", "SON"))
  )

seasonal_min <- diag_monthly %>%
  group_by(SSP, Year, Season) %>%
  summarise(
    seasonal_mineralised_pg = sum(mineralised_pg_monthly, na.rm = TRUE),
    seasonal_bioavailable_pg = sum(total_inorg_pg_monthly, na.rm = TRUE),
    mean_k_t = mean(k_T_pool_weighted_mean, na.rm = TRUE),
    .groups = "drop"
  )


write_csv(
  seasonal_min,
  file.path(out_dir, "seasonal_mineralisation_bioavailable.csv")
)

p_seasonal_min <- ggplot(
  seasonal_min,
  aes(x = Year, y = seasonal_mineralised_pg, color = SSP)
) +
  geom_line(linewidth = 0.5, alpha = 0.8) +
  geom_vline(xintercept = 2015, linetype = "dashed") +
  facet_wrap(~ Season, ncol = 2, scales = "free_y") +
  scale_color_manual(values = ssp_colors, drop = FALSE) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(
    x = "Year",
    y = "Seasonal mineralised N [Pg N season⁻¹]",
    color = "SSP scenario"
  )

ggsave(
  file.path(out_dir, "seasonal_mineralised_N.png"),
  p_seasonal_min,
  width = 9,
  height = 6,
  dpi = 300
)

p_seasonal_share <- seasonal_min %>%
  group_by(SSP, Year) %>%
  mutate(
    annual_mineralised_pg = sum(seasonal_mineralised_pg, na.rm = TRUE),
    seasonal_share = seasonal_mineralised_pg / annual_mineralised_pg
  ) %>%
  ungroup() %>%
  ggplot(aes(x = Year, y = seasonal_share, color = SSP)) +
  geom_line(linewidth = 0.5, alpha = 0.8) +
  geom_vline(xintercept = 2015, linetype = "dashed") +
  facet_wrap(~ Season, ncol = 2) +
  scale_color_manual(values = ssp_colors, drop = FALSE) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(
    x = "Year",
    y = "Seasonal share of annual mineralisation [-]",
    color = "SSP scenario"
  )

ggsave(
  file.path(out_dir, "seasonal_share_mineralised_N.png"),
  p_seasonal_share,
  width = 9,
  height = 6,
  dpi = 300
)

# Mean seasonal cycle of mineralised N
# Uses diag_monthly created from the diagnostics CSVs
monthly_clim_both <- diag_monthly %>%
  filter(Year >= 2050, Year <= 2070) %>%
  group_by(SSP, Month) %>%
  summarise(
    Mineralised = mean(mineralised_pg_monthly, na.rm = TRUE),
    Total_inorg = mean(total_inorg_pg_monthly, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(
    cols = c(Mineralised, Total_inorg),
    names_to = "Variable",
    values_to = "Pg_N"
  )



p_clim_both <- ggplot(
  monthly_clim_both,
  aes(x = Month, y = Pg_N, color = Variable, linetype = Variable)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ SSP, ncol = 2) +
  scale_x_continuous(
    breaks = 1:12,
    labels = month.abb
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  ) +
  labs(
    title = "Mean seasonal cycle of mineralised and bioavailable N",
    x = "Month",
    y = "N flux [Pg N month⁻¹]",
    color = NULL,
    linetype = NULL
  )

print(p_clim_both)

ggsave(
  file.path(out_dir, "mean_seasonal_cycle_mineralised_bioavailable_N.png"),
  p_clim_both,
  width = 9,
  height = 6,
  dpi = 300
)

# k_env_pool_weighted_mean, k_T_pool_weighted_mean, k_env_unweighted_monthly, k_T_mean_monthly

monthly_kT <- diag_monthly %>%
  filter(
    (Year >= 1880 & Year <= 1900) |
      (Year >= 2000 & Year <= 2020) |
      (Year >= 2080 & Year <= 2099)
  ) %>%
  mutate(
    Period = case_when(
      Year >= 1880 & Year <= 1900 ~ "1880–1900",
      Year >= 2000 & Year <= 2020 ~ "2000–2020",
      Year >= 2080 & Year <= 2099 ~ "2080–2099"
    ),
    Period = factor(
      Period,
      levels = c("1880–1900", "2000–2020", "2080–2099")
    )
  ) %>%
  group_by(SSP, Period, Month) %>%
  summarise(
    kT = mean(k_T_mean_monthly, na.rm = TRUE),
    .groups = "drop"
  )

period_cols <- c(
  "1880–1900" = "steelblue",
  "2000–2020" = "orange",
  "2080–2099" = "firebrick"
)

p_clim_kT <- ggplot(
  monthly_kT,
  aes(
    x = Month,
    y = kT,
    colour = Period,
    group = Period
  )
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  facet_wrap(~ SSP, ncol = 2) +
  scale_colour_manual(values = period_cols) +
  scale_x_continuous(
    breaks = 1:12,
    labels = month.abb
  ) +
  theme_bw() +
  labs(
    x = "Month",
    y = "k_T_mean_monthly",
    colour = "Period"
  ) +
  theme(
    legend.position = "bottom"
  )

print(p_clim_kT)



library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

# -----------------------------
# Settings
# -----------------------------

base_dir <- "monthly_mineralised/latest_2perc_baseline"
out_dir <- file.path(base_dir, "yearly_nc_analysis")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ssps <- c("126", "245", "370", "585")
years <- 1850:2099

common_extent <- ext(-180, 180, 60, 90)

ssp_labels <- c(
  "126" = "SSP1-2.6",
  "245" = "SSP2-4.5",
  "370" = "SSP3-7.0",
  "585" = "SSP5-8.5"
)

ssp_colors <- c(
  "SSP1-2.6" = "blue",
  "SSP2-4.5" = "orange",
  "SSP3-7.0" = "#D73027",
  "SSP5-8.5" = "#7B3294"
)

# -----------------------------
# Variables to analyse
# -----------------------------

vars <- c(
  "thawed_N_annual_change" = "Thawed N annual change",
  "organic_N_pool_remaining" = "Organic N pool remaining",
  "new_thawed_total_N" = "New thawed total N",
  "net_organic_N_pool_change" = "Net organic N pool change",
  "new_thawed_organic_N" = "New thawed organic N", 
  "rapid_bioavailable_N" = "Rapid Bioavailable N"
)

# -----------------------------
# Find yearly NetCDF files
# -----------------------------

find_nc_file <- function(ssp, var_name) {
  
  files <- list.files(
    base_dir,
    pattern = "\\.nc$",
    full.names = TRUE
  )
  
  matches <- files[
    grepl(var_name, basename(files), fixed = TRUE) &
      grepl(ssp, basename(files), fixed = TRUE) &
      grepl("yearly", basename(files), fixed = TRUE)
  ]
  
  if (length(matches) == 0) {
    warning("No file found for ", var_name, " SSP ", ssp)
    return(NA_character_)
  }
  
  if (length(matches) > 1) {
    message("Multiple files found for ", var_name, " SSP ", ssp, ". Using:")
    message(matches[1])
  }
  
  matches[1]
}

# -----------------------------
# Summarise one yearly raster
# -----------------------------

summarise_yearly_raster <- function(nc_file, ssp, var_name, var_label) {
  
  if (is.na(nc_file) || !file.exists(nc_file)) {
    return(NULL)
  }
  
  cat("Reading:", nc_file, "\n")
  
  r <- rast(nc_file)
  r <- crop(r, common_extent)
  
  n <- nlyr(r)
  file_years <- years[seq_len(min(n, length(years)))]
  r <- r[[seq_along(file_years)]]
  
  area_r <- cellSize(r[[1]], unit = "m")
  area_r <- mask(area_r, r[[1]])
  
  if (var_name == "max_ALD") {
    
    total_value <- rep(NA_real_, nlyr(r))
    
    mean_value <- global(
      r,
      "mean",
      na.rm = TRUE
    )[, 1]
    
    total_unit <- NA_character_
    mean_unit <- "m"
    
  } else {
    
    total_value <- global(
      r * area_r,
      "sum",
      na.rm = TRUE
    )[, 1] / 1e12
    
    mean_value <- global(
      r,
      "mean",
      na.rm = TRUE
    )[, 1]
    
    total_unit <- "Pg N"
    mean_unit <- "kg N m-2"
  }
  
  tibble(
    Year = file_years,
    SSP_raw = ssp,
    SSP = ssp_labels[[ssp]],
    Variable_raw = var_name,
    Variable = var_label,
    total_pg = total_value,
    mean_value = mean_value,
    total_unit = total_unit,
    mean_unit = mean_unit
  )
}

# -----------------------------
# Read and summarise all files
# -----------------------------

yearly_df <- bind_rows(lapply(ssps, function(ssp) {
  
  bind_rows(lapply(names(vars), function(v) {
    
    f <- find_nc_file(ssp, v)
    
    summarise_yearly_raster(
      nc_file = f,
      ssp = ssp,
      var_name = v,
      var_label = vars[[v]]
    )
  }))
}))

yearly_df <- yearly_df %>%
  mutate(
    SSP = factor(
      SSP,
      levels = c("SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5")
    ),
    Variable = factor(Variable, levels = unname(vars))
  )

write_csv(
  yearly_df,
  file.path(out_dir, "yearly_nc_summary_all_variables.csv")
)

# -----------------------------
# Plot function
# -----------------------------

plot_yearly_var <- function(var_label, y_col = "total_pg") {
  
  df_plot <- yearly_df %>%
    filter(Variable == var_label)
  
  if (nrow(df_plot) == 0) {
    warning("No data for ", var_label)
    return(NULL)
  }
  
  if (y_col == "total_pg") {
    y_lab <- paste0(var_label, " [Pg N]")
  } else {
    unit_label <- unique(na.omit(df_plot$mean_unit))[1]
    y_lab <- paste0(var_label, " [", unit_label, "]")
  }
  
  p <- ggplot(df_plot, aes(x = Year, y = .data[[y_col]], color = SSP)) +
    geom_line(linewidth = 0.7, na.rm = TRUE) +
    geom_vline(xintercept = 2015, linetype = "dashed") +
    scale_color_manual(values = ssp_colors, drop = FALSE) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(
      x = "Year",
      y = y_lab,
      color = "SSP scenario"
    )
  
  safe_name <- gsub("[^A-Za-z0-9]+", "_", var_label)
  
  ggsave(
    file.path(out_dir, paste0("yearly_", safe_name, "_", y_col, ".png")),
    p,
    width = 8,
    height = 5,
    dpi = 300
  )
  
  p
}

# -----------------------------
# Make plots
# -----------------------------

plots_total_pg <- lapply(
  setdiff(unname(vars), "Max ALD"),
  plot_yearly_var,
  y_col = "total_pg"
)

plots_mean <- lapply(
  unname(vars),
  plot_yearly_var,
  y_col = "mean_value"
)

# -----------------------------
# Extra diagnostics
# -----------------------------

comparison_df <- yearly_df %>%
  filter(Variable != "Max ALD") %>%
  select(Year, SSP, Variable, total_pg) %>%
  pivot_wider(names_from = Variable, values_from = total_pg) %>%
  mutate(
    organic_fraction_new_thaw =
      `New thawed organic N` / `New thawed total N`,
    
    net_pool_vs_input =
      `Net organic N pool change` / `New thawed organic N`,
    
    mineralisation_loss_estimate =
      `New thawed organic N` - `Net organic N pool change`
  )

write_csv(
  comparison_df,
  file.path(out_dir, "yearly_nc_variable_comparison.csv")
)

p_org_fraction <- ggplot(
  comparison_df,
  aes(x = Year, y = organic_fraction_new_thaw, color = SSP)
) +
  geom_line(linewidth = 0.7, na.rm = TRUE) +
  geom_vline(xintercept = 2015, linetype = "dashed") +
  scale_color_manual(values = ssp_colors, drop = FALSE) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(
    x = "Year",
    y = "New thawed organic N / new thawed total N [-]",
    color = "SSP scenario"
  )

ggsave(
  file.path(out_dir, "organic_fraction_new_thaw.png"),
  p_org_fraction,
  width = 8,
  height = 5,
  dpi = 300
)

p_loss_estimate <- ggplot(
  comparison_df,
  aes(x = Year, y = mineralisation_loss_estimate, color = SSP)
) +
  geom_line(linewidth = 0.7, na.rm = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 2015, linetype = "dashed") +
  scale_color_manual(values = ssp_colors, drop = FALSE) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(
    x = "Year",
    y = "New thawed organic N - net organic pool change [Pg N yr⁻¹]",
    color = "SSP scenario"
  )

ggsave(
  file.path(out_dir, "estimated_mineralisation_loss_from_pool.png"),
  p_loss_estimate,
  width = 8,
  height = 5,
  dpi = 300
)

print(p_org_fraction)
print(p_loss_estimate)

# -----------------------------
# Print plots
# -----------------------------

print(p_seasonal_min)
print(p_seasonal_share)
print(p_org_fraction)
print(p_loss_estimate)


# -----------------------------
# Cumulative net organic N pool change
# and cumulative thawed N yearly increment
# -----------------------------
levels(yearly_df$Variable)

cum_df <- yearly_df %>%
  filter(
    Variable %in% c(
      "Net organic N pool change",
      "Thawed N annual change"
    )
  ) %>%
  arrange(SSP, Variable, Year) %>%
  group_by(SSP, Variable) %>%
  mutate(
    cumulative_pg = cumsum(total_pg)
  ) %>%
  ungroup()

write_csv(
  cum_df,
  file.path(out_dir, "cumulative_net_organic_and_thawed_increment.csv")
)

p_cum <- ggplot(
  cum_df,
  aes(x = Year, y = cumulative_pg, color = SSP)
) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  geom_vline(xintercept = 2015, linetype = "dashed") +
  facet_wrap(~ Variable, scales = "free_y", ncol = 1) +
  scale_color_manual(values = ssp_colors, drop = FALSE) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 11)
  ) +
  labs(
    x = "Year",
    y = "Cumulative N [Pg N]",
    color = "SSP scenario"
  )

ggsave(
  file.path(out_dir, "cumulative_net_organic_pool_change_and_thawed_increment.png"),
  p_cum,
  width = 8,
  height = 7,
  dpi = 300
)

print(p_cum)

cum_anom_df <- cum_df %>%
  group_by(SSP, Variable) %>%
  mutate(
    ref_mean = mean(cumulative_pg[Year >= 2000 & Year <= 2020], na.rm = TRUE),
    cumulative_anom_pg = cumulative_pg - ref_mean
  ) %>%
  ungroup()

p_cum_anom <- ggplot(
  cum_anom_df,
  aes(x = Year, y = cumulative_anom_pg, color = SSP)
) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 2015, linetype = "dashed") +
  facet_wrap(~ Variable, scales = "free_y", ncol = 1) +
  scale_color_manual(values = ssp_colors, drop = FALSE) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(
    x = "Year",
    y = "Cumulative N anomaly [Pg N]",
    color = "SSP scenario"
  )

ggsave(
  file.path(out_dir, "cumulative_anomaly_net_organic_pool_change_and_thawed_increment.png"),
  p_cum_anom,
  width = 8,
  height = 7,
  dpi = 300
)

print(p_cum_anom)



names(yearly_df)

str(yearly_df)

head(yearly_df)



#### cumulative thawed organic N vs total N 


# ============================================================
# Cumulative newly thawed organic N
# All SSPs, 1850–2099
#
# Input file example:
# arctic_new_thawed_organic_N_yearly_585_mean_1850_2099_w_temp.nc
#
# Assumptions:
# - Raster units are kg N m^-2
# - Each NetCDF layer represents one year
# - The first layer represents 1850
# - 1850 is an initialization/baseline year and is set to zero
# ============================================================

#=============================================


# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

library(terra)
library(dplyr)
library(purrr)
library(readr)
library(ggplot2)

terraOptions(
  memfrac = 0.4,
  progress = 1
)


# ------------------------------------------------------------
# 2. Settings
# ------------------------------------------------------------

base_dir <- "monthly_mineralised/arctic_final_mean"

out_dir <- file.path(
  base_dir,
  "cumulative_new_thawed_organic_N_analysis"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

ssps <- c("126", "245", "370", "585")

ssp_labels <- c(
  "126" = "SSP1-2.6",
  "245" = "SSP2-4.5",
  "370" = "SSP3-7.0",
  "585" = "SSP5-8.5"
)

ssp_colors <- c(
  "SSP1-2.6" = "blue",
  "SSP2-4.5" = "orange",
  "SSP3-7.0" = "#D73027",
  "SSP5-8.5" = "#7B3294"
)

start_year <- 1850
end_year   <- 2099


# ------------------------------------------------------------
# 3. Find one SSP NetCDF
# ------------------------------------------------------------

find_new_thawed_organic_file <- function(ssp) {
  
  file_pattern <- paste0(
    "^arctic_new_thawed_organic_N_yearly_",
    ssp,
    "_mean_1850_2099_w_temp\\.nc$"
  )
  
  candidates <- list.files(
    path = base_dir,
    pattern = file_pattern,
    recursive = TRUE,
    full.names = TRUE
  )
  
  if (length(candidates) == 0) {
    
    stop(
      "Could not find newly thawed organic N file for SSP ",
      ssp,
      " inside:\n",
      normalizePath(base_dir, mustWork = FALSE)
    )
  }
  
  if (length(candidates) > 1) {
    
    warning(
      "Multiple matching files found for SSP ",
      ssp,
      ". Using:\n",
      candidates[1]
    )
  }
  
  candidates[1]
}


# ------------------------------------------------------------
# 4. Extract annual pan-Arctic Pg N
# ------------------------------------------------------------
extract_annual_new_thawed_organic <- function(ssp) {
  
  nc_file <- find_new_thawed_organic_file(ssp)
  
  message("\nProcessing SSP ", ssp)
  message("Reading: ", nc_file)
  
  r <- rast(nc_file)
  
  expected_n_years <- end_year - start_year + 1
  
  if (nlyr(r) != expected_n_years) {
    stop(
      "Incorrect number of layers for SSP ",
      ssp,
      ". Expected ",
      expected_n_years,
      ", found ",
      nlyr(r)
    )
  }
  
  years <- start_year:end_year
  
  area_m2 <- cellSize(
    r[[1]],
    unit = "m",
    mask = TRUE
  )
  
  annual_pg <- numeric(nlyr(r))
  
  for (i in seq_len(nlyr(r))) {
    
    total_kg <- global(
      r[[i]] * area_m2,
      fun = "sum",
      na.rm = TRUE
    )[1, 1]
    
    annual_pg[i] <- total_kg / 1e12
  }
  
  tibble(
    SSP_raw = ssp,
    SSP = unname(ssp_labels[ssp]),
    Year = years,
    annual_new_thawed_organic_N_pg = annual_pg
  ) %>%
    mutate(
      cumulative_new_thawed_organic_N_pg =
        cumsum(annual_new_thawed_organic_N_pg)
    )
}

# ------------------------------------------------------------
# 5. Process all SSPs
# ------------------------------------------------------------

cumulative_thawed_organic_df <- map_dfr(
  ssps,
  extract_annual_new_thawed_organic
) %>%
  mutate(
    SSP = factor(
      SSP,
      levels = unname(ssp_labels[ssps])
    )
  ) %>%
  arrange(SSP, Year)


# ------------------------------------------------------------
# 6. Check original and cumulative starting values
# ------------------------------------------------------------

starting_values <- cumulative_thawed_organic_df %>%
  group_by(SSP) %>%
  slice_min(
    order_by = Year,
    n = 1,
    with_ties = FALSE
  ) %>%
  select(
    SSP,
    Year,
    annual_new_thawed_organic_N_pg,
    cumulative_new_thawed_organic_N_pg
  )

print(starting_values)


# ------------------------------------------------------------
# 7. Save time series
# ------------------------------------------------------------

write_csv(
  cumulative_thawed_organic_df,
  file.path(
    out_dir,
    "cumulative_new_thawed_organic_N_all_SSPs_1850_2099.csv"
  )
)


# ------------------------------------------------------------
# 8. Plot
# ------------------------------------------------------------
p_cumulative_thawed_organic <- ggplot(
  cumulative_thawed_organic_df,
  aes(
    x = Year,
    y = cumulative_new_thawed_organic_N_pg,
    colour = SSP
  )
) +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = ssp_colors) +
  labs(
    title = "Cumulative newly thawed organic nitrogen",
    subtitle = "Pan-Arctic land area north of 60°N",
    x = "Year",
    y = "Cumulative newly thawed organic N [Pg N]",
    colour = NULL
  ) +
  theme_bw()
print(p_cumulative_thawed_organic)


# ------------------------------------------------------------
# 9. Save plot
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    out_dir,
    "cumulative_new_thawed_organic_N_all_SSPs_1850_2099.png"
  ),
  plot = p_cumulative_thawed_organic,
  width = 9,
  height = 5.5,
  dpi = 300
)




# ============================================================
# 1. Read annual  new_thawed_total_N for all SSPs
# ============================================================

cumulative_new_thawed_total_N <- map_dfr(ssps, function(ssp) {
  
  nc_file <- find_nc_file(
    ssp = ssp,
    var_name = "new_thawed_total_N"
  )
  
  summarise_yearly_raster(
    nc_file = nc_file,
    ssp = ssp,
    var_name = "new_thawed_total_N",
    var_label = vars[["new_thawed_total_N"]]
  )
}) %>%
  arrange(SSP_raw, Year) %>%
  group_by(SSP_raw, SSP) %>%
  mutate(
    annual_new_thawed_total_N = total_pg,
    
    # Cumulative newly new_thawed_total_N since 1850
    cumulative_new_thawed_total_N = cumsum(
      replace_na(annual_new_thawed_total_N, 0)
    )
  ) %>%
  ungroup()


# Inspect result
print(cumulative_new_thawed_total_N)


# ============================================================
# 3. Plot all SSPs together
# ============================================================

p_cumulative_new_thawed_total_N <- ggplot(
  cumulative_new_thawed_total_N,
  aes(
    x = Year,
    y = cumulative_new_thawed_total_N,
    colour = SSP
  )
) +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = ssp_colors) +
  scale_x_continuous(
    breaks = seq(1850, 2100, by = 25),
    limits = c(1850, 2099),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    labels = scales::label_number(
      accuracy = 0.01,
      big.mark = ","
    ),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "Cumulative newly thawed total nitrogen",
    subtitle = "Pan-Arctic land area north of 60°N",
    x = "Year",
    y = "Cumulative newly thawed total N [Pg N]",
    colour = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

print(p_cumulative_new_thawed_total_N)


# ============================================================
# 4. Save plot
# ============================================================

ggsave(
  filename = file.path(
    out_dir,
    "cumulative_newly_thawed_organic_N_all_SSPs_1850_2099.png"
  ),
  plot = p_cumulative_thawed_organic,
  width = 9,
  height = 5.5,
  dpi = 300
)




### to compare rapid bioavailable, mineralised and total bioavailable N: 
unique(yearly_df$Variable)

#### actual saved rapid bioavailable
rapid_inorganic_df <- yearly_df %>%
  filter(Variable == "Rapid Bioavailable N") %>%
  mutate(
    Variable = "Rapid inorganic N"
  )

print(names(rapid_inorganic_df))
print(head(rapid_inorganic_df))

p_rapid_inorganic <- ggplot(
  rapid_inorganic_df,
  aes(
    x = Year,
    y = total_pg,
    colour = SSP,
    group = SSP
  )
) +
  geom_line(linewidth = 1) +
  theme_bw() +
  labs(
    title = "Annual rapid inorganic nitrogen",
    x = "Year",
    y = expression("Rapid inorganic N (Pg N yr"^{-1}*")"),
    colour = "Scenario"
  )

print(p_rapid_inorganic)


library(dplyr)
library(ggplot2)
unique(diag_monthly$Variable)
library(dplyr)
library(tidyr)
library(ggplot2)

# Annual mineralised + bioavailable from monthly diagnostics
annual_diag <- diag_monthly %>%
  group_by(SSP, Year) %>%
  summarise(
    `Mineralised N` = sum(mineralised_pg_monthly, na.rm = TRUE),
    `Bioavailable N` = sum(bioavailable_pg_monthly, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(`Mineralised N`, `Bioavailable N`),
    names_to = "Variable",
    values_to = "total_pg"
  )

# Rapid inorganic N from yearly thaw increment
rapid_inorganic_df <- yearly_df %>%
  filter(Variable == "Rapid Bioavailable N")

# Combine all three
compare_df <- bind_rows(
  annual_diag,
  rapid_inorganic_df
)

# Plot
ggplot(compare_df,
       aes(x = Year, y = total_pg, colour = Variable)) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~ SSP) +
  theme_bw() +
  labs(
    x = "Year",
    y = expression("Annual N flux (Pg N yr"^{-1}*")"),
    colour = NULL
  )
#
check_df <- compare_df %>%
  select(Year, SSP, Variable, total_pg) %>%
  tidyr::pivot_wider(names_from = Variable,
                     values_from = total_pg) %>%
  mutate(
    Difference = `Bioavailable N` - `Rapid Bioavailable N`
  )



cum_future_df <- annual_diag %>%
  filter(Year >= 2015) %>%
  arrange(SSP, Variable, Year) %>%
  group_by(SSP, Variable) %>%
  mutate(
    cumulative_pg = cumsum(total_pg)
  ) %>%
  ungroup()

ggplot(cum_future_df,
       aes(x = Year, y = cumulative_pg, colour = SSP)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ Variable, scales = "free_y") +
  theme_bw() +
  labs(
    x = "Year",
    y = "Cumulative N since 2015 (Pg N)",
    colour = "SSP scenario"
  )

cum_anomaly_df <- annual_diag %>%
  arrange(SSP, Variable, Year) %>%
  group_by(SSP, Variable) %>%
  mutate(
    cumulative_pg = cumsum(total_pg),
    cumulative_anomaly_pg = cumulative_pg - cumulative_pg[Year == 2015]
  ) %>%
  ungroup()

ggplot(cum_anomaly_df %>% filter(Year >= 2015),
       aes(x = Year, y = cumulative_anomaly_pg, colour = SSP)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ Variable, scales = "free_y") +
  theme_bw() +
  labs(
    x = "Year",
    y = "Cumulative anomaly since 2015 (Pg N)",
    colour = "SSP scenario"
  )

library(dplyr)
library(tidyr)
library(ggplot2)

bio_check_df <- diag_monthly %>%
  group_by(SSP, Year) %>%
  summarise(
    annual_mineralised_pg = sum(mineralised_pg_monthly, na.rm = TRUE),
    annual_bioavailable_pg = sum(bioavailable_pg_monthly, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    yearly_df %>%
      filter(Variable == "Thawed N yearly increment") %>%
      select(SSP, Year, thawed_increment_pg = total_pg),
    by = c("SSP", "Year")
  ) %>%
  mutate(
    expected_rapid_inorg_pg = thawed_increment_pg * 0.1136,
    expected_bioavailable_pg = annual_mineralised_pg + expected_rapid_inorg_pg,
    difference_pg = annual_bioavailable_pg - expected_bioavailable_pg
  )

summary(bio_check_df$difference_pg)

bio_check_df %>%
  summarise(
    max_abs_difference_pg = max(abs(difference_pg), na.rm = TRUE),
    mean_abs_difference_pg = mean(abs(difference_pg), na.rm = TRUE)
  )

ggplot(bio_check_df,
       aes(x = Year, y = difference_pg, colour = SSP)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 1) +
  theme_bw() +
  labs(
    x = "Year",
    y = "Bioavailable check difference (Pg N)",
    colour = "SSP scenario"
  )



library(dplyr)
library(tidyr)
library(ggplot2)

component_df <- diag_monthly %>%
  group_by(SSP, Year) %>%
  summarise(
    `Mineralised N` = sum(mineralised_pg_monthly, na.rm = TRUE),
    `Bioavailable N` = sum(bioavailable_pg_monthly, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  left_join(
    yearly_df %>%
      filter(Variable == "Thawed N annual change") %>%
      select(SSP, Year, thawed_N_annual_change_pg = total_pg),
    by = c("SSP", "Year")
  ) %>%
  mutate(
    `Rapid inorganic N` = thawed_N_annual_change_pg * 0.1136,
    `Mineralised fraction of bioavailable N` =
      `Mineralised N` / `Bioavailable N`
  )

component_long <- component_df %>%
  select(
    SSP,
    Year,
    `Mineralised N`,
    `Rapid inorganic N`,
    `Bioavailable N`
  ) %>%
  pivot_longer(
    cols = c(`Mineralised N`, `Rapid inorganic N`, `Bioavailable N`),
    names_to = "Component",
    values_to = "PgN"
  )

ggplot(component_long,
       aes(x = Year, y = PgN, colour = Component)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ SSP) +
  theme_bw() +
  labs(
    x = "Year",
    y = expression("Annual N flux (Pg N yr"^{-1}*")"),
    colour = NULL
  )

ggplot(component_df,
       aes(x = Year,
           y = `Mineralised fraction of bioavailable N`,
           colour = SSP)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  labs(
    x = "Year",
    y = "Mineralised / bioavailable N",
    colour = "SSP scenario"
  )


component_long <- component_df %>%
  select(
    SSP,
    Year,
    `Rapid inorganic N`,
    `Mineralised N`
  ) %>%
  pivot_longer(
    cols = c(`Rapid inorganic N`, `Mineralised N`),
    names_to = "Component",
    values_to = "PgN"
  )

ggplot(component_long,
       aes(x = Year,
           y = PgN,
           fill = Component)) +
  geom_area(position = "stack") +
  facet_wrap(~ SSP) +
  theme_bw() +
  labs(
    x = "Year",
    y = expression("Annual bioavailable N (Pg N yr"^{-1}*")"),
    fill = NULL
  )


library(dplyr)
library(tidyr)
library(ggplot2)

# Annual mineralised N
annual_components <- diag_monthly %>%
  group_by(SSP, Year) %>%
  summarise(
    mineralised_pg = sum(mineralised_pg_monthly, na.rm = TRUE),
    bioavailable_pg = sum(bioavailable_pg_monthly, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    yearly_df %>%
      filter(Variable == "Thawed N annual change") %>%
      select(SSP, Year, thawed_N_annual_change_pg = total_pg),
    by = c("SSP", "Year")
  ) %>%
  mutate(
    rapid_inorganic_pg = thawed_N_annual_change_pg * 0.1136
  )

cumulative_components <- annual_components %>%
  arrange(SSP, Year) %>%
  group_by(SSP) %>%
  mutate(
    cumulative_mineralised = cumsum(mineralised_pg),
    cumulative_rapid_inorganic = cumsum(rapid_inorganic_pg)
  ) %>%
  ungroup()

plot_df <- cumulative_components %>%
  select(
    SSP,
    Year,
    cumulative_mineralised,
    cumulative_rapid_inorganic
  ) %>%
  pivot_longer(
    cols = starts_with("cumulative"),
    names_to = "Component",
    values_to = "PgN"
  ) %>%
  mutate(
    Component = recode(
      Component,
      cumulative_mineralised = "Mineralised N",
      cumulative_rapid_inorganic = "Rapid inorganic N"
    )
  )

ggplot(plot_df,
       aes(x = Year,
           y = PgN,
           colour = SSP)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ Component, scales = "free_y") +
  theme_bw() +
  labs(
    x = "Year",
    y = "Cumulative N (Pg N)",
    colour = "SSP scenario"
  )

future_components <- annual_components %>%
  filter(Year >= 2015) %>%
  arrange(SSP, Year) %>%
  group_by(SSP) %>%
  mutate(
    cumulative_mineralised = cumsum(mineralised_pg),
    cumulative_rapid_inorganic = cumsum(rapid_inorganic_pg)
  ) %>%
  ungroup()


annual_thaw <- yearly_df %>%
  filter(Variable == "Thawed N annual change") %>%
  arrange(SSP, Year) %>%
  group_by(SSP) %>%
  mutate(
    cumulative_thaw = cumsum(total_pg)
  )

ggplot(annual_thaw,
       aes(Year, cumulative_thaw, colour = SSP)) +
  geom_line(linewidth = 1) +
  theme_bw() +
  labs(
    y = "Cumulative thawed total N (Pg N)"
  )

cumulative_thawed_df <- yearly_df %>%
  
  filter(Variable == "Thawed N annual change") %>%
  
  arrange(SSP, Year) %>%
  
  group_by(SSP) %>%
  
  mutate(
    
    cumulative_thawed_pg = cumsum(total_pg)
    
  ) %>%
  
  ungroup()

ggplot(cumulative_thawed_df,
       aes(x = Year,
           y = cumulative_thawed_pg,
           colour = SSP)) +
  geom_line(linewidth = 1.2) +
  geom_vline(xintercept = 2015, linetype = "dashed") +
  theme_bw(base_size = 16) +
  labs(
    x = "Year",
    y = "Cumulative thawed permafrost N (Pg N)",
    colour = "SSP scenario"
  )


cumulative_thawed_df <- yearly_df %>%
  filter(Variable == "Thawed N annual change") %>%
  arrange(SSP, Year) %>%
  group_by(SSP) %>%
  mutate(
    cumulative_thawed_pg = cumsum(total_pg),
    reference_2000_2020 = mean(
      cumulative_thawed_pg[Year >= 2000 & Year <= 2020],
      na.rm = TRUE
    ),
    cumulative_thawed_anomaly = cumulative_thawed_pg - reference_2000_2020
  ) %>%
  ungroup()

ggplot(
  cumulative_thawed_df,
  aes(
    x = Year,
    y = cumulative_thawed_anomaly,
    colour = SSP
  )
) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_line(linewidth = 1.2) +
  ylim(-7, 30)+
  scale_colour_manual(values = ssp_colors) +
  geom_vline(xintercept = 2015, linetype = "dashed") +
  theme_bw(base_size = 16) +
  labs(
    x = "Year",
    y = "Cumulative thawed permafrost N anomaly (Pg N)",
    colour = "SSP scenario"
  )











######### actual mineralisation rate
library(dplyr)
library(ggplot2)

# ============================================================
# 0. Setup
# ============================================================

ssp_list <- c("126", "245", "370", "585")   # adjust to match your actual folder naming

out_dir <- "monthly_mineralised/test_arctic_long"

# ============================================================
# 1. Load + compute realized turnover rate for each SSP
# ============================================================

annual_rate_all <- lapply(ssp_list, function(ssp) {
  
  diagnostic_file <- file.path(
    out_dir,
    paste0("arctic_monthly_diagnostics_w_temp", ssp, "_1850_2100.csv")
  )
  
  cat("Loading:", diagnostic_file, "\n")
  
  diagnostic_df <- read.csv(diagnostic_file, stringsAsFactors = FALSE)
  diagnostic_df$Date <- as.Date(diagnostic_df$Date)
  
  annual_rate_df <- diagnostic_df %>%
    group_by(Year) %>%
    summarise(
      mineralised_pg_yearly = sum(mineralised_pg_monthly, na.rm = TRUE),
      organic_pool_pg = first(organic_pool_remaining_pg),
      .groups = "drop"
    ) %>%
    mutate(
      realized_annual_turnover_fraction = ifelse(
        organic_pool_pg > 0.01,   # arbitrary floor — tune this
        mineralised_pg_yearly / organic_pool_pg,
        NA_real_
      ),
      SSP = ssp
    )
  annual_rate_df
  
}) %>% bind_rows()

# ============================================================
# 2. Recode SSP labels for nicer legend
# ============================================================

ssp_labels <- c(
  "126" = "SSP1-2.6",
  "245" = "SSP2-4.5",
  "370" = "SSP3-7.0",
  "585" = "SSP5-8.5"
)

annual_rate_all <- annual_rate_all %>%
  mutate(SSP = recode(SSP, !!!ssp_labels))

ssp_col<- c(
  "SSP1-2.6" = "blue",
  "SSP2-4.5" = "orange",
  "SSP3-7.0" = "red",
  "SSP5-8.5" = "purple"
)
# ============================================================
# 3. Save env csv
# ============================================================

write.csv(
  annual_rate_all,
  file.path(out_dir, "arctic_realized_annual_turnover_all_ssp_1850_2099.csv"),
  row.names = FALSE
)

# ============================================================
# 4. Plot all SSPs together
# ============================================================
  
p_turnover_all <- ggplot(
  annual_rate_all,
  aes(x = Year,
      y = realized_annual_turnover_fraction * 100,
      color = SSP)
) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = ssp_col) +
  theme_bw() +
  coord_cartesian(xlim = c(1950, 2100), ylim = c(0, 0.3)) +
  labs(
    title = "Actual annual mineralisation rate",
    subtitle = "Mineralised N / remaining organic N pool, by SSP scenario",
    x = "Year",
    y = "% of organic pool mineralised per year",
    color = "SSP scenario"
  )

print(p_turnover_all)

ggsave(
  filename = file.path(out_dir, "realized_annual_turnover_all_ssp.png"),
  plot = p_turnover_all,
  width = 8, height = 5.5, dpi = 300
)


effective_k =
  annual_mineralised /
  organic_pool_remaining







### rapid bioavailable
library(terra)
library(dplyr)
library(ggplot2)

# -----------------------------
# Setup
# -----------------------------

base_dir <- "monthly_mineralised/test_arctic_long"
out_dir <- file.path(base_dir, "rapid_inorganic_analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ssps <- c("126", "245", "370", "585")

ssp_labels <- c(
  "126" = "SSP1-2.6",
  "245" = "SSP2-4.5",
  "370" = "SSP3-7.0",
  "585" = "SSP5-8.5"
)

ssp_colors <- c(
  "SSP1-2.6" = "blue",
  "SSP2-4.5" = "orange",
  "SSP3-7.0" = "#D73027",
  "SSP5-8.5" = "#7B3294"
)

# -----------------------------
# Read and spatially sum files
# -----------------------------

rapid_list <- lapply(ssps, function(ssp) {
  
  file <- file.path(
    base_dir,
    paste0(
      "arctic_rapid_bioavailable_N_yearly_",
      ssp,
      "_mean_1850_2099_w_temp.nc"
    )
  )
  
  if (!file.exists(file)) {
    stop("Missing file: ", file)
  }
  
  rapid <- rast(file)
  
  area_rast <- cellSize(rapid[[1]], unit = "m")
  
  rapid_pg <- global(
    rapid * area_rast,
    "sum",
    na.rm = TRUE
  )[, 1] / 1e12
  
  years <- as.integer(format(time(rapid), "%Y"))
  
  # Fallback if the NetCDF time axis is missing
  if (length(years) != nlyr(rapid) || anyNA(years)) {
    years <- 1850:(1850 + nlyr(rapid) - 1)
  }
  
  data.frame(
    Year = years,
    SSP = as.character(ssp_labels[[ssp]]),
    Rapid_inorganic_Pg = rapid_pg
  )
})

rapid_df <- bind_rows(rapid_list) %>%
  arrange(SSP, Year) %>%
  group_by(SSP) %>%
  mutate(
    Cumulative_rapid_inorganic_Pg = cumsum(Rapid_inorganic_Pg),
    smoothed_20yr = zoo::rollmean(Rapid_inorganic_Pg, k = 20, fill = NA, align = "center")
  ) %>%
  ungroup()

rapid_df$SSP <- factor(
  rapid_df$SSP,
  levels = unname(ssp_labels)
)

# -----------------------------
# Annual plot
# -----------------------------

p_annual <- ggplot(
  rapid_df,
  aes(
    x = Year,
    y = Rapid_inorganic_Pg,
    colour = SSP
  )
) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = ssp_colors) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(
    title = "Annual rapid inorganic nitrogen",
    x = "Year",
    y = expression("Rapid inorganic N (Pg N yr"^{-1}*")"),
    colour = "Scenario"
  )

print(p_annual)



# -----------------------------
# Cumulative plot
# -----------------------------

p_cumulative <- ggplot(
  rapid_df,
  aes(
    x = Year,
    y = Cumulative_rapid_inorganic_Pg,
    colour = SSP
  )
) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = ssp_colors) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(
    title = "Cumulative rapidly bioavailable N",
    x = "Year",
    y = "[Pg N]",
    colour = "Scenario"
  )

print(p_cumulative)

rapid_df <- bind_rows(rapid_list) %>%
  arrange(SSP, Year) %>%
  group_by(SSP) %>%
  mutate(
    Cumulative_rapid_inorganic_Pg = cumsum(Rapid_inorganic_Pg),
    smoothed_20yr = zoo::rollmean(
      Rapid_inorganic_Pg,
      k = 20,
      fill = NA,
      align = "center"
    )
  ) %>%
  mutate(
    reference_2000_2014 = mean(
      Cumulative_rapid_inorganic_Pg[
        Year >= 2000 & Year <= 2014
      ],
      na.rm = TRUE
    ),
    Cumulative_relative_2000_2014 =
      Cumulative_rapid_inorganic_Pg - reference_2000_2014
  ) %>%
  ungroup()

p_cumulative <- ggplot(
  rapid_df,
  aes(
    x = Year,
    y = Cumulative_relative_2000_2014,
    colour = SSP
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey50"
  ) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = ssp_colors) +
  theme_bw() +
  theme(
    legend.position = "bottom"
  ) +
  labs(
    title = "Cumulative rapidly bioavailable N anomaly",
    x = "Year",
    y = expression(
      "Pg N"
    ),
    colour = "Scenario"
  )

print(p_cumulative)

# -----------------------------
# Save outputs
# -----------------------------

ggsave(
  file.path(out_dir, "annual_rapid_inorganic_N_all_SSPs.png"),
  p_annual,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(out_dir, "cumulative_rapid_inorganic_N_all_SSPs.png"),
  p_cumulative,
  width = 8,
  height = 5,
  dpi = 300
)

write.csv(
  rapid_df,
  file.path(out_dir, "rapid_inorganic_N_all_SSPs.csv"),
  row.names = FALSE
)





# ------------------------------------------------------------------#













##### with or without temp

library(terra)
library(dplyr)
library(ggplot2)

ssp <- "370"

dir_temp_8  <- "monthly_mineralised/whole_region_mean_wtemp_8deg/_w_temp_w_sm_yearly_nc_370"
dir_temp_12 <- "monthly_mineralised/whole_region_mean_wtemp_12deg/_w_temp_w_sm_yearly_nc_370"
dir_no_temp <- "monthly_mineralised/whole_region_mean_notemp/_w_temp_w_sm_yearly_nc_370"

years <- 1850:2099

read_monthly_series <- function(folder, label){
  
  monthly_list <- vector("list", length(years))
  
  for(i in seq_along(years)){
    
    yr <- years[i]
    
    f <- file.path(
      folder,
      paste0(
        "region_bioavailable_N_monthly_",
        ssp,
        "_",
        yr,
        "_w_temp_sm.nc"
      )
    )
    
    r <- rast(f)
    
    area <- cellSize(r[[1]], unit="m")
    
    monthly_pg <- global(
      r * area,
      "sum",
      na.rm=TRUE
    )[,1] / 1e12
    
    monthly_list[[i]] <- data.frame(
      Date = seq(
        as.Date(paste0(yr,"-01-15")),
        by="month",
        length.out=12
      ),
      Year = yr,
      Month = 1:12,
      Bioavailable = monthly_pg
    )
  }
  
  bind_rows(monthly_list) %>%
    arrange(Date) %>%
    mutate(
      cumulative_bioavailable = cumsum(Bioavailable),
      Run = label
    )
}

bio_8  <- read_monthly_series(dir_temp_8,  "Tref = 8°C")
bio_12 <- read_monthly_series(dir_temp_12, "Tref = 12°C")
bio_no <- read_monthly_series(dir_no_temp, "No temperature")

bio_all <- bind_rows(bio_8, bio_12, bio_no)
bio_temp <- bind_rows(bio_no, bio_8)
bio_temp <- bind_rows(bio_no, bio_12)

ggplot(bio_all,
       aes(Date, Bioavailable, colour=Run)) +
  geom_line() +
  theme_bw() +
  labs(
    y="Bioavailable N (Pg month⁻¹)",
    x=NULL
  )

ggplot(bio_all,
       aes(Date, cumulative_bioavailable, colour=Run)) +
  geom_line(linewidth=1) +
  theme_bw() +
  labs(
    y="Cumulative bioavailable N (Pg)",
    x=NULL
  )

comparison <- bio_all %>%
  select(Date, Run, cumulative_bioavailable) %>%
  tidyr::pivot_wider(
    names_from = Run,
    values_from = cumulative_bioavailable
  ) %>%
  mutate(
    diff_8 = `Tref = 8°C` - `No temperature`,
    diff_12 = `Tref = 12°C` - `No temperature`
  )

ggplot(comparison) +
  geom_line(aes(Date, diff_8, colour="8°C")) +
  geom_line(aes(Date, diff_12, colour="12°C")) +
  theme_bw() +
  labs(
    y="Additional cumulative bioavailable N (Pg)",
    colour=""
  )











#############################################################################


#### only 585
library(readr)
library(dplyr)
library(ggplot2)
library(zoo)

# ============================================================
# Settings
# ============================================================

base_dir <- "monthly_mineralised/60degN_test18_region"

out_dir <- file.path(base_dir, "diagnostic_csv_analysis_60N")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ssp <- "126"

# ============================================================
# Read diagnostic CSV
# ============================================================

diag_monthly <- read_csv(
  file.path(
    base_dir,
    paste0("arctic_monthly_diagnostics_w_temp", ssp, "_1850_2100.csv")
  ),
  show_col_types = FALSE
)

print(names(diag_monthly))
print(head(diag_monthly))
print(tail(diag_monthly))
# ============================================================
# Annual summaries
# ============================================================

diag_annual <- diag_monthly %>%
  group_by(Year) %>%
  summarise(
    mineralised_pg_yr = sum(mineralised_pg_monthly, na.rm = TRUE),
    bioavailable_pg_yr = sum(bioavailable_pg_monthly, na.rm = TRUE),
    mean_k_t = mean(k_T_mean_monthly, na.rm = TRUE),
    mean_k_t_weighted = mean(k_T_pool_weighted_mean, na.rm = TRUE),
    organic_pool_remaining_pg = mean(organic_pool_remaining_pg, na.rm = TRUE),
    thawed_N_annual_change_pg = mean(thawed_N_annual_change_pg, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Year) %>%
  mutate(
    cumulative_mineralised_pg = cumsum(mineralised_pg_yr),
    cumulative_bioavailable_pg = cumsum(bioavailable_pg_yr),
    
    #ref_cum_mineralised = mean(
    #  cumulative_mineralised_pg[Year >= 2000 & Year <= 2020],
    #  na.rm = TRUE
    #),
    
    #
    
    #cumulative_mineralised_anom_pg =
    #  cumulative_mineralised_pg #- ref_cum_mineralised,
    
    #cumulative_bioavailable_anom_pg =
     # cumulative_bioavailable_pg #- ref_cum_bioavailable,
    
    mineralised_20yr =
      rollapply(mineralised_pg_yr, 20,
                mean, align = "center", fill = NA),
    
    bioavailable_20yr =
      rollapply(bioavailable_pg_yr, 20,
                mean, align = "center", fill = NA),
    
    
  )

# ============================================================
# Plot function
# ============================================================

plot_diag <- function(df, yvar, ylab, filename) {
  
  p <- ggplot(df, aes(x = Year, y = .data[[yvar]])) +
    geom_line(
      linewidth = 0.9,
      colour = "#7B3294",
      na.rm = TRUE
    ) +
    geom_vline(
      xintercept = 2015,
      linetype = "dashed"
    ) +
    theme_minimal() +
    labs(
      x = "Year",
      y = ylab
    )
  
  ggsave(
    file.path(out_dir, filename),
    p,
    width = 8,
    height = 5,
    dpi = 300
  )
  
  p
}

# ============================================================
# Plots
# ============================================================

p1 <- plot_diag(
  diag_annual,
  "mineralised_pg_yr",
  "Annual mineralised N [Pg N yr⁻¹]",
  "annual_mineralised_N_585.png"
)

# Annual + 20-year mean
p2 <- ggplot(diag_annual, aes(x = Year)) +
  geom_line(
    aes(y = bioavailable_pg_yr),
    colour = "#7B3294",
    alpha = 0.35,
    linewidth = 0.5
  ) +
  geom_line(
    aes(y = bioavailable_20yr),
    colour = "#7B3294",
    linewidth = 1.2
  ) +
  geom_vline(
    xintercept = 2015,
    linetype = "dashed"
  ) +
  theme_minimal() +
  labs(
    x = "Year",
    y = "Annual bioavailable N [Pg N yr⁻¹]"
  )

ggsave(
  file.path(out_dir, "annual_bioavailable_N_with_20yr_mean_585.png"),
  p2,
  width = 8,
  height = 5,
  dpi = 300
)

p3 <- plot_diag(
  diag_annual,
  "cumulative_mineralised_anom_pg",
  "Cumulative mineralised N anomaly [Pg N]",
  "cumulative_mineralised_N_anomaly_585.png"
)

p4 <- plot_diag(
  diag_annual,
  "cumulative_bioavailable_anom_pg",
  "Cumulative bioavailable N anomaly [Pg N]",
  "cumulative_bioavailable_N_anomaly_585.png"
)

p5 <- plot_diag(
  diag_annual,
  "mean_k_t",
  "Mean kT",
  "mean_kT_585.png"
)

p6 <- plot_diag(
  diag_annual,
  "organic_pool_remaining_pg",
  "Organic N pool remaining [Pg N]",
  "organic_pool_remaining_585.png"
)

p7 <- plot_diag(
  diag_annual,
  "thawed_N_annual_change_pg",
  "Newly thawed N annual change [Pg N yr⁻¹]",
  "thawed_N_annual_change_585.png"
)

print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)
print(p7)



library(terra)
library(dplyr)
library(ggplot2)

out_dir <- "monthly_mineralised/mean_2perc_baseline"
ssp_list <- c("126", "245", "370", "585")   # adjust to match your actual SSP codes
years_to_plot <- 1850:2099

# -----------------------------
# 1. Read k_env yearly files and compute area-weighted spatial mean
#    both monthly and annually, for each SSP
# -----------------------------

read_k_env_ssp <- function(ssp, years) {
  
  out_dir_yearly <- file.path(out_dir, paste0("yearly_nc_", ssp))
  
  monthly_rows <- vector("list", length(years))
  
  for (i in seq_along(years)) {
    yr <- years[i]
    
    f <- file.path(out_dir_yearly, paste0("arctic_k_env_monthly_w_temp_sm_", ssp, "_", yr, ".nc"))
    
    if (!file.exists(f)) {
      warning("Missing file: ", f)
      next
    }
    
    r <- rast(f)   # 12 layers, Jan-Dec
    area_rast <- cellSize(r[[1]], unit = "m")
    
    monthly_means <- terra::global(r, "mean", weights = area_rast, na.rm = TRUE)[, 1]
    
    monthly_rows[[i]] <- data.frame(
      SSP = ssp,
      Year = yr,
      Month = 1:12,
      Date = seq(as.Date(paste0(yr, "-01-16")), by = "month", length.out = 12),
      k_env_monthly = monthly_means
    )
  }
  
  bind_rows(monthly_rows)
}

k_env_monthly_all <- bind_rows(lapply(ssp_list, read_k_env_ssp, years = years_to_plot))

# -----------------------------
# 2. Annual mean (average across the 12 months within each year)
# -----------------------------

k_env_annual_all <- k_env_monthly_all %>%
  group_by(SSP, Year) %>%
  summarise(
    k_env_annual = mean(k_env_monthly, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    SSP = recode(
      SSP,
      "126" = "SSP1-2.6",
      "245" = "SSP2-4.5",
      "370" = "SSP3-7.0",
      "585" = "SSP5-8.5"
    )
  )
# -----------------------------
# 3. Define the three time windows and label the monthly data
# -----------------------------

time_windows <- list(
  "1880-1900" = 1880:1900,
  "2000-2020" = 2000:2020,
  "2080-2099" = 2080:2099
)

assign_window <- function(yr) {
  for (w in names(time_windows)) {
    if (yr %in% time_windows[[w]]) return(w)
  }
  NA_character_
}

k_env_monthly_windows <- k_env_monthly_all %>%
  mutate(Window = sapply(Year, assign_window)) %>%
  filter(!is.na(Window)) %>%
  mutate(Window = factor(Window, levels = names(time_windows)))

# -----------------------------
# 3b. Average across years WITHIN each window, per month, per SSP
#     -> one 12-month seasonal cycle per SSP per window
# -----------------------------

k_env_seasonal <- k_env_monthly_windows %>%
  group_by(SSP, Window, Month) %>%
  summarise(k_env_mean = mean(k_env_monthly, na.rm = TRUE), .groups = "drop")
# -----------------------------
# 4. Plot styling
# -----------------------------

ssp_labels <- c(
  "126" = "SSP1-2.6",
  "245" = "SSP2-4.5",
  "370" = "SSP3-7.0",
  "585" = "SSP5-8.5"
)

window_colors <- c(
  "1880-1900" = "steelblue",
  "2000-2020" = "darkorange",
  "2080-2099" = "firebrick"
)


ssp_colors <- c(
  "SSP1-2.6" = "blue",
  "SSP2-4.5" = "orange",
  "SSP3-7.0" = "#D73027",
  "SSP5-8.5" = "#7B3294"
)



# -----------------------------
# 5. Plot: one line per time window, faceted by SSP
# -----------------------------

p_seasonal <- ggplot(k_env_seasonal, aes(x = Month, y = k_env_mean, colour = Window)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  facet_wrap(~ SSP, nrow = 2, labeller = as_labeller(ssp_labels)) +
  scale_colour_manual(values = window_colors, name = "Time period") +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Mean seasonal cycle of k_env by SSP scenario and time period",
    subtitle = "Spatial mean, weighted by the soil layer thickness",
    x = "Month",
    y = expression(k[env]~"[-]")
  )

print(p_seasonal)

# -----------------------------
# 6. Plot: annual mean time series, full 1850-2099, coloured by SSP
# -----------------------------

p_annual <- ggplot(k_env_annual_all, aes(x = Year, y = k_env_annual, colour = SSP)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = ssp_colors, labels = ssp_labels, name = "SSP scenario") +
  theme_bw() +
  labs(
    title = "Annual mean k_env by SSP scenario",
    x = "Year",
    y = expression(k[env]~"[-]")
  )

print(p_annual)


#### with std deviation


library(terra)
library(dplyr)
library(ggplot2)

# ============================================================
# SETTINGS
# ============================================================

out_dir_mean <- "monthly_mineralised/mean_2perc_baseline"
out_dir_sd   <- "monthly_mineralised/std_2perc_baseline"

ssp_list <- c("126", "245", "370", "585")
years_to_plot <- 1850:2099


# ============================================================
# 1. READ k_env MEAN + STD FILES
# ============================================================
read_k_env_ssp <- function(ssp, years) {
  
  out_dir_yearly_mean <- file.path(
    out_dir_mean,
    paste0("yearly_nc_", ssp)
  )
  
  out_dir_yearly_sd <- file.path(
    out_dir_sd,
    paste0("yearly_nc_", ssp)
  )
  
  monthly_rows <- vector("list", length(years))
  
  for (i in seq_along(years)) {
    
    yr <- years[i]
    
    # -----------------------------
    # Mean file
    # -----------------------------
    
    f_mean <- file.path(
      out_dir_yearly_mean,
      paste0(
        "arctic_k_env_monthly_w_temp_sm_",
        ssp, "_", yr, ".nc"
      )
    )
    
    # -----------------------------
    # SD file
    # -----------------------------
    
    f_sd <- file.path(
      out_dir_yearly_sd,
      paste0(
        "arctic_k_env_monthly_w_temp_sm_",
        ssp, "_", yr, ".nc"
      )
    )
    
    if (!file.exists(f_mean)) {
      warning("Missing mean file: ", f_mean)
      next
    }
    
    if (!file.exists(f_sd)) {
      warning("Missing SD file: ", f_sd)
      next
    }
    
    # -----------------------------
    # Read rasters
    # -----------------------------
    
    r_mean <- rast(f_mean)
    r_sd   <- rast(f_sd)
    
    # -----------------------------
    # Check dimensions
    # -----------------------------
    
    if (nlyr(r_mean) != 12) {
      warning(
        "Mean file does not have 12 layers: ",
        f_mean
      )
    }
    
    if (nlyr(r_sd) != 12) {
      warning(
        "SD file does not have 12 layers: ",
        f_sd
      )
    }
    
    # -----------------------------
    # Area weighting
    # -----------------------------
    
    area_rast <- cellSize(
      r_mean[[1]],
      unit = "m"
    )
    
    # -----------------------------
    # Monthly spatial mean
    # -----------------------------
    
    monthly_means <- terra::global(
      r_mean,
      "mean",
      weights = area_rast,
      na.rm = TRUE
    )[, 1]
    
    # -----------------------------
    # Monthly spatial SD
    #
    # IMPORTANT:
    # This follows your diagnostic
    # aggregation approach:
    #
    # mean quantity -> mean()
    # SD quantity   -> mean()
    # -----------------------------
    
    monthly_sds <- terra::global(
      r_sd,
      "mean",
      weights = area_rast,
      na.rm = TRUE
    )[, 1]
    
    # -----------------------------
    # Store
    # -----------------------------
    
    monthly_rows[[i]] <- data.frame(
      SSP = ssp,
      Year = yr,
      Month = 1:12,
      Date = seq(
        as.Date(paste0(yr, "-01-16")),
        by = "month",
        length.out = 12
      ),
      k_env_monthly = monthly_means,
      k_env_monthly_sd = monthly_sds
    )
  }
  
  bind_rows(monthly_rows)
}
# ============================================================
# 2. READ ALL SSPs
# ============================================================

k_env_monthly_all <- bind_rows(
  lapply(
    ssp_list,
    read_k_env_ssp,
    years = years_to_plot
  )
)


# ============================================================
# 3. ANNUAL MEAN + ANNUAL SD
# ============================================================

ssp_labels <- c(
  "126" = "SSP1-2.6",
  "245" = "SSP2-4.5",
  "370" = "SSP3-7.0",
  "585" = "SSP5-8.5"
)

k_env_monthly_all <- bind_rows(
  lapply(
    ssp_list,
    read_k_env_ssp,
    years = years_to_plot
  )
) %>%
  mutate(
    SSP = recode(
      SSP,
      !!!ssp_labels
    ),
    SSP = factor(
      SSP,
      levels = c(
        "SSP1-2.6",
        "SSP2-4.5",
        "SSP3-7.0",
        "SSP5-8.5"
      )
    )
  )

k_env_annual_all <- k_env_monthly_all %>%
  group_by(SSP, Year) %>%
  summarise(
    
    # Mean: mean of monthly values
    k_env_annual = mean(
      k_env_monthly,
      na.rm = TRUE
    ),
    
    # SD: mean of monthly SD values
    k_env_annual_sd = mean(
      k_env_monthly_sd,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


# ============================================================
# 4. DEFINE TIME WINDOWS
# ============================================================

time_windows <- list(
  "1880-1900" = 1880:1900,
  "2000-2020" = 2000:2020,
  "2080-2099" = 2080:2099
)


assign_window <- function(yr) {
  
  for (w in names(time_windows)) {
    
    if (yr %in% time_windows[[w]]) {
      return(w)
    }
    
  }
  
  NA_character_
}


k_env_monthly_windows <- k_env_monthly_all %>%
  mutate(
    Window = sapply(
      Year,
      assign_window
    )
  ) %>%
  filter(!is.na(Window)) %>%
  mutate(
    Window = factor(
      Window,
      levels = names(time_windows)
    )
  )


# ============================================================
# 5. SEASONAL CYCLE
# ============================================================

k_env_seasonal <- k_env_monthly_windows %>%
  group_by(SSP, Window, Month) %>%
  summarise(
    
    # Mean across years
    k_env_mean = mean(
      k_env_monthly,
      na.rm = TRUE
    ),
    
    # SD across years, following your
    # existing diagnostic convention:
    # mean of the monthly SD values
    k_env_sd = mean(
      k_env_monthly_sd,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

# ============================================================
# 6. PLOT STYLING
# ============================================================

ssp_labels <- c(
  "126" = "SSP1-2.6",
  "245" = "SSP2-4.5",
  "370" = "SSP3-7.0",
  "585" = "SSP5-8.5"
)

window_colors <- c(
  "1880-1900" = "steelblue",
  "2000-2020" = "darkorange",
  "2080-2099" = "firebrick"
)


# ============================================================
# 7. SEASONAL CYCLE WITH SD
# ============================================================

p_seasonal <- ggplot(
  k_env_seasonal,
  aes(
    x = Month,
    y = k_env_mean,
    colour = Window,
    fill = Window
  )
) +
  
  # SD ribbon
  geom_ribbon(
    aes(
      ymin = k_env_mean - k_env_sd,
      ymax = k_env_mean + k_env_sd
    ),
    alpha = 0.15,
    colour = NA
  ) +
  
  # Mean line
  geom_line(
    linewidth = 1
  ) +
  
  geom_point(
    size = 1.5
  ) +
  
  facet_wrap(
    ~ SSP,
    nrow = 2,
    labeller = as_labeller(ssp_labels)
  ) +
  
  scale_colour_manual(
    values = window_colors,
    name = "Time period"
  ) +
  
  scale_fill_manual(
    values = window_colors,
    guide = "none"
  ) +
  
  scale_x_continuous(
    breaks = 1:12,
    labels = month.abb
  ) +
  
  theme_bw() +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  
  labs(
    title = "Mean seasonal cycle of k_env by SSP scenario and time period",
    subtitle = "Spatial mean ± spatial SD",
    x = "Month",
    y = expression(k[env] ~ "[-]")
  )

print(p_seasonal)


# ============================================================
# 8. ANNUAL MEAN TIME SERIES WITH SD
# ============================================================

p_annual <- ggplot(
  k_env_annual_all,
  aes(
    x = Year,
    y = k_env_annual,
    colour = SSP,
    fill = SSP
  )
) +
  
  # SD ribbon
  geom_ribbon(
    aes(
      ymin = k_env_annual - k_env_annual_sd,
      ymax = k_env_annual + k_env_annual_sd
    ),
    alpha = 0.15,
    colour = NA
  ) +
  
  # Mean line
  geom_line(
    linewidth = 1
  ) +
  
  scale_colour_manual(
    values = ssp_colors,
    labels = ssp_labels,
    name = "SSP scenario"
  ) +
  
  scale_fill_manual(
    values = ssp_colors,
    guide = "none"
  ) +
  
  theme_bw() +
  
  labs(
    title = "Annual mean k_env by SSP scenario",
    subtitle = "Spatial mean ± spatial SD",
    x = "Year",
    y = expression(k[env] ~ "[-]")
  )

print(p_annual)
# -----------------------------
# 7. Save
# -----------------------------

ggsave(file.path(out_dir, "k_env_seasonal_by_period_ssp.png"), p_seasonal, width = 12, height = 5)
ggsave(file.path(out_dir, "k_env_annual_by_ssp.png"), p_annual, width = 9, height = 5)






library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(zoo)

# =========================================================
# Paths
# =========================================================

mean_dir <- "monthly_mineralised/mean_2perc_baseline"
std_dir  <- "monthly_mineralised/std_2perc_baseline"   # <-- ADJUST if your std folder lives elsewhere

out_dir <- file.path(mean_dir, "diagnostic_csv_analysis_60N")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ssps <- c("126", "245", "370", "585")

ssp_labels <- c(
  "126" = "SSP1-2.6",
  "245" = "SSP2-4.5",
  "370" = "SSP3-7.0",
  "585" = "SSP5-8.5"
)

ssp_colors <- c(
  "SSP1-2.6" = "blue",
  "SSP2-4.5" = "orange",
  "SSP3-7.0" = "#D73027",
  "SSP5-8.5" = "#7B3294"
)

# Columns we expect to carry an SD counterpart (edit list if needed)
sd_value_cols <- c(
  "mineralised_pg_monthly",
  "total_inorg_pg_monthly",
  "k_T_mean_monthly",
  "k_T_pool_weighted_mean",
  "k_env_pool_weighted_mean",
  "organic_pool_remaining_pg",
  "thawed_N_annual_change_pg"
)

# =========================================================
# Readers
# =========================================================

read_diag_csv <- function(dir, ssp) {
  
  f <- file.path(
    dir,
    paste0("arctic_monthly_diagnostics_w_temp", ssp, "_1850_2100.csv")
  )
  
  if (!file.exists(f)) {
    warning("Missing file: ", f)
    return(NULL)
  }
  
  read_csv(f, show_col_types = FALSE) %>%
    mutate(
      SSP_raw = ssp,
      SSP = ssp_labels[ssp]
    )
}

# --- Mean (existing) diagnostics ---
diag_monthly_mean <- bind_rows(lapply(ssps, function(s) read_diag_csv(mean_dir, s))) %>%
  mutate(
    SSP = factor(SSP, levels = c("SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5"))
  )

# --- Std diagnostics ---
diag_monthly_sd_raw <- bind_rows(lapply(ssps, function(s) read_diag_csv(std_dir, s))) %>%
  mutate(
    SSP = factor(SSP, levels = c("SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5"))
  )

# Rename the value columns in the SD table to *_sd so they don't collide on join
sd_cols_present <- intersect(sd_value_cols, names(diag_monthly_sd_raw))

diag_monthly_sd <- diag_monthly_sd_raw %>%
  select(SSP_raw, SSP, Year, any_of(c("Month", sd_cols_present))) %>%
  rename_with(~ paste0(.x, "_sd"), all_of(sd_cols_present))

# Join key: Year (+ Month if present in both)
join_keys <- intersect(c("SSP_raw", "SSP", "Year", "Month"), names(diag_monthly_mean))
join_keys <- intersect(join_keys, names(diag_monthly_sd))

diag_monthly <- diag_monthly_mean %>%
  left_join(diag_monthly_sd, by = join_keys)

# Check column names
print(names(diag_monthly))
print(head(diag_monthly))

# =========================================================
# Annual summaries (mean + SD)
# =========================================================
# Matching the reference approach exactly: the SD is aggregated with
# the SAME function used for the corresponding mean quantity (sum -> sum,
# mean -> mean). No quadrature/error-propagation anywhere.

diag_annual <- diag_monthly %>%
  group_by(SSP, Year) %>%
  summarise(
    mineralised_pg_yr = sum(mineralised_pg_monthly, na.rm = TRUE),
    total_inorg_pg_yr = sum(total_inorg_pg_monthly, na.rm = TRUE),
    mean_k_t = mean(k_T_mean_monthly, na.rm = TRUE),
    mean_k_t_weighted = mean(k_T_pool_weighted_mean, na.rm = TRUE),
    k_env_weighted_mean = mean(k_env_pool_weighted_mean, na.rm = TRUE),
    organic_pool_remaining_pg = mean(organic_pool_remaining_pg, na.rm = TRUE),
    thawed_N_annual_change_pg = mean(thawed_N_annual_change_pg, na.rm = TRUE),
    
    # ---- SD: same aggregation function as its mean counterpart above ----
    mineralised_pg_yr_sd = sum(mineralised_pg_monthly_sd, na.rm = TRUE),
    total_inorg_pg_yr_sd = sum(total_inorg_pg_monthly_sd, na.rm = TRUE),
    mean_k_t_sd = mean(k_T_mean_monthly_sd, na.rm = TRUE),
    mean_k_t_weighted_sd = mean(k_T_pool_weighted_mean_sd, na.rm = TRUE),
    k_env_weighted_mean_sd = mean(k_env_pool_weighted_mean_sd, na.rm = TRUE),
    organic_pool_remaining_pg_sd = mean(organic_pool_remaining_pg_sd, na.rm = TRUE),
    thawed_N_annual_change_pg_sd = mean(thawed_N_annual_change_pg_sd, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  group_by(SSP) %>%
  arrange(Year, .by_group = TRUE) %>%
  mutate(
    # Cumulative annual thawed-N change
    cumulative_thawed_N_pg = cumsum(replace_na(thawed_N_annual_change_pg, 0)),
    
    cumulative_mineralised_pg = cumsum(mineralised_pg_yr),
    cumulative_inorg_pg = cumsum(total_inorg_pg_yr),
    
    # Cumulative SD: same function (cumsum) as used for the cumulative mean
    cumulative_mineralised_pg_sd = cumsum(mineralised_pg_yr_sd),
    cumulative_inorg_pg_sd = cumsum(total_inorg_pg_yr_sd),
    
    mineralised_20yr = zoo::rollapply(mineralised_pg_yr, 20, mean, align = "center", fill = NA),
    inorg_20yr = zoo::rollapply(total_inorg_pg_yr, 20, mean, align = "center", fill = NA),
    
    mineralised_20yr_sd = zoo::rollapply(mineralised_pg_yr_sd, 20, mean, align = "center", fill = NA),
    inorg_20yr_sd = zoo::rollapply(total_inorg_pg_yr_sd, 20, mean, align = "center", fill = NA)
  ) %>%
  ungroup()

# =========================================================
# Anomaly relative to present-day (2000-2020) reference period
# =========================================================
# Single-stage: subtract the 2000-2020 reference-period mean directly
# from the raw series. Applied identically to Mean and SD (SD is just
# walked through the same subtraction -- no quadrature, no baseline
# re-centering step).

compute_anomaly_pair <- function(df, value_col, sd_col,
                                 ref_range = c(2000, 2020)) {
  
  ref <- df %>% filter(Year >= ref_range[1], Year <= ref_range[2])
  
  ref_mean <- mean(ref[[value_col]], na.rm = TRUE)
  ref_sd   <- mean(ref[[sd_col]], na.rm = TRUE)
  
  df %>%
    mutate(
      "{value_col}_anom" := .data[[value_col]] - ref_mean,
      "{sd_col}_anom"     := .data[[sd_col]] - ref_sd
    )
}

diag_annual <- diag_annual %>%
  group_by(SSP) %>%
  group_modify(~ compute_anomaly_pair(.x, "mineralised_pg_yr", "mineralised_pg_yr_sd")) %>%
  group_modify(~ compute_anomaly_pair(.x, "total_inorg_pg_yr", "total_inorg_pg_yr_sd")) %>%
  group_modify(~ compute_anomaly_pair(.x, "cumulative_mineralised_pg", "cumulative_mineralised_pg_sd")) %>%
  group_modify(~ compute_anomaly_pair(.x, "cumulative_inorg_pg", "cumulative_inorg_pg_sd")) %>%
  ungroup() %>%
  rename(
    annual_mineralised_anom_pg     = mineralised_pg_yr_anom,
    annual_mineralised_anom_pg_sd  = mineralised_pg_yr_sd_anom,
    annual_inorg_anom_pg           = total_inorg_pg_yr_anom,
    annual_inorg_anom_pg_sd        = total_inorg_pg_yr_sd_anom,
    cumulative_mineralised_anom_pg    = cumulative_mineralised_pg_anom,
    cumulative_mineralised_anom_pg_sd = cumulative_mineralised_pg_sd_anom,
    cumulative_inorg_anom_pg          = cumulative_inorg_pg_anom,
    cumulative_inorg_anom_pg_sd       = cumulative_inorg_pg_sd_anom
  )

# Rolling 20-year mean of the anomalies (plain mean, same as the reference)
diag_annual <- diag_annual %>%
  group_by(SSP) %>%
  arrange(Year, .by_group = TRUE) %>%
  mutate(
    cumulative_mineralised_anom_20yr =
      zoo::rollapply(cumulative_mineralised_anom_pg, 20, mean, align = "center", fill = NA),
    cumulative_mineralised_anom_20yr_sd =
      zoo::rollapply(cumulative_mineralised_anom_pg_sd, 20, mean, align = "center", fill = NA),
    
    cumulative_inorg_anom_20yr =
      zoo::rollapply(cumulative_inorg_anom_pg, 20, mean, align = "center", fill = NA),
    cumulative_inorg_anom_20yr_sd =
      zoo::rollapply(cumulative_inorg_anom_pg_sd, 20, mean, align = "center", fill = NA)
  ) %>%
  ungroup()

# =========================================================
# Plot function (mean line + optional SD ribbon)
# =========================================================

plot_diag <- function(df, yvar, ylab, title, filename, sd_var = NULL) {
  
  p <- ggplot(df, aes(x = Year, y = .data[[yvar]], color = SSP))
  
  if (!is.null(sd_var) && sd_var %in% names(df)) {
    p <- p +
      geom_ribbon(
        aes(
          ymin = .data[[yvar]] - .data[[sd_var]],
          ymax = .data[[yvar]] + .data[[sd_var]],
          fill = SSP
        ),
        color = NA,
        alpha = 0.18,
        na.rm = TRUE
      ) +
      scale_fill_manual(values = ssp_colors, drop = FALSE, guide = "none")
  }
  
  p <- p +
    geom_line(linewidth = 0.8, na.rm = TRUE) +
    geom_vline(xintercept = 2015, linetype = "dashed", color = "black") +
    scale_color_manual(values = ssp_colors, drop = FALSE) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12)
    ) +
    labs(
      x = "Year",
      y = ylab,
      color = "SSP scenario"
      # title = title
    )
  
  ggsave(
    file.path(out_dir, filename),
    p,
    width = 8,
    height = 5,
    dpi = 300
  )
  
  return(p)
}

# =========================================================
# Make plots (with SD ribbons where available)
# =========================================================

p1 <- plot_diag(
  diag_annual,
  "annual_mineralised_anom_pg",
  "Annual mineralised N [Pg N yr\u207b\u00b9]",
  "Mineralised N, annual",
  "annual_mineralised_N.png",
  sd_var = "annual_mineralised_anom_pg_sd"
)

p2 <- plot_diag(
  diag_annual,
  "cumulative_mineralised_anom_pg",
  "Cumulative mineralised N anomaly [Pg N]",
  "Cumulative mineralised N anomaly (rel. 2000-2020)",
  "cumulative_mineralised_N_anomaly.png",
  sd_var = "cumulative_mineralised_anom_pg_sd"
)

p3 <- plot_diag(
  diag_annual,
  "annual_inorg_anom_pg",
  "Annual total inorganic N [Pg N yr\u207b\u00b9]",
  "Total inorganic N, annual",
  "annual_total_inorg_N.png",
  sd_var = "annual_inorg_anom_pg_sd"
)

p4 <- plot_diag(
  diag_annual,
  "cumulative_inorg_anom_pg",
  "Cumulative inorganic N anomaly [Pg N]",
  "Cumulative inorganic N anomaly (rel. 2000-2020)",
  "cumulative_inorg_N_anomaly.png",
  sd_var = "cumulative_inorg_anom_pg_sd"
)


ALD <- ALD + theme(legend.position = "none")

combined_plot <- ALD | plot_total_thawed +
  plot_annotation(tag_levels = "a", tag_suffix = ")")&
  
  theme(legend.position = "right") &
  
  plot_layout(guides = "collect")



combined_plot
ggsave("Fig1.png", combined_plot,
       width = 15,
       height = 6,
       units = "cm", 
       #scale = 1.4,
       dpi = 500,  #70
       #device = cairo_pdf
)






### mineralisation rates comparison between time periods


library(terra)
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)

# ============================================================
# Settings
# ============================================================

base_dir <- "monthly_mineralised/mean_2perc_baseline"
out_dir  <- file.path(base_dir, "period_summary_inorganic_N_gridded")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ssps <- c("126", "245", "370", "585")

ssp_labels <- c(
  "126" = "SSP1-2.6",
  "245" = "SSP2-4.5",
  "370" = "SSP3-7.0",
  "585" = "SSP5-8.5"
)

ssp_colors <- c(
  "SSP1-2.6" = "blue",
  "SSP2-4.5" = "orange",
  "SSP3-7.0" = "#D73027",
  "SSP5-8.5" = "#7B3294"
)

periods <- list(
  "Present-day (2000-2020)"     = 2000:2020,
  "Mid-century (2040-2060)"     = 2040:2060,
  "End-of-century (2080-2099)"  = 2080:2099
)

period_levels <- names(periods)

# ============================================================
# 1. Function: read one year's monthly total_inorg raster,
#    sum to an annual kg N m-2 yr-1 layer
# ============================================================

read_annual_total_inorg <- function(ssp, yr) {
  
  f <- file.path(
    base_dir,
    paste0("yearly_nc_", ssp),
    paste0(
      "arctic_total_inorg_N_monthly_w_temp_sm_",
      ssp, "_", yr, ".nc"
    )
  )
  
  if (!file.exists(f)) {
    warning("Missing file: ", f)
    return(NULL)
  }
  
  r <- rast(f)   # 12 monthly layers, kg N m-2 month total
  
  if (nlyr(r) != 12) {
    warning("Expected 12 layers, found ", nlyr(r), " in ", f)
  }
  
  annual_kg_m2 <- app(r, sum, na.rm = TRUE)   # kg N m-2 yr-1
  names(annual_kg_m2) <- as.character(yr)
  
  annual_kg_m2
}

# ============================================================
# 2. For each SSP and period, build the mean annual raster
#    across the years in that period, convert to g N m-2 yr-1
# ============================================================

build_period_raster <- function(ssp, period_years) {
  
  yearly_rasters <- lapply(period_years, function(yr) {
    read_annual_total_inorg(ssp, yr)
  })
  
  yearly_rasters <- Filter(Negate(is.null), yearly_rasters)
  
  if (length(yearly_rasters) == 0) {
    warning("No years found for SSP ", ssp, " in this period.")
    return(NULL)
  }
  
  yearly_stack <- rast(yearly_rasters)
  
  # Mean across years within the period, in kg N m-2 yr-1
  mean_annual_kg_m2 <- app(yearly_stack, mean, na.rm = TRUE)
  
  # Convert kg N m-2 yr-1 -> g N m-2 yr-1
  mean_annual_g_m2 <- mean_annual_kg_m2 * 1000
  
  mean_annual_g_m2
}

# ============================================================
# 3. Loop over SSPs and periods: build rasters + scalar summaries
# ============================================================

period_raster_list <- list()
period_summary_rows <- list()

for (ssp in ssps) {
  
  for (period_name in period_levels) {
    
    yrs <- periods[[period_name]]
    
    cat("Processing SSP", ssp, "-", period_name, "\n")
    
    r_period <- build_period_raster(ssp, yrs)
    
    if (is.null(r_period)) next
    
    key <- paste(ssp, period_name, sep = "_")
    period_raster_list[[key]] <- r_period
    
    # Area-weighted spatial mean, min, max -> scalar summary
    area_rast <- cellSize(r_period, unit = "m")
    
    spatial_mean <- global(
      r_period, "mean", weights = area_rast, na.rm = TRUE
    )[1, 1]
    
    spatial_vals <- global(r_period, "range", na.rm = TRUE)
    
    period_summary_rows[[key]] <- data.frame(
      SSP_raw = ssp,
      SSP = ssp_labels[[ssp]],
      Period = period_name,
      mean_gN_m2_yr = spatial_mean,
      min_gN_m2_yr  = spatial_vals[1, 1],
      max_gN_m2_yr  = spatial_vals[1, 2]
    )
  }
}

period_summary_area <- bind_rows(period_summary_rows) %>%
  mutate(
    SSP = factor(SSP, levels = unname(ssp_labels)),
    Period = factor(Period, levels = period_levels)
  ) %>%
  arrange(SSP, Period)

print(period_summary_area)

write_csv(
  period_summary_area,
  file.path(out_dir, "inorganic_N_flux_period_summary_gN_m2_yr_gridded.csv")
)

# ============================================================
# 4. Save the period-mean rasters (for spatial maps later)
# ============================================================

for (key in names(period_raster_list)) {
  
  writeRaster(
    period_raster_list[[key]],
    file.path(out_dir, paste0("mean_total_inorg_gN_m2_yr_", key, ".nc")),
    overwrite = TRUE,
    filetype = "NetCDF"
  )
}

# ============================================================
# 5. Bar plot: scalar summary by period and SSP
# ============================================================

p_period_bar <- ggplot(
  period_summary_area,
  aes(x = Period, y = mean_gN_m2_yr, fill = SSP)
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = ssp_colors, drop = FALSE) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 20, hjust = 1)
  ) +
  labs(
    title = "Mean total inorganic N flux by period and scenario",
    subtitle = "Area-weighted spatial mean across 60-90°N land area",
    x = NULL,
    y = expression("Total inorganic N flux (g N m"^{-2}*" yr"^{-1}*")"),
    fill = "SSP scenario"
  )

print(p_period_bar)

ggsave(
  file.path(out_dir, "inorganic_N_flux_period_summary_gN_m2_yr_barplot.png"),
  p_period_bar,
  width = 9,
  height = 5.5,
  dpi = 300
)

# ============================================================
# 6. Example: plot the spatial map for one SSP/period
# ============================================================

example_key <- paste("585", "End-of-century (2080-2099)", sep = "_")

if (example_key %in% names(period_raster_list)) {
  plot(
    period_raster_list[[example_key]],
    main = "Mean total inorganic N flux, SSP5-8.5, 2080-2099 (g N m-2 yr-1)"
  )
}

cat("Done. Outputs saved in:\n", out_dir, "\n")




