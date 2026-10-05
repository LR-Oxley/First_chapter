# ============================================================
# Pan-Arctic time series from the .nc outputs of the
# mineralisation model, with uncertainty from perturbed-input runs:
#   mean  run: model driven by the mean thawed N and ALD
#   plus  run: model driven by mean + 1 SD
#   minus run: model driven by mean - 1 SD
# Every output is processed identically in all three runs (pan-Arctic
# totals, cumulative sums, anomaly, running mean); the uncertainty is
#   sd(t) = |output_plus(t) - output_minus(t)| / 2
#
# A) Yearly files (one layer per year, 1850-2099), in each run folder:
#   arctic_rapid_inorg_N_yearly_<ssp>_<run>_1850_2099_w_temp.nc            flux  (kg N m-2 yr-1)
#   arctic_new_thawed_total_N_yearly_<ssp>_<run>_1850_2099_w_temp.nc       flux  (kg N m-2 yr-1)
#   arctic_new_thawed_organic_N_yearly_<ssp>_<run>_1850_2099_w_temp.nc     flux  (kg N m-2 yr-1)
#   arctic_organic_N_pool_remaining_yearly_<ssp>_<run>_1850_2099_w_temp.nc stock (kg N m-2)
#   arctic_permafrost_thawed_total_N_<ssp>_<run>_1850_2099_w_temp.nc       stock (kg N m-2)
#
# B) Monthly files (12 layers Jan-Dec, one file per year), in
#    <run folder>/yearly_nc_<ssp>/:
#   arctic_mineralised_N_monthly_w_temp_sm_<ssp>_<year>.nc   flux (kg N m-2 month-1)
#   arctic_total_inorg_N_monthly_w_temp_sm_<ssp>_<year>.nc   flux (kg N m-2 month-1)
#   arctic_k_env_monthly_w_temp_sm_<ssp>_<year>.nc           rate modifier (-)
#   arctic_k_T_temperature_only_monthly_<ssp>_<year>.nc      rate modifier (-)
#
#   fluxes:          12 months summed -> kg N m-2 yr-1 -> pan-Arctic Pg N yr-1
#   rate modifiers:  area-weighted pan-Arctic mean per month, then averaged
#                    over the months in k_months (default: all 12)
#
# Steps:
#   1. settings
#   2. read yearly files  -> pan-Arctic value per year and run (cached)
#   3. read monthly files -> pan-Arctic value per year and run (cached, slow)
#   4. cumulative sums of the fluxes (per run)
#   5. anomaly (optional) + 20-year running mean (per run),
#      then uncertainty = |plus - minus| / 2
#   6. plot function
#   7. one figure per variable
# ============================================================

library(terra)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)


# ------------------------------------------------------------
# 1. Settings
# ------------------------------------------------------------
# run folders; set plus/minus to NULL to plot without uncertainty
run_dirs <- list(
  mean  = "monthly_mineralised/mean_2perc_baseline",
  plus  = "monthly_mineralised/plus_2perc_baseline",
  minus = "monthly_mineralised/minus_2perc_baseline"
)
out_dir <- file.path(run_dirs$mean, "diagnostic_nc_analysis_60N")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ext_sub    <- ext(-179.95, 179.95, 60, 90)
start_year <- 1850          # used only if the yearly files carry no time information

ssps <- c("126", "245", "370", "585")
ssp_labels <- c("126" = "SSP1-2.6", "245" = "SSP2-4.5",
                "370" = "SSP3-7.0", "585" = "SSP5-8.5")
ssp_colors <- c("SSP1-2.6" = "blue", "SSP2-4.5" = "orange",
                "SSP3-7.0" = "#D73027", "SSP5-8.5" = "#7B3294")

# TRUE  = plot the anomaly relative to ref_start-ref_end
# FALSE = plot the absolute values
plot_anomaly <- TRUE
ref_start <- 2000
ref_end   <- 2020

# TRUE = reuse the cached pan-Arctic CSVs instead of re-reading the rasters
# (set to FALSE after re-running the model or changing k_months)
read_yearly        <- TRUE
read_monthly       <- FALSE      # slow: 250 files per variable, SSP and run
use_cache          <- TRUE
cache_file         <- file.path(out_dir, "pan_arctic_yearly_from_nc_runs.csv")
cache_file_monthly <- file.path(out_dir, "pan_arctic_from_monthly_nc_runs.csv")

# A) yearly files: file name (first %s = SSP code, second %s = run), kind, label
yearly_vars <- tribble(
  ~key,                 ~file_fmt,                                                             ~kind,   ~label,
  "rapid_inorg",        "arctic_rapid_inorg_N_yearly_%s_%s_1850_2099_w_temp.nc",              "flux",  "Pre-thaw inorganic N",
  "new_thawed_total",   "arctic_new_thawed_total_N_yearly_%s_%s_1850_2099_w_temp.nc",         "flux",  "Newly thawed total N",
  "new_thawed_organic", "arctic_new_thawed_organic_N_yearly_%s_%s_1850_2099_w_temp.nc",       "flux",  "Newly thawed organic N",
  "organic_pool",       "arctic_organic_N_pool_remaining_yearly_%s_%s_1850_2099_w_temp.nc",   "stock", "Organic N pool",
  "thawed_total",       "arctic_permafrost_thawed_total_N_%s_%s_1850_2099_w_temp.nc",         "stock", "Thawed permafrost N"
)

# B) monthly files: name part of arctic_<name>_<ssp>_<year>.nc
#    kind "flux" = summed over months; kind "rate" = averaged over k_months
monthly_vars <- tribble(
  ~key,          ~name,                              ~kind,  ~label,
  "mineralised", "mineralised_N_monthly_w_temp_sm",  "flux", "Mineralised N",
  "total_inorg", "total_inorg_N_monthly_w_temp_sm",  "flux", "Total inorganic N",
  "k_env",       "k_env_monthly_w_temp_sm",          "rate", "Environmental modifier k_env",
  "k_T",         "k_T_temperature_only_monthly",     "rate", "Temperature modifier k_T"
)
monthly_file_fmt <- "arctic_%s_%s_%d.nc"     # name, SSP code, year
monthly_years    <- 1850:2099

# months used to average the rate modifiers (1:12 = whole year, 6:8 = summer)
k_months <- 1:12

# plotted as 20-year running mean only (no yearly noise)
smooth_only <- c("rapid_inorg")    # e.g. c("rapid_inorg", "k_env", "k_T")

# also plot the cumulative sum of each flux (Pg N)
plot_cumulative <- TRUE

# figure size
fig_width_cm  <- 8.5
fig_height_cm <- 6
fig_dpi       <- 500


# ------------------------------------------------------------
# 2. Yearly files -> pan-Arctic totals, per run
#
# kg N m-2 (yr-1) x cell area (m2) = kg N (yr-1); / 1e12 = Pg N (yr-1)
# ------------------------------------------------------------
read_pan_arctic <- function(f) {
  
  r <- crop(rast(f), ext_sub)
  area <- cellSize(r[[1]], unit = "m")
  
  pg <- global(r * area, "sum", na.rm = TRUE)[, 1] / 1e12
  
  yrs <- suppressWarnings(as.integer(format(time(r), "%Y")))
  if (length(yrs) != nlyr(r) || any(is.na(yrs))) {
    yrs <- start_year + seq_len(nlyr(r)) - 1
  }
  
  tibble(Year = yrs, value = pg)
}

pan_yearly <- NULL

if (read_yearly) {
  
  if (use_cache && file.exists(cache_file)) {
    
    cat("Reading cached yearly pan-Arctic totals:", cache_file, "\n")
    pan_yearly <- read_csv(cache_file, show_col_types = FALSE)
    
  } else {
    
    rows <- list()
    
    for (run in names(run_dirs)) {
      if (is.null(run_dirs[[run]])) next
      
      for (s in ssps) {
        for (k in seq_len(nrow(yearly_vars))) {
          
          f <- file.path(run_dirs[[run]], sprintf(yearly_vars$file_fmt[k], s, run))
          if (!file.exists(f)) {
            warning("Missing file: ", f)
            next
          }
          cat("Reading", yearly_vars$key[k], "SSP", s, "run", run, "\n")
          
          tab <- read_pan_arctic(f)
          tab$SSP <- ssp_labels[[s]]
          tab$var <- yearly_vars$key[k]
          tab$run <- run
          rows[[length(rows) + 1]] <- tab
        }
      }
    }
    
    pan_yearly <- bind_rows(rows)
    write_csv(pan_yearly, cache_file)
  }
}


# ------------------------------------------------------------
# 3. Monthly files -> one pan-Arctic value per year, per run
#
# how = "sum":  sum of the 12 months per cell (kg N m-2 yr-1),
#               x cell area, summed over the domain -> Pg N yr-1
# how = "mean": area-weighted domain mean for each month in k_months,
#               then averaged over those months (dimensionless)
# ------------------------------------------------------------
read_pan_arctic_monthly <- function(folder, name, s, how) {
  
  out <- vector("list", length(monthly_years))
  area <- NULL
  n_missing <- 0
  
  for (i in seq_along(monthly_years)) {
    
    yr <- monthly_years[i]
    f  <- file.path(folder, sprintf(monthly_file_fmt, name, s, yr))
    
    if (!file.exists(f)) {
      n_missing <- n_missing + 1
      out[[i]] <- tibble(Year = yr, value = NA_real_)
      next
    }
    
    r <- crop(rast(f), ext_sub)
    if (is.null(area)) area <- cellSize(r[[1]], unit = "m")
    
    if (how == "sum") {
      annual <- sum(r, na.rm = TRUE)
      v <- global(annual * area, "sum", na.rm = TRUE)[1, 1] / 1e12
    } else {
      layers <- intersect(k_months, seq_len(nlyr(r)))   # layers are Jan-Dec
      v <- mean(global(r[[layers]], "mean", weights = area, na.rm = TRUE)[, 1], na.rm = TRUE)
    }
    
    out[[i]] <- tibble(Year = yr, value = v)
  }
  
  if (n_missing > 0) {
    warning(n_missing, " monthly files missing for ", name, " SSP ", s, " in ", folder)
  }
  bind_rows(out)
}


pan_monthly <- NULL

if (read_monthly) {
  
  if (use_cache && file.exists(cache_file_monthly)) {
    
    cat("Reading cached monthly pan-Arctic totals:", cache_file_monthly, "\n")
    pan_monthly <- read_csv(cache_file_monthly, show_col_types = FALSE)
    
  } else {
    
    rows <- list()
    
    for (run in names(run_dirs)) {
      if (is.null(run_dirs[[run]])) next
      
      for (s in ssps) {
        
        folder <- file.path(run_dirs[[run]], paste0("yearly_nc_", s))
        if (!dir.exists(folder)) {
          warning("Missing folder: ", folder)
          next
        }
        
        for (k in seq_len(nrow(monthly_vars))) {
          
          how <- if (monthly_vars$kind[k] == "flux") "sum" else "mean"
          cat("Reading monthly", monthly_vars$key[k], "SSP", s, "run", run, "\n")
          
          tab <- read_pan_arctic_monthly(folder, monthly_vars$name[k], s, how)
          tab$SSP <- ssp_labels[[s]]
          tab$var <- monthly_vars$key[k]
          tab$run <- run
          rows[[length(rows) + 1]] <- tab
        }
      }
    }
    
    pan_monthly <- bind_rows(rows)
    write_csv(pan_monthly, cache_file_monthly)
  }
}

pan <- bind_rows(pan_yearly)#, pan_monthly)

all_vars <- bind_rows(
  yearly_vars  %>% select(key, kind, label),
  #monthly_vars %>% select(key, kind, label)
)

annual <- pan %>%
  left_join(all_vars %>% rename(var = key), by = "var")


# ------------------------------------------------------------
# 4. Cumulative sums of the fluxes (Pg N), separately per run
# ------------------------------------------------------------
if (plot_cumulative) {
  cumulative <- annual %>%
    filter(kind == "flux") %>%
    group_by(run, SSP, var) %>%
    arrange(Year, .by_group = TRUE) %>%
    mutate(value = cumsum(replace_na(value, 0)),
           var   = paste0(var, "_cumulative"),
           label = paste("", tolower(label)),
           kind  = "cumulative") %>%
    ungroup()
  annual <- bind_rows(annual, cumulative)
}


# ------------------------------------------------------------
# Thawed organic N pool: 2080-2099, present day and increase,
# with uncertainty = |plus - minus| / 2
# ------------------------------------------------------------
organic <- annual %>%
  filter(var == "organic_pool") %>%
  group_by(SSP, run) %>%
  summarise(
    pool_2000_2020 = mean(value[Year >= 2000 & Year <= 2020], na.rm = TRUE),
    pool_2080_2099 = mean(value[Year >= 2080 & Year <= 2099], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(increase = pool_2080_2099 - pool_2000_2020) %>%
  pivot_longer(c(pool_2000_2020, pool_2080_2099, increase),
               names_to = "quantity", values_to = "value") %>%
  pivot_wider(names_from = run, values_from = value) %>%
  mutate(sd   = abs(plus - minus) / 2,
         text = sprintf("%.1f\\,$\\pm$\\,%.1f", mean, sd))

print(organic %>%
        select(quantity, SSP, text) %>%
        pivot_wider(names_from = SSP, values_from = text),
      width = Inf)
# ------------------------------------------------------------
# 5. Anomaly (optional) and 20-year running mean, per run;
#    then uncertainty = |plus - minus| / 2
# ------------------------------------------------------------
per_run <- annual %>%
  group_by(run, SSP, var) %>%
  arrange(Year, .by_group = TRUE) %>%
  mutate(
    ref_mean = mean(value[Year >= ref_start & Year <= ref_end], na.rm = TRUE),
    y        = if (plot_anomaly) value - ref_mean else value,
    y_20     = zoo::rollapply(y, 20, mean, align = "center", fill = NA)
  ) %>%
  ungroup()

series <- per_run %>%
  select(run, SSP, var, kind, label, Year, y, y_20) %>%
  pivot_wider(names_from = run, values_from = c(y, y_20))

# runs that were not found -> no uncertainty
for (col in c("y_mean", "y_plus", "y_minus", "y_20_mean", "y_20_plus", "y_20_minus")) {
  if (!col %in% names(series)) series[[col]] <- NA_real_
}

series <- series %>%
  mutate(y_sd    = abs(y_plus    - y_minus)    / 2,
         y_20_sd = abs(y_20_plus - y_20_minus) / 2) %>%
  rename(y = y_mean, y_20 = y_20_mean) %>%
  mutate(SSP = factor(SSP, levels = unname(ssp_labels)))

write_csv(series, file.path(out_dir,
                            paste0("pan_arctic_series_", ifelse(plot_anomaly, "anomaly", "absolute"), ".csv")))


# ------------------------------------------------------------
# 6. Plot function
# Historical (1850-2014) = black + grey ribbon
# Future (2015-2099)      = SSP colours + coloured ribbons
# ------------------------------------------------------------
plot_diag <- function(df,
                      yvar,
                      ylab,
                      title,
                      filename,
                      sd_var = NULL,
                      running_var = NULL,
                      running_sd_var = NULL,
                      show_raw = TRUE) {
  
  df <- df %>% arrange(SSP, Year)
  
  df_hist   <- df %>% filter(Year <= 2014)
  df_future <- df %>% filter(Year >= 2015)
  
  # SSP used to represent the shared historical period
  historical_ssp <- unique(df$SSP)[1]
  df_hist <- df_hist %>% filter(SSP == historical_ssp)
  
  p <- ggplot(df, aes(x = Year, group = SSP))
  
  # raw values
  if (show_raw) {
    p <- p +
      geom_line(data = df_hist, aes(y = .data[[yvar]]),
                color = "black", linewidth = 0.2, alpha = 0.45, na.rm = TRUE) +
      geom_line(data = df_future, aes(y = .data[[yvar]], color = SSP),
                linewidth = 0.2, alpha = 0.45, na.rm = TRUE)
  }
  
  # raw SD ribbon (optional)
  if (show_raw && !is.null(sd_var) && sd_var %in% names(df)) {
    p <- p +
      geom_ribbon(data = df_hist,
                  aes(ymin = .data[[yvar]] - .data[[sd_var]],
                      ymax = .data[[yvar]] + .data[[sd_var]]),
                  fill = "grey", color = NA, alpha = 0.03) +
      geom_ribbon(data = df_future,
                  aes(ymin = .data[[yvar]] - .data[[sd_var]],
                      ymax = .data[[yvar]] + .data[[sd_var]],
                      fill = SSP),
                  color = NA, alpha = 0.03)
  }
  
  # running mean
  if (!is.null(running_var) && running_var %in% names(df)) {
    p <- p +
      geom_line(data = df_hist, aes(y = .data[[running_var]]),
                color = "black", linewidth = 0.2, na.rm = TRUE) +
      geom_line(data = df_future, aes(y = .data[[running_var]], color = SSP),
                linewidth = 0.2, na.rm = TRUE)
  }
  
  # running-mean SD ribbon
  if (!is.null(running_var) && !is.null(running_sd_var) &&
      running_var %in% names(df) && running_sd_var %in% names(df)) {
    p <- p +
      geom_ribbon(data = df_hist,
                  aes(ymin = .data[[running_var]] - .data[[running_sd_var]],
                      ymax = .data[[running_var]] + .data[[running_sd_var]]),
                  fill = "grey", color = NA, alpha = 0.5) +
      geom_ribbon(data = df_future,
                  aes(ymin = .data[[running_var]] - .data[[running_sd_var]],
                      ymax = .data[[running_var]] + .data[[running_sd_var]],
                      fill = SSP),
                  color = NA, alpha = 0.2)
  }
  
  p +
    geom_vline(xintercept = 2015, linetype = "dashed", color = "black") +
    scale_color_manual(values = ssp_colors, drop = FALSE) +
    scale_fill_manual(values = ssp_colors, drop = FALSE, guide = "none") +
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.text    = element_text(size = 4),
          axis.title   = element_text(size = 4),
          legend.text  = element_text(size = 4),
          legend.title = element_text(size = 4)) +
    labs(x = "Year", y = ylab, color = "", title = title)
}


# ------------------------------------------------------------
# 7. One figure per variable
# ------------------------------------------------------------
plots <- list()

for (v in unique(series$var)) {
  
  d    <- series %>% filter(var == v)
  kind <- d$kind[1]
  lab  <- d$label[1]
  
  ylab <- switch(kind,
                 flux = bquote(.(lab) ~ "[Pg N" ~ yr^-1 * "]"),
                 rate = bquote(.(lab) ~ "[-]"),
                 bquote(.(lab) ~ "[Pg N]"))            # stock, cumulative
  
  tag <- paste0(v, "_", ifelse(plot_anomaly, "anomaly", "absolute"))
  
  # no ribbons if the plus/minus runs are missing
  has_sd <- any(is.finite(d$y_sd))
  
  p <- plot_diag(d,
                 yvar           = "y",
                 ylab           = ylab,
                 title          = "",
                 filename       = paste0(tag, ".png"),
                 sd_var         = if (has_sd) "y_sd" else NULL,
                 running_var    = "y_20",
                 running_sd_var = if (has_sd) "y_20_sd" else NULL,
                 show_raw       = !(v %in% smooth_only))
  
  if (!plot_anomaly) p <- p + geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3)
  
  ggsave(file.path(out_dir, paste0(tag, ".png")), p,
         width = fig_width_cm, height = fig_height_cm, units = "cm", dpi = fig_dpi)
  ggsave(file.path(out_dir, paste0(tag, ".pdf")), p,
         width = fig_width_cm, height = fig_height_cm, units = "cm", device = cairo_pdf)
  
  plots[[v]] <- p
  print(p)
}

cat("\nDone. Figures and tables in", out_dir, "\n")











# ============================================================
# Pan-Arctic time series from the diagnostics CSVs of the
# mineralisation model (fast: no rasters are read).
#
# Uncertainty from perturbed-input runs:
#   mean  run: model driven by the mean thawed N and ALD
#   plus  run: model driven by mean + 1 SD (thawed N and ALD)
#   minus run: model driven by mean - 1 SD (thawed N and ALD)
# Every output is processed identically in all three runs (annual sums,
# cumulative sums, anomaly, running mean); the uncertainty is
#   sd(t) = |output_plus(t) - output_minus(t)| / 2
#
# Input, per SSP and run:
#   <run dir>/arctic_monthly_diagnostics_w_temp<ssp>_1850_2100.csv
#   one row per month, pan-Arctic values already computed by the model:
#     mineralised_pg_monthly, total_inorg_pg_monthly,
#     inorg_rapid_available              Pg N month-1
#     organic_pool_remaining_pg,
#     thawed_N_annual_change_pg          yearly values, repeated each month
#     k_T_*, k_env_*                     area-weighted mean modifiers (-)
#
# Steps:
#   1. settings
#   2. read the CSVs (mean, plus, minus run)
#   3. monthly -> annual (sum / yearly value / mean over k_months)
#   4. cumulative sums of the fluxes (per run)
#   5. anomaly (optional) + 20-year running mean (per run),
#      then uncertainty = |plus - minus| / 2
#   6. plot function
#   7. one figure per variable
#   8. seasonal cycle per period (mean run)
#   9. pan-Arctic mean per m2 (annual and cumulative mineralised N,
#      pre-thaw inorganic N, total inorganic N)
# ============================================================

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)


# ------------------------------------------------------------
# 1. Settings
# ------------------------------------------------------------
# run folders; set plus/minus to NULL to plot without uncertainty
run_dirs <- list(
  mean  = "monthly_mineralised/mean_2perc_baseline",
  plus  = "monthly_mineralised/plus_2perc_baseline",
  minus = "monthly_mineralised/minus_2perc_baseline"
)
out_dir <- file.path(run_dirs$mean, "diagnostic_csv_analysis_60N")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

csv_fmt <- "arctic_monthly_diagnostics_w_temp%s_1850_2100.csv"   # %s = SSP code

ssps <- c("126", "245", "370", "585")
ssp_labels <- c("126" = "SSP1-2.6", "245" = "SSP2-4.5",
                "370" = "SSP3-7.0", "585" = "SSP5-8.5")
ssp_colors <- c("SSP1-2.6" = "blue", "SSP2-4.5" = "orange",
                "SSP3-7.0" = "#D73027", "SSP5-8.5" = "#7B3294")

last_year <- 2099

# TRUE  = plot the anomaly relative to ref_start-ref_end
# FALSE = plot the absolute values
plot_anomaly <- TRUE
ref_start <- 2000
ref_end   <- 2020

# variables to analyse
#   col:    column in the CSV
#   kind:   flux (Pg N yr-1), stock (Pg N) or rate (dimensionless)
#   how:    sum     = sum of the 12 months
#           first   = yearly value (repeated in every month of the CSV)
#           kmonths = mean over the months in k_months
csv_vars <- tribble(
  ~col,                         ~kind,   ~how,      ~label,
  "mineralised_pg_monthly",     "flux",  "sum",     "Mineralised N",
  "total_inorg_pg_monthly",     "flux",  "sum",     "Total inorganic N",
  "inorg_rapid_available",      "flux",  "sum",     "Pre-thaw inorganic N",
  "thawed_N_annual_change_pg",  "flux",  "first",   "Annual change in thawed N",
  "organic_pool_remaining_pg",  "stock", "first",   "Organic N pool",
  "k_env_pool_weighted_mean",   "rate",  "kmonths", "Environmental modifier k_env",
  "k_T_pool_weighted_mean",     "rate",  "kmonths", "Temperature modifier k_T",
  "k_T_mean_monthly",           "rate",  "kmonths", "k_T (depth-weighted mean)"
)

# months used to average the rate modifiers (1:12 = whole year, 6:8 = summer)
k_months <- 1:12

# plotted as 20-year running mean only (no yearly noise)
smooth_only <- c("inorg_rapid_available")

# also plot the cumulative sum of each flux (Pg N)
plot_cumulative <- TRUE

# TRUE = annual fluxes plotted in Tg N yr-1 (cumulative amounts and stocks stay in Pg N)
flux_in_Tg <- TRUE

# figure size
fig_width_cm  <- 8.5
fig_height_cm <- 6
fig_dpi       <- 500

# seasonal cycle figures (Step 8): mean seasonal cycle per period and SSP
#   scale: multiply the CSV values by this (1000 = Pg -> Tg)
seasonal_vars <- list(
  k_env_depth_unweighted_mean = list(
    ylab  = expression(k[env] ~ "[-]"),
    scale = 1),
  mineralised_pg_monthly = list(
    ylab  = expression("Mineralised N [Tg N month"^-1 * "]"),
    scale = 1000)
)
seasonal_periods <- tribble(
  ~period,                      ~start, ~end,
  "Preindustrial (1850-1900)",  1850,   1900,
  "Present day (2000-2020)",    2000,   2020,
  "End of century (2080-2099)", 2080,   2099
)
period_colours <- c("steelblue", "#F39C12", "firebrick")   # same order as seasonal_periods
seasonal_width_cm  <- 12
seasonal_height_cm <- 7

# Pan-Arctic MEAN per m2 (Step 9): total / domain area
#   annual fluxes  -> g N m-2 yr-1
#   cumulative     -> g N m-2
per_area_vars <- c("mineralised_pg_monthly",
                   "inorg_rapid_available",
                   "total_inorg_pg_monthly",
                   "mineralised_pg_monthly_cumulative",
                   "total_inorg_pg_monthly_cumulative")
# domain = model land mask (cells with data in the thawed-N input file);
# or set domain_area_m2 to a number to skip reading the raster
domain_file    <- file.path("total_thawed_extended", "arctic_total_thawed_126_60N_mean.nc")
domain_area_m2 <- NULL


# ------------------------------------------------------------
# 2. Read the CSVs
# ------------------------------------------------------------
read_diag <- function(dir, s) {
  f <- file.path(dir, sprintf(csv_fmt, s))
  if (!file.exists(f)) {
    warning("Missing file: ", f)
    return(NULL)
  }
  read_csv(f, show_col_types = FALSE) %>%
    filter(Year <= last_year)
}


# ------------------------------------------------------------
# 3. Monthly rows -> one value per year and variable
# ------------------------------------------------------------
to_annual <- function(df) {
  
  if (is.null(df)) return(NULL)
  out <- list()
  
  for (k in seq_len(nrow(csv_vars))) {
    
    col <- csv_vars$col[k]
    how <- csv_vars$how[k]
    
    if (!col %in% names(df)) {
      warning("Column not in CSV: ", col)
      next
    }
    
    a <- df %>%
      group_by(Year) %>%
      summarise(value = switch(how,
                               sum     = sum(.data[[col]], na.rm = TRUE),
                               first   = first(.data[[col]]),
                               kmonths = mean(.data[[col]][Month %in% k_months], na.rm = TRUE)),
                .groups = "drop") %>%
      mutate(var = col)
    
    out[[length(out) + 1]] <- a
  }
  bind_rows(out)
}

rows <- list()
monthly_rows <- list()    # raw monthly rows of the mean run, for the seasonal cycle

for (run in names(run_dirs)) {
  
  if (is.null(run_dirs[[run]])) next
  
  for (s in ssps) {
    
    cat("Reading", run, "run, SSP", s, "\n")
    raw <- read_diag(run_dirs[[run]], s)
    if (is.null(raw)) next
    
    if (run == "mean") {
      monthly_rows[[length(monthly_rows) + 1]] <- raw %>% mutate(SSP = ssp_labels[[s]])
    }
    
    tab <- to_annual(raw)
    tab$SSP <- ssp_labels[[s]]
    tab$run <- run
    rows[[length(rows) + 1]] <- tab
  }
}

annual <- bind_rows(rows) %>%
  left_join(csv_vars %>% select(var = col, kind, label), by = "var")

write_csv(annual, file.path(out_dir, "pan_arctic_annual_from_csv.csv"))


# ------------------------------------------------------------
# 4. Cumulative sums of the fluxes (Pg N), separately per run
# ------------------------------------------------------------
if (plot_cumulative) {
  cumulative <- annual %>%
    filter(kind == "flux") %>%
    group_by(run, SSP, var) %>%
    arrange(Year, .by_group = TRUE) %>%
    mutate(value = cumsum(replace_na(value, 0)),
           var   = paste0(var, "_cumulative"),
           label = paste("Cumulative", tolower(label)),
           kind  = "cumulative") %>%
    ungroup()
  annual <- bind_rows(annual, cumulative)
}


# ------------------------------------------------------------
# 5. Anomaly (optional) and 20-year running mean, per run;
#    then uncertainty = |plus - minus| / 2
# ------------------------------------------------------------
per_run <- annual %>%
  group_by(run, SSP, var) %>%
  arrange(Year, .by_group = TRUE) %>%
  mutate(
    ref_mean = mean(value[Year >= ref_start & Year <= ref_end], na.rm = TRUE),
    y        = if (plot_anomaly) value - ref_mean else value,
    y_20     = zoo::rollapply(y, 20, mean, align = "center", fill = NA)
  ) %>%
  ungroup()

series <- per_run %>%
  select(run, SSP, var, kind, label, Year, y, y_20) %>%
  pivot_wider(names_from = run, values_from = c(y, y_20))

# runs that were not found -> no uncertainty
for (col in c("y_plus", "y_minus", "y_20_plus", "y_20_minus")) {
  if (!col %in% names(series)) series[[col]] <- NA_real_
}

series <- series %>%
  mutate(y_sd    = abs(y_plus    - y_minus)    / 2,
         y_20_sd = abs(y_20_plus - y_20_minus) / 2) %>%
  rename(y = y_mean, y_20 = y_20_mean) %>%
  mutate(SSP = factor(SSP, levels = unname(ssp_labels)))

write_csv(series, file.path(out_dir,
                            paste0("pan_arctic_series_", ifelse(plot_anomaly, "anomaly", "absolute"), ".csv")))


# ------------------------------------------------------------
# 6. Plot function
# Historical (1850-2014) = black + grey ribbon
# Future (2015-2099)      = SSP colours + coloured ribbons
# ------------------------------------------------------------
plot_diag <- function(df,
                      yvar,
                      ylab,
                      title,
                      filename,
                      sd_var = NULL,
                      running_var = NULL,
                      running_sd_var = NULL,
                      show_raw = TRUE) {
  
  df <- df %>% arrange(SSP, Year)
  
  df_hist   <- df %>% filter(Year <= 2014)
  df_future <- df %>% filter(Year >= 2015)
  
  # SSP used to represent the shared historical period
  historical_ssp <- unique(df$SSP)[1]
  df_hist <- df_hist %>% filter(SSP == historical_ssp)
  
  p <- ggplot(df, aes(x = Year, group = SSP))
  
  # raw values
  if (show_raw) {
    p <- p +
      geom_line(data = df_hist, aes(y = .data[[yvar]]),
                color = "black", linewidth = 0.2, alpha = 0.45, na.rm = TRUE) +
      geom_line(data = df_future, aes(y = .data[[yvar]], color = SSP),
                linewidth = 0.2, alpha = 0.45, na.rm = TRUE)
  }
  
  # raw SD ribbon (optional)
  if (show_raw && !is.null(sd_var) && sd_var %in% names(df)) {
    p <- p +
      geom_ribbon(data = df_hist,
                  aes(ymin = .data[[yvar]] - .data[[sd_var]],
                      ymax = .data[[yvar]] + .data[[sd_var]]),
                  fill = "grey", color = NA, alpha = 0.03) +
      geom_ribbon(data = df_future,
                  aes(ymin = .data[[yvar]] - .data[[sd_var]],
                      ymax = .data[[yvar]] + .data[[sd_var]],
                      fill = SSP),
                  color = NA, alpha = 0.03)
  }
  
  # running mean
  if (!is.null(running_var) && running_var %in% names(df)) {
    p <- p +
      geom_line(data = df_hist, aes(y = .data[[running_var]]),
                color = "black", linewidth = 0.2, na.rm = TRUE) +
      geom_line(data = df_future, aes(y = .data[[running_var]], color = SSP),
                linewidth = 0.2, na.rm = TRUE)
  }
  
  # running-mean SD ribbon
  if (!is.null(running_var) && !is.null(running_sd_var) &&
      running_var %in% names(df) && running_sd_var %in% names(df)) {
    p <- p +
      geom_ribbon(data = df_hist,
                  aes(ymin = .data[[running_var]] - .data[[running_sd_var]],
                      ymax = .data[[running_var]] + .data[[running_sd_var]]),
                  fill = "grey", color = NA, alpha = 0.5) +
      geom_ribbon(data = df_future,
                  aes(ymin = .data[[running_var]] - .data[[running_sd_var]],
                      ymax = .data[[running_var]] + .data[[running_sd_var]],
                      fill = SSP),
                  color = NA, alpha = 0.2)
  }
  
  p +
    geom_vline(xintercept = 2015, linetype = "dashed", color = "black") +
    scale_color_manual(values = ssp_colors, drop = FALSE) +
    scale_fill_manual(values = ssp_colors, drop = FALSE, guide = "none") +
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.text    = element_text(size = 4),
          axis.title   = element_text(size = 4),
          legend.text  = element_text(size = 4),
          legend.title = element_text(size = 4)) +
    labs(x = "Year", y = ylab, color = "", title = title)
}


# ------------------------------------------------------------
# 7. One figure per variable
# ------------------------------------------------------------
plots <- list()

for (v in unique(series$var)) {
  
  d    <- series %>% filter(var == v)
  kind <- d$kind[1]
  lab  <- d$label[1]
  
  # annual fluxes in Tg N yr-1 (cumulative and stocks stay in Pg N)
  if (kind == "flux" && flux_in_Tg) {
    d <- d %>% mutate(across(any_of(c("y", "y_sd", "y_20", "y_20_sd",
                                      "y_plus", "y_minus", "y_20_plus", "y_20_minus")),
                             ~ .x * 1000))
  }
  flux_unit <- if (flux_in_Tg) "[Tg N" else "[Pg N"
  
  ylab <- switch(kind,
                 flux = bquote(.(lab) ~ .(flux_unit) ~ yr^-1 * "]"),
                 rate = bquote(.(lab) ~ "[-]"),
                 bquote(.(lab) ~ "[Pg N]"))            # stock, cumulative
  
  tag <- paste0(v, "_", ifelse(plot_anomaly, "anomaly", "absolute"))
  
  # no ribbons if the plus/minus runs are missing
  has_sd <- any(is.finite(d$y_sd))
  
  p <- plot_diag(d,
                 yvar           = "y",
                 ylab           = ylab,
                 title          = "",
                 filename       = paste0(tag, ".png"),
                 sd_var         = if (has_sd) "y_sd" else NULL,
                 running_var    = "y_20",
                 running_sd_var = if (has_sd) "y_20_sd" else NULL,
                 show_raw       = !(v %in% smooth_only))
  
  if (!plot_anomaly) p <- p + geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3)
  
  ggsave(file.path(out_dir, paste0(tag, ".png")), p,
         width = fig_width_cm, height = fig_height_cm, units = "cm", dpi = fig_dpi)
  ggsave(file.path(out_dir, paste0(tag, ".pdf")), p,
         width = fig_width_cm, height = fig_height_cm, units = "cm", device = cairo_pdf)
  
  plots[[v]] <- p
  print(p)
}


# ------------------------------------------------------------
# 8. Seasonal cycle: mean value per calendar month, per period and SSP
#    (mean run only; one panel per SSP, one line per period)
# ------------------------------------------------------------
monthly_all <- bind_rows(monthly_rows) %>%
  mutate(SSP = factor(SSP, levels = unname(ssp_labels)))

# assign each row to a period (rows outside all periods are dropped)
monthly_all$period <- NA_character_
for (k in seq_len(nrow(seasonal_periods))) {
  in_period <- monthly_all$Year >= seasonal_periods$start[k] &
    monthly_all$Year <= seasonal_periods$end[k]
  monthly_all$period[in_period] <- seasonal_periods$period[k]
}

seasonal_cols <- intersect(names(seasonal_vars), names(monthly_all))
missing_cols  <- setdiff(names(seasonal_vars), names(monthly_all))
if (length(missing_cols) > 0) warning("Columns not in CSV: ", paste(missing_cols, collapse = ", "))

seasonal <- monthly_all %>%
  filter(!is.na(period)) %>%
  group_by(SSP, period, Month) %>%
  summarise(across(all_of(seasonal_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(period = factor(period, levels = seasonal_periods$period))

write_csv(seasonal, file.path(out_dir, "seasonal_cycle_by_period.csv"))

for (col in seasonal_cols) {
  
  d <- seasonal %>%
    transmute(SSP, period, Month, value = .data[[col]] * seasonal_vars[[col]]$scale)
  
  p_season <- ggplot(d, aes(x = Month, y = value, colour = period, group = period)) +
    geom_line(linewidth = 0.4) +
    geom_point(size = 0.8) +
    facet_wrap(~SSP, ncol = 2) +
    scale_x_continuous(breaks = 1:12, labels = month.abb) +
    scale_colour_manual(values = setNames(period_colours, seasonal_periods$period)) +
    labs(x = "Month", y = seasonal_vars[[col]]$ylab, colour = "") +
    theme_bw(base_size = 7) +
    theme(axis.text.x      = element_text(angle = 45, hjust = 1),
          strip.background = element_rect(fill = "grey85"),
          legend.position  = "bottom")
  
  ggsave(file.path(out_dir, paste0("seasonal_", col, ".png")), p_season,
         width = seasonal_width_cm, height = seasonal_height_cm, units = "cm", dpi = fig_dpi)
  ggsave(file.path(out_dir, paste0("seasonal_", col, ".pdf")), p_season,
         width = seasonal_width_cm, height = seasonal_height_cm, units = "cm", device = cairo_pdf)
  print(p_season)
}


# ------------------------------------------------------------
# 9. Pan-Arctic mean per m2 (total / domain area)
#    Pg N -> g N: x 1e15; / domain area (m2)
#    (terra is called with terra:: so it does not mask dplyr::select)
# ------------------------------------------------------------
if (is.null(domain_area_m2)) {
  r_dom <- terra::crop(terra::rast(domain_file), terra::ext(-179.95, 179.95, 60, 90))
  land  <- !is.na(r_dom[[terra::nlyr(r_dom)]])            # same land mask as the model
  cell_area <- terra::cellSize(r_dom[[1]], unit = "m")
  domain_area_m2 <- terra::global(terra::mask(cell_area, land, maskvalues = 0),
                                  "sum", na.rm = TRUE)[1, 1]
}
cat("Domain area used for pan-Arctic means [10^6 km2]:", round(domain_area_m2 / 1e12, 3), "\n")

pg_to_gm2 <- 1e15 / domain_area_m2

series_area <- series %>%
  filter(var %in% per_area_vars) %>%
  mutate(across(any_of(c("y", "y_sd", "y_20", "y_20_sd",
                         "y_plus", "y_minus", "y_20_plus", "y_20_minus")),
                ~ .x * pg_to_gm2))

write_csv(series_area, file.path(out_dir,
                                 paste0("pan_arctic_mean_per_m2_", ifelse(plot_anomaly, "anomaly", "absolute"), ".csv")))

for (v in unique(series_area$var)) {
  
  d    <- series_area %>% filter(var == v)
  kind <- d$kind[1]
  lab  <- d$label[1]
  
  ylab <- if (kind == "flux") {
    bquote(.(lab) ~ "[g N" ~ m^-2 ~ yr^-1 * "]")
  } else {
    bquote(.(lab) ~ "[g N" ~ m^-2 * "]")                # cumulative
  }
  
  tag <- paste0(v, "_per_m2_", ifelse(plot_anomaly, "anomaly", "absolute"))
  has_sd <- any(is.finite(d$y_sd))
  
  p <- plot_diag(d,
                 yvar           = "y",
                 ylab           = ylab,
                 title          = "",
                 filename       = paste0(tag, ".png"),
                 sd_var         = if (has_sd) "y_sd" else NULL,
                 running_var    = "y_20",
                 running_sd_var = if (has_sd) "y_20_sd" else NULL,
                 show_raw       = !(v %in% smooth_only))
  
  if (!plot_anomaly) p <- p + geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3)
  
  ggsave(file.path(out_dir, paste0(tag, ".png")), p,
         width = fig_width_cm, height = fig_height_cm, units = "cm", dpi = fig_dpi)
  ggsave(file.path(out_dir, paste0(tag, ".pdf")), p,
         width = fig_width_cm, height = fig_height_cm, units = "cm", device = cairo_pdf)
  
  plots[[paste0(v, "_per_m2")]] <- p
  print(p)
}


# ------------------------------------------------------------
# 10. Numbers for the results text (Tg N, except the organic pool in Pg N)
#     Every quantity is computed separately for the mean, plus and minus
#     run; uncertainty = |plus - minus| / 2.
# ------------------------------------------------------------
pd <- c(2000, 2020)     # present day
pind <- c(1850, 1900)   # preindustrial
fu <- c(2080, 2099)     # end of century

in_p <- function(y, p) y >= p[1] & y <= p[2]

# one number per SSP and run from an annual series (values converted to Tg)
calc <- function(varname, fun, scale = 1000) {
  annual %>%
    filter(var == varname) %>%
    group_by(SSP, run) %>%
    arrange(Year, .by_group = TRUE) %>%
    summarise(value = fun(Year, value * scale), .groups = "drop")
}

quantities <- list(
  # pre-thaw inorganic N, annual (Tg N yr-1)
  prethaw_1850_1970_below_today =
    calc("inorg_rapid_available", function(y, v) mean(v[in_p(y, pd)]) - mean(v[in_p(y, c(1850, 1970))])),
  prethaw_2080_2099_vs_today =
    calc("inorg_rapid_available", function(y, v) mean(v[in_p(y, fu)]) - mean(v[in_p(y, pd)])),
  # pre-thaw inorganic N, cumulative at 2099 relative to 2000-2020 (Tg N)
  prethaw_cumulative_2099_vs_today =
    calc("inorg_rapid_available_cumulative", function(y, v) v[y == 2099] - mean(v[in_p(y, pd)])),
  # thawed organic N pool, mean 2080-2099 (Pg N)
  organic_pool_2080_2099_PgN =
    calc("organic_pool_remaining_pg", function(y, v) mean(v[in_p(y, fu)]), scale = 1),
  # mineralisation, annual (Tg N yr-1)
  mineralised_historical_increase =
    calc("mineralised_pg_monthly", function(y, v) mean(v[in_p(y, pd)]) - mean(v[in_p(y, pind)])),
  mineralised_2080_2099_vs_today =
    calc("mineralised_pg_monthly", function(y, v) mean(v[in_p(y, fu)]) - mean(v[in_p(y, pd)])),
  # total inorganic N, cumulative (Tg N)
  total_released_1850_to_2020 =
    calc("total_inorg_pg_monthly", function(y, v) sum(v[y <= 2020])),
  total_cumulative_2099_vs_today =
    calc("total_inorg_pg_monthly_cumulative", function(y, v) v[y == 2099] - mean(v[in_p(y, pd)])),
  mineralised_cumulative_2099_vs_today =
    calc("mineralised_pg_monthly_cumulative", function(y, v) v[y == 2099] - mean(v[in_p(y, pd)]))
)

numbers <- bind_rows(quantities, .id = "quantity") %>%
  pivot_wider(names_from = run, values_from = value) %>%
  mutate(sd   = if (all(c("plus", "minus") %in% names(.))) abs(plus - minus) / 2 else NA_real_,
         text = ifelse(is.na(sd),
                       sprintf("%.2f", mean),
                       sprintf("%.2f\\,$\\pm$\\,%.2f", mean, sd)))

# peak of the pre-thaw release (20-year running mean of the anomaly, after 2015)
prethaw_runs <- annual %>%
  filter(var == "inorg_rapid_available") %>%
  group_by(SSP, run) %>%
  arrange(Year, .by_group = TRUE) %>%
  mutate(anom_Tg = (value - mean(value[in_p(Year, pd)])) * 1000,
         anom_20 = zoo::rollapply(anom_Tg, 20, mean, align = "center", fill = NA)) %>%
  ungroup()

peak_year <- prethaw_runs %>%
  filter(run == "mean", Year >= 2015, is.finite(anom_20)) %>%
  group_by(SSP) %>%
  slice_max(anom_20, n = 1, with_ties = FALSE) %>%
  select(SSP, peak_year = Year)

prethaw_peak <- prethaw_runs %>%
  inner_join(peak_year, by = c("SSP", "Year" = "peak_year")) %>%
  select(SSP, run, Year, anom_20) %>%
  pivot_wider(names_from = run, values_from = anom_20) %>%
  mutate(sd   = if (all(c("plus", "minus") %in% names(.))) abs(plus - minus) / 2 else NA_real_,
         text = sprintf("%.2f\\,$\\pm$\\,%.2f (peak year %d)", mean, sd, Year))

# share of mineralisation in the cumulative release, and SSP5-8.5 / SSP1-2.6 ratio
cum_mean <- numbers %>%
  filter(quantity %in% c("total_cumulative_2099_vs_today", "mineralised_cumulative_2099_vs_today")) %>%
  select(SSP, quantity, mean) %>%
  pivot_wider(names_from = quantity, values_from = mean) %>%
  mutate(mineralisation_share_pct = 100 * mineralised_cumulative_2099_vs_today /
           total_cumulative_2099_vs_today)

ratio_high_low <- with(cum_mean,
                       total_cumulative_2099_vs_today[SSP == "SSP5-8.5"] / total_cumulative_2099_vs_today[SSP == "SSP1-2.6"])

write_csv(numbers,      file.path(out_dir, "numbers_for_text.csv"))
write_csv(prethaw_peak, file.path(out_dir, "numbers_prethaw_peak.csv"))

cat("\n=== Numbers for the text (mean \u00b1 |plus - minus|/2) ===\n")
print(numbers %>% select(quantity, SSP, text) %>%
        pivot_wider(names_from = SSP, values_from = text), width = Inf)
cat("\n=== Peak of pre-thaw inorganic N release (Tg N yr-1 above 2000-2020, 20-yr mean) ===\n")
print(prethaw_peak %>% select(SSP, text), width = Inf)
cat("\n=== Share of mineralisation in the cumulative release by 2099 (%) ===\n")
print(cum_mean %>% select(SSP, mineralisation_share_pct))
cat("\nRatio SSP5-8.5 / SSP1-2.6 (cumulative release by 2099):", round(ratio_high_low, 2), "\n")

cat("\nDone. Figures and tables in", out_dir, "\n")