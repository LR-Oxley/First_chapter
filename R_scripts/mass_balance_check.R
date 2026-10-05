library(terra)
library(dplyr)

# -----------------------------
# 0. Point to your saved outputs
# -----------------------------
out_dir <- "monthly_mineralised/mean_2perc_baseline"
ssp <- "585"   # set to whichever SSP you want to check

diagnostic_csv <- file.path(
  out_dir,
  paste0("arctic_monthly_diagnostics_w_temp", ssp, "_1850_2100.csv")
)

new_thawed_organic_N_nc <- file.path(
  out_dir,
  paste0("arctic_new_thawed_organic_N_yearly_", ssp, "_mean_1850_2099_w_temp.nc")
)

# -----------------------------
# 1. Load monthly diagnostics (has mineralised_pg_monthly AND
#    organic_pool_remaining_pg already joined in per year)
# -----------------------------
diagnostic_df <- read.csv(diagnostic_csv)

# yearly mineralised N (Pg), summed over the 12 months
mineralised_pg_yearly <- diagnostic_df %>%
  group_by(Year) %>%
  summarise(mineralised_pg_yearly = sum(mineralised_pg_monthly, na.rm = TRUE)) %>%
  arrange(Year)

# organic_pool_remaining_pg is constant within a year (it was joined onto every
# month), so just take one value per year
organic_pool_pg_yearly <- diagnostic_df %>%
  group_by(Year) %>%
  summarise(organic_pool_remaining_pg = first(organic_pool_remaining_pg)) %>%
  arrange(Year)

# -----------------------------
# 2. Load new_thawed_organic_N raster and compute yearly Pg N
# -----------------------------
new_thawed_organic_N <- rast(new_thawed_organic_N_nc)
area_rast <- cellSize(new_thawed_organic_N[[1]], unit = "m")

new_thawed_organic_N_pg_yearly <- global(
  new_thawed_organic_N * area_rast,
  "sum",
  na.rm = TRUE
)[, 1] / 1e12

years_from_nc <- as.integer(format(time(new_thawed_organic_N), "%Y"))

new_thawed_df <- data.frame(
  Year = years_from_nc,
  new_thawed_organic_N_pg = new_thawed_organic_N_pg_yearly
)

# -----------------------------
# Load new_thawed_total_N and compute yearly Pg N
# -----------------------------
new_thawed_total_N_nc <- file.path(
  out_dir,
  paste0("arctic_new_thawed_total_N_yearly_", ssp, "_mean_1850_2099_w_temp.nc")
)

new_thawed_total_N <- rast(new_thawed_total_N_nc)

new_thawed_total_N_pg_yearly <- global(
  new_thawed_total_N * area_rast,
  "sum",
  na.rm = TRUE
)[, 1] / 1e12

new_thawed_total_df <- data.frame(
  Year = as.integer(format(time(new_thawed_total_N), "%Y")),
  new_thawed_total_N_pg = new_thawed_total_N_pg_yearly
)

# -----------------------------
# 2. Join into mass_balance_df and add cumulative sum
# -----------------------------
mass_balance_df <- mass_balance_df %>%
  left_join(new_thawed_total_df, by = "Year") %>%
  arrange(Year) %>%
  mutate(
    cum_new_thawed_total_N = cumsum(new_thawed_total_N_pg)
  )

print(tail(mass_balance_df))

# -----------------------------
# 3. Assemble mass_balance_df exactly as before
# -----------------------------
mass_balance_df <- new_thawed_df %>%
  left_join(mineralised_pg_yearly, by = "Year") %>%
  left_join(organic_pool_pg_yearly, by = "Year") %>%
  arrange(Year) %>%
  mutate(
    cum_new_thawed_organic_N = cumsum(new_thawed_organic_N_pg),
    cum_mineralised_N        = cumsum(mineralised_pg_yearly),
    predicted_pool           = cum_new_thawed_organic_N - cum_mineralised_N,
    residual                 = predicted_pool - organic_pool_remaining_pg,
    relative_residual        = residual / pmax(abs(predicted_pool), 1e-9)
  )

print(tail(mass_balance_df))

library(ggplot2)
library(tidyr)

mass_balance_long <- mass_balance_df %>%
  select(Year, cum_new_thawed_total_N, cum_new_thawed_organic_N,
         cum_mineralised_N, predicted_pool, organic_pool_remaining_pg) %>%
  pivot_longer(cols = -Year, names_to = "term", values_to = "Pg_N")


p_mass_balance <- ggplot(mass_balance_long, aes(x = Year, y = Pg_N)) +
  geom_line(data = ~ filter(.x, term == "organic_pool_remaining_pg"),
            aes(colour = term, linetype = term), linewidth = 0.2) +
  geom_line(data = ~ filter(.x, term == "predicted_pool"),
            aes(colour = term, linetype = term), linewidth = 0.2) +
  geom_line(data = ~ filter(.x, term == "cum_mineralised_N"),
            aes(colour = term, linetype = term), linewidth = 0.2) +
  geom_line(data = ~ filter(.x, term == "cum_new_thawed_organic_N"),
            aes(colour = term, linetype = term), linewidth = 0.2) +
  geom_line(data = ~ filter(.x, term == "cum_new_thawed_total_N"),
            aes(colour = term, linetype = term), linewidth = 0.2) +   # drawn last = on top
  theme_bw() +
  labs(
    title = "N pool mass balance check",
    y = "Pg N", colour = NULL, linetype = NULL
  ) +
  scale_colour_manual(
    values = c(
      cum_new_thawed_total_N    = "black",
      cum_new_thawed_organic_N  = "darkorange",
      cum_mineralised_N         = "firebrick",
      predicted_pool            = "grey60",
      organic_pool_remaining_pg = "blue"
    ),
    breaks = c("cum_new_thawed_total_N", "cum_new_thawed_organic_N",
               "cum_mineralised_N", "predicted_pool", "organic_pool_remaining_pg"),
    labels = c(
      cum_new_thawed_total_N    = "Cumulative new thawed total N",
      cum_new_thawed_organic_N  = "Cumulative new thawed organic N",
      cum_mineralised_N         = "Cumulative mineralised N",
      predicted_pool            = "Predicted pool (thawed organic N \u2212 mineralised)",
      organic_pool_remaining_pg = "Actual organic pool remaining"
    )
  ) +
  scale_linetype_manual(
    values = c(
      cum_new_thawed_total_N    = "solid",
      cum_new_thawed_organic_N  = "solid",
      cum_mineralised_N         = "solid",
      predicted_pool            = "dashed",
      organic_pool_remaining_pg = "solid"
    ),
    guide = "none"
  ) +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 8),
    legend.text = element_text(size = 8)
  ) +
  guides(colour = guide_legend(nrow = 3))

print(p_mass_balance)


p_mass_balance <- p_mass_balance+
  theme(
    title = element_text(size =4),
    legend.margin = margin(t = 0),
    legend.box.margin = margin(0.5, 0.5, 0.5, 0.5),
    plot.margin = margin(1, 1, 1, 1),
    legend.position = "bottom", 
    legend.key.size = unit(0.3, "cm"),      # size of each colour swatch
    legend.text = element_text(size =5),     # label text
    legend.title = element_text(size = 5),    # legend title text
    
    # --- shrink the gap between axis titles and the panel ---
    axis.title.x = element_text(margin = margin(t = 2)),   # was default ~half a line
    axis.title.y = element_text(margin = margin(r = 2)),
    
    # --- pull tick labels closer to the axis (less padding) ---
    axis.text.x = element_text(margin = margin(t = 1)),
    axis.text.y = element_text(margin = margin(r = 1)),
    
    # --- optionally shrink the fonts a touch, which shrinks their box too ---
    axis.title = element_text(size = 4),
    axis.text  = element_text(size = 4)
  )

ggsave(
  "p_mass_balance.pdf",
  p_mass_balance,
  width = 8.5,
  height = 6,
  units = "cm", 
  dpi = 500
)




p_residual <- ggplot(mass_balance_df, aes(x = Year, y = residual)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 1) +
  theme_bw() +
  labs(title = "Mass balance residual", y = "Residual (Pg N)")

print(p_residual)


thawed_permafrost_N_nc <- file.path(
  out_dir, paste0("arctic_permafrost_thawed_total_N_", ssp, "_mean_1850_2099_w_temp.nc")
)
thawed_permafrost_N <- rast(thawed_permafrost_N_nc)

f_org <- 1 - 0.01  # f_inorg_rapid = 0.01 from your params

initial_pool_pg <- global(
  (thawed_permafrost_N[[1]] * f_org) * area_rast,
  "sum", na.rm = TRUE
)[1, 1] / 1e12

mass_balance_df <- mass_balance_df %>%
  mutate(predicted_pool_incl_init = initial_pool_pg + cum_new_thawed_organic_N - cum_mineralised_N,
         residual_incl_init = predicted_pool_incl_init - organic_pool_remaining_pg)


depth_check <- read.csv(file.path(out_dir, paste0("depth_distribution_conservation_check_", ssp, "_1850_2099.csv")))

ggplot(depth_check, aes(Year, relative_error)) +
  geom_line() +
  geom_vline(xintercept = 2005, linetype = "dashed", colour = "red") +
  theme_bw() +
  labs(title = "Depth-redistribution relative error over time")

mass_balance_df$relative_residual_final <- mass_balance_df$residual / mass_balance_df$organic_pool_remaining_pg

ggplot(mass_balance_df, aes(Year, relative_residual_final)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line() +
  theme_bw() +
  labs(title = "Residual as fraction of organic pool", y = "Relative residual")


library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

annual <- readr::read_csv(file.path("monthly_mineralised/diagnostic_csv_analysis_60N",
                                    "pan_arctic_annual_from_csv.csv"))

f_pre <- 0.01   # pre-thaw inorganic fraction

mb <- annual %>%
  filter(run == "mean", SSP == "SSP5-8.5",
         var %in% c("thawed_N_annual_change_pg", "mineralised_pg_monthly",
                    "organic_pool_remaining_pg")) %>%
  select(Year, var, value) %>%
  pivot_wider(names_from = var, values_from = value) %>%
  arrange(Year) %>%
  mutate(
    cum_thawed_total = cumsum(thawed_N_annual_change_pg),
    cum_thawed_org   = cumsum(thawed_N_annual_change_pg * (1 - f_pre)),
    cum_mineralised  = cumsum(mineralised_pg_monthly),
    actual           = organic_pool_remaining_pg,
    # anchor the prediction to the first year's simulated pool
    predicted        = cum_thawed_org - cum_mineralised +
      (first(actual) - (first(cum_thawed_org) - first(cum_mineralised))),
    residual         = actual - predicted,
    rel_error_pct    = 100 * residual / pmax(abs(actual), 1e-6)
  )

# --- panel a: pools ---
pools <- mb %>%
  select(Year, `Cumulative newly thawed total N` = cum_thawed_total,
         `Cumulative newly thawed organic N` = cum_thawed_org,
         `Cumulative mineralised N` = cum_mineralised,
         `Simulated organic pool` = actual,
         `Expected organic pool (thawed organic ??? mineralised)` = predicted) %>%
  pivot_longer(-Year, names_to = "series", values_to = "PgN")

line_types <- c("Cumulative newly thawed total N"                      = "solid",
                "Cumulative newly thawed organic N"                    = "dotted",
                "Cumulative mineralised N"                             = "solid",
                "Simulated organic pool"                               = "solid",
                "Expected organic pool (thawed organic ??? mineralised)" = "dashed")
line_cols  <- c("Cumulative newly thawed total N"                      = "black",
                "Cumulative newly thawed organic N"                    = "#E69F00",
                "Cumulative mineralised N"                             = "#D55E00",
                "Simulated organic pool"                               = "#0072B2",
                "Expected organic pool (thawed organic ??? mineralised)" = "grey50")

p_a <- ggplot(pools, aes(Year, PgN, colour = series, linetype = series)) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = line_cols) +
  scale_linetype_manual(values = line_types) +
  labs(x = NULL, y = "Pg N", colour = NULL, linetype = NULL, tag = "a") +
  theme_minimal(base_size = 8) +
  theme(legend.position = "right", legend.key.width = unit(0.8, "cm"))

# --- panel b: residual (should be ~0) ---
p_b <- ggplot(mb, aes(Year, residual * 1000)) +        # Pg -> Tg
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_line(colour = "#0072B2", linewidth = 0.6) +
  labs(x = "Year", y = "Simulated ??? expected\npool [Tg N]", tag = "b") +
  theme_minimal(base_size = 8)

p_mb <- p_a / p_b + plot_layout(heights = c(3, 1))
print(p_mb)
ggsave("s_mass_balance.pdf", p_mb, width = 17, height = 11, units = "cm")

cat("Max |residual|:", signif(max(abs(mb$residual)) * 1000, 3), "Tg N;",
    "max relative error:", signif(max(abs(mb$rel_error_pct[mb$Year > 1900])), 3), "%\n")




library(terra)
pool <- rast(file.path("monthly_mineralised/mean_2perc_baseline",
                       "arctic_organic_N_pool_remaining_yearly_585_mean_1850_2099_w_temp.nc"))

n_valid <- global(!is.na(pool), "sum")[, 1]
par(mar = c(4, 4, 1, 1))
plot(1850:2099, n_valid, type = "l", xlab = "Year", ylab = "Cells with valid pool")

# years in which cells drop out
yrs <- 1850:2099
data.frame(Year = yrs[-1], cells_lost = -diff(n_valid))[diff(n_valid) < 0, ]
