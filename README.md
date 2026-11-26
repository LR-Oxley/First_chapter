Code for paper: " --- "

# Nitrogen mobilised by permafrost thaw in the pan-Arctic, 1850–2099

R code for the manuscript on nitrogen (N) mobilised by permafrost thaw in the pan-Arctic (60–90° N).
It estimates how much previously frozen soil N enters the active layer as permafrost thaws and how much of it becomes inorganic N.
The calculation runs from 1850 to 2099 for five CMIP6 models and four SSP scenarios.

The code only tracks N that was frozen before thaw, i.e. N stored below the pre-industrial (1850–1900) active layer.
It does not cover N cycling in soil that was already part of the active layer before 1850.

---

## Workflow

The scripts run in this order. Each one reads the output of the one before it.

```
ALD_models_plot.R  ->  N_thaw_models.R  ->  4D_min_models.R  ->  analysis scripts
   (z*)                 (N_AL)               (N_thaw, N_min, N_inorg)  (figures, numbers)
```

| Step | Script | What it does | Main output |
|---|---|---|---|
| 1 | `ALD_models_plot.R` | Bias-corrects the active-layer depth z\* of each CMIP6 model against ESA Permafrost_cci (2000–2020) and limits it to 0–5 m. Also makes the z\* figures, maps and tables. | `output_data_filled/ALD_corr_<ssp>_60N_<model>.nc` |
| 2 | `N_thaw_models.R` | Scales a land-cover-specific N depth profile to the soil N stock and integrates it from the surface to z\*. The result is the total N in the active layer, N_AL. | `output_data_filled/arctic_total_thawed_<ssp>_60N_<model>.nc` |
| 3 | `4D_min_models.R` | Computes the newly thawed N relative to 1850–1900 and splits it into pre-thaw inorganic and organic N. It then mineralises the organic pool every month, by depth layer, using soil temperature and moisture. | `<out_base>/<model>_2perc_baseline/arctic_monthly_diagnostics_w_temp<ssp>_1850_2100.csv` and yearly/monthly NetCDF maps |
| 4a | `analysis_csv_models.R` | Builds pan-Arctic time series, cumulative sums, seasonal cycles, per-m² values and the numbers used in the text, all as multi-model mean [min, max]. | figures and CSV tables in `diagnostic_csv_analysis_models_new/` |
| 4b | `N_mass_balance_models.R` | Checks the N mass balance (frozen, thawed organic, released as inorganic N) for every model run. | figures and CSV tables in `N_mass_balance_models/` |

---

## Usage

Steps 2 and 3 run once per SSP and model:

```bash
# step 1: all models and SSPs in one run
Rscript ALD_models_plot.R

# steps 2 and 3: one SSP and one model per call
Rscript N_thaw_models.R 585 Nor
Rscript 4D_min_models.R 585 Nor

# step 4: after all runs have finished
Rscript analysis_csv_models.R
Rscript N_mass_balance_models.R
```

Arguments:

- `ssp`: `126`, `245`, `370`, `585`
- `model`: `Nor` (NorESM2-LM), `MRI` (MRI-ESM2-0), `MIROC6`, `FGOALS` (FGOALS-f3-L), `E3SM` (E3SM-1-0; only SSP3-7.0 and SSP5-8.5)

Before running, set the output directory `out_base` at the top of `4D_min_models.R`, and `base_dir` in the analysis scripts, to a path on your system.

To run every model and scenario:

```bash
for ssp in 126 245 370 585; do
  for model in Nor MRI MIROC6 FGOALS E3SM; do
    if [ "$model" = "E3SM" ] && { [ "$ssp" = "126" ] || [ "$ssp" = "245" ]; }; then continue; fi
    Rscript N_thaw_models.R $ssp $model
    Rscript 4D_min_models.R $ssp $model
  done
done
```

---

## Input data

The input data are not included in this repository. The scripts expect them in `input_data/`:

| File | Content | Source |
|---|---|---|
| `ALD_<ssp>/<model>_regrid_filled.nc` | Active-layer depth of each CMIP6 model, regridded and gap-filled | CMIP6 |
| `ALD_<ssp>/<model>_regrid.nc` | Same, without gap-filling (used only to evaluate the models) | CMIP6 |
| `ESA_ALD_files/ESA_ALD_present_day.nc` | Observed active-layer depth, 2000–2020 | ESA Permafrost_cci |
| `TN_30deg_corr.nc` | Soil organic N stock, 0–3 m (kg N m⁻²) | Palmtag et al. (2022) |
| `LC_60deg.nc` | Land-cover classes (IGBP) | MODIS |
| `mean_verttemp_ssp<ssp>_60deg.nc` | Monthly soil temperature by depth (K) | CMIP6 |
| `sm_<ssp>_60deg.nc` | Relative soil moisture | CMIP6 |

The plotting scripts also need `R_scripts/plot_style.R`, which sets the colours, theme and saving options.

---

## Model in short

Symbols follow the Methods of the manuscript. The code uses the same names, with an underscore in place of the comma (N_thaw,org → `N_thaw_org`).

1. **N in the active layer.** For each land-cover class *i*:

   N_AL = 1 / (1 − f_inorg) · ∫₀^z\* n_org,conc(i, z) dz

   where n_org,conc(i, z) = a_i + b_i · exp(−k_i z), scaled to the soil N stock between 0 and 3 m. Wetlands have a uniform profile.

2. **Thawed N.**
   - N_thaw,cum(t) = N_AL(t) − mean N_AL over 1850–1900
   - N_thaw(t) = N_thaw,cum(t) − N_thaw,cum(t − 1)

3. **Split.**
   - Pre-thaw inorganic N: N_thaw,inorg = f_inorg · N_thaw (f_inorg = 0.01), released evenly over June–August.
   - The rest goes into an organic pool, N_thaw,org(x, y, d). It is spread over soil layers down to 5 m in proportion to the N profile.

4. **Mineralisation.** Every month, in every layer:

   N_min = k_base / 12 · ξ_env · N_thaw,org, with ξ_env = ξ_T · ξ_θ

   - k_base = 0.02 yr⁻¹
   - ξ_T is a peaked Arrhenius function of soil temperature (E_a = 80 kJ mol⁻¹, E_d = 200 kJ mol⁻¹, T_opt = 28 °C), normalised to 1 at 25 °C and set to 0 in frozen layers.
   - ξ_θ depends on relative soil moisture (optimum 0.65).

5. **Total mobilised inorganic N.** N_inorg = N_thaw,inorg + Σ_d N_min

---

## Diagnostics CSV

`arctic_monthly_diagnostics_w_temp<ssp>_1850_2100.csv` has one row per month and is the main input of the analysis scripts:

| Column | Symbol | Unit |
|--------|--------|------|
| `N_min_pg_monthly` | N_min | Pg N month⁻¹ |
| `N_inorg_pg_monthly` | N_inorg | Pg N month⁻¹ |
| `N_thaw_inorg_pg_monthly` | N_thaw,inorg | Pg N month⁻¹ |
| `N_thaw_pg` | N_thaw (annual, repeated in each month) | Pg N yr⁻¹ |
| `N_thaw_org_pg` | N_thaw,org at the end of the year | Pg N |
| `xi_T_mean_monthly` | ξ_T (depth-weighted mean) | – |
| `xi_T_pool_weighted_mean` | ξ_T, weighted by the organic pool | – |
| `xi_theta_mean_monthly` | ξ_θ | – |
| `xi_env_pool_weighted_mean` | ξ_env, weighted by the organic pool | – |
| `xi_env_depth_unweighted_mean` | ξ_env, unweighted mean over depth | – |

CSVs from earlier versions used other column names (`mineralised_*`, `total_inorg_*`, `pre_thaw_inorg_available`, `organic_pool_remaining_pg`, `k_*`). The analysis scripts rename these automatically when they read the file.

---

## Requirements

- R (≥ 4.2)
- Packages: `terra`, `sf`, `dplyr`, `tidyr`, `readr`, `ggplot2`, `zoo`, `patchwork`, `rnaturalearth`

```r
install.packages(c("terra", "sf", "dplyr", "tidyr", "readr", "ggplot2",
                   "zoo", "patchwork", "rnaturalearth"))
```

A single run of `4D_min_models.R` handles monthly, depth-resolved rasters for 250 years. It takes several hours and needs a lot of memory, so it is best run on a cluster, one job per SSP and model.

---

## Citation

If you use this code, please cite the manuscript (reference to be added on publication).

## Licence

<!-- add a licence, e.g. MIT or CC BY 4.0 -->
