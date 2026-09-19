# np-crisis (round-2 revision package)

Reproducibility package for "International Reserves, Crisis Risk, and Risk
Tolerance" (revised version). This is a **clean,
self-contained** repository: every table and figure in the revised paper,
its Supplementary Material, and the response letter can be regenerated from
the files here.

> DRAFT NOTE: this file is the round-2 replacement for `README.md`.
> Rename it to `README.md` when committing the round-2 package.

## What changed in round 2

1. **Reserve-coverage predictors excluded from the backtest.** The revised
   analysis file is `NP_Crisis_r2_noReserves.Rmd`; it removes `RES_BMGS`,
   `RES_BM`, `RES_GDP`, `RES_CH_GDP`, and `RES_D_S` from the predictor set,
   consistent with Section 2.2 of the paper. The risk-tolerance calibration
   was already reserve-free, so its outputs are unchanged relative to the
   previous round.
2. **New simulation study** (Section S.6 of the Supplementary Material) in
   `simulations_r2/`.
3. **Two calibration models documented.** In
   `NP_Crisis_r2_noReserves.Rmd`, `i.model = 1` (signal extraction)
   produces Figure 3 of the main text; `i.model = 5` (AdaBoost) produces
   the robustness series discussed in Section 5.2 and shown in Figure R.1
   of the response letter. The AdaBoost option requires
   `np-crisis-helper-funcs-r2.R`, which generalizes one branch of
   `reserves_match.mse.alpha()`.
4. `rebuild_tables_r2_noReserves.R` extracts the NP-ROC AUC and AUPC
   summary tables including the full-feature column (d = 66 after the
   exclusion), writing `AUC_test_corrected.csv` and
   `pAUC_test_corrected.csv`.

## Contents

- `NP_Crisis.Rmd`: round-1 analysis (kept for the record; superseded).
- `NP_Crisis_r2_noReserves.Rmd`: round-2 main analysis. **Knit this file.**
- `np-crisis-helper-funcs.R` / `np-crisis-helper-funcs-r2.R`: helper
  functions (the `-r2` version is sourced by the round-2 Rmd).
- `rebuild_tables_r2_noReserves.R`: builds the corrected Table 3 inputs.
- `simulations_r2/`: the simulation study (see its own README); standalone,
  needs only base R plus the `parallel` package, no project data.
- `data/`: input datasets.
- `outputs/`: round-1 outputs (kept for the record).
- `hpc_boot/`: self-contained SLURM workflow (NYU Torch) that computed the
  500-replication bootstrap band for the AdaBoost robustness recalibration;
  the merged band is `outputs_r2_noReserves/rds/alpha_calibration_bootstrap_ci.rds`
  (the AdaBoost run; the signal-extraction equivalents are under
  `outputs_r2_noReserves/rds/se_based_run1/`).
- `outputs_r2_noReserves/`: round-2 outputs, including the `rds/` cache
  that makes reproduction fast by default, and `se_based_run1/` with the
  signal-extraction calibration outputs backed up separately.

## Mapping from paper objects to files

| Paper object | Produced by | Output |
|---|---|---|
| Table 3 (NP-ROC AUC, AUPC, classic AUC) | `NP_Crisis_r2_noReserves.Rmd` + `rebuild_tables_r2_noReserves.R` | `outputs_r2_noReserves/*_corrected.csv`, `ROC_AUC_test_90 CI.csv` |
| Figure 2 (classic and NP-ROC curves) | `NP_Crisis_r2_noReserves.Rmd` | `outputs_r2_noReserves/roc_classic_nproc_models.png` |
| Figure 3 (revealed risk tolerance) | `NP_Crisis_r2_noReserves.Rmd` with `i.model = 1` | `outputs_r2_noReserves/rds/se_based_run1/` (identical to round 1) |
| Section 5.2 robustness series / letter Figure R.1 | `NP_Crisis_r2_noReserves.Rmd` with `i.model = 5` | `outputs_r2_noReserves/rds/alpha_calibration_*.rds` |
| Supplement S.6 (simulation study, Figures and tables) | `simulations_r2/run-all.R` | `simulations_r2/outputs/` |

## Reproduce the results

From the repo root:

1. Restore the R package environment (recommended):
   - `R -q -f scripts/restore-renv.R`
2. Render the round-2 report:
   - `Rscript -e "rmarkdown::render('NP_Crisis_r2_noReserves.Rmd')"`
3. Rebuild the corrected Table 3 inputs:
   - `Rscript rebuild_tables_r2_noReserves.R`
4. Run the simulation study (about 3 minutes in full mode on a laptop):
   - `cd simulations_r2 && Rscript run-all.R full`

## Notes

- By default `use_cached_rds <- TRUE`, so the shipped cache in
  `outputs_r2_noReserves/rds/` is reused and rendering is fast. Set it to
  `FALSE` to recompute everything; the full backtest (five classifiers over
  the alpha and screening grids) plus the 500-replication bootstrap takes
  several hours.
- The shipped cache was generated with the reserve-coverage exclusion in
  place; do not mix it with the round-1 cache in `outputs/rds/`.
- To reproduce the AdaBoost robustness calibration, set `i.model <- 5` in
  the calibration chunk and move the `alpha_calibration_*.rds` files out of
  `outputs_r2_noReserves/rds/` first (cached calibration results would
  otherwise be reused). The signal-extraction versions are preserved in
  `outputs_r2_noReserves/rds/se_based_run1/`.
- The umbrella threshold requires at least `log(delta)/log(1-alpha)`
  left-out crisis observations, so the smallest feasible missed-crisis
  bound at this sample size is about 0.06; grids therefore start at 0.10.
- `rmarkdown::render()` requires Pandoc; `R --vanilla -q -f scripts/knit.R`
  knits without it.
