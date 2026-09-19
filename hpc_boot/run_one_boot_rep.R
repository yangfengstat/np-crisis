# =============================================================================
# run_one_boot_rep.R -- one bootstrap replication of the AdaBoost risk-
# tolerance calibration, for a SLURM array task.
#
# Usage:  Rscript run_one_boot_rep.R <rep_id>
#         (rep_id falls back to $SLURM_ARRAY_TASK_ID)
# Env:    QUICK=1 runs a reduced grid (2 alphas x 2 years) as a plumbing test.
#
# Inputs: boot_inputs.rds (frozen locally from NP_Crisis_r2_noReserves.Rmd:
#         data, calibration table, grids, i.model = 5) + helper_funcs_hpc.R.
# Output: reps/rep_<id>.rds with (Country.Name, Year.cutoff, Alpha, rep).
#
# Seeding: set.seed(123000 + rep_id) per replication (documented deviation
# from the in-Rmd stream seeding; statistically equivalent for the bootstrap).
# =============================================================================
source("helper_funcs_hpc.R")

b <- suppressWarnings(as.integer(commandArgs(trailingOnly = TRUE)[1]))
if (is.na(b)) b <- suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")))
stopifnot(!is.na(b))

inp <- readRDS("boot_inputs.rds")
list2env(inp, envir = globalenv())  # firstvar must be a global for reserves_match.mse.alpha

if (Sys.getenv("QUICK") == "1") {
  alpha.choice_cal <- c(0.5, 0.9)
  year.cutoff_choice <- c(2006, 2007)
}

set.seed(123000 + b)
boot_countries <- sample(country_boot_pool, size = length(country_boot_pool), replace = TRUE)
boot_list <- lapply(boot_countries, function(cty) current_calibration[current_calibration$Country_Name == cty, ])
current_boot <- do.call(rbind, boot_list)
current_boot$outcome <- current_boot$SuddenStop_GI

t0 <- Sys.time()
m <- reserves_match.mse.alpha(alpha.choice_cal, current_boot, year.cutoff_choice,
                              test.size, Country.list, CrossSectional_Calibration,
                              model.choice_cal, i.model, delta = delta_0)
m$Reserves.error.squared <- m$Reserves.error^2
DT <- data.table::data.table(m)
ann <- DT[, .SD[which.min(Reserves.error.squared)], by = c("Country.Name", "Year.cutoff")]
out <- data.frame(Country.Name = ann$Country.Name, Year.cutoff = ann$Year.cutoff,
                  Alpha = ann$Alpha, rep = b)

dir.create("reps", showWarnings = FALSE)
saveRDS(out, sprintf("reps/rep_%03d.rds", b))
cat(sprintf("rep %d done: %d rows in %.1f min\n", b, nrow(out),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
