# =============================================================================
# merge_boot_reps.R -- combine rep_*.rds into the bootstrap CI objects,
# replicating the aggregation of the "np-crisis alpha bootstrap" chunk in
# NP_Crisis_r2_noReserves.Rmd (pool all rows, group by Year.cutoff,
# 90 percent quantile band + sd, then the Year.cutoff + 2 shift).
#
# Usage: Rscript merge_boot_reps.R
# Outputs: alpha_calibration_bootstrap_ci.rds / .csv
# =============================================================================
files <- list.files("reps", pattern = "^rep_.*\\.rds$", full.names = TRUE)
cat(length(files), "replication files found\n")
stopifnot(length(files) > 0)
df <- do.call(rbind, lapply(files, readRDS))
cat(nrow(df), "rows;", length(unique(df$rep)), "unique replications\n")

bootstrap_level <- 0.90
yrs <- sort(unique(df$Year.cutoff))
alpha_ci <- do.call(rbind, lapply(yrs, function(y) {
  z <- df$Alpha[df$Year.cutoff == y]
  data.frame(Year.cutoff = y,
             Alpha_low  = as.numeric(quantile(z, probs = (1 - bootstrap_level) / 2, na.rm = TRUE)),
             Alpha_high = as.numeric(quantile(z, probs = 1 - (1 - bootstrap_level) / 2, na.rm = TRUE)),
             Alpha_se   = sd(z, na.rm = TRUE))
}))
alpha_ci$Year.cutoff <- alpha_ci$Year.cutoff + 2   # as in the Rmd

saveRDS(alpha_ci, "alpha_calibration_bootstrap_ci.rds")
write.csv(alpha_ci, "alpha_calibration_bootstrap_ci.csv", row.names = FALSE)
cat("wrote alpha_calibration_bootstrap_ci.rds / .csv\n")
print(alpha_ci, digits = 3)
