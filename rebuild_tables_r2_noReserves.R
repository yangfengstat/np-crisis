# Rebuild NP-ROC AUC and AUPC summary tables for the no-reserves rerun,
# including the full feature set (66 predictors), which the hardcoded
# size grid c(10,...,70,71) in the Rmd missed. New file; touches nothing else.
library(zoo)
res <- readRDS("outputs_r2_noReserves/rds/nproc_results.rds")
res <- as.data.frame(res, stringsAsFactors = FALSE)
model.choice <- c("SE", "logistic", "svm", "randomforest", "ada")
test <- res[, c("missed.crises", "screen.size", paste0("test_", model.choice))]
names(test) <- c("missed.crises", "screen.size", model.choice)
sizes <- sort(unique(test$screen.size))
cat("screen sizes present:", sizes, "\n")
pauc_np <- function(x, y, max_fpr = 0.25) {
  keep <- is.finite(x) & is.finite(y); x <- x[keep]; y <- y[keep]
  if (length(x) < 2) return(NA_real_)
  id <- order(x); x <- x[id]; y <- y[id]
  if (max_fpr <= min(x)) return(0)
  if (max_fpr < max(x)) {
    y_max <- stats::approx(x, y, xout = max_fpr, ties = "ordered")$y
    keep <- x < max_fpr; x <- c(x[keep], max_fpr); y <- c(y[keep], y_max)
  } else { keep <- x <= max_fpr; x <- x[keep]; y <- y[keep] }
  if (length(x) < 2) return(0)
  sum(diff(x) * zoo::rollmean(y, 2))
}
auc <- pauc <- as.data.frame(matrix(NA_real_, nrow = length(model.choice), ncol = length(sizes)))
colnames(auc) <- colnames(pauc) <- paste0("features_", sizes)
for (s in sizes) for (i in seq_along(model.choice)) {
  d <- test[test$screen.size == s, ]
  x <- d$missed.crises; y <- d[, i + 2]; id <- order(x)
  auc[i, which(sizes == s)] <- sum(diff(x[id]) * zoo::rollmean(y[id], 2))
  pauc[i, which(sizes == s)] <- pauc_np(x, y)
}
auc <- cbind(model = model.choice, auc); pauc <- cbind(model = model.choice, pauc)
write.csv(auc, "outputs_r2_noReserves/AUC_test_corrected.csv", row.names = FALSE)
write.csv(pauc, "outputs_r2_noReserves/pAUC_test_corrected.csv", row.names = FALSE)
cat("\nNP-ROC AUC (corrected, incl. full set):\n"); print(auc, digits = 3)
cat("\nAUPC (corrected, incl. full set):\n"); print(pauc, digits = 3)
