# =============================================================================
# Panel A: population-level missed-crisis (type II) error control.
# Shows: NP umbrella violation rate <= delta; plug-in quantile ~ 1/2;
# symmetric rule provides no control at all.
# =============================================================================
source("sim-helpers.R")

if (!exists("CFG_A")) {
  CFG_A <- list(B = 200, ns = c(2000), alphas = c(0.10, 0.15, 0.20, 0.30, 0.40, 0.50),
                delta = 0.1, seed = 20260609, cores = max(1, parallel::detectCores() - 2))
}

pars <- sim_params()
all_res <- list()
for (n in CFG_A$ns) {
  res <- panelA_run(CFG_A$B, n, CFG_A$alphas, CFG_A$delta, pars,
                    seed = CFG_A$seed + n, cores = CFG_A$cores)
  res$n <- n
  all_res[[as.character(n)]] <- res
}
resA <- do.call(rbind, all_res)
summA <- do.call(rbind, lapply(split(resA, resA$n), function(rr) {
  s <- panelA_summarize(rr, CFG_A$delta); s$n <- rr$n[1]; s
}))

saveRDS(resA, "outputs/panelA_raw.rds")
write.csv(summA, "outputs/panelA_summary.csv", row.names = FALSE)

pdf("outputs/panelA_error_control.pdf", width = 10, height = 5)
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.5, 1))
n0 <- CFG_A$ns[length(CFG_A$ns)]
s0 <- summA[summA$n == n0, ]
meths <- c("NP umbrella", "plug-in quantile", "symmetric")
ltys <- c(1, 2, 3); pchs <- c(19, 17, 15)
plot(NA, xlim = range(CFG_A$alphas), ylim = c(0, 1), xlab = expression(alpha),
     ylab = expression("population missed-crisis rate " * R[1]),
     main = paste0("Realized type II error (n = ", n0, ")"))
abline(0, 1, col = "grey60")
for (i in seq_along(meths)) {
  si <- s0[s0$method == meths[i], ]
  lines(si$alpha, si$R1, lty = ltys[i]); points(si$alpha, si$R1, pch = pchs[i])
}
legend("topleft", meths, lty = ltys, pch = pchs, bty = "n", cex = 0.85)
plot(NA, xlim = range(CFG_A$alphas), ylim = c(0, 1), xlab = expression(alpha),
     ylab = expression("violation rate  P(" * R[1] > alpha * ")"),
     main = paste0("Violation rate vs nominal delta = ", CFG_A$delta))
abline(h = CFG_A$delta, col = "grey60")
for (i in seq_along(meths)) {
  si <- s0[s0$method == meths[i], ]
  lines(si$alpha, si$viol, lty = ltys[i]); points(si$alpha, si$viol, pch = pchs[i])
}
dev.off()

cat("Panel A done. Violation-rate summary (n =", n0, "):\n")
print(summA[summA$n == n0, c("alpha", "method", "R1", "viol", "infeasible")],
      row.names = FALSE, digits = 3)
