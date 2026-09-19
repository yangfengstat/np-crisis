# Round-2 simulation study (new files; nothing in the parent folder is modified)

Implements the Monte Carlo study reported in Section S.6 of the Supplementary
Material: a numerical assessment of the proposed framework and of the
reliability of the policy tool under known ground truth.

## Files

- `sim-helpers.R`: DGP, closed-form population errors, NP umbrella
  (mirrors `npc_aux.core()` in `../np-crisis-helper-funcs.R`), comparator
  rules, insurance-model welfare functions (mirror `rho.star()` and
  Supplementary Material S.2), and the four panel engines.
- `sim-panelA.R`: population-level missed-crisis control; violation rates
  vs the nominal delta for NP umbrella / plug-in quantile / symmetric rule.
- `sim-panelB.R`: B1 exact Prop-1 asymmetry (downward deviations cost
  2.3x-5.9x upward ones); B1b Eq-12 weight diagnostic (with a consistent
  estimator the class priors dominate, so omega1/omega0 is NOT generically
  above one; flagged for the authors); B2 pipeline welfare cost under low
  and elevated true risk (optimal alpha falls as risk rises; symmetric rule
  dominated by about 10x); B3 analytic crossover at alpha_bar = Phi(-Delta/2).
- `sim-panelC.R`: recovery of a known time-varying risk tolerance alpha*_t
  through the two-stage pipeline (validates "revealed risk tolerance").
- `sim-panelD.R`: robustness: time-varying r_t (known / ignored),
  calibration jitter, and misspecified (QDA) DGP.
- `run-all.R`: master driver. `Rscript run-all.R quick` or `... full`.

## Provisional design choices (flagged for sign-off in the spec)

- Signal strength: Bayes AUC 0.70, so Delta = 0.742 and alpha_bar = 0.355.
- Base rate pi0 = 0.042 (82/2014 in the data), n_train = 2000, d = 30 (s = 5).
- Reserve observation noise tau = 3 percent.
- True alpha*_t path: two humps (0.5 -> 0.9 -> 0.45 -> 0.9 -> 0.4) over
  T = 25 periods, mirroring the shape of Figure 3 in the paper.
- Insurance calibration: benchmark midpoints; p from delta.op and pi.bar
  exactly as in `rho.star()`: p = 1 - delta.op/((1-pi.bar)(pi.bar+delta.op)).

## Notes and known deviations from the spec

- The umbrella threshold requires m >= log(delta)/log(1-alpha) left-out
  crisis observations; with n = 2000 and pi0 = 4.2 percent (about 42 left
  out), alpha = 0.05 is infeasible, so the alpha grids start at 0.10. The
  same constraint binds in the real data (82 crises).
- Panel C reports Monte Carlo bands rather than the spec's bootstrap-CI
  coverage (cheaper; bootstrap coverage can be added once the design is
  approved).
- Comparators share the SAME fitted scoring function and differ only in the
  threshold rule, isolating the umbrella's contribution.
