# AdaBoost bootstrap on NYU HPC (Torch or Greene)

Computes the 500-replication bootstrap band for the AdaBoost robustness
recalibration of the risk tolerance (Section 5.2 robustness / response-letter
Figure R.1). Fully self-contained: all data and configuration are frozen in
`boot_inputs.rds`, so no Stata files or project setup are needed on the
cluster. Each replication is a single-core job of about 45-60 minutes; as a
500-task array the whole band completes in roughly one queue-limited hour.

## Bundle contents

- `boot_inputs.rds`: frozen inputs (data, calibration table, grids, i.model = 5).
- `helper_funcs_hpc.R`: calibration functions with a minimal library block.
- `run_one_boot_rep.R`: one replication per SLURM array task.
- `merge_boot_reps.R`: combines `reps/rep_*.rds` into the CI objects.
- `submit_boot.sbatch`: the array job (adjust the module/account lines).

## Steps

1. Copy the bundle to the cluster:

```bash
scp -r hpc_boot <netid>@login.torch.hpc.nyu.edu:/scratch/<netid>/np-boot
```

2. On the login node, install the five required packages once (any recent R
   module works; adjust the module name after checking `module avail r`):

```bash
module load r  # check the exact name with: module avail r
Rscript -e 'install.packages(c("nproc","ada","imputeMissings","purrr","data.table"), repos="https://cloud.r-project.org")'
```

3. Plumbing test (about 1 minute; runs replication 1 on a reduced grid):

```bash
cd /scratch/<netid>/np-boot
QUICK=1 Rscript run_one_boot_rep.R 1 && rm reps/rep_001.rds
```

4. Submit the array:

```bash
mkdir -p logs reps
sbatch submit_boot.sbatch
```

5. When all tasks finish (`squeue -u <netid>` empty; expect 500 files in
   `reps/`), merge and copy the result back:

```bash
Rscript merge_boot_reps.R
```

```bash
scp <netid>@login.torch.hpc.nyu.edu:/scratch/<netid>/np-boot/alpha_calibration_bootstrap_ci.* .
```

6. Hand `alpha_calibration_bootstrap_ci.rds` back to Claude (drop it in
   `np-crisis-code/outputs_r2_noReserves/rds/`): the ada band then gets folded
   into the response-letter Figure R.1 and the repo package snapshot.

## Notes

- Seeding is per replication (`set.seed(123000 + rep_id)`), a documented
  deviation from the in-Rmd stream seeding; statistically equivalent for
  bootstrap purposes.
- Missing tasks are fine to resubmit individually, e.g.
  `sbatch --array=17,243 submit_boot.sbatch`; the merge script uses whatever
  `rep_*.rds` files exist.
- The aggregation in `merge_boot_reps.R` replicates the paper pipeline
  exactly (pool rows, group by year, 90 percent band + sd, year shift), so
  the result is directly comparable to the signal-extraction band shipped
  with the paper.
