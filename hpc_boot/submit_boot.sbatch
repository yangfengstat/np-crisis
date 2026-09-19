#!/bin/bash
# SLURM array job: 500 independent bootstrap replications, one per task.
# Each task is single-core, about 45-60 minutes.
# Adjust the module line (run: module avail r) and add your account or
# partition flags if your allocation requires them, e.g.:
#   #SBATCH --account=<your_account>
#   #SBATCH --partition=<partition>
#SBATCH --job-name=np-boot
#SBATCH --account=torch_pr_365_general
#SBATCH --array=1-500
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB
#SBATCH --time=02:00:00
#SBATCH --output=logs/rep_%a.out

export OMP_NUM_THREADS=1
module purge
module load r/4.5.1

export R_LIBS_USER="$HOME/R/np-boot-lib"
cd "$SLURM_SUBMIT_DIR"
mkdir -p logs reps
Rscript run_one_boot_rep.R "$SLURM_ARRAY_TASK_ID"
