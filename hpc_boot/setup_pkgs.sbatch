#!/bin/bash
# One-time setup job: installs the five R packages into $HOME/R/np-boot-lib
# and runs a QUICK plumbing replication. Submit this FIRST, check
# logs/setup.out for "rep 999 done: 66 rows", then submit submit_boot.sbatch.
#SBATCH --job-name=np-setup
#SBATCH --account=torch_pr_365_general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8GB
#SBATCH --time=00:40:00
#SBATCH --output=logs/setup.out

export OMP_NUM_THREADS=1
module purge
module load r/4.5.1
export R_LIBS_USER="$HOME/R/np-boot-lib"
mkdir -p "$R_LIBS_USER"
Rscript -e 'install.packages(c("nproc","ada","imputeMissings","purrr","data.table"), repos="https://cloud.r-project.org")'
cd "$SLURM_SUBMIT_DIR"
QUICK=1 Rscript run_one_boot_rep.R 999
