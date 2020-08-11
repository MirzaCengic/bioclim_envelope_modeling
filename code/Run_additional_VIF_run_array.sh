#!/bin/bash

#SBATCH --partition=milkun
#SBATCH --time=48:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=vif_ens
#SBATCH --mem-per-cpu=12G
#SBATCH --array=1-300
# Standard out and Standard Error output files with the job number in the name.
#SBATCH -o "/Logs/log_VIFs_%a.out"
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=m.cengic@science.ru.nl

export TMPDIR=/scratch/mdls
mkdir -p $TMPDIR

srun /opt/R-3.6.3/bin/R --vanilla --no-save --args ${SLURM_ARRAY_TASK_ID} < /Run_additional_VIF10_run_models.R
