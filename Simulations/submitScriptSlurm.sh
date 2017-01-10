#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=16:00:00
#SBATCH --mem=5GB
trap "echo recieved SIGUSR1;" SIGUSR1;
R CMD BATCH --no-save --no-restore -- jobScript.R jobScript.Rout.$SLURM_JOBID
