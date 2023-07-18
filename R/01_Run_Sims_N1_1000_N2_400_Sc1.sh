#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH --mem=2g
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH --mail-type=END,FAIL,REQUEUE,STAGE_OUT
#SBATCH --mail-user=EMAIL_ADDRESS
module load r/4.1.3

#define variables for number of sims and subdirectory name where results are stored
N_SIMS=2000
SUBDIR_NAME=N1_1000_N2_400_SC1

sbatch --output=/dev/null --error=/dev/null --time=2:00:00 --array=1-$N_SIMS --job-name=01_N1_1000_N2_400_SC1 --wait R CMD BATCH --no-save --no-restore "--args $SUBDIR_NAME" 01_Analyze_SS_MS.R 01_Analyze_SS_MS_N1_1000_N2_400_SC1.Rout






