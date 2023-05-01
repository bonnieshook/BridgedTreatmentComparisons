#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH --mem=2g
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH --mail-type=END,FAIL,REQUEUE,STAGE_OUT
#SBATCH --mail-user=EMAIL_ADDRESS
module load r/4.1.3

#set number of sims once
N_SIMS=2000

##### Enter sizes of non-focal (N1) and target (N2) populations, scenario (SCEN), and subdirectory name where results will be stored (SUBDIR_NAME)
N1=1000
N2=400

SCEN=1
SUBDIR_NAME=N1_1000_N2_400_SC1
sbatch --output=/dev/null --error=/dev/null --time=12:00:00 --job-name=00_DGM_090822_N1_1000_N2_400_Sc1 --wait R CMD BATCH --no-save --no-restore "--args $N_SIMS $SCEN $N1 $N2 $SUBDIR_NAME" 00_DGM_090822.R 00_DGM_090822_N1_1000_N2_400_Sc1.Rout








