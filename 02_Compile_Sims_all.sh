#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH --mem=1g
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH --mail-type=END,FAIL,REQUEUE,STAGE_OUT
#SBATCH --mail-user=ADDRESS@email.unc.edu
module load r/4.1.3

#set number of sims once
N_SIMS=2000


##### N1=1000, N2=400 SIMS
SUBDIR_NAME=N1_1000_N2_400_SC1
sbatch --output=/dev/null --error=/dev/null --time=2:00:00 --job-name=02_Compile_N1_1000_N2_400_Sc1 --wait R CMD BATCH --no-save --no-restore "--args $N_SIMS $SUBDIR_NAME" 02_compile_sims.r 02_Compile_N1_1000_N2_400_Sc1.Rout

SUBDIR_NAME=N1_1000_N2_400_SC2
sbatch --output=/dev/null --error=/dev/null --time=2:00:00 --job-name=02_Compile_N1_1000_N2_400_Sc2 --wait R CMD BATCH --no-save --no-restore "--args $N_SIMS $SUBDIR_NAME" 02_compile_sims.r 02_Compile_N1_1000_N2_400_Sc2.Rout

SUBDIR_NAME=N1_1000_N2_400_SC3
sbatch --output=/dev/null --error=/dev/null --time=2:00:00 --job-name=02_Compile_N1_1000_N2_400_Sc3 --wait R CMD BATCH --no-save --no-restore "--args $N_SIMS $SUBDIR_NAME" 02_compile_sims.r 02_Compile_N1_1000_N2_400_Sc3.Rout

SUBDIR_NAME=N1_1000_N2_400_SC4
sbatch --output=/dev/null --error=/dev/null --time=2:00:00 --job-name=02_Compile_N1_1000_N2_400_Sc4 --wait R CMD BATCH --no-save --no-restore "--args $N_SIMS $SUBDIR_NAME" 02_compile_sims.r 02_Compile_N1_1000_N2_400_Sc4.Rout

SUBDIR_NAME=N1_1000_N2_400_SC5
sbatch --output=/dev/null --error=/dev/null --time=2:00:00 --job-name=02_Compile_N1_1000_N2_400_Sc5 --wait R CMD BATCH --no-save --no-restore "--args $N_SIMS $SUBDIR_NAME" 02_compile_sims.r 02_Compile_N1_1000_N2_400_Sc5.Rout

##### N1=1000, N2=2000 SIMS
SUBDIR_NAME=N1_1000_N2_2000_SC1
sbatch --output=/dev/null --error=/dev/null --time=2:00:00 --job-name=02_Compile_N1_1000_N2_2000_Sc1 --wait R CMD BATCH --no-save --no-restore "--args $N_SIMS $SUBDIR_NAME" 02_compile_sims.r 02_Compile_N1_1000_N2_2000_Sc1.Rout

SUBDIR_NAME=N1_1000_N2_2000_SC2
sbatch --output=/dev/null --error=/dev/null --time=2:00:00 --job-name=02_Compile_N1_1000_N2_2000_Sc2 --wait R CMD BATCH --no-save --no-restore "--args $N_SIMS $SUBDIR_NAME" 02_compile_sims.r 02_Compile_N1_1000_N2_2000_Sc2.Rout

SUBDIR_NAME=N1_1000_N2_2000_SC3
sbatch --output=/dev/null --error=/dev/null --time=2:00:00 --job-name=02_Compile_N1_1000_N2_2000_Sc3 --wait R CMD BATCH --no-save --no-restore "--args $N_SIMS $SUBDIR_NAME" 02_compile_sims.r 02_Compile_N1_1000_N2_2000_Sc3.Rout

SUBDIR_NAME=N1_1000_N2_2000_SC4
sbatch --output=/dev/null --error=/dev/null --time=2:00:00 --job-name=02_Compile_N1_1000_N2_2000_Sc4 --wait R CMD BATCH --no-save --no-restore "--args $N_SIMS $SUBDIR_NAME" 02_compile_sims.r 02_Compile_N1_1000_N2_2000_Sc4.Rout

SUBDIR_NAME=N1_1000_N2_2000_SC5
sbatch --output=/dev/null --error=/dev/null --time=2:00:00 --job-name=02_Compile_N1_1000_N2_2000_Sc5 --wait R CMD BATCH --no-save --no-restore "--args $N_SIMS $SUBDIR_NAME" 02_compile_sims.r 02_Compile_N1_1000_N2_2000_Sc5.Rout


##### N1=400, N2=1000 SIMS
SUBDIR_NAME=N1_400_N2_1000_SC1
sbatch --output=/dev/null --error=/dev/null --time=2:00:00 --job-name=02_Compile_N1_400_N2_1000_Sc1 --wait R CMD BATCH --no-save --no-restore "--args $N_SIMS $SUBDIR_NAME" 02_compile_sims.r 02_Compile_N1_400_N2_1000_Sc1.Rout

SUBDIR_NAME=N1_400_N2_1000_SC2
sbatch --output=/dev/null --error=/dev/null --time=2:00:00 --job-name=02_Compile_N1_400_N2_1000_Sc2 --wait R CMD BATCH --no-save --no-restore "--args $N_SIMS $SUBDIR_NAME" 02_compile_sims.r 02_Compile_N1_400_N2_1000_Sc2.Rout

SUBDIR_NAME=N1_400_N2_1000_SC3
sbatch --output=/dev/null --error=/dev/null --time=2:00:00 --job-name=02_Compile_N1_400_N2_1000_Sc3 --wait R CMD BATCH --no-save --no-restore "--args $N_SIMS $SUBDIR_NAME" 02_compile_sims.r 02_Compile_N1_400_N2_1000_Sc3.Rout

SUBDIR_NAME=N1_400_N2_1000_SC4
sbatch --output=/dev/null --error=/dev/null --time=2:00:00 --job-name=02_Compile_N1_400_N2_1000_Sc4 --wait R CMD BATCH --no-save --no-restore "--args $N_SIMS $SUBDIR_NAME" 02_compile_sims.r 02_Compile_N1_400_N2_1000_Sc4.Rout

SUBDIR_NAME=N1_400_N2_1000_SC5
sbatch --output=/dev/null --error=/dev/null --time=2:00:00 --job-name=02_Compile_N1_400_N2_1000_Sc5 --wait R CMD BATCH --no-save --no-restore "--args $N_SIMS $SUBDIR_NAME" 02_compile_sims.r 02_Compile_N1_400_N2_1000_Sc5.Rout





