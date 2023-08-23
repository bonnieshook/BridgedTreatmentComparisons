# Fusing Trial Data for Treatment Comparisons: Single versus Multi-Span Bridging

---

## R code for simulations 

These programs implement the simulation study described in Shook-Sa et al, "Fusing Trial Data for Treatment Comparisons: Single versus Multi-Span Bridging" in R

The simulations are implemented in the following order:
1. Obtain the true values of the ATE for each scenario empirically using the "00_DGM_091222_truth.R" program
2. Generate data for the simulation study for each sample size and scenario of interest using the 00_DGM_090822.R program, called with the (example) shell script: 00_DGM_shell_n1_1000_n2_400_sc1.sh. Note a separate shell script
   of this form was used for each sample size and scenario of interest.
3. Analyze the data from each generated sample with the SS, MS, and naive estimators presented in the manuscript using the 01_Analyze_SS_MS.R program, which calls the estimators from the 00_Estimators_01.31.23.R program.
   This code is run in parallel using the (example) shell script 01_Run_Sims_N1_1000_N2_400_Sc1.sh. For each scenario of interest, this program outputs a separate CSV file with the results for one iteration of the simulation.
   Note a separate shell script of this form was used for each sample size and scenario of interest.
4. Combine the results from all iterations of the simulation using the 02_compile_sims.R program, called with the shell script 02_Compile_Sims_all.sh. This program creates a final dataset with one row for each iteration of the 
   simulations.
5. Calculate summary measures for each scenario and sample size considered using 03_combine_sim_table_01.31.23.R.

Developed by: Bonnie Shook-Sa
