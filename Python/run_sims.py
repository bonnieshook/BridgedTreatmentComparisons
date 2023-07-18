#####################################################################################################################
# Fusing Trial Data for Treatment Comparisons: Single versus Multi-Span Bridging
#   Runs the simulation replication, then prints and saves all results to tables
#
# Paul Zivich (2023/04/27)
#####################################################################################################################

# Importing dependencies
import warnings
import numpy as np
import pandas as pd
from delicatessen import MEstimator

from dgm import generate_data, calculate_truth
from efuncs import ee_naive_ms, ee_naive_ss, ee_bridge_ss, ee_bridge_ms
from metrics import calculate_metrics

########################################################################
# Setting up simulation meta information
np.random.seed(7777777)                                     # Random seed
scenarios = [1, 2, 3, 4, 5]                                 # Scenarios to consider
sims = 2000                                                 # Number of iteration per scenario
n_pairs = (400, 1000), (1000, 400), (2000, 1000)            # Pairs of sample sizes to consider
pd.set_option('display.max_columns', None)                  # Have all columns displayed in prints to Console
warnings.filterwarnings("ignore", category=RuntimeWarning)  # Ignore RuntimeWarnings (divide by zero in root-finding)


########################################################################
# Defining estimating functions to pass to delicatessen

def psi_naive_ms(theta):
    return ee_naive_ms(theta=theta,
                       y=y, a=a, r=r, m=m)


def psi_naive_ss(theta):
    return ee_naive_ss(theta=theta,
                       y=y, a=a, r=r, m=m)


def psi_bridge_ss(theta):
    return ee_bridge_ss(theta=theta,
                        y=y, a=a, r=r, m=m,
                        W=W, V3=V, V1=V)


def psi_bridge_ms(theta):
    return ee_bridge_ms(theta=theta,
                        y=y, a=a, r=r, m=m,
                        W=W, V3=V, V2a=V, V2b=V, V1=V)


########################################################################
# Running simulation for outcome type, sample sizes, and scenarios
for outcome_type in [True, False]:                                          # Go through continuous & binary outcomes
    continuous = outcome_type                                               # Store whether we are looking at continuous
    if continuous:                                                          # Starting values for root-finding
        ss_init = [60., 300., 240.]                                         # ... putting continuous starting values
        ms_init = [60., 0., 300., 270., 270., 240.]                         # ... in the general vicinity
    else:
        ss_init = [0.5, 0.9, 0.4]                                           # ... putting binary starting values
        ms_init = [0.5, 0.0, 0.9, 0.65, 0.65, 0.4]                          # ... in the general vicinity

    # Go through each pair of N's provided
    for ns in n_pairs:                                                      # For each outcome type, look at N combos
        n1, n0 = ns                                                         # Extract n2, n1 for simulation
        scenario_result = []                                                # Storage for all scenario results

        for scenario in scenarios:                                          # Go through each scenario (1-5)
            # Setting up for the corresponding scenario
            truth = calculate_truth(n=20000000, scenario=scenario)          # Computing truth with 20 mil observations
            if continuous:                                                  # Extract depending on outcome type
                truth = truth[0]
            else:
                truth = truth[1]

            # Storage for scenario-specific results
            nss_est, nss_var, nss_cov = [], [], []                          # Naive SS results per iteration
            nms_est, nms_var, nms_cov = [], [], []                          # Naive MS results per iteration
            nms_diag_est, nms_diag_cov = [], []                             # Naive MS Diagnostic results per iteration
            bss_est, bss_var, bss_cov = [], [], []                          # Bridge SS results per iteration
            bms_est, bms_var, bms_cov = [], [], []                          # Bridge MS results per iteration
            bms_diag_est, bms_diag_cov = [], []                             # Bridge MS diagnostic results per iteration

            # Running simulations for that scenario
            for i in range(sims):                                           # Run for each iteration
                d = generate_data(n1=n1, n0=n0, scenario=scenario)          # Generate trial data for scenario
                r = np.asarray(d['S'])                                      # Sampling indicator
                a = np.asarray(d['A'])                                      # Treatment arm
                W = np.asarray(d[['C', 'X1', 'X2']])                        # Sampling model design matrix
                V = np.asarray(d[['C', 'X1']])                              # Missing model design matrix
                m = np.asarray(1 - d['M'])                                  # Missing indicator
                if continuous:                                              # Outcome depends on outcome type
                    y = np.asarray(d['Y'])
                else:
                    y = np.asarray(d['B'])

                # Naive -- Single-Span
                estr = MEstimator(psi_naive_ss, init=ss_init)               # M-estimator for Naive SS
                estr.estimate(solver="lm", maxiter=20000)                   # Computing parameters
                nss_est.append(estr.theta[0] - truth)                       # Storing bias
                nss_var.append(estr.variance[0, 0])                         # Storing variance
                ci = estr.confidence_intervals()                            # Extract confidence intervals
                if ci[0, 0] < truth < ci[0, 1]:                             # Store confidence interval results
                    nss_cov.append(1)
                else:
                    nss_cov.append(0)

                # Naive -- Multi-Span
                estr = MEstimator(psi_naive_ms, init=ms_init)               # M-estimator for Naive MS
                estr.estimate(solver="lm", maxiter=20000)                   # Computing parameters
                nms_est.append(estr.theta[0] - truth)                       # Storing bias
                nms_var.append(estr.variance[0, 0])                         # Storing variance
                ci = estr.confidence_intervals()                            # Extract confidence intervals
                if ci[0, 0] < truth < ci[0, 1]:                             # Store confidence intervals
                    nms_cov.append(1)
                else:
                    nms_cov.append(0)
                nms_diag_est.append(estr.theta[1] - 0)                      # Storing diagnostic
                if ci[1, 0] < 0 < ci[1, 1]:                                 # Extract diagnostic confidence intervals
                    nms_diag_cov.append(1)                                  # Store diagnostic confidence intervals
                else:
                    nms_diag_cov.append(0)

                # Bridge -- Single-Span
                inits = ss_init + [
                         0., 1., 0.,
                         1.5, -0.5,
                         1.5, -0.5
                         ]
                estr = MEstimator(psi_bridge_ss, init=inits)                # M-estimator for Bridge SS
                try:                                                        # Catching any issues...
                    estr.estimate(solver="lm", maxiter=20000)               # Compute parameters
                    bss_est.append(estr.theta[0] - truth)                   # Storing bias
                    bss_var.append(estr.variance[0, 0])                     # Storing variance
                    ci = estr.confidence_intervals()                        # Extract confidence intervals
                    if ci[0, 0] < truth < ci[0, 1]:                         # Store confidence intervals
                        bss_cov.append(1)
                    else:
                        bss_cov.append(0)
                except RuntimeError:                                        # Fail to root-find, output NaN
                    bss_est.append(np.nan)
                    bss_var.append(np.nan)
                    bss_cov.append(np.nan)

                # Bridge -- Multi-Span
                inits = ms_init + [
                         0., 1., 0.,
                         1.5, -0.5,
                         1.5, -0.5,
                         1.5, -0.5,
                         1.5, -0.5
                         ]
                estr = MEstimator(psi_bridge_ms, init=inits)                # M-estimator for Bridge MS
                try:                                                        # Catching any issues...
                    estr.estimate(solver="lm", maxiter=20000)               # Compute parameters
                    bms_est.append(estr.theta[0] - truth)                   # Storing bias
                    bms_var.append(estr.variance[0, 0])                     # Storing variance
                    ci = estr.confidence_intervals()                        # Extract confidence intervals
                    if ci[0, 0] < truth < ci[0, 1]:                         # Store confidence intervals
                        bms_cov.append(1)
                    else:
                        bms_cov.append(0)
                    bms_diag_est.append(estr.theta[1] - 0)                  # Storing diagnostic
                    if ci[1, 0] < 0 < ci[1, 1]:                             # Extract diagnostic confidence intervals
                        bms_diag_cov.append(1)                              # Store diagnostic confidence intervals
                    else:
                        bms_diag_cov.append(0)
                except RuntimeError:                                        # Fail to root-find, output NaN
                    bms_est.append(np.nan)
                    bms_var.append(np.nan)
                    bms_cov.append(np.nan)
                    bms_diag_est.append(np.nan)
                    bms_diag_cov.append(np.nan)

            # Calculating metrics by estimator for each scenario
            nms = calculate_metrics(scenario=scenario, estimator='Naive MS',
                                    bias=nms_est, se=np.sqrt(nms_var), coverage=nms_cov,
                                    bias_diagnostic=nms_diag_est, coverage_diagnostic=nms_diag_cov)
            scenario_result.append(nms)
            nss = calculate_metrics(scenario=scenario, estimator='Naive SS',
                                    bias=nss_est, se=np.sqrt(nss_var), coverage=nss_cov,
                                    bias_diagnostic=None, coverage_diagnostic=None)
            scenario_result.append(nss)
            bms = calculate_metrics(scenario=scenario, estimator='Bridge MS',
                                    bias=bms_est, se=np.sqrt(bms_var), coverage=bms_cov,
                                    bias_diagnostic=bms_diag_est, coverage_diagnostic=bms_diag_cov)
            scenario_result.append(bms)
            bss = calculate_metrics(scenario=scenario, estimator='Bridge SS',
                                    bias=bss_est, se=np.sqrt(bss_var), coverage=bss_cov,
                                    bias_diagnostic=None, coverage_diagnostic=None)
            scenario_result.append(bss)

        # Converting results into a dataframe
        cols = ["Scenario", "Estimator", "Bias", "ASE", "ESE", "SER", "RMSE", "C", "Diag", "DiagC"]
        results = pd.DataFrame(scenario_result, columns=cols)
        results = results.set_index("Scenario")
        print("==================================")
        print("N_1:       ", n1)
        print("N_0:       ", n0)
        print("ITERATIONS:", sims)
        print("CONTINUOUS:", continuous)
        print("==================================")
        if not continuous:
            results["Bias"] *= 100
            results["ASE"] *= 100
            results["ESE"] *= 100
            results["RMSE"] *= 100
            results["Diag"] *= 100
        print(results.round(2))
        print("\n")

        # Saving results for outcome type and sample size as .csv
        if continuous:
            ctype = "c"
        else:
            ctype = "b"
        file_name = ctype+"_n"+str(n1)+"n"+str(n0)+".csv"
        results.to_csv("results/" + file_name, index=True)


# Output
#
# ==================================
# N_1:        400
# N_0:        1000
# ITERATIONS: 2000
# CONTINUOUS: True
# ==================================
#           Estimator   Bias   ASE   ESE   SER   RMSE     C   Diag  DiagC
# Scenario
# 1          Naive MS   0.16  5.03  5.07  0.99   5.07  0.95  -0.13   0.94
# 1          Naive SS   0.03  3.55  3.52  1.01   3.52  0.95    NaN    NaN
# 1         Bridge MS   0.16  5.02  5.07  0.99   5.07  0.95  -0.10   0.95
# 1         Bridge SS  -0.01  2.17  2.18  0.99   2.18  0.95    NaN    NaN
# 2          Naive MS  17.75  5.66  5.61  1.01  18.61  0.12  20.35   0.00
# 2          Naive SS  38.10  4.31  4.35  0.99  38.35  0.00    NaN    NaN
# 2         Bridge MS  -0.16  6.01  5.92  1.02   5.92  0.95  -0.01   0.95
# 2         Bridge SS  -0.19  3.66  3.66  1.00   3.66  0.94    NaN    NaN
# 3          Naive MS  37.13  5.49  5.31  1.03  37.51  0.00   1.01   0.94
# 3          Naive SS  38.14  4.30  4.26  1.01  38.38  0.00    NaN    NaN
# 3         Bridge MS  19.06  5.76  5.73  1.00  19.91  0.09 -19.04   0.00
# 3         Bridge SS   0.06  3.66  3.67  1.00   3.67  0.95    NaN    NaN
# 4          Naive MS  47.53  5.65  5.47  1.03  47.84  0.00  20.57   0.00
# 4          Naive SS  68.10  4.30  4.24  1.02  68.23  0.00    NaN    NaN
# 4         Bridge MS  29.65  6.01  5.89  1.02  30.23  0.00   0.19   0.95
# 4         Bridge SS  29.91  3.66  3.67  1.00  30.13  0.00    NaN    NaN
# 5          Naive MS  57.68  5.66  5.60  1.01  57.95  0.00  10.27   0.21
# 5          Naive SS  67.95  4.31  4.26  1.01  68.09  0.00    NaN    NaN
# 5         Bridge MS  39.87  6.00  5.97  1.01  40.31  0.00  -9.98   0.14
# 5         Bridge SS  29.83  3.66  3.71  0.99  30.06  0.00    NaN    NaN
#
# ==================================
# N_1:        1000
# N_0:        400
# ITERATIONS: 2000
# CONTINUOUS: True
# ==================================
#           Estimator   Bias   ASE   ESE   SER   RMSE     C   Diag  DiagC
# Scenario
# 1          Naive MS  -0.03  5.02  5.04  1.00   5.04  0.95   0.05   0.94
# 1          Naive SS   0.02  3.55  3.61  0.98   3.61  0.95    NaN    NaN
# 1         Bridge MS  -0.03  5.02  5.04  1.00   5.04  0.95   0.02   0.95
# 1         Bridge SS   0.01  2.20  2.24  0.98   2.24  0.94    NaN    NaN
# 2          Naive MS  17.68  6.34  6.26  1.01  18.75  0.20  20.40   0.00
# 2          Naive SS  38.08  5.11  5.13  1.00  38.42  0.00    NaN    NaN
# 2         Bridge MS  -0.09  6.56  6.54  1.00   6.54  0.95   0.12   0.96
# 2         Bridge SS   0.06  3.62  3.63  1.00   3.63  0.95    NaN    NaN
# 3          Naive MS  36.92  6.28  6.31  1.00  37.45  0.00   1.09   0.94
# 3          Naive SS  38.01  5.12  5.16  0.99  38.36  0.00    NaN    NaN
# 3         Bridge MS  18.75  6.45  6.60  0.98  19.88  0.18 -18.83   0.00
# 3         Bridge SS  -0.07  3.62  3.71  0.97   3.71  0.94    NaN    NaN
# 4          Naive MS  47.96  6.32  6.41  0.99  48.39  0.00  20.26   0.00
# 4          Naive SS  68.22  5.10  5.35  0.95  68.43  0.00    NaN    NaN
# 4         Bridge MS  30.16  6.55  6.63  0.99  30.88  0.00   0.03   0.94
# 4         Bridge SS  29.96  3.64  3.59  1.01  30.18  0.00    NaN    NaN
# 5          Naive MS  57.90  6.33  6.51  0.97  58.26  0.00  10.26   0.21
# 5          Naive SS  68.16  5.11  5.16  0.99  68.35  0.00    NaN    NaN
# 5         Bridge MS  40.10  6.55  6.56  1.00  40.64  0.00  -9.89   0.23
# 5         Bridge SS  30.10  3.63  3.73  0.97  30.33  0.00    NaN    NaN
#
# ==================================
# N_1:        2000
# N_0:        1000
# ITERATIONS: 2000
# CONTINUOUS: True
# ==================================
#           Estimator   Bias   ASE   ESE   SER   RMSE     C   Diag  DiagC
# Scenario
# 1          Naive MS   0.05  3.29  3.26  1.01   3.26  0.96  -0.03   0.96
# 1          Naive SS   0.03  2.33  2.29  1.02   2.29  0.96    NaN    NaN
# 1         Bridge MS   0.05  3.29  3.25  1.01   3.25  0.96  -0.04   0.96
# 1         Bridge SS   0.02  1.42  1.37  1.04   1.37  0.96    NaN    NaN
# 2          Naive MS  17.85  4.10  4.06  1.01  18.31  0.01  20.31   0.00
# 2          Naive SS  38.16  3.29  3.30  1.00  38.30  0.00    NaN    NaN
# 2         Bridge MS   0.08  4.26  4.22  1.01   4.23  0.95   0.03   0.95
# 2         Bridge SS   0.08  2.36  2.42  0.98   2.42  0.94    NaN    NaN
# 3          Naive MS  37.21  4.06  4.07  1.00  37.43  0.00   1.01   0.92
# 3          Naive SS  38.22  3.29  3.36  0.98  38.37  0.00    NaN    NaN
# 3         Bridge MS  19.15  4.19  4.11  1.02  19.58  0.00 -19.08   0.00
# 3         Bridge SS  -0.04  2.36  2.29  1.03   2.29  0.96    NaN    NaN
# 4          Naive MS  47.88  4.10  4.05  1.01  48.05  0.00  20.27   0.00
# 4          Naive SS  68.15  3.29  3.25  1.01  68.23  0.00    NaN    NaN
# 4         Bridge MS  30.09  4.26  4.26  1.00  30.39  0.00  -0.04   0.95
# 4         Bridge SS  29.96  2.36  2.35  1.00  30.05  0.00    NaN    NaN
# 5          Naive MS  57.85  4.10  4.17  0.98  58.00  0.00  10.37   0.01
# 5          Naive SS  68.21  3.29  3.32  0.99  68.29  0.00    NaN    NaN
# 5         Bridge MS  40.03  4.27  4.39  0.97  40.27  0.00 -10.01   0.02
# 5         Bridge SS  29.93  2.37  2.44  0.97  30.02  0.00    NaN    NaN
#
# ==================================
# N_1:        400
# N_0:        1000
# ITERATIONS: 2000
# CONTINUOUS: False
# ==================================
#           Estimator   Bias   ASE   ESE   SER   RMSE     C   Diag  DiagC
# Scenario
# 1          Naive MS  -0.03  5.39  5.39  1.00   5.39  0.96  -0.03   0.96
# 1          Naive SS  -0.07  3.35  3.37  0.99   3.37  0.95    NaN    NaN
# 1         Bridge MS  -0.04  5.41  5.40  1.00   5.40  0.95   0.09   0.95
# 1         Bridge SS   0.09  3.04  3.02  1.01   3.02  0.95    NaN    NaN
# 2          Naive MS  -4.37  5.46  5.48  1.00   7.01  0.88  13.65   0.05
# 2          Naive SS   9.27  3.88  3.93  0.99  10.07  0.35    NaN    NaN
# 2         Bridge MS   0.16  6.47  6.45  1.00   6.45  0.95  -0.16   0.95
# 2         Bridge SS  -0.04  3.81  3.88  0.98   3.88  0.94    NaN    NaN
# 3          Naive MS  11.83  5.43  5.54  0.98  13.06  0.41  -2.45   0.88
# 3          Naive SS   9.38  3.88  3.90  0.99  10.16  0.33    NaN    NaN
# 3         Bridge MS  19.89  6.40  6.56  0.98  20.95  0.13 -19.78   0.00
# 3         Bridge SS   0.10  3.81  3.90  0.98   3.90  0.94    NaN    NaN
# 4          Naive MS   6.82  5.23  5.29  0.99   8.63  0.74  13.74   0.04
# 4          Naive SS  20.56  3.55  3.50  1.01  20.86  0.00    NaN    NaN
# 4         Bridge MS  16.46  6.11  6.14  0.99  17.57  0.23  -0.03   0.95
# 4         Bridge SS  16.44  3.56  3.50  1.02  16.81  0.00    NaN    NaN
# 5          Naive MS  13.51  5.31  5.12  1.04  14.44  0.28   7.33   0.55
# 5          Naive SS  20.84  3.54  3.51  1.01  21.13  0.00    NaN    NaN
# 5         Bridge MS  25.73  6.14  6.00  1.02  26.42  0.01  -9.07   0.40
# 5         Bridge SS  16.63  3.56  3.57  1.00  17.00  0.00    NaN    NaN
#
# ==================================
# N_1:        1000
# N_0:        400
# ITERATIONS: 2000
# CONTINUOUS: False
# ==================================
#           Estimator   Bias   ASE   ESE   SER   RMSE     C   Diag  DiagC
# Scenario
# 1          Naive MS   0.04  5.83  5.93  0.98   5.93  0.95   0.05   0.94
# 1          Naive SS   0.09  4.02  4.14  0.97   4.14  0.94    NaN    NaN
# 1         Bridge MS   0.07  5.84  5.92  0.99   5.92  0.95   0.00   0.94
# 1         Bridge SS   0.03  3.31  3.38  0.98   3.38  0.94    NaN    NaN
# 2          Naive MS  -4.46  4.63  4.59  1.01   6.40  0.84  13.78   0.03
# 2          Naive SS   9.32  3.03  3.02  1.00   9.80  0.15    NaN    NaN
# 2         Bridge MS  -0.09  7.10  7.17  0.99   7.17  0.94   0.11   0.94
# 2         Bridge SS   0.04  4.30  4.44  0.97   4.44  0.94    NaN    NaN
# 3          Naive MS  11.62  4.93  4.84  1.02  12.59  0.35  -2.36   0.92
# 3          Naive SS   9.25  3.03  3.01  1.01   9.73  0.16    NaN    NaN
# 3         Bridge MS  19.71  7.10  7.13  1.00  20.96  0.23 -19.75   0.02
# 3         Bridge SS  -0.07  4.33  4.30  1.01   4.30  0.94    NaN    NaN
# 4          Naive MS   6.84  4.31  4.31  1.00   8.09  0.64  13.79   0.03
# 4          Naive SS  20.62  2.50  2.49  1.00  20.77  0.00    NaN    NaN
# 4         Bridge MS  16.31  6.39  6.50  0.98  17.56  0.29   0.18   0.95
# 4         Bridge SS  16.52  3.48  3.57  0.97  16.90  0.01    NaN    NaN
# 5          Naive MS  13.39  4.56  4.55  1.00  14.14  0.16   7.47   0.51
# 5          Naive SS  20.86  2.49  2.54  0.98  21.02  0.00    NaN    NaN
# 5         Bridge MS  25.65  6.44  6.60  0.98  26.49  0.04  -8.92   0.53
# 5         Bridge SS  16.72  3.50  3.66  0.96  17.11  0.01    NaN    NaN
#
# ==================================
# N_1:        2000
# N_0:        1000
# ITERATIONS: 2000
# CONTINUOUS: False
# ==================================
#           Estimator   Bias   ASE   ESE   SER   RMSE     C   Diag  DiagC
# Scenario
# 1          Naive MS  -0.07  3.79  3.91  0.97   3.91  0.95   0.05   0.95
# 1          Naive SS  -0.02  2.59  2.57  1.01   2.57  0.95    NaN    NaN
# 1         Bridge MS  -0.07  3.79  3.92  0.97   3.92  0.95   0.01   0.95
# 1         Bridge SS  -0.02  2.14  2.15  1.00   2.15  0.95    NaN    NaN
# 2          Naive MS  -4.36  3.10  3.07  1.01   5.33  0.72  13.72   0.00
# 2          Naive SS   9.36  2.06  2.00  1.03   9.57  0.01    NaN    NaN
# 2         Bridge MS   0.10  4.62  4.59  1.01   4.59  0.95  -0.05   0.95
# 2         Bridge SS   0.05  2.78  2.72  1.02   2.72  0.95    NaN    NaN
# 3          Naive MS  11.69  3.27  3.28  1.00  12.14  0.06  -2.38   0.85
# 3          Naive SS   9.31  2.06  2.10  0.98   9.54  0.01    NaN    NaN
# 3         Bridge MS  19.72  4.60  4.66  0.99  20.26  0.01 -19.65   0.00
# 3         Bridge SS   0.10  2.78  2.82  0.99   2.82  0.95    NaN    NaN
# 4          Naive MS   6.82  2.89  2.89  1.00   7.41  0.35  13.81   0.00
# 4          Naive SS  20.62  1.73  1.70  1.02  20.69  0.00    NaN    NaN
# 4         Bridge MS  16.29  4.19  4.22  0.99  16.83  0.03   0.11   0.95
# 4         Bridge SS  16.43  2.31  2.36  0.98  16.60  0.00    NaN    NaN
# 5          Naive MS  13.38  3.05  3.03  1.01  13.72  0.01   7.39   0.17
# 5          Naive SS  20.78  1.72  1.70  1.01  20.85  0.00    NaN    NaN
# 5         Bridge MS  25.63  4.21  4.21  1.00  25.97  0.00  -8.99   0.20
# 5         Bridge SS  16.65  2.31  2.30  1.00  16.80  0.00    NaN    NaN
