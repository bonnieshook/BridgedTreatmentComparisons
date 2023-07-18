############################################################################################################
# Fusing Trial Data for Treatment Comparisons: Single versus Multi-Span Bridging
#   Helper function(s) to compute simulation performance metrics
#
# Paul Zivich (2023/04/27)
############################################################################################################

# Importing dependencies
import numpy as np


def calculate_metrics(scenario, estimator, bias, se, coverage, bias_diagnostic=None, coverage_diagnostic=None):
    """Calculate the performance metrics for the simulation experiment

    Parameters
    ----------
    scenario : int, str
        Label for the scenario
    estimator : str
        Label for the input estimator
    bias : ndarray
        Array of bias estimates
    se : ndarray
        Array of standard error estimates
    coverage : ndarray
        Array of whether confidence interval covers the truth
    bias_diagnostic : ndarray, None
        Optional array of bias estimates for the diagnostic
    coverage_diagnostic : ndarray, None
        Optional array of whether confidence interval covers zero for the diagnostic

    Returns
    -------
    list
    """
    b = np.nanmean(bias)             # Compute bias
    ase = np.nanmean(se)             # Compute average standard error
    ese = np.nanstd(bias, ddof=1)    # Compute empirical standard error
    ser = ase /ese                   # Compute standard error ratio
    c = np.nanmean(coverage)         # Compute confidence interval coverage
    rmse = np.sqrt(b**2 + ese**2)    # Compute root mean squared error

    # Return metrics based on whether diagnostic is also calculated
    if bias_diagnostic is None:
        return scenario, estimator, b, ase, ese, ser, rmse, c, np.nan, np.nan
    else:
        bd = np.nanmean(bias_diagnostic)         # Compute difference from zero for diagnostic
        cd = np.nanmean(coverage_diagnostic)     # Compute coverage of zero for diagnostics
        return scenario, estimator, b, ase, ese, ser, rmse, c, bd, cd
