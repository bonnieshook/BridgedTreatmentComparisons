############################################################################################################
# Fusing Trial Data for Treatment Comparisons: Single versus Multi-Span Bridging
#   Estimating functions for use in simulation experiments and applied example
#
# Paul Zivich (2023/04/27)
############################################################################################################

# Importing dependencies
import numpy as np
from delicatessen.estimating_equations import ee_regression
from delicatessen.utilities import inverse_logit


def nuisance_mean(theta, outcome, subset):
    """Function to calculate the complete-case mean, conditional to some subset of the data. This function is used to
    simplify the implementations of the Naive estimators.
    """
    return subset * (np.nan_to_num(outcome, copy=True, nan=0) - theta)


def ee_naive_ss(theta, y, a, r, m):
    """Stacked estimating functions for the Naive Single-Span estimator.

    Parameters
    ----------
    theta : ndarray, list
        Parameters to estimate
    y : ndarray, list
        Outcome (allows for binary or continuous)
    a : ndarray, list
        Assigned treatment arm
    r : ndarray, list
        Population / sample indicator, with 1 being the target and 0 being the non-focal
    m : ndarray, list
        Missing indicator for y

    Returns
    -------
    ndarray :
        Estimating functions for the ATE, mean in A=3 of R=2 trial, mean in A=1 of R=1 trial
    """
    # Estimating function for the risk difference
    riskdiff = np.ones(y.shape[0])*(theta[1] - theta[2]) - theta[0]

    # Estimating function for the mean in each trial arm
    yhat_a2s1 = nuisance_mean(theta=theta[1], outcome=y, subset=r*(a == 3)*m)
    yhat_a0s0 = nuisance_mean(theta=theta[2], outcome=y, subset=(1-r)*(a == 1)*m)

    # Returning stacked estimating functions
    return np.vstack((riskdiff,
                      yhat_a2s1, yhat_a0s0,))


def ee_naive_ms(theta, y, a, r, m):
    """Stacked estimating functions for the Naive Multi-Span estimator.

    Parameters
    ----------
    theta : ndarray, list
        Parameters to estimate
    y : ndarray, list
        Outcome (allows for binary or continuous)
    a : ndarray, list
        Assigned treatment arm
    r : ndarray, list
        Population / sample indicator, with 1 being the target and 0 being the non-focal
    m : ndarray, list
        Missing indicator for y

    Returns
    -------
    ndarray :
        Estimating functions for the ATE, diagnostic, mean in A=3 of R=2 trial, mean in A=2 of R=2 trial, mean in A=2
        of R=1 trial, mean in A=1 of R=1 trial
    """
    # Estimating function for the risk difference
    riskdiff = np.ones(y.shape[0])*((theta[2] - theta[3]) + (theta[4] - theta[5])) - theta[0]

    # Estimating function for the diagnostic
    diag = np.ones(y.shape[0])*(theta[3] - theta[4]) - theta[1]

    # Estimating functions for the mean in each trial arm
    yhat_a2s1 = nuisance_mean(theta=theta[2], outcome=y, subset=r*(a == 3)*m)
    yhat_a1s1 = nuisance_mean(theta=theta[3], outcome=y, subset=r*(a == 2)*m)
    yhat_a1s0 = nuisance_mean(theta=theta[4], outcome=y, subset=(1-r)*(a == 2)*m)
    yhat_a0s0 = nuisance_mean(theta=theta[5], outcome=y, subset=(1-r)*(a == 1)*m)

    # Returning stacked estimating functions
    return np.vstack((riskdiff, diag,
                      yhat_a2s1, yhat_a1s1, yhat_a1s0, yhat_a0s0,))


def ee_bridge_ss(theta, y, a, r, m, W, V3, V1):
    """Stacked estimating functions for the Bridge Single-Span estimator. This estimator assumes that randomization for
    both trials was 1:1.

    Parameters
    ----------
    theta : ndarray, list
        Parameters to estimate
    y : ndarray, list
        Outcome (allows for binary or continuous)
    a : ndarray, list
        Assigned treatment arm
    r : ndarray, list
        Population / sample indicator, with 1 being the target and 0 being the non-focal
    m : ndarray, list
        Missing indicator for y
    W : ndarray, list
        Design matrix for the selection / sampling model
    V3 : ndarray, list
        Design matrix for the missing model, restricted to A=3
    V1 : ndarray, list
        Design matrix for the missing model, restricted to A=1

    Returns
    -------
    ndarray :
        Estimating functions for the ATE, mean in A=3 of R=2 trial, mean in A=1 of R=1 trial, sampling nuisance model
        parameters, missing nuisance model for A=3 parameters, missing nuisance model for A=1 parameters
    """
    # Extracting parameters as named variables to ease coding
    w_len = W.shape[1]                               # Getting shape of W design matrix
    v3_len = V3.shape[1]                             # Getting shape of V3 design matrix
    interest_param = theta[0]                        # Average treatment effect
    r_a2_s1, r_a0_s0 = theta[1:3]                    # Mean in different arms of trials
    beta = theta[3:3+w_len]                          # Parameters for sampling model
    gamma_a2s1 = theta[3+w_len:3+w_len+v3_len]       # Parameters for missing model in R=2
    gamma_a0s0 = theta[3+w_len+v3_len:]              # Parameters for missing model in R=1
    iptw = 1 / 0.5                                   # Here, the propensity score is known

    # Calculate the sampling model, Pr(R | X; beta), nuisance parameters
    pr_s_nuisance = ee_regression(theta=beta,        # Using built-in regression
                                  X=W, y=1-r,        # ... given variables
                                  model='logistic')  # ... with logistic model
    pr_s_nuisance = pr_s_nuisance * (a != 2)         # Shared arm does not contribute to model
    pr_s1 = inverse_logit(np.dot(W, beta))           # Predicted probability of being R=2
    odds = pr_s1 / (1 - pr_s1)                       # Predicted odds of being R=2
    iosw = r*1 + (1-r)/odds                          # Inverse Odds of Sampling Weights

    # Calculate the missing models, Pr(M | X, R, A; gamma), nuisance parameters
    pr_m_nuisance_a2s1 = ee_regression(theta=gamma_a2s1,             # Using built-in regression
                                       X=V3, y=m,                    # ... given variables
                                       model='logistic')             # ... with logistic model
    pr_m_nuisance_a2s1 = pr_m_nuisance_a2s1 * (a == 3)               # Only A=3 arm contributes

    pr_m_nuisance_a0s0 = ee_regression(theta=gamma_a0s0,             # Using built-in regression
                                       X=V1, y=m,                    # ... given variables
                                       model='logistic')             # ... with logistic model
    pr_m_nuisance_a0s0 = pr_m_nuisance_a0s0 * (a == 1)               # Only single arm contributes

    pr_m = (inverse_logit(np.dot(V3, gamma_a2s1))*r*(a == 3)         # Compute all probability in R=2
            + 1*(a == 2)                                             # ... prevent NaN errors
            + inverse_logit(np.dot(V1, gamma_a0s0))*(1-r)*(a == 1))  # ... probability in R=1
    ipmw = np.where(m, 1/pr_m, 0)                                    # Inverse probability of missing weight

    # Estimating functions for the mean in each arm for estimator
    ipw = iptw * ipmw * iosw                                         # Overall inverse probability weight
    y_nan = np.nan_to_num(y, copy=True, nan=0.)                      # Convert missing y to 0 to prevent NaN errors
    yhat_a3s1 = r * m * (a == 3) * ipw * (y_nan - r_a2_s1)           # Risk in R=2, A=3 reweighted
    yhat_a1s0 = (1-r) * m * (a == 1) * ipw * (y_nan - r_a0_s0)       # Risk in R=1, A=1 reweighted

    # Estimating function for average treatment effect
    ate_brdg = np.ones(y.shape[0])*(r_a2_s1 - r_a0_s0) - interest_param

    # Returning stacked estimating functions
    return np.vstack((ate_brdg,
                      yhat_a3s1, yhat_a1s0,
                      pr_s_nuisance,
                      pr_m_nuisance_a2s1, pr_m_nuisance_a0s0
                      ))


def ee_bridge_ms(theta, y, a, r, m, W, V3, V2a, V2b, V1):
    """Stacked estimating functions for the Bridge Multi-Span estimator. This estimator assumes that randomization for
    both trials was 1:1.

    Parameters
    ----------
    theta : ndarray, list
        Parameters to estimate
    y : ndarray, list
        Outcome (allows for binary or continuous)
    a : ndarray, list
        Assigned treatment arm
    r : ndarray, list
        Population / sample indicator, with 1 being the target and 0 being the non-focal
    m : ndarray, list
        Missing indicator for y
    W : ndarray, list
        Design matrix for the selection / sampling model
    V3 : ndarray, list
        Design matrix for the missing model, restricted to A=3
    V2a : ndarray, list
        Design matrix for the missing model, restricted to A=2,R=2
    V2b : ndarray, list
        Design matrix for the missing model, restricted to A=2,R=1
    V1 : ndarray, list
        Design matrix for the missing model, restricted to A=1

    Returns
    -------
    ndarray :
        Estimating functions for the ATE, diagnostic, mean in A=3 of R=2 trial, mean in A=2 of R=2 trial, mean in A=2
        of R=1 trial, mean in A=1 of R=1 trial, sampling nuisance model parameters, missing nuisance model for A=3
        parameters, missing nuisance model for A=2,R=2 parameters, missing nuisance model for A=2,R=1 parameters,
        missing nuisance model for A=1 parameters
    """
    # Extracting parameters as named variables to ease coding
    w_len = W.shape[1]                                                          # Getting shape of W design matrix
    v3_len, v2a_len, v2b_len = V3.shape[1], V2a.shape[1], V2b.shape[1]          # Getting shapes of V design matrices
    interest_param = theta[0]                                                   # Average treatment effect
    diagnostic = theta[1]                                                       # Diagnostic
    r_a2_s1, r_a1_s1, r_a1_s0, r_a0_s0 = theta[2:6]                             # Mean in different arms of trials
    beta = theta[6:6+w_len]                                                     # Parameters for sampling model
    gamma_a2s1 = theta[6+w_len:6+w_len+v3_len]                                  # Parameters missing model in R=2,A=3
    gamma_a1s1 = theta[6+w_len+v3_len:6+w_len+v3_len+v2a_len]                   # Parameters missing model in R=2,A=2
    gamma_a1s0 = theta[6+w_len+v3_len+v2a_len:6+w_len+v3_len+v2a_len+v2b_len]   # Parameters missing model in R=1,A=2
    gamma_a0s0 = theta[6+w_len+v3_len+v2a_len+v2b_len:]                         # Parameters missing model in R=1,A=1
    iptw = 1 / 0.5                                                              # Here, the propensity score is known

    # Calculate the sampling model, Pr(R | X; beta), nuisance parameters
    pr_s_nuisance = ee_regression(theta=beta,               # Using built-in regression
                                  X=W, y=1-r,               # ... given variables
                                  model='logistic')         # ... with logistic model
    pr_s1 = inverse_logit(np.dot(W, beta))                  # Predicted probability of R=2
    odds = pr_s1 / (1 - pr_s1)                              # Predicted odds of R=2
    iosw = r*1 + (1-r)/odds                                 # Inverse Odds of Sampling Weights

    # Calculate the missing models, Pr(M | X, A, R; gamma), nuisance parameters
    pr_m_nuisance_a2s1 = ee_regression(theta=gamma_a2s1,             # Using built-in regression
                                       X=V3, y=m,                    # ... given variables
                                       model='logistic')             # ... with logistic model
    pr_m_nuisance_a2s1 = pr_m_nuisance_a2s1 * (a == 3)               # Only A=3 arm contributes

    pr_m_nuisance_a1s1 = ee_regression(theta=gamma_a1s1,             # Using built-in regression
                                       X=V2a, y=m,                   # ... given variables
                                       model='logistic')             # ... with logistic model
    pr_m_nuisance_a1s1 = pr_m_nuisance_a1s1 * r * (a == 2)           # Only A=2,R=2 arm contributes

    pr_m_nuisance_a1s0 = ee_regression(theta=gamma_a1s0,             # Using built-in regression
                                       X=V2b, y=m,                   # ... given variables
                                       model='logistic')             # ... with logistic model
    pr_m_nuisance_a1s0 = pr_m_nuisance_a1s0 * (1-r) * (a == 2)       # Only A=2,R=1 arm contributes

    pr_m_nuisance_a0s0 = ee_regression(theta=gamma_a0s0,             # Using built-in regression
                                       X=V1, y=m,                    # ... given variables
                                       model='logistic')             # ... with logistic model
    pr_m_nuisance_a0s0 = pr_m_nuisance_a0s0 * (1-r) * (a == 1)       # Only A=1 arm contributes

    pr_m = (inverse_logit(np.dot(V3, gamma_a2s1))*r*(a == 3)         # Compute all probability in A=3
            + inverse_logit(np.dot(V2a, gamma_a1s1))*r*(a == 2)      # ... probability in A=2,R=2
            + inverse_logit(np.dot(V2b, gamma_a1s0))*(1-r)*(a == 2)  # ... probability in A=2,R=1
            + inverse_logit(np.dot(V1, gamma_a0s0))*(1-r)*(a == 1))  # ... probability in A=1
    ipmw = np.where(m, 1/pr_m, 0)                                    # Inverse probability of missing weight

    # Estimating functions for the mean in each arm for estimator
    ipw = iptw * ipmw * iosw                                         # Overall inverse probability weight
    y_nan = np.nan_to_num(y, copy=True, nan=0.)                      # Convert missing y to 0 to prevent errors
    yhat_a2s1 = r * m * (a == 3) * ipw * (y_nan - r_a2_s1)           # Risk in R=2, A=3 reweighted
    yhat_a1s1 = r * m * (a == 2) * ipw * (y_nan - r_a1_s1)           # Risk in R=2, A=2 reweighted
    yhat_a1s0 = (1-r) * m * (a == 2) * ipw * (y_nan - r_a1_s0)       # Risk in R=1, A=2 reweighted
    yhat_a0s0 = (1-r) * m * (a == 1) * ipw * (y_nan - r_a0_s0)       # Risk in R=1, A=1 reweighted

    # Estimating function for risk difference
    ate_brdg = np.ones(y.shape[0])*((r_a2_s1 - r_a1_s1) + (r_a1_s0 - r_a0_s0)) - interest_param

    # Estimating function for diagnostic
    ate_diag = np.ones(y.shape[0])*(r_a1_s1 - r_a1_s0) - diagnostic

    # Returning stacked estimating functions
    return np.vstack((ate_brdg,
                      ate_diag,
                      yhat_a2s1, yhat_a1s1, yhat_a1s0, yhat_a0s0,
                      pr_s_nuisance,
                      pr_m_nuisance_a2s1, pr_m_nuisance_a1s1, pr_m_nuisance_a1s0, pr_m_nuisance_a0s0
                      ))
