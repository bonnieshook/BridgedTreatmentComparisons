############################################################################################################
# Fusing Trial Data for Treatment Comparisons: Single versus Multi-Span Bridging
#   Functions for data generating mechanisms for simulation experiments
#
# Paul Zivich (2023/04/27)
############################################################################################################

# Importing dependencies
import numpy as np
import pandas as pd
from scipy.stats import logistic


def generate_data(n1, n0, scenario):
    """Generate the data set for the specified trial sizes and scenario.

    Parameters
    ----------
    n1 : int
        Number of observations to generate for the trial in the target population
    n0 : int
        Number of observations to generate for the trial in the secondary population
    scenario : int
        Scenario for the data generating mechanism, 1-5 are valid options.

    Returns
    -------
    DataFrame
    """
    # Generating target population
    d1 = _target_pop_(n=n1, scenario=scenario)    # Generates trial in target population
    d0 = _second_pop_(n=n0, scenario=scenario)    # Generates trial in secondary (non-focal) population
    d = pd.concat([d1, d0], ignore_index=True)    # Puts data sets together

    # Data processing
    d['Y'] = (np.where(d['A'] == 1, 1, 0)*d['Ya1']       # Apply causal consistency formula
              + np.where(d['A'] == 2, 1, 0)*d['Ya2']     # ... for continuous outcome
              + np.where(d['A'] == 3, 1, 0)*d['Ya3'])
    d['B'] = (np.where(d['A'] == 1, 1, 0)*d['Ba1']       # Apply causal consistency formula
              + np.where(d['A'] == 2, 1, 0)*d['Ba2']     # ... for binary outcome
              + np.where(d['A'] == 3, 1, 0)*d['Ba3'])
    d['Y'] = np.where(d['M'] == 1, np.nan, d['Y'])       # Setting as NaN if missing
    d['B'] = np.where(d['M'] == 1, np.nan, d['B'])       # Setting as NaN if missing
    d['C'] = 1                                           # Intercept term for design matrices

    # Returning setup data
    return d


def calculate_truth(n, scenario):
    """Estimate the true value empirically by comparing potential outcomes of a simulation

    Parameters
    ----------
    n : int
        Number of observations to generate for the target population
    scenario : int
        Scenario for the data generating mechanism, 1-5 are valid options.

    Returns
    -------
    float, float
    """
    d = _target_pop_(n=n, scenario=scenario)    # Generate observations in target population
    return [np.mean(d['Ya3'] - d['Ya1']),       # Calculate ATE using potential outcomes
            np.mean(d['Ba3'] - d['Ba1'])]       # ... for both outcome types and return


def _target_pop_(n, scenario):
    """Internal function called to generate data for the target population according to the requested scenario. See
    Table 2 of the paper.
    """
    # Covariate 1 (baseline IDU)
    if scenario == 1:
        p1 = 0.25
    elif scenario in [2, 3, 4, 5]:
        p1 = 0.2
    else:
        raise ValueError("Invalid scenario")
    x1 = np.random.binomial(n=1, p=p1, size=n)

    # Covariate 2 (baseline CD4 cell count)
    if scenario == 1:
        x2 = np.random.normal(175 - 10*x1, scale=30, size=n)
    else:
        x2 = np.random.normal(175 + 10 - 20*x1, scale=30, size=n)
    x2 = np.where(x2 < 0, 0, x2)                                    # Bounding CD4 at zero

    # Potential outcomes under each treatment by scenario
    if scenario == 1:
        ya1 = np.random.normal(50 - 5*x1 + 1.1*x2, scale=20, size=n)
        ya2 = np.random.normal(80 - 5*x1 + 1.1*x2, scale=20, size=n)
        ya3 = np.random.normal(110 - 5*x1 + 1.1*x2, scale=20, size=n)
    elif scenario == 2:
        ya1 = np.random.normal(35 - 80*x1 + 1.0*x2, scale=20, size=n)
        ya2 = np.random.normal(30 - 10*x1 + 1.1*x2, scale=20, size=n)
        ya3 = np.random.normal(40 + 20*x1 + 1.2*x2, scale=20, size=n)
    elif scenario == 3:
        ya1 = np.random.normal(35 - 80*x1 + 1.0*x2, scale=20, size=n)
        ya2 = np.random.normal(40 + 10*x1 + 1.0*x2, scale=20, size=n)
        ya3 = np.random.normal(40 + 20*x1 + 1.2*x2, scale=20, size=n)
    elif scenario == 4:
        ya1 = np.random.normal(45 - 80*x1 + 1.0*x2 - 0, scale=20, size=n)
        ya2 = np.random.normal(30 - 10*x1 + 1.1*x2, scale=20, size=n)
        ya3 = np.random.normal(30 + 20*x1 + 1.2*x2 + 20, scale=20, size=n)
    else:
        ya1 = np.random.normal(45 - 80*x1 + 1.0*x2 - 0, scale=20, size=n)
        ya2 = np.random.normal(30 - 10*x1 + 1.1*x2 + 0, scale=20, size=n)
        ya3 = np.random.normal(30 + 20*x1 + 1.2*x2 + 20, scale=20, size=n)

    # Bounding CD4 at zero
    ya1 = np.where(ya1 < 0, 0, ya1)
    ya2 = np.where(ya2 < 0, 0, ya2)
    ya3 = np.where(ya3 < 0, 0, ya3)

    # Binary outcome version from CD4 (CD4 > 250)
    ba1 = np.where(ya1 > 250, 1, 0)
    ba2 = np.where(ya2 > 250, 1, 0)
    ba3 = np.where(ya3 > 250, 1, 0)

    # Missing data mechanism
    if scenario == 1:
        m = np.random.binomial(n=1, p=0.15, size=n)
    else:
        m = np.random.binomial(n=1, p=logistic.cdf(-2.1 + 0.5*x1), size=n)

    # Treatment assignment mechanism
    a = np.random.binomial(n=1, p=0.5, size=n) + 2

    # Creating as a pandas data set
    d = pd.DataFrame()
    d['X1'] = x1
    d['X2'] = x2
    d['A'] = a
    d['Ya1'] = ya1
    d['Ya2'] = ya2
    d['Ya3'] = ya3
    d['Ba1'] = ba1
    d['Ba2'] = ba2
    d['Ba3'] = ba3
    d['M'] = m
    d['S'] = 1
    return d


def _second_pop_(n, scenario):
    """Internal function called to generate data for the secondary population according to the requested scenario. See
    Table 2 of the paper.
    """
    # Covariate 1 (baseline IDU)
    if scenario == 1:
        p1 = 0.25
    elif scenario in [2, 3, 4, 5]:
        p1 = 0.5
    else:
        raise ValueError("Invalid scenario")
    x1 = np.random.binomial(n=1, p=p1, size=n)

    # Covariate 2 (baseline CD4 cell count)
    if scenario == 1:
        x2 = np.random.normal(175 - 10*x1, scale=30, size=n)
    else:
        x2 = np.random.normal(175 - 20*x1, scale=30, size=n)
    x2 = np.where(x2 < 0, 0, x2)                               # Bounding CD4 by zero

    # Potential outcomes under each treatment by scenario
    if scenario == 1:
        ya1 = np.random.normal(50 - 5*x1 + 1.1*x2, scale=20, size=n)
        ya2 = np.random.normal(80 - 5*x1 + 1.1*x2, scale=20, size=n)
        ya3 = np.random.normal(110 - 5*x1 + 1.1*x2, scale=20, size=n)
    elif scenario == 2:
        ya1 = np.random.normal(35 - 80*x1 + 1.0*x2, scale=20, size=n)
        ya2 = np.random.normal(30 - 10*x1 + 1.1*x2, scale=20, size=n)
        ya3 = np.random.normal(40 + 20*x1 + 1.2*x2, scale=20, size=n)
    elif scenario == 3:
        ya1 = np.random.normal(35 - 80*x1 + 1.0*x2, scale=20, size=n)
        ya2 = np.random.normal(45 - 10*x1 + 1.1*x2, scale=20, size=n)
        ya3 = np.random.normal(40 + 20*x1 + 1.2*x2, scale=20, size=n)
    elif scenario == 4:
        ya1 = np.random.normal(45 - 80*x1 + 1.0*x2 - 30, scale=20, size=n)
        ya2 = np.random.normal(30 - 10*x1 + 1.1*x2, scale=20, size=n)
        ya3 = np.random.normal(30 + 20*x1 + 1.2*x2 + 0, scale=20, size=n)
    else:
        ya1 = np.random.normal(45 - 80*x1 + 1.0*x2 - 30, scale=20, size=n)
        ya2 = np.random.normal(30 - 10*x1 + 1.1*x2 + 10, scale=20, size=n)
        ya3 = np.random.normal(30 + 20*x1 + 1.2*x2 + 0, scale=20, size=n)

    # Bounding CD4 by zero
    ya1 = np.where(ya1 < 0, 0, ya1)
    ya2 = np.where(ya2 < 0, 0, ya2)
    ya3 = np.where(ya3 < 0, 0, ya3)

    # Generating binary outcome (CD4 > 250)
    ba1 = np.where(ya1 > 250, 1, 0)
    ba2 = np.where(ya2 > 250, 1, 0)
    ba3 = np.where(ya3 > 250, 1, 0)

    # Missing data mechanism
    if scenario == 1:
        m = np.random.binomial(n=1, p=0.15, size=n)
    else:
        m = np.random.binomial(n=1, p=logistic.cdf(-2.0 + 0.5*x1), size=n)

    # Assigned treatment
    a = np.random.binomial(n=1, p=0.5, size=n) + 1

    # Creating as a pandas data set
    d = pd.DataFrame()
    d['X1'] = x1
    d['X2'] = x2
    d['A'] = a
    d['Ya1'] = ya1
    d['Ya2'] = ya2
    d['Ya3'] = ya3
    d['Ba1'] = ba1
    d['Ba2'] = ba2
    d['Ba3'] = ba3
    d['M'] = m
    d['S'] = 0
    return d
