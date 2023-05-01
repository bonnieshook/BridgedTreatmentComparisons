############################################################################################################
# Fusing Trial Data for Treatment Comparisons: Single versus Multi-Span Bridging
#   Example replication in Python
#
# Paul Zivich (2023/04/27)
############################################################################################################

############################################################################
# Importing dependencies

import numpy as np
import pandas as pd
from delicatessen import MEstimator

from efuncs import ee_bridge_ss, ee_bridge_ms

############################################################################
# Estimating Functions for delicatessen


def psi_singlespan(theta):
    return ee_bridge_ss(theta=theta,
                        y=y, a=a, r=r, m=m,
                        W=W, V3=V3, V1=V1)


def psi_multispan(theta):
    return ee_bridge_ms(theta=theta,
                        y=y, a=a, r=r, m=m,
                        W=W, V3=V3, V2a=V22, V2b=V21, V1=V1)


############################################################################
# Design Matrices

w_cols = ['intercept', 'male', 'black', 'hispanic', 'idu', 'karnof_cat2', 'karnof_cat3',
          'age_ms_SP1', 'age_ms_SP2', 'age_ms_SP3', 'age_ms_SP4',
          'CD4BL_ms_SP1', 'CD4BL_ms_SP2', 'CD4BL_ms_SP3', 'CD4BL_ms_SP4']
w_len = len(w_cols)
v_cols = ['intercept', 'male', 'black', 'hispanic', 'idu', 'karnof_cat2', 'karnof_cat3',
          'age_ms_SP1', 'age_ms_SP2', 'age_ms_SP3', 'age_ms_SP4',
          'CD4BL_ms_SP1', 'CD4BL_ms_SP2', 'CD4BL_ms_SP3', 'CD4BL_ms_SP4']
v_len = len(v_cols)
vs_cols = ['intercept', 'male', 'black', 'hispanic', 'idu', 'karnof_cat23',
          'age_ms_SP1', 'age_ms_SP2', 'age_ms_SP3', 'age_ms_SP4',
          'CD4BL_ms_SP1', 'CD4BL_ms_SP2', 'CD4BL_ms_SP3', 'CD4BL_ms_SP4']
vs_len = len(vs_cols)

############################################################################
# Unrestricted

########################################
# Loading Pre-Processed Data

d = pd.read_csv("data/actg_unrestricted.csv")
d['intercept'] = 1
d['karnof_cat2'] = np.where(d['karnof_cat'] == 2, 1, np.nan)
d['karnof_cat2'] = np.where(d['karnof_cat'] != 2, 0, d['karnof_cat2'])
d['karnof_cat3'] = np.where(d['karnof_cat'] == 3, 1, np.nan)
d['karnof_cat3'] = np.where(d['karnof_cat'] != 3, 0, d['karnof_cat3'])
d['karnof_cat23'] = np.where(d['karnof_cat'] <= 2, 1, np.nan)
d['karnof_cat23'] = np.where(d['karnof_cat'] > 2, 0, d['karnof_cat23'])

W = np.asarray(d[w_cols])
V3 = np.asarray(d[vs_cols])
V21 = np.asarray(d[v_cols])
V22 = np.asarray(d[v_cols])
V1 = np.asarray(d[v_cols])
a = np.asarray(d['TRT'])
r = np.asarray(1 - d['ACTG320'])
m = np.where(d['CD4WK8'].isna(), 0, 1)

V = np.asarray(d[v_cols])


########################################
# Single-Span Estimator

# Continuous
y = np.asarray(d['CD4WK8'])
inits = [0., 175., 175., ] + [0., ]*w_len + [0., ]*vs_len + [0., ]*v_len
estrc_ss = MEstimator(psi_singlespan, init=inits)
estrc_ss.estimate(solver="lm", maxiter=20000)
cic_ss = estrc_ss.confidence_intervals()

# Binary
y = np.where(d['CD4WK8'] > 250, 1, 0)
y = np.where(d['CD4WK8'].isna(), np.nan, y)
inits = [0., 0.5, 0.5, ] + list(estrc_ss.theta[3:])
estrb_ss = MEstimator(psi_singlespan, init=inits)
estrb_ss.estimate(solver="lm", maxiter=20000)
cib_ss = estrb_ss.confidence_intervals()

########################################
# Multi-Span Estimator

# Continuous
y = np.asarray(d['CD4WK8'])
inits = [0., 0., 175., 175., 175., 175., ] + [0., ]*w_len + [0., ]*vs_len + [0., ]*v_len*3
estrc_ms = MEstimator(psi_multispan, init=inits)
estrc_ms.estimate(solver="lm", tolerance=1e-12, maxiter=20000)
cic_ms = estrc_ms.confidence_intervals()

# Binary
y = np.where(d['CD4WK8'] > 250, 1, 0)
y = np.where(d['CD4WK8'].isna(), np.nan, y)
inits = [0., 0., 0.5, 0.5, 0.5, 0.5, ] + list(estrc_ms.theta[6:])
estrb_ms = MEstimator(psi_multispan, init=inits)
estrb_ms.estimate(solver="lm", tolerance=1e-12, maxiter=20000)
cib_ms = estrb_ms.confidence_intervals()

########################################
# Results

print("=========================")
print("CD4 - Unrestricted")
print("=========================")
print("Continuous")
print("-------------------------")
print("SS:   ", np.round(estrc_ss.theta[0], 1), np.round(cic_ss[0, :], 1))
print("MS:   ", np.round(estrc_ms.theta[0], 1), np.round(cic_ms[0, :], 1))
print("Diag: ", np.round(estrc_ms.theta[1], 1), np.round(cic_ms[1, :], 1))
print("-------------------------")
print("Binary")
print("-------------------------")
print("SS:   ", np.round(estrb_ss.theta[0], 2), np.round(cib_ss[0, :], 2))
print("MS:   ", np.round(estrb_ms.theta[0], 2), np.round(cib_ms[0, :], 2))
print("Diag: ", np.round(estrb_ms.theta[1], 2), np.round(cib_ms[1, :], 2))
print("=========================")


############################################################################
# CD4 < 400

########################################
# Loading Pre-Processed Data

d = pd.read_csv("data/actg_400.csv")
d['karnof_cat23'] = np.where(d['karnof_cat3'] == 1, 1, 0)
d['karnof_cat23'] = np.where(d['karnof_cat3'].isna(), np.nan, d['karnof_cat23'])
d['intercept'] = 1


W = np.asarray(d[w_cols])
V3 = np.asarray(d[vs_cols])
V21 = np.asarray(d[v_cols])
V22 = np.asarray(d[v_cols])
V1 = np.asarray(d[v_cols])
a = np.asarray(d['TRT'])
r = np.asarray(1 - d['ACTG320'])
m = np.where(d['CD4WK8'].isna(), 0, 1)

########################################
# Single-Span Estimator

# Continuous
y = np.asarray(d['CD4WK8'])
inits = [0., 175., 175., ] + [0., ]*w_len + [0., ]*vs_len + [0., ]*v_len
estrc_ss = MEstimator(psi_singlespan, init=inits)
estrc_ss.estimate(solver="lm", maxiter=20000)
cic_ss = estrc_ss.confidence_intervals()

# Binary
y = np.where(d['CD4WK8'] > 250, 1, 0)
y = np.where(d['CD4WK8'].isna(), np.nan, y)
inits = [0., 0.5, 0.5, ] + list(estrc_ss.theta[3:])
estrb_ss = MEstimator(psi_singlespan, init=inits)
estrb_ss.estimate(solver="lm", maxiter=20000)
cib_ss = estrb_ss.confidence_intervals()

########################################
# Multi-Span Estimator

# Continuous
y = np.asarray(d['CD4WK8'])
inits = [0., 0., 175., 175., 175., 175., ] + [0., ]*w_len + [0., ]*vs_len + [0., ]*v_len*3
estrc_ms = MEstimator(psi_multispan, init=inits)
estrc_ms.estimate(solver="lm", tolerance=1e-12, maxiter=20000)
cic_ms = estrc_ms.confidence_intervals()

# Binary
y = np.where(d['CD4WK8'] > 250, 1, 0)
y = np.where(d['CD4WK8'].isna(), np.nan, y)
inits = [0., 0., 0.5, 0.5, 0.5, 0.5, ] + list(estrc_ms.theta[6:])
estrb_ms = MEstimator(psi_multispan, init=inits)
estrb_ms.estimate(solver="lm", tolerance=1e-12, maxiter=20000)
cib_ms = estrb_ms.confidence_intervals()

########################################
# Results

print("=========================")
print("CD4 < 400")
print("=========================")
print("Continuous")
print("-------------------------")
print("SS:   ", np.round(estrc_ss.theta[0], 1), np.round(cic_ss[0, :], 1))
print("MS:   ", np.round(estrc_ms.theta[0], 1), np.round(cic_ms[0, :], 1))
print("Diag: ", np.round(estrc_ms.theta[1], 1), np.round(cic_ms[1, :], 1))
print("-------------------------")
print("Binary")
print("-------------------------")
print("SS:   ", np.round(estrb_ss.theta[0], 2), np.round(cib_ss[0, :], 2))
print("MS:   ", np.round(estrb_ms.theta[0], 2), np.round(cib_ms[0, :], 2))
print("Diag: ", np.round(estrb_ms.theta[1], 2), np.round(cib_ms[1, :], 2))
print("=========================")


############################################################################
# CD4 < 300

########################################
# Loading Pre-Processed Data

d = pd.read_csv("data/actg_300.csv")
d['intercept'] = 1
d['karnof_cat2'] = np.where(d['karnof_cat'] == 2, 1, np.nan)
d['karnof_cat2'] = np.where(d['karnof_cat'] != 2, 0, d['karnof_cat2'])
d['karnof_cat3'] = np.where(d['karnof_cat'] == 3, 1, np.nan)
d['karnof_cat3'] = np.where(d['karnof_cat'] != 3, 0, d['karnof_cat3'])
d['karnof_cat23'] = np.where(d['karnof_cat'] <= 2, 1, np.nan)
d['karnof_cat23'] = np.where(d['karnof_cat'] > 2, 0, d['karnof_cat23'])

W = np.asarray(d[w_cols])
V3 = np.asarray(d[vs_cols])
V21 = np.asarray(d[v_cols])
V22 = np.asarray(d[v_cols])
V1 = np.asarray(d[v_cols])
a = np.asarray(d['TRT'])
r = np.asarray(1 - d['ACTG320'])
m = np.where(d['CD4WK8'].isna(), 0, 1)

########################################
# Single-Span Estimator

# Continuous
y = np.asarray(d['CD4WK8'])
inits = [0., 175., 175., ] + [0., ]*w_len + [0., ]*vs_len + [0., ]*v_len
estrc_ss = MEstimator(psi_singlespan, init=inits)
estrc_ss.estimate(solver="lm", maxiter=20000)
cic_ss = estrc_ss.confidence_intervals()

# Binary
y = np.where(d['CD4WK8'] > 250, 1, 0)
y = np.where(d['CD4WK8'].isna(), np.nan, y)
inits = [0., 0.5, 0.5, ] + list(estrc_ss.theta[3:])
estrb_ss = MEstimator(psi_singlespan, init=inits)
estrb_ss.estimate(solver="lm", maxiter=20000)
cib_ss = estrb_ss.confidence_intervals()

########################################
# Multi-Span Estimator

# Continuous
y = np.asarray(d['CD4WK8'])
inits = [0., 0., 175., 175., 175., 175., ] + [0., ]*w_len + [0., ]*vs_len + [0., ]*v_len*3
estrc_ms = MEstimator(psi_multispan, init=inits)
estrc_ms.estimate(solver="lm", tolerance=1e-12, maxiter=20000)
cic_ms = estrc_ms.confidence_intervals()

# Binary
y = np.where(d['CD4WK8'] > 250, 1, 0)
y = np.where(d['CD4WK8'].isna(), np.nan, y)
inits = [0., 0., 0.5, 0.5, 0.5, 0.5, ] + list(estrc_ms.theta[6:])
estrb_ms = MEstimator(psi_multispan, init=inits)
estrb_ms.estimate(solver="lm", tolerance=1e-12, maxiter=20000)
cib_ms = estrb_ms.confidence_intervals()

########################################
# Results

print("=========================")
print("CD4 < 300")
print("=========================")
print("Continuous")
print("-------------------------")
print("SS:   ", np.round(estrc_ss.theta[0], 1), np.round(cic_ss[0, :], 1))
print("MS:   ", np.round(estrc_ms.theta[0], 1), np.round(cic_ms[0, :], 1))
print("Diag: ", np.round(estrc_ms.theta[1], 1), np.round(cic_ms[1, :], 1))
print("-------------------------")
print("Binary")
print("-------------------------")
print("SS:   ", np.round(estrb_ss.theta[0], 2), np.round(cib_ss[0, :], 2))
print("MS:   ", np.round(estrb_ms.theta[0], 2), np.round(cib_ms[0, :], 2))
print("Diag: ", np.round(estrb_ms.theta[1], 2), np.round(cib_ms[1, :], 2))
print("=========================")


############################################################################
# CD4 < 200

########################################
# Loading Pre-Processed Data

d = pd.read_csv("data/actg_200.csv")
d['intercept'] = 1
d['karnof_cat2'] = np.where(d['karnof_cat'] == 2, 1, np.nan)
d['karnof_cat2'] = np.where(d['karnof_cat'] != 2, 0, d['karnof_cat2'])
d['karnof_cat3'] = np.where(d['karnof_cat'] == 3, 1, np.nan)
d['karnof_cat3'] = np.where(d['karnof_cat'] != 3, 0, d['karnof_cat3'])
d['karnof_cat23'] = np.where(d['karnof_cat'] <= 2, 1, np.nan)
d['karnof_cat23'] = np.where(d['karnof_cat'] > 2, 0, d['karnof_cat23'])

W = np.asarray(d[w_cols])
V3 = np.asarray(d[vs_cols])
V21 = np.asarray(d[v_cols])
V22 = np.asarray(d[v_cols])
V1 = np.asarray(d[v_cols])
a = np.asarray(d['TRT'])
r = np.asarray(1 - d['ACTG320'])
m = np.where(d['CD4WK8'].isna(), 0, 1)

########################################
# Single-Span Estimator

# Continuous
y = np.asarray(d['CD4WK8'])
inits = [0., 175., 175., ] + [0., ]*w_len + [0., ]*vs_len + [0., ]*v_len
estrc_ss = MEstimator(psi_singlespan, init=inits)
estrc_ss.estimate(solver="lm", maxiter=20000)
cic_ss = estrc_ss.confidence_intervals()

# Binary
y = np.where(d['CD4WK8'] > 250, 1, 0)
y = np.where(d['CD4WK8'].isna(), np.nan, y)
inits = [0., 0.5, 0.5, ] + list(estrc_ss.theta[3:])
estrb_ss = MEstimator(psi_singlespan, init=inits)
estrb_ss.estimate(solver="lm", maxiter=20000)
cib_ss = estrb_ss.confidence_intervals()

########################################
# Multi-Span Estimator

# Continuous
y = np.asarray(d['CD4WK8'])
inits = [0., 0., 175., 175., 175., 175., ] + [0., ]*w_len + [0., ]*vs_len + [0., ]*v_len*3
estrc_ms = MEstimator(psi_multispan, init=inits)
estrc_ms.estimate(solver="lm", tolerance=1e-12, maxiter=20000)
cic_ms = estrc_ms.confidence_intervals()

# Binary
y = np.where(d['CD4WK8'] > 250, 1, 0)
y = np.where(d['CD4WK8'].isna(), np.nan, y)
inits = [0., 0., 0.5, 0.5, 0.5, 0.5, ] + list(estrc_ms.theta[6:])
estrb_ms = MEstimator(psi_multispan, init=inits)
estrb_ms.estimate(solver="lm", tolerance=1e-12, maxiter=20000)
cib_ms = estrb_ms.confidence_intervals()

########################################
# Results

print("=========================")
print("CD4 < 200")
print("=========================")
print("Continuous")
print("-------------------------")
print("SS:   ", np.round(estrc_ss.theta[0], 1), np.round(cic_ss[0, :], 1))
print("MS:   ", np.round(estrc_ms.theta[0], 1), np.round(cic_ms[0, :], 1))
print("Diag: ", np.round(estrc_ms.theta[1], 1), np.round(cic_ms[1, :], 1))
print("-------------------------")
print("Binary")
print("-------------------------")
print("SS:   ", np.round(estrb_ss.theta[0], 2), np.round(cib_ss[0, :], 2))
print("MS:   ", np.round(estrb_ms.theta[0], 2), np.round(cib_ms[0, :], 2))
print("Diag: ", np.round(estrb_ms.theta[1], 2), np.round(cib_ms[1, :], 2))
print("=========================")

############################################################################
# Output

# =========================
# CD4 - Unrestricted
# =========================
# Continuous
# -------------------------
# SS:    99.6 [ 75.5 123.6]
# MS:    30.9 [-13.   74.9]
# Diag:  64.8 [35.3 94.3]
# -------------------------
# Binary
# -------------------------
# SS:    0.23 [0.15 0.31]
# MS:    0.05 [-0.09  0.19]
# Diag:  0.16 [0.07 0.26]
# =========================
# =========================
# CD4 < 400
# =========================
# Continuous
# -------------------------
# SS:    69.0 [51.6 86.5]
# MS:    27.6 [-10.4  65.6]
# Diag:  37.7 [12.4 63. ]
# -------------------------
# Binary
# -------------------------
# SS:    0.18 [0.1  0.26]
# MS:    0.04 [-0.09  0.18]
# Diag:  0.12 [0.03 0.21]
# =========================
# =========================
# CD4 < 300
# =========================
# Continuous
# -------------------------
# SS:    56.9 [41.9 71.9]
# MS:    25.8 [-6.9 58.4]
# Diag:  26.7 [ 4.7 48.6]
# -------------------------
# Binary
# -------------------------
# SS:    0.15 [0.07 0.22]
# MS:    0.02 [-0.1   0.14]
# Diag:  0.12 [0.03 0.2 ]
# =========================
# =========================
# CD4 < 200
# =========================
# Continuous
# -------------------------
# SS:    63.2 [48.  78.4]
# MS:    57.4 [28.7 86. ]
# Diag:  8.9 [-10.4  28.3]
# -------------------------
# Binary
# -------------------------
# SS:    0.14 [0.06 0.21]
# MS:    0.08 [-0.03  0.18]
# Diag:  0.07 [-0.01  0.14]
# =========================
