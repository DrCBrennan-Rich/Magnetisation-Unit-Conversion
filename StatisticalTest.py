# -*- coding: utf-8 -*-
"""
@author: cbrph
"""
import numpy as np
from scipy import stats

def Analyse_Measurements(data, alpha=0.05):

    data = np.asarray(data)

    if len(data) < 3:
        raise ValueError("Need at least 3 measurements.")

    mean = np.mean(data)
    StandardDeviation = np.std(data, ddof=1)          # sample standard deviation
    StandardErrorMean = StandardDeviation/np.sqrt(len(data))

    W, p = stats.shapiro(data)

    normal = p > alpha

    print(f"Number of measurements : {len(data)}")
    print(f"Mean                   : {mean:.6g}")
    print(f"Standard deviation     : {StandardDeviation:.6g}")
    print(f"Standard error         : {StandardErrorMean:.6g}")
    print(f"Skewness               : {stats.skew(data):.3f}")
    print(f"Kurtosis               : {stats.kurtosis(data):.3f}")
    print()
    print("Shapiro-Wilk normality test")
    print(f"W = {W:.4f}")
    print(f"p = {p:.4f}")

    if normal:
        print("\n✓ No evidence against a normal distribution.")
    else:
        print("\n⚠ Data may not be normally distributed.")

    return {
        "mean": mean,
        "StandardDeviation": StandardDeviation,
        "StandardErrorMean": StandardErrorMean,
        "W": W,
        "p": p,
        "normal": normal
    }

#List of data points should be inserted here
Measurements = [17.70, 17.428, 17.982, 17.365, 17.937, 17.326, 17.847]

Results = Analyse_Measurements(Measurements)
