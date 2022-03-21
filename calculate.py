#!/usr/bin/python
from __future__ import division, print_function
from cmath import exp
import numpy as np
import pandas as pd
import numpy.linalg as linalg
from math import sqrt
from scipy.stats import norm
from collections import OrderedDict
import math

def calculate(ggr_df, h1, h2, intercept, N2, m):
    N1 = len(ggr_df)
    Ns = np.sum(ggr_df['ovp'])
    y = ggr_df['Phenotype'] * ggr_df['gz']
    x1 = ggr_df['ggg'] + (N2 - Ns) * ggr_df['grg']
    w = (h1 / m * ggr_df['gg'] + (1 - h1)) * (h2 / m * ggr_df['grrg'] + (1 - h2) * ggr_df['ggg'] + (N2 - Ns) * ggr_df['grg']) / N2
    w = np.array(1 / w)
    x0 = ggr_df['gg'] * ggr_df['ovp']
    if Ns > 0:
        X = np.vstack((x0, x1)).T
        xwx = (X.T * w).dot(X)
        xwy = (X.T * w).dot(y)
        beta = linalg.inv(xwx).dot(xwy)
        rhomn = beta[1]
        rhoen = beta[0]

    else:
        rhomn = np.sum(x1 * w * y) / np.sum(x1 * w * x1)
        rhoen = 0

    w = (rhomn * x1 + rhoen * x0) ** 2 + (h1 / m * ggr_df['gg'] + (1 - h1)) * (h2 / m * ggr_df['grrg'] + (1 - h2) * ggr_df['ggg'] + (N2 - Ns) * ggr_df['grg']) / N2
    w = np.array(1 / w)
    y = y - rhoen * x0
    ywy = np.sum(y * w * y)
    xwx = 1 / np.sum(x1 * w * x1)
    xwy = np.sum(x1 * w * y)
    rho = xwx * xwy * m * np.sqrt(N2)
    sigma2 = (ywy - xwy ** 2 * xwx) / (N1 - 1)
    se_rho = np.sqrt(sigma2 * xwx) * m * np.sqrt(N2)

    if intercept != 0.0:
        T = np.sum(w * x1 * ggr_df['gg']) / np.sum(w * x1 * x1)
        rho = rho - T * m * intercept * np.sqrt(N2 / N1)

    out = pd.DataFrame(OrderedDict(
        [
            ('rho', [rho]),
            ('se', [se_rho]),
            ('pvalue', [norm.sf(abs(rho / se_rho)) * 2]),
            ('corr', [rho / sqrt(h1 * h2)]),
            ('h1', [h1]),
            ('h2', [h2]),
            ('m', [m]),
            ('N1', [N1]),
            ('N2', [N2]),
            ('Intercept', [intercept])
        ]
    ))
    return out
