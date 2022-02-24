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

def calculate(ggr_df, h1, h2, unknown, N2, m):
    N1 = len(ggr_df)
    Ns = np.sum(ggr_df['ovp'])
    y = ggr_df['Phenotype'] * ggr_df['gz']
    x1 = ggr_df['ggg'] + (N2 - Ns) * ggr_df['grg']
    w = (h1 / m * ggr_df['gg'] + (1 - h1)) * (h2 / m * ggr_df['grrg'] + (1 - h2) * ggr_df['ggg'] + (N2 - Ns) * ggr_df['grg']) / N2
    w = np.array(1 / w)
    x0 = ggr_df['gg'] * ggr_df['ovp']
    if Ns > 0:
        #x0 = ggr_df['ovp']
        X = np.vstack((x0, x1)).T
        xwx = (X.T * w).dot(X)
        xwy = (X.T * w).dot(y)
        beta = linalg.inv(xwx).dot(xwy)
        rhomn = beta[1]
        rhoen = beta[0]
        # w = (rhomn * x1 + rhoen * x0) ** 2 + (h1 / m * ggr_df['gg'] + (1 - h1)) * (h2 / m * ggr_df['grrg'] + (1 - h2) * ggr_df['ggg'] + (N2 - Ns) * ggr_df['grg']) / N2
        # w = 1 / w
        # w = np.array(w)
        # y = y - rhoen * x0
        # xwx = linalg.inv((X.T * w).dot(X))
        # xwy = (X.T * w).dot(y)
        # beta = xwx.dot(xwy)
        # rho = beta[1] * m * np.sqrt(N2)
        # sigma2 = (ywy - xwy.T.dot(xwx).dot(xwy)) / (N1 - 2)
        # se_rho = np.sqrt(sigma2 * xwx[1][1]) * m * np.sqrt(N2)

    else:
        rhomn = np.sum(x1 * w * y) / np.sum(x1 * w * x1)
        rhoen = 0
        if unknown:
            w0 = (h1 / m * ggr_df['gg'] + (1 - h1)) * (h2 / m * ggr_df['grrg'] + (1 - h2) * ggr_df['ggg'] + (N2 - Ns) * ggr_df['grg']) / N2
            rhomn, rhoen, x0 = _ns_regression(y, x1, ggr_df['gg'], w0, rhomn, N1)

        # ywy = np.sum(y * w * y)
        # xwx = 1 / np.sum(x1 * w * x1)
        # xwy = np.sum(x1 * w * y)
        # rho = xwx * xwy * m * np.sqrt(N2)
        # sigma2 = (ywy - xwy ** 2 * xwx) / (N1 - 1)
        # se_rho = np.sqrt(sigma2 * xwx) * m * np.sqrt(N2)

    w = (rhomn * x1 + rhoen * x0) ** 2 + (h1 / m * ggr_df['gg'] + (1 - h1)) * (h2 / m * ggr_df['grrg'] + (1 - h2) * ggr_df['ggg'] + (N2 - Ns) * ggr_df['grg']) / N2
    w = np.array(1 / w)
    y = y - rhoen * x0
    ywy = np.sum(y * w * y)
    xwx = 1 / np.sum(x1 * w * x1)
    xwy = np.sum(x1 * w * y)
    rho = xwx * xwy * m * np.sqrt(N2)
    sigma2 = (ywy - xwy ** 2 * xwx) / (N1 - 1)
    se_rho = np.sqrt(sigma2 * xwx) * m * np.sqrt(N2)

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
            ('N2', [N2])
        ]
    ))
    return out


def _ns_regression(y, x1, gg, w0, rhomn, N1):
    rhomn0 = 0
    rhoen = 0
    ovp = np.random.binomial(1, 0.25, N1)
    while abs(rhomn - rhomn0) < 0.001:
        rhomn0 = rhomn
        x0 = gg * ovp
        w = (rhomn * x1 + rhoen * x0) ** 2 + w0
        w = np.array(1 / w)

        X = np.vstack((x0, x1)).T
        xwx = linalg.inv((X.T * w).dot(X))
        xwy = (X.T * w).dot(y)
        beta = xwx.dot(xwy)
        rhomn = beta[1]
        rhoen = beta[0]
        ywy = np.sum(y * w * y)
        sigma2 = (ywy - xwy.T.dot(xwx).dot(xwy)) / (N1 - 2)

        y0 = rhomn * x1
        y1 = y0 + rhoen * gg
        nominator = np.exp(-((y - y1) * np.sqrt(w / sigma2)) ** 2 / 2)
        denominator = nominator + np.exp(-((y - y0) * np.sqrt(w / sigma2)) ** 2 / 2)
        ovp = nominator / denominator

        if np.sum(ovp) > N1 / 2:
            ovp = 1 - ovp
            rhomn = 0
            rhoen = 0
    return rhomn, rhoen, gg * ovp
