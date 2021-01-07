#!/usr/bin/python
from __future__ import division, print_function
import numpy as np
import pandas as pd
import numpy.linalg as linalg
from math import sqrt
from sklearn import linear_model
from scipy.stats import norm
from collections import OrderedDict

def calculate(ggr_df, N1, N2, h1, h2):
    m = len(ggr_df)
    y = np.array(ggr_df['y'])
    ggr = np.array(ggr_df['ggr'])
    l2 = np.array(ggr_df['l2'])
    w1 = N1 * (N2 * h2 * l2 / m + 1)
    w2 = (np.sum(y) / np.sum(ggr) * ggr) ** 2

    w = 1 / (w1 + w2)
    lm = linear_model.LinearRegression().fit(pd.DataFrame(ggr), pd.DataFrame(y), sample_weight=w)
    intercept = lm.intercept_[0]
    rho_comb = lm.coef_[0]
    w2 = (intercept + rho_comb * ggr) ** 2
    w = 1 / (w1 + w2)

    # Calculate Jackknife variance estimate

    nblock = 200
    q_block = np.empty(200)

    X = np.vstack([np.ones(len(ggr)), ggr]).T
    tot_X = X.T.dot((X.T * w).T)
    tot_Y = X.T.dot(y * w)

    for j, (b_y, b_x, b_w) in enumerate(zip(np.array_split(y, nblock), np.array_split(ggr, nblock), np.array_split(w, nblock))):
        cur_X = np.vstack([np.ones(len(b_x)), b_x]).T
        cur_tot_X = tot_X - cur_X.T.dot((cur_X.T * b_w).T)
        cur_tot_Y = tot_Y - cur_X.T.dot(b_y * b_w)
        coeff = linalg.inv(cur_tot_X).dot(cur_tot_Y)[1]
        q_block[j] = coeff * m / sqrt(N2)

    rho = np.mean(q_block)
    se_rho = sqrt((nblock - 1) * np.sum((q_block - rho) ** 2) / nblock)

    out = pd.DataFrame(OrderedDict(
        [
            ('rho', rho),
            ('se', se_rho),
            ('pvalue', norm.sf(abs(rho / se_rho)) * 2),
            ('corr', rho / sqrt(h1 * h2)),
            ('h1', h1),
            ('h2', h2),
            ('m', m),
            ('N1', N1),
            ('N2', N2)
        ]
    ))
    return out
