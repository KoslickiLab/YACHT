import numpy as np
import cvxpy as cp
import pandas as pd
import csv
import compute_weight as cw
from scipy.sparse import load_npz
import argparse
import warnings


def solve_lp(A, y, w):
    """
    Runs the linear program for nonnegative quantile regression with weight w on the equation Ax = y.
    :param A: matrix (reference database)
    :param y: vector (sample kmer counts)
    :param w: False positive weight vector
    :return: vector x (estimated organism kmer counts)
    """
    K, N = np.shape(A)
    x = cp.Variable(N)
    u = cp.Variable(K)
    v = cp.Variable(K)
    tau = 1 / (w + 1)
    objective = cp.Minimize(
        tau @ u + (1 - tau) @ v
    )
    constraints = [
        x >= 0,
        u >= 0,
        v >= 0,
        u - v + (A @ x) == y,
    ]
    prob = cp.Problem(objective, constraints)
    result = prob.solve(solver=cp.SCIPY, verbose=False)
    recov_y = A @ x.value
    resid = y - (A @ x.value)
    return x.value, resid