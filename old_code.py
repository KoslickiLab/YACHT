def recover_abundance_from_vectors(A, y, w):
    """
    Runs the linear program for quantile regression with weight w on the equation Ax = y.
    :param A: matrix (reference database)
    :param y: vector (sample kmer counts)
    :param w: False positive weight
    :return: vector x (estimated organism counts)
    """
    K, N = np.shape(A)
    x = cp.Variable(N)
    u = cp.Variable(K)
    v = cp.Variable(K)
    tau = 1 / (w + 1)
    ones_K = np.ones(K)
    objective = cp.Minimize(
        tau * (ones_K @ u) + (1 - tau) * (ones_K @ v)
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