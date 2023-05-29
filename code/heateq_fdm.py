from math import cos, exp, pi, sin

import numpy as np
import numpy.linalg as lin
import scipy.sparse as sp


def heateq_fdm(a, b, T, N, M, theta):
    h = (b - a) / (N + 1)
    k = T / M
    e = np.ones(N, 1)
    I = sp.eye(N)
    G = sp.diags((-e[0 : N - 1], 2 * e, -e[0 : N - 1]), (-1, 0, 1))
    G = h ** (-2) * G
    B = I + k * theta * G
    C = I - k * (1 - theta) * G
    x = np.linspace(a + h, b - h, N)
    u0 = x * sin(pi * x)
    f = -(1 - np.pi**2) * x * sin(np.pi * x) - 2 * np.pi * cos(np.pi * x)
    u = np.zeros((N, M + 1))
    err = np.zeros(M + 1)
    u[:, 0] = u0
    for j in range(M):
        F = f * (theta * np.exp(-(j + 1) * k) + (1 - theta) * np.exp(-j * k))
        RHS = np.matmul(C, u[:, j]) + k * F
        u[:, j + 1] = lin.solve(B, RHS)
        err[j + 1] = lin.norm(u[:, j + 1] - exp(-k * (j + 1)) * u0)
    return np.sqrt(h) * max(err)
