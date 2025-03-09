import numpy as np
from scipy.linalg import hessenberg

# Given matrix A
A = np.array([
    [4.0, 1.0, 0.0, 0.0, 0.0],
    [1.0, 4.0, 1.0, 0.0, 0.0],
    [0.0, 1.0, 4.0, 1.0, 0.0],
    [0.0, 0.0, 1.0, 4.0, 1.0],
    [0.0, 0.0, 0.0, 1.0, 4.0]
], dtype=float)

# Compute the tridiagonal form using Hessenberg decomposition (since it reduces to tridiagonal for symmetric matrices)
T, _ = hessenberg(A, calc_q=True)

print (T)