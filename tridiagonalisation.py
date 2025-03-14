import numpy as np

def householder_tridiagonalize(matrix):
    if not isinstance(matrix, np.ndarray) or len(matrix.shape) != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("Input should be a square matrix represented as a NumPy array.")

    n = matrix.shape[0]
    tridiagonal = matrix.astype(float)

    for i in range(n - 2):
        x = tridiagonal[i + 1:, i]
        norm_x = np.linalg.norm(x)

        if norm_x == 0:
            continue

        u = x.copy()
        u[0] += np.sign(u[0]) * norm_x
        u /= np.linalg.norm(u)

        H = np.eye(n - i - 1) - 2.0 * np.outer(u, u)

        tridiagonal[i + 1:, i:] = H @ tridiagonal[i + 1:, i:]
        tridiagonal[:, i + 1:] = tridiagonal[:, i + 1:] @ H

    return tridiagonal

# Example usage
matrix = np.array([
    [4.0, 1.0, 2.0, 0.0, 0.0],
    [1.0, 3.0, 1.0, 2.0, 0.0],
    [2.0, 1.0, 4.0, 1.0, 0.0],
    [0.0, 2.0, 1.0, 5.0, 1.0],
    [0.0, 0.0, 0.0, 1.0, 3.0]
], dtype=float)

tridiagonal_matrix = householder_tridiagonalize(matrix)
print("Tridiagonalized Matrix:")
print(tridiagonal_matrix)
