import numpy as np

def givens_rotation(a, b):
    """Compute the Givens rotation parameters c and s."""
    norm = np.sqrt(a**2 + b**2)
    if norm < 1e-10:
        return 1.0, 0.0
    else:
        c = a / norm
        s = -b / norm
        return c, s

def apply_givens(A, c, s, i, j):
    """Apply the Givens rotation to rows and columns i and j of matrix A."""
    # Apply to rows
    row_i = c * A[i, :] - s * A[j, :]
    row_j = s * A[i, :] + c * A[j, :]
    A[i, :], A[j, :] = row_i, row_j

    # Apply to columns to maintain symmetry
    col_i = c * A[:, i] - s * A[:, j]
    col_j = s * A[:, i] + c * A[:, j]
    A[:, i], A[:, j] = col_i, col_j

def tridiagonalize_band(A, bandwidth):
    """Tridiagonalize a symmetric band matrix using Givens rotations."""
    n = A.shape[0]
    for j in range(n - 2):  # Iterate through columns
        for i in range(min(j + bandwidth, n - 1), j + 1, -1):  # Work from the bottom row upward
            if abs(A[i, j]) > 1e-10:  # If significant, apply Givens rotation
                c, s = givens_rotation(A[i - 1, j], A[i, j])
                apply_givens(A, c, s, i - 1, i)

    # Extract diagonal and subdiagonal elements
    d = np.diag(A)
    e = np.diag(A, k=1)
    return d, e

# Example symmetric band matrix
A = np.array([
    [4.0, 1.0, 2.0, 0.0, 0.0],
    [1.0, 3.0, 1.0, 2.0, 0.0],
    [2.0, 1.0, 4.0, 1.0, 0.0],
    [0.0, 2.0, 1.0, 5.0, 1.0],
    [0.0, 0.0, 0.0, 1.0, 3.0]
], dtype=float)

# Bandwidth (number of diagonals above the main diagonal)
bandwidth = 2

# Perform tridiagonalization
d, e = tridiagonalize_band(A, bandwidth)

# Output the final tridiagonal matrix
print("Diagonal (d):", d)
print("Subdiagonal (e):", e)

# Output the modified matrix
print("\nModified Matrix A:")
print(A)