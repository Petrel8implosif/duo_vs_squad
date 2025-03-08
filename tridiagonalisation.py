import numpy as np

def givens_rotation(a, b):
    """Calculate the Givens rotation matrix coefficients c and s."""
    r = np.hypot(a, b)
    c = a / r
    s = -b / r
    return c, s

def apply_givens_rotation(A, i, j, n, c, s):
    """Apply Givens rotation to rows and columns of A."""
    for k in range(n):
        # Rotate rows i and j
        temp_i = c * A[i, k] - s * A[j, k]
        temp_j = s * A[i, k] + c * A[j, k]
        A[i, k] = temp_i
        A[j, k] = temp_j

    for k in range(n):
        # Rotate columns i and j
        temp_i = c * A[k, i] - s * A[k, j]
        temp_j = s * A[k, i] + c * A[k, j]
        A[k, i] = temp_i
        A[k, j] = temp_j

def tridiagonalize(matrix):
    """Tridiagonalize a symmetric matrix using Givens rotations."""
    n = matrix.shape[0]
    A = matrix.copy()
    d = np.zeros(n)
    e = np.zeros(n - 1)

    for j in range(n - 1):
        i = j + 1
        a = A[j, j]
        b = A[i, j]
        
        if np.abs(b) < 1e-12:
            continue  # Skip if the subdiagonal element is effectively zero

        # Compute Givens rotation coefficients
        c, s = givens_rotation(a, b)

        # Apply the Givens rotation
        apply_givens_rotation(A, j, i, n, c, s)

        # Store the diagonal and subdiagonal elements
        d[j] = A[j, j]
        e[j] = A[i, j]

    d[-1] = A[-1, -1]  # The last diagonal element

    return d, e

# Example usage
matrix = np.array([[4, 1, 2, 0, 0],
                   [1, 3, 1, 2, 0],
                   [2, 1, 4, 1, 0],
                   [0, 2, 1, 5, 1],
                   [0, 0, 0, 1, 3]], dtype=float)

d, e = tridiagonalize(matrix)

print("Diagonal (d):", d)
print("Subdiagonal (e):", e)
