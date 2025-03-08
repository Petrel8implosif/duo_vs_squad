#include <stdio.h>
#include <math.h>

/**
 * Tridiagonalizes a symmetric band matrix using similarity transformations (Givens rotations).
 *
 * @param A Row-major stored band matrix of size n × (k + 1)
 * @param n Size of the matrix
 * @param k Number of nonzero subdiagonals
 * @param d Output array of size n (diagonal elements of the tridiagonal matrix)
 * @param e Output array of size n-1 (subdiagonal elements)
 */
void tridiagonalize_band(double *A, int n, int k, double *d, double *e) {
    for (int j = 0; j < n - 1; j++) {
        int i = j + 1;

        double a = A[j * (k + 1)];    // Diagonal element
        double b = A[i * (k + 1) + j]; // Subdiagonal element

        if (fabs(b) < 1e-12) continue; // Avoid unnecessary computation

        // Compute Givens rotation coefficients
        double norm = sqrt(a * a + b * b);
        double c = a / norm;
        double s = -b / norm;

        // Apply Givens rotation to zero out subdiagonal
        for (int col = 0; col < n; col++) {
            double temp_j = c * A[j * n + col] - s * A[i * n + col];
            double temp_i = s * A[j * n + col] + c * A[i * n + col];
            A[j * n + col] = temp_j;
            A[i * n + col] = temp_i;
        }

        for (int row = 0; row < n; row++) {
            double temp_j = c * A[row * n + j] - s * A[row * n + i];
            double temp_i = s * A[row * n + j] + c * A[row * n + i];
            A[row * n + j] = temp_j;
            A[row * n + i] = temp_i;
        }

        // Store the diagonal and subdiagonal
        d[j] = A[j * n + j];
        e[j] = A[i * n + j];
    }

    // Store the last diagonal element
    d[n - 1] = A[(n - 1) * n + (n - 1)];
}

int main() {
    int n = 5; // Matrix size
    int k = 2; // Bandwidth

    // Symmetric band matrix stored in row-major format (n × n)
    double A[5][5] = {
        {4, 1, 2, 0, 0},
        {1, 3, 1, 2, 0},
        {2, 1, 4, 1, 0},
        {0, 2, 1, 5, 1},
        {0, 0, 0, 1, 3}
    };

    double d[5], e[4]; // Output arrays

    tridiagonalize_band((double *)A, n, k, d, e);

    printf("\nFinal Tridiagonal Matrix:\n");
    printf("Diagonal (d): ");
    for (int i = 0; i < n; i++) {
        printf("%lf ", d[i]);
    }
    printf("\nSubdiagonal (e): ");
    for (int i = 0; i < n - 1; i++) {
        printf("%lf ", e[i]);
    }
    printf("\n");

    return 0;
}
