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
    for (int i = 0; i < n; i++) {
        for (int j = 2; j < k + 1; j++) { // Start from the 2nd superdiagonal
            double elim = A[i * (k + 1) + j];
            if (fabs(elim) > 1e-10) {
                // Perform Givens rotation
                double a = A[i * (k + 1) + j - 1];
                double b = elim;
                double c, s;
                if (fabs(b) > fabs(a)) {
                    double tau = -a / b;
                    s = 1 / sqrt(1 + tau * tau);
                    c = s * tau;
                } else {
                    double tau = -b / a;
                    c = 1 / sqrt(1 + tau * tau);
                    s = c * tau;
                }

                // Update diagonal elements during Givens rotation
                A[i * (k + 1) + j - 1] = c * a - s * b;
                A[i * (k + 1) + j] = 0.0;

                for (int l = j + 1; l < k + 1; l++) {
                    double temp = c * A[i * (k + 1) + l - 1] - s * A[i * (k + 1) + l];
                    A[i * (k + 1) + l] = s * A[i * (k + 1) + l - 1] + c * A[i * (k + 1) + l];
                    A[i * (k + 1) + l - 1] = temp;
                }

                // Update columns below the diagonal
                for (int l = i + 1; l < n; l++) {
                    double temp = c * A[l * (k + 1) + j - 1] - s * A[l * (k + 1) + j];
                    A[l * (k + 1) + j] = s * A[l * (k + 1) + j - 1] + c * A[l * (k + 1) + j];
                    A[l * (k + 1) + j - 1] = temp;
                }
            }
        }

        // Extract and update diagonal and subdiagonal elements
        d[i] = A[i * (k + 1)];
        if (i < n - 1) {
            e[i] = A[i * (k + 1) + 1];
        }
    }
}



int main() {
    int n = 5; // Matrix size
    int k = 2; // Bandwidth

    // Symmetric band matrix stored in row-major format (n × (k + 1))
    double A[] = {
        4, 1, 2,  // Row 0: diagonal, 1st superdiagonal, 2nd superdiagonal
        3, 1, 2,  // Row 1: diagonal, 1st superdiagonal, 2nd superdiagonal
        4, 1, 0,  // Row 2: diagonal, 1st superdiagonal, 2nd superdiagonal
        5, 1, 0,  // Row 3: diagonal, 1st superdiagonal, 2nd superdiagonal
        3, 0, 0   // Row 4: diagonal, (rest are padding)
    };

    double d[5], e[4]; // Output arrays

    tridiagonalize_band(A, n, k, d, e);

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

    printf("\nModified Matrix A:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k + 1; j++) {
            printf("%lf ", A[i * (k + 1) + j]);
        }
        printf("\n");
    }

    return 0;
}
