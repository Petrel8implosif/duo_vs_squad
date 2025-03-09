#include <stdio.h>
#include <math.h>

#define N 5  // Matrix size

// Function to compute Givens rotation
void givens_rotation(double a, double b, double *c, double *s) {
    double norm = sqrt(a * a + b * b);
    if (norm < 1e-10) {
        *c = 1.0;
        *s = 0.0;
    } else {
        *c = a / norm;
        *s = -b / norm;
    }
}

// Function to apply Givens rotation to a matrix
void apply_givens(double A[N][N], double c, double s, int i, int j) {
    for (int k = 0; k < N; k++) {
        double temp = c * A[i][k] - s * A[j][k];
        A[j][k] = s * A[i][k] + c * A[j][k];
        A[i][k] = temp;
    }

    for (int k = 0; k < N; k++) {
        double temp = c * A[k][i] - s * A[k][j];
        A[k][j] = s * A[k][i] + c * A[k][j];
        A[k][i] = temp;
    }
}

// Function to tridiagonalize a symmetric matrix using Givens rotations
void tridiagonalize(double A[N][N], double d[N], double e[N - 1]) {
    for (int j = 0; j < N - 2; j++) {  // Iterate through columns
        for (int i = N - 1; i > j + 1; i--) {  // Work from the bottom row upward
            if (fabs(A[i][j]) > 1e-10) { // If significant, apply Givens rotation
                double c, s;
                givens_rotation(A[i - 1][j], A[i][j], &c, &s);

                // Apply Givens rotation to zero out A[i][j]
                apply_givens(A, c, s, i - 1, i);
            }
        }
    }

    // Extract diagonal and subdiagonal elements
    for (int i = 0; i < N; i++) {
        d[i] = A[i][i];  // Diagonal elements
        if (i < N - 1) {
            e[i] = A[i][i + 1];  // First subdiagonal elements
        }
    }
}

// Main function
int main() {
    // Example symmetric matrix stored in full N x N format
    double A[N][N] = {
        {4.0, 1.0, 2.0, 0.0, 0.0},
        {1.0, 3.0, 1.0, 2.0, 0.0},
        {2.0, 1.0, 4.0, 1.0, 0.0},
        {0.0, 2.0, 1.0, 5.0, 1.0},
        {0.0, 0.0, 0.0, 1.0, 3.0}
    };

    double d[N], e[N - 1];  // Output arrays for tridiagonal matrix

    // Perform tridiagonalization
    tridiagonalize(A, d, e);

    // Output the final tridiagonal matrix
    printf("\nFinal Tridiagonal Matrix:\n");
    printf("Diagonal (d): ");
    for (int i = 0; i < N; i++) {
        printf("%lf ", d[i]);
    }
    printf("\nSubdiagonal (e): ");
    for (int i = 0; i < N - 1; i++) {
        printf("%lf ", e[i]);
    }
    printf("\n");

    // Output modified matrix
    printf("\nModified Matrix A:\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%lf ", A[i][j]);
        }
        printf("\n");
    }

    return 0;
}