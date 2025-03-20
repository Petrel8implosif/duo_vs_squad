#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define TOL 1e-12

// Adds a column of zeros on the left of the matrix
double *add_column_of_zeroes(double *A_compact, int n, int m) {
    double *A_exp = (double *)malloc(n * (m + 1) * sizeof(double));
    if (A_exp == NULL) {
        return NULL;
    }
    for (int i = 0; i < n; i++) {
        A_exp[i * (m + 1)] = 0.0; // Add zero at the beginning of each row
        for (int j = 0; j < m; j++) {
            A_exp[i * (m + 1) + (j + 1)] = A_compact[i * m + j];
        }
    }
    return A_exp;
}

void print_matrix(double *A, int n, int m) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            printf("%.6f ", A[i * m + j]);
        }
        printf("\n");
    }
}

void givens_rotation(double a, double b, double *c, double *s) {
    double r = hypot(a, b);
    if (r < TOL) {
        *c = 1.0;
        *s = 0.0;
    } else {
        *c = a / r;
        *s = -b / r;
    }
}

// Apply Givens rotation to rows and columns of a full matrix
void apply_givens(double *A, int n, double c, double s, int i, int j) {
    // Apply to rows i and j
    for (int col = 0; col < n; col++) {
        double temp_i = c * A[i * n + col] - s * A[j * n + col];
        double temp_j = s * A[i * n + col] + c * A[j * n + col];
        A[i * n + col] = temp_i;
        A[j * n + col] = temp_j;
    }
    // Apply to columns i and j
    for (int row = 0; row < n; row++) {
        double temp_i = c * A[row * n + i] - s * A[row * n + j];
        double temp_j = s * A[row * n + i] + c * A[row * n + j];
        A[row * n + i] = temp_i;
        A[row * n + j] = temp_j;
    }
}

void tridiagonalize_band(double *A_compact, int n, int m, double *d, double *e) {
    double *A_full = add_column_of_zeroes(A_compact, n, m);
    if (A_full == NULL) return;

    // Tridiagonalize using Givens rotations on the full matrix
    for (int j = 0; j < n - 2; j++) {
        for (int i = n - 1; i > j + 1; i--) {
            if (fabs(A_full[i * n + j]) > TOL) {
                double c, s;
                givens_rotation(A_full[(i-1)*n + j], A_full[i*n + j], &c, &s);
                apply_givens(A_full, n, c, s, i-1, i);
            }
        }
    }

    print_matrix(A_full, n, n);

    // Extract diagonal and subdiagonal
    for (int i = 0; i < n; i++) {
        d[i] = A_full[i * n + i];
        if (i < n - 1) {
            e[i] = A_full[i * n + (i + 1)];
        }
    }

    free(A_full);
}



int main() {
    int n = 5, m = 2;
    // Correct compact storage representing the symmetric band matrix
    double A_compact[] = {

            0.0, 0.0, 4.0,  // Row 0
            0.0, 1.0, 3.0,  // Row 1
            2.0, 1.0, 4.0,  // Row 2
            2.0, 1.0, 5.0,  // Row 3
            0.0, 1.0, 3.0   // Row 4

    };

    double d[n], e[n-1];
    tridiagonalize_band(A_compact, n, m, d, e);

    printf("Diagonal (d): ");
    for (int i = 0; i < n; i++) printf("%.6f ", d[i]);
    printf("\nSubdiagonal (e): ");
    for (int i = 0; i < n-1; i++) printf("%.6f ", e[i]);
    printf("\n");

    return 0;
}