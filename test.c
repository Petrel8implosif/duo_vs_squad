#include <stdio.h>
#include <math.h>

// Function to compute Givens rotation coefficients
void compute_givens(double a, double b, double *c, double *s) {
    if (b == 0) {
        *c = 1.0;
        *s = 0.0;
    } else {
        double r = sqrt(a * a + b * b);
        *c = a / r;
        *s = -b / r;
    }
}

// Function to apply Givens rotation to a band matrix
void apply_givens_band(double *A, int n, int k, int i, int j, double c, double s) {
    if (i >= n || j >= n || i >= j) return;

    // Update row elements within band structure
    for (int col = 0; col < k; col++) {
        if (i + col < n && j + col < n) { 
            double temp = c * A[i * k + col] - s * A[j * k + col];
            A[j * k + col] = s * A[i * k + col] + c * A[j * k + col];
            A[i * k + col] = temp;
        }
    }

    // Update column elements within band structure
    for (int row = i; row < j; row++) {
        if (row + k < n) {
            double temp = c * A[row * k + (k - 1)] - s * A[(row + 1) * k + (k - 2)];
            A[(row + 1) * k + (k - 2)] = s * A[row * k + (k - 1)] + c * A[(row + 1) * k + (k - 2)];
            A[row * k + (k - 1)] = temp;
        }
    }
}

// Function to tridiagonalize a symmetric band matrix stored in compact form
void tridiagonalize_band(double *A, int n, int k, double *d, double *e) {
    for (int i = 0; i < n - 2; i++) {
        for (int j = k - 1; j > 0; j--) {
            if (fabs(A[i * k + j]) > 1e-10) { 
                double c, s;
                compute_givens(A[i * k], A[i * k + j], &c, &s);
                apply_givens_band(A, n, k, i, i + j, c, s);
            }
        }
    }

    // Extract the diagonal and subdiagonal elements
    for (int i = 0; i < n; i++) {
        d[i] = A[i * k]; // Diagonal elements
        if (i < n - 1) {
            e[i] = A[i * k + 1]; // Subdiagonal elements
        }
    }
}

// Example usage
int main() {
    int n = 5; // Matrix size
    int k = 3; // Bandwidth

    // Example symmetric band matrix in compact form
    double A[] = {
        4.0, 1.0, 2.0, // format: diagonal, subdiagonal, sub-subdiagonal
        3.0, 1.0, 2.0,
        4.0, 1.0, 0.0,
        5.0, 1.0, 0.0,
        3.0, 0.0, 0.0
    };

    double d[n], e[n - 1]; // Arrays for diagonal and subdiagonal elements

    tridiagonalize_band(A, n, k, d, e);

    // Output results
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