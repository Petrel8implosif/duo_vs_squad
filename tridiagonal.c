#include <stdio.h>
#include <math.h>

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

// Function to apply Givens rotation to a row-major stored matrix
void apply_givens(double *A, int n, int k, double c, double s, int i, int j) {
    for (int col = 0; col < n; col++) {
        double temp = c * A[i * n + col] - s * A[j * n + col];
        A[j * n + col] = s * A[i * n + col] + c * A[j * n + col];
        A[i * n + col] = temp;
    }

    for (int row = 0; row < n; row++) {
        double temp = c * A[row * n + i] - s * A[row * n + j];
        A[row * n + j] = s * A[row * n + i] + c * A[row * n + j];
        A[row * n + i] = temp;
    }
}

// Function to tridiagonalize a symmetric band matrix
void tridiagonalize_full(double *A, int n, int k, double *d, double *e) {
    for (int j = 0; j < n - 2; j++) {
        for (int i = n - 1; i > j + 1; i--) {
            if (fabs(A[i * n + j]) > 1e-10) {
                double c, s;
                givens_rotation(A[(i - 1) * n + j], A[i * n + j], &c, &s);
                apply_givens(A, n, k, c, s, i - 1, i);
            }
        }
    }
    
    for (int i = 0; i < n; i++) {
        d[i] = A[i * n + i];
        if (i < n - 1) {
            e[i] = A[i * n + i + 1];
        }
    }
    e[4] = 0.0;

}

// Main function
int main() {
    int n = 5;  // Matrix size
    int k = 2;  // Bandwidth (assuming a tridiagonal matrix)
    double A[5 * 5] = {
        4.0, 1.0, 2.0, 0.0, 0.0,
        1.0, 3.0, 1.0, 2.0, 0.0,
        2.0, 1.0, 4.0, 1.0, 0.0,
        0.0, 2.0, 1.0, 5.0, 1.0,
        0.0, 0.0, 0.0, 1.0, 3.0
    };
    double d[5], e[5];

    tridiagonalize_full(A, n, k, d, e);

    printf("\nFinal Tridiagonal Matrix:\n");
    printf("Diagonal (d): ");
    for (int i = 0; i < n; i++) {
        printf("%lf ", d[i]);
    }
    printf("\nSubdiagonal (e): ");
    for (int i = 0; i < n ; i++) {
        printf("%lf ", e[i]);
    }
    printf("\n");

    return 0;
}
