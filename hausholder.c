#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Function to compute the norm of a vector
double vector_norm(double *v, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += v[i] * v[i];
    }
    return sqrt(sum);
}

// Function to perform Householder transformation
void tridiagonalize_full(double *A, int n, int k, double *d, double *e) {
    for (int k = 0; k < n - 2; k++) {
        double *v = (double *)calloc(n, sizeof(double));
        double norm_x;

        // Compute the Householder vector
        for (int i = k + 1; i < n; i++) {
            v[i] = A[i *n + k];
        }
        norm_x = vector_norm(&v[k + 1], n - k - 1);
        if (v[k + 1] >= 0) {
            v[k + 1] += norm_x;
        } else {
            v[k + 1] -= norm_x;
        }

        double vtv = 0.0;
        for (int i = k + 1; i < n; i++) {
            vtv += v[i] * v[i];
        }
        if (vtv < 1e-10) { // Avoid division by zero
            free(v);
            continue;
        }
        double beta = 2.0 / vtv;

        // Compute p = beta * A * v
        double *p = (double *)calloc(n, sizeof(double));
        for (int i = 0; i < n; i++) {
            for (int j = k + 1; j < n; j++) {
                p[i] += A[i *n + j] * v[j];
            }
            p[i] *= beta;
        }

        // Compute t = beta * (v ⋅ p) / 2
        double t = 0.0;
        for (int j = k + 1; j < n; j++) {
            t += v[j] * p[j];
        }
        t *= beta;
        t /= 2.0;

        // Compute q = p - t * v
        double *q = (double *)malloc(n * sizeof(double));
        for (int i = 0; i < n; i++) {
            q[i] = p[i] - t * v[i];
        }

        // Update A: A = A - v*q^T - q*v^T
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                A[i *n + j] -= v[i] * q[j] + q[i] * v[j];
            }
        }

        

        free(v);
        free(p);
        free(q);
    }
    // update d and e with the diagonal and subdiagonal elements respectively
    for (int k = 0; k < n; k++) {
    d[k] = A[k *n + k];
    e[k] = A[(k + 1) *n + k];
    }
}

// Helper function to print a A
void print_A(double **A, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%10.4f ", A[i][j]);
        }
        printf("\n");
    }
}

int main() {
    int n = 5;  // A size
    int k = 2;  // Bandwidth (assuming a tridiagonal A)
    double A[5 * 5] = {
        4.0, 1.0, 2.0, 0.0, 0.0,
        1.0, 3.0, 1.0, 2.0, 0.0,
        2.0, 1.0, 4.0, 1.0, 0.0,
        0.0, 2.0, 1.0, 5.0, 1.0,
        0.0, 0.0, 0.0, 1.0, 3.0
    };
    double d[5], e[5];

    tridiagonalize_full(A, n, k, d, e);

    printf("\nFinal Tridiagonal A:\n");
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

