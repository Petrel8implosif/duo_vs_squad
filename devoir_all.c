#include "devoir_1.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

int iter = 0;
#define N 5
#define EPSILON 1e-12
#define SQUARE(x) ((x) * (x))
#define idxQ(i, j, m) ((i) + (j) * (m))
#define idxR(i, j) ((i) + (j) * n)

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

        for (int i = k + 1; i < n; i++) {
            v[i] = A[i * n + k];
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
        if (vtv < 1e-14) {
            free(v);
            continue;
        }
        double beta = 2.0 / vtv;

        double *p = (double *)calloc(n, sizeof(double));
        for (int i = 0; i < n; i++) {
            for (int j = k + 1; j < n; j++) {
                p[i] += A[i * n + j] * v[j];
            }
            p[i] *= beta;
        }

        double t = 0.0;
        for (int j = k + 1; j < n; j++) {
            t += v[j] * p[j];
        }
        t *= beta;
        t /= 2.0;

        double *q = (double *)malloc(n * sizeof(double));
        for (int i = 0; i < n; i++) {
            q[i] = p[i] - t * v[i];
        }

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                A[i * n + j] -= v[i] * q[j] + q[i] * v[j];
            }
        }

        free(v);
        free(p);
        free(q);
    }
    for (int i = 0; i < n; i++) {
        d[i] = A[i * n + i];
    }
    for (int i = 0; i < n - 1; i++) {
        e[i] = A[(i + 1) * n + i];
    }
}

int step_qr_tridiag(double *d, double *e, int m, double eps) {
    if (m <= 1) return m;
    double delta = (d[m-2] - d[m-1]) / 2.0;
    double bsq = e[m-2] * e[m-2];
    double mu = d[m-1] - ((delta >= 0) ? 1 : -1) * bsq / (fabs(delta) + sqrt(delta * delta + bsq));
    
    for (int k = 0; k < m; k++) 
        d[k] -= mu;
    
    for (int k = 1; k < m; k++) {
        if (fabs(e[k]) < eps * (fabs(d[k-1]) + fabs(d[k]))) {
            e[k] = 0.0;
        }
    }
    for (int k = m - 1; k >= 1; k--) {
        if (fabs(e[k]) > eps * (fabs(d[k-1]) + fabs(d[k]))) 
            break;
        m--;
    }
    for (int k = 0; k < m; k++) 
        d[k] += mu;
    return m;
}

int qr_eigs_band(double *A, int n, int k, double eps, int max_iter, double *d) {
    double *dr = (double *)calloc(n, sizeof(double));
    double *e = (double *)calloc(n-1, sizeof(double));
    tridiagonalize_full(A, n, k, dr, e);

    int m = n;
    for (int i = 0; i < n; i++) {
        m = step_qr_tridiag(dr, e, m, eps);
        if (m <= 0) break;
        d[n - 1 - i] = dr[m - 1];
    }

    free(dr);
    free(e);
    return (iter > max_iter) ? -1 : iter;
}

int main() {
    double A[N * N] = {4, 1, 2, 3, 1, 1, 5, 2, 2, 2, 2, 1, 6, 3, 3, 3, 2, 3, 7, 4, 4, 3, 3, 4, 8};
    double d[N];
    qr_eigs_band(A, N, 2, EPSILON, 100, d);
    printf("Eigenvalues:\n");
    for (int i = 0; i < N; i++) {
        printf("%lf\n", d[i]);
    }
    return 0;
}
