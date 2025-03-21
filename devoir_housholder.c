#include "devoir_1.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int iter = 0;
#define N 5
#define EPSILON 1e-12
#define SQUARE(x) ((x) * (x))
#define idxQ(i, j, m) ((i) + (j) * (m))
#define idxR(i, j) ((i) + (j) * n)

// Function to compute the sign of a number
int sign(double x) {
    return (x > 0) - (x < 0);
}

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
    for (int i = 0; i < n-1; i++) {
        e[i+1] = A[(i + 1) * n + i];
    }
    //print e and d
    printf("Diagonal (d): ");
    for (int i = 0; i < n; i++) {
        printf("%lf ", d[i]);
    }
    printf("\nSubdiagonal (e): ");
    for (int i = 0; i < n; i++) {
        printf("%lf ", e[i]);
    }
    printf("\n");
}

int step_qr_tridiag(double *d, double *e, int m, double eps) {
    // Compute wilkinson shift
    double const delta = (d[m-2] - d[m-1]) / 2.;
    double const bsq = e[m-2]*e[m-2];
    double const mu = d[m-1] - sign(delta) * bsq / (fabs(delta) + sqrt(delta*delta + bsq));

    for (int k = 0; k < m; k++) 
        d[k] -= mu;

    double inv_hyp;
    double *c = malloc(m * sizeof(double)); // cosines
    double *s = malloc(m * sizeof(double)); // sines

    // Compute QR decomposition using Givens rotations
    double ekp = e[1], ek, dk, ck, sk;
    for (int k = 1; k < m; k++) {
        dk = d[k-1];

        inv_hyp = 1. / sqrt(dk*dk + ekp*ekp);
        c[k]   = ck = dk  * inv_hyp;
        s[k]   = sk = ekp * inv_hyp; 

        // A = G^TA
        d[k-1] =  ck * dk + sk * ekp;
        ek = e[k];
        e[k]   =  ck * ek + sk * d[k];
        if (k < m-1) {
            ekp = e[k+1];
            e[k+1] = ck;
        }

        d[k] = -sk * ek + ck * d[k];
    }

    // A = AG
    for (int k = 1; k < m; k++) {
        dk = d[k-1];
        ek = e[k];
        ck = c[k];
        sk = s[k];

        d[k-1] = ck * dk + sk * ek;
        e[k]   = sk * d[k];
        d[k]  = ck;
    }

    // Adding back the shift
    for (int k = 0; k < m; k++) 
        d[k] += mu;

    // Zeroing out small subdiagonal elements
    for (int k = 1; k < m; k++)
        if (fabs(e[k]) < EPSILON) e[k] = 0.;

    // Update active dimension
    for (int k = m-1; k >= 1; k--) {
        if (fabs(e[k]) > eps * (fabs(d[k-1]) + fabs(d[k]))) 
            break;
        m--;
    }

    free(c);
    free(s);

    return m;
}

// Corrected QR eigensolver
int qr_eigs_full(double *A, int n, int k, double eps, int max_iter, double *d) {
    double *e = malloc((n) * sizeof(double));
    double *diag = malloc(n * sizeof(double));
    
    tridiagonalize_full(A, n, k, diag, e);
    
    int iter = 0;
    int m = n;
    
    while (m > 1 && iter < max_iter) {
        int new_m = step_qr_tridiag(diag, e, m, eps);
        if (new_m == m) break;
        m = new_m;
        iter++;
    }
    
    // Copy and sort eigenvalues
    memcpy(d, diag, n*sizeof(double));
    for (int i = 0; i < n; i++) {
        for (int j = i+1; j < n; j++) {
            if (d[i] < d[j]) {
                double tmp = d[i];
                d[i] = d[j];
                d[j] = tmp;
            }
        }
    }
    
    free(e);
    free(diag);
    return (iter >= max_iter) ? -1 : iter;
}

int main() {
    double A[N*N] = {
        1.0, 2.0, 0.0, 0.0, 0.0,
        2.0, 3.0, 4.0, 0.0, 0.0,
        0.0, 4.0, 5.0, 6.0, 0.0,
        0.0, 0.0, 6.0, 7.0, 8.0,
        0.0, 0.0, 0.0, 8.0, 9.0
    };
    double d[N];
    int result = qr_eigs_full(A, N, 2, EPSILON, 100, d);
    
    if (result == -1) {
        printf("Convergence failed\n");
    } else {
        printf("Eigenvalues:\n");
        for (int i = 0; i < N; i++) {
            printf("%.6f\n", d[i]);
        }
    }
    return 0;
}
