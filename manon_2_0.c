#include "devoir_1.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#define N 5
#define SQUARE(x) ((x) * (x))
#define idxQ(i, j, m) ((i) + (j) * (m))
#define idxR(i, j) ((i) + (j) * n)


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
    e[k+1] = A[(k + 1) *n + k];
    }
    e[0] = 0.0;
}

double sign(double x) {
    return (x > 0) - (x < 0);
}

int isz(double x) {
    return fabs(x) < 1e-10;
}

int step_qr_tridiag(double *d, double *e, int m, double eps) {
    // Compute wilkinson shift
    double const delta = (d[m-2] - d[m-1]) / 2.0;
    double const bsq = e[m-2]*e[m-2];
    double const mu = d[m-1] - sign(delta) * bsq / (fabs(delta) + sqrt(delta * delta + bsq));

    // Subtracting the shift
    for (int k = 0; k < m; k++) 
        d[k] -= mu;

    double inv_hyp;
    double *c = malloc(m * sizeof(double)); // cosines
    double *s = malloc(m * sizeof(double)); // sines

    // Compute QR decomposition using Givens rotations
    double ekp = e[1], ek, dk, ck, sk;
    for (int k = 1; k < m; k++) {
        dk = d[k-1];

        inv_hyp = 1.0 / sqrt(dk * dk + ekp * ekp);
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
        if (isz(e[k])) e[k] = 0.0;

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

int qr_eigs_full(double *A, int n, int k, double eps, int max_iter, double *d){
    double *e = (double *)calloc(n, sizeof(double));
    tridiagonalize_full(A, n, k, d, e);
    int m = n;
    int i = 0;
    int m_temp = m;
    for (i = 0; i < max_iter && m > 1; i++) {
        m = step_qr_tridiag(d, e, m, eps);
        if (i == 100) {
            break;
        }
    }
    if (i == max_iter) {
        return -1;
    }
    free(e);
    return i;
}

double *create_matrix(int nx, int ny, double lx, double ly, int storage) {
    int lda, k;
    int size = nx * ny;
    double dx2 = SQUARE(lx / (nx + 1));
    double dy2 = SQUARE(ly / (ny + 1));
    double alpha, beta, gamma;
    double *L;

    k = nx;
    alpha = 1. / dx2;
    beta = 1. / dy2;
    gamma = 2 * (alpha + beta);

    if (storage == 2) {
        lda = k + 1;
        L = (double *)calloc(size * lda, sizeof(double));
        for (int l = 0; l < size; l++) {
            L[l * lda + k - k] = -beta;
            if (l % k != 0)
                L[l * lda + k - 1] = -alpha;
            L[l * lda + k - 0] = +gamma;
        }
    } else if (storage == 1) {
        lda = 2 * k + 1;
        L = (double *)calloc(size * lda, sizeof(double));
        for (int l = 0; l < size; l++) {
            L[l * lda + k - k] = -beta;
            if (l % k != 0)
                L[l * lda + k - 1] = -alpha;
            L[l * lda + k + 0] = +gamma;
            if (l % k != k - 1)
                L[l * lda + k + 1] = -alpha;
            L[l * lda + k + k] = -beta;
        }
    } else {
        lda = size;
        L = (double *)calloc(size * lda, sizeof(double));
        for (int idx, i = 0; i < ny; i++) {
            for (int j = 0; j < nx; j++) {
                idx = i * k + j;
                L[idx * lda + idx] = gamma;
                if (0 < i)
                    L[idx * lda + idx - k] = -beta;
                if (i < ny - 1)
                    L[idx * lda + idx + k] = -beta;
                if (0 < j)
                    L[idx * lda + idx - 1] = -alpha;
                if (j < nx - 1)
                    L[idx * lda + idx + 1] = -alpha;
            }
        }
    }
    return L;
}

int main() {
    /*oui();
    return 0;*/
    double lx = 10.0;
    double ly = 10.0;
    int nx = 2;
    int ny = 2;
    double *A;
    double *d = (double *)calloc(nx*ny, sizeof(double));
    A = create_matrix(nx, ny, lx, ly, 0);
    FILE *file = fopen("A_devoir.txt", "w");
    if (file != NULL) {
        for (int i = 0; i < nx*ny; i++) {
            for (int j = 0; j < nx*ny; j++) {
                fprintf(file, "%f ", A[i * nx*ny + j]);
            }
            fprintf(file, "\n");
        }
        fclose(file);
    } else {
        printf("Error opening file!\n");
    }
    int n = nx * ny;

    int k = 0;
        for (int row = 0; row < n; row++) {
            for (int col = 0; col < n; col++) {
                if (A[row * n + col] != 0) {
                    int band_width = abs(row - col);
                    if (band_width > k) {
                        k = band_width;
                    }
                }
            }
        }
    double eps = 1e-12;
    int max_iter = 10000;
    printf("k = %d\n", k);
    int iteration = qr_eigs_full(A, n, k, eps, max_iter,d);
    printf("Number of iterations");
    printf("Number of iterations: %d\n", iteration);
    printf("Eigenvalues:\n");
    for (int i = 0; i < n; i++) {
        printf("%f\n", d[i]);
    }
    free(d);
    return 0;
}