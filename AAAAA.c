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
    printf("Tridiagonalized matrix A:\n");
    for (int row = 0; row < n; row++) {
        for (int col = 0; col < n; col++) {
            printf("%f ", A[row * n + col]);
        }
        printf("\n");
    }
    
    for (int i = 0; i < n-1; i++) {
        d[i] = A[i * n + i];
        e[i] = A[i * n + i + 1];
    }
    d[n-1] = A[(n-1) * n + (n-1)];

}

int step_qr_tridiag(double *d, double *e, int m, double eps){
    double t1 = d[m-2];
    double t2 = d[m-1];
    double tn = e[m-2];
    double delta = (t1 - t2)/2.0;
    double mu = t2 - copysign(tn*tn/(fabs(delta) + hypot(delta, tn)), delta);

    // Application du shift
    for(int i = 0; i < m; i++) 
        d[i] -= mu;

    double *cos = (double *)calloc(m-1, sizeof(double));
    double *sin = (double *)calloc(m-1, sizeof(double));
    double temp = e[0];
    double c,s;

    for (int i = 0; i < m - 1; i++) {
        // gauche
        double a = d[i];
        double b = e[i];
        double hypot = sqrt(a * a + b * b);
        if(hypot < 1e-10){
            cos[i] = 1.0;
            sin[i] = 0.0;
        }
        else{
            cos[i] = a / hypot;
            sin[i] = -b / hypot;
        }
        c = cos[i];
        s = sin[i];
        d[i] = hypot;
        e[i] = c * temp - s * d[i+1];
        d[i+1] = s * temp + c * d[i+1];
        if (i < m - 2){
            temp = c * e[i+1];
        }
    }
    double *A = (double *)calloc(m * m, sizeof(double));
    for (int i = 0; i < m; i++) {
        A[i * m + i] = d[i];
        if (i < m - 1) {
            A[i * m + (i + 1)] = e[i];
            A[(i + 1) * m + i] = e[i];
        }
    }
    printf("Gauche A:\n");
    for (int row = 0; row < m; row++) {
        for (int col = 0; col < m; col++) {
            printf("%f ", A[row * m + col]);
        }
        printf("\n");
    }
    free(A);

    for(int i = 0; i < m-1; i++){
        //droite
        c = cos[i];
        s = sin[i];
        d[i] = c * d[i] - s * e[i];
        temp = e[i];
        e[i] = c * e[i] - s * d[i+1];
        d[i+1] = s * temp + c * d[i+1];

    }
    double *B = (double *)calloc(m * m, sizeof(double));
    for (int i = 0; i < m; i++) {
        B[i * m + i] = d[i];
        if (i < m - 1) {
            B[i * m + (i + 1)] = e[i];
            B[(i + 1) * m + i] = e[i];
        }
    }
    printf("Droite A:\n");
    for (int row = 0; row < m; row++) {
        for (int col = 0; col < m; col++) {
            printf("%f ", B[row * m + col]);
        }
        printf("\n");
    }
    free(B);
    
    free(cos);
    free(sin);
    for(int i = 0; i < m; i++){
        d[i] += mu;
    }
    printf("d: ");
    for (int i = 0; i < m; i++) {
        printf("%f ", d[i]);
    }
    printf("\n");

    printf("e: ");
    for (int i = 0; i < m - 1; i++) {
        printf("%f ", e[i]);
    }
    printf("\n");
  

    if(fabs(e[m-2]) > eps * (fabs(d[m-2]) + fabs(d[m-1]))){
        return m;
    }
    return m - 1;
}

int qr_eigs_band(double *A, int n, int k, double eps, int max_iter, double *d){
    double *e = (double *)calloc(n-1, sizeof(double));
    tridiagonalize_full(A, n, k, d, e);
    int m = n;
    int i = 0;
    int m_temp = m;
    for (i = 0; i < max_iter && m > 1; i++) {
        m = step_qr_tridiag(d, e, m, eps);
        if(i > 1000){
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
    int nx = 4;
    int ny = 4;
    double *A;
    double *dr = (double *)calloc(nx*ny, sizeof(double));
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
    int iteration = qr_eigs_band(A, n, k, eps, max_iter,dr);
    printf("Number of iterations");
    printf("Number of iterations: %d\n", iteration);
    printf("Eigenvalues:\n");
    for (int i = 0; i < n; i++) {
        printf("%f\n", dr[i]);
    }
    free(dr);
    return 0;
}