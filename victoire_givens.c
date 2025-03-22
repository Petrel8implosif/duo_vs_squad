#include "devoir_1.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#define SQUARE(x) ((x) * (x))
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
    
    for (int i = 0; i < n-1; i++) {
        d[i] = A[i * n + i];
        
        e[i+1] = A[i * n + i + 1];
    }
    d[n-1] = A[(n-1) * n + (n-1)];

}


int step_qr_tridiag(double *d, double *e, int m, double eps){
    double *A = (double *)calloc(m * m, sizeof(double));
    for (int i = 0; i < m; i++) {
        A[i * m + i] = d[i];
        if (i < m - 1) {
            A[i * m + (i + 1)] = e[i];
            A[(i + 1) * m + i] = e[i];
        }
    }
    double t_nn = A[m * m-1];  // calcul du shift de wilkinson
    double t_n1n1 = A[(m-2) * m + (m-2)];
    double t_n1n = A[(m-1) * m + (m-2)];

    double dd = (t_n1n1 - t_nn) / 2.0;
    if (dd == 0) {
        double mu = t_nn - fabs(t_n1n);
    }
    double sign_d = (dd > 0) - (dd < 0); // Equivalent à sign(d)
    double mu = t_nn + dd - sign_d * sqrt(dd * dd + t_n1n * t_n1n);

 
    

    for (int i = 0; i < m; i++) {
        A[i * m + i] -= mu;
    }

    double *cos = (double *)calloc(m-1, sizeof(double));
    double *sin = (double *)calloc(m-1, sizeof(double));

    for (int i = 0; i < m - 1; i++) {
        double c = 1.0;
        double s = 0.0;
        double a = A[i * m + i];
        double b = A[(i+1) * m + i];
        double hypot = sqrt(a * a + b * b);
        if(hypot < 1e-10){
            c = 1.0;
            s = 0.0;
        }
        else{
            c = a / hypot;
            s = -b / hypot;
        }
        cos[i] = c;
        sin[i] = s;


        for (int k = i; k < i+3 && k < m; k++) {  // gauche 
            if(k == m || k == i+2) {
            double R_ik = A[i * m + k];
            double R_jk = A[(i+1) * m + k];
            A[i * m + k] = c * R_ik - s * R_jk;
            A[(i+1) * m + k] = c * R_jk;
            }
            else {
            double R_ik = A[i * m + k];
            double R_jk = A[(i+1) * m + k];
            A[i * m + k] = c * R_ik - s * R_jk;
            A[(i+1) * m + k] = s * R_ik + c * R_jk;
            }
        }
    }

        double c1 = cos[0];
        double s1 = sin[0];
    
        double R_ik = A[0];
        double R_jk = A[1];
        A[0] = c1 * R_ik - s1 * R_jk;
        A[1] = s1 * R_ik + c1 * R_jk;
        
        double rR_ik = A[m];
        double rR_jk = A[m + 1];
        A[m] = c1 * rR_ik - s1 * rR_jk;
        A[m + 1] = c1 * rR_jk;
    



    for(int i = 0; i < m-2; i++){
        //droite
            double c = 1.0;
            double s = 0.0;
            c = cos[i+1];
            s = sin[i+1];
            for (int k = i; k < i+3; k++) {  
                    double R_ik = A[k * m + i+1];
                    double R_jk = A[k * m + i+2];
                    A[k * m + i+1] = c * R_ik - s * R_jk;
                    A[k * m + i+2] = s * R_ik + c * R_jk;
        }
    }
    
    free(cos);
    free(sin);
    for(int i = 0; i < m; i++){
        A[i * m + i] += mu;
    }
    for (int i = 0; i < m; i++) {
        d[i] = A[i * m + i];
    }
    for (int i = 0; i < m - 1; i++) {
        e[i] = A[i * m + (i+1)];
    }

    if(fabs(A[(m-1) * m + (m-2)]) > eps * (fabs(A[(m-1) * m + (m-1)]) + fabs(A[(m-2) * m + (m-2)]))){
        free(A);
        return m;
    }
    free(A);
    return m - 1;
}

int qr_eigs_full(double *A, int n, int k, double eps, int max_iter, double *d){
    double *dr = (double *)calloc(n, sizeof(double));
    double *e = (double *)calloc(n-1, sizeof(double));
    tridiagonalize_full(A, n, k, dr, e);
    int m = n;
    int i = 0;
    int m_temp = m;
    for (i = 0; i < max_iter && m > 1; i++) {
        m = step_qr_tridiag(dr, e, m, eps);
        if (m != m_temp) {
            d[m] = dr[m];
        }
        m_temp = m;
    }
    d[0] = dr[0];
    if (i == max_iter) {
        return -1;
    }
    free(dr);
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
    int nx = 10;
    int ny = 10;
    double *E;
    double *d = (double *)calloc(nx*ny, sizeof(double));
    E = create_matrix(nx, ny, lx, ly, 0);
    FILE *file = fopen("A_devoir.txt", "w");
    if (file != NULL) {
        for (int i = 0; i < nx*ny; i++) {
            for (int j = 0; j < nx*ny; j++) {
                fprintf(file, "%f ", E[i * nx*ny + j]);
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
                if (E[row * n + col] != 0) {
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
    int iteration = qr_eigs_full(E, n, k, eps, max_iter,d);
    printf("Number of iterations");
    printf("Number of iterations: %d\n", iteration);
    printf("Eigenvalues:\n");
    for (int i = 0; i < n; i++) {
        printf("%f\n", d[i]);
    }
    free(d);
    free(E);
    return 0;
}