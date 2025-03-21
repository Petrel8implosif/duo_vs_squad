#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cblas.h>
#include <complex.h>
#define N 8
#define EPSILON 1e-12
#define SQUARE(x) ((x) * (x))
#define idxQ(i, j, m) ((i) + (j) * (m))
#define idxR(i, j) ((i) + (j) * n)
// version full avec A mais il faudrait la construire à partir de d et e.

int step_qr_tridiag(double *d, double *e, int m, double eps){
    double *A = (double *)calloc(m * m, sizeof(double));
    for (int i = 0; i < m; i++) {
        A[i * m + i] = d[i];
        if (i < m - 1) {
            A[i * m + (i + 1)] = e[i];
            A[(i + 1) * m + i] = e[i];
        }
    }
    printf("Initial matrix A:\n");
    for (int row = 0; row < m; row++) {
        for (int col = 0; col < m; col++) {
            printf("%f ", A[row * m + col]);
        }
        printf("\n");
    }
    double t_nn = A[m * m-1];
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
        double c, s;
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
        double c,s;
        c = cos[i+1];
        s = sin[i+1];
        for (int k = i; k < i+3; k++) {  
                double R_ik = A[k * m + i+1];
                double R_jk = A[k * m + i+2];
                A[k * m + i+1] = c * R_ik - s * R_jk;
                A[k * m + i+2] = s * R_ik + c * R_jk;
    }
    free(cos);
    free(sin);

    }
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
        return m;
    }
    free(A);
    return m - 1;
}


int spectral_decomposition(double *d, double*e, int n, double eps, double* resultat) {
    int i;
    int m = n;
    for (i =0; i < 1000 && m>1; i++) {
        m = step_qr_tridiag(d, e, m, eps);
    }
    return i;
}

int main() {
    double *A = (double *)malloc(N * N * sizeof(double));
    double initial_A[N * N] = {
        1.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        5.0, 9.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 2.0, 8.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 8.0, 5.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 6.0, 3.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 3.0, 7.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 1.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 8.0
    };
    memcpy(A, initial_A, N * N * sizeof(double));
    FILE *file_origin = fopen("A_devoir.txt", "w");
    if (file_origin != NULL) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                fprintf(file_origin, "%f ", A[i * N + j]);
            }
            fprintf(file_origin, "\n");
        }
        fclose(file_origin);
    } else {
        printf("Error opening file for writing.\n");
    }
    double* d = (double *)malloc(N * sizeof(double));
    double* e = (double *)malloc((N - 1) * sizeof(double));
    for (int i = 0; i < N; i++) {
        d[i] = A[i * N + i];
    }
    for (int i = 0; i < N - 1; i++) {
        e[i] = A[i * N + i + 1];
    }

    printf("GOOOOO");
    int iter = spectral_decomposition(d, e, N, EPSILON, d);

    printf("Matrix A after QR decomposition:\n");

    printf("Eigenvalues:\n");
    for (int i = 0; i < N; i++) {
        printf("%f\n", d[i]);
    }
    free(d);
    free(e);
    printf("Number of iterations: %d\n", iter);
    return 0;
}