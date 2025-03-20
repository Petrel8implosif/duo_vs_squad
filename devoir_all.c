#include "devoir_1.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
int iter = 0;
#define N 5
#define EPSILON 1e-9
#define SQUARE(x) ((x) * (x))
#define idxQ(i, j, m) ((i) + (j) * (m))
#define idxR(i, j) ((i) + (j) * n)


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//-------------
// Givens rotation: we want c and s so that
// [ c  -s ] [ a ] = [ r ]
// [ s   c ] [ b ]   [ 0 ]
// Note: here we choose s = -b/norm so that the element is zeroed.
void givens_rotation(double a, double b, double *c, double *s) {
    double norm = sqrt(a * a + b * b);
    if (norm < 1e-10) {
        *c = 1.0;
        *s = 0.0;
    } else {
        *c = a / norm;
        *s = -b / norm; // Corrected so that s zeroes the element.
    }
}

//-------------
// Apply a Givens rotation (c,s) to rows i and j of an n x n matrix A.
void apply_givens(double *A, int n, double c, double s, int i, int j) {
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

//-------------
// This function is not modified per the directions.
void tridiagonalize_full(double *A, int n, int k, double *d, double *e) {
    for (int j = 0; j < n - 2; j++) {
        for (int i = n - 1; i > j + 1; i--) {
            if (fabs(A[i * n + j]) > 1e-12) {
                double c, s;
                givens_rotation(A[(i - 1) * n + j], A[i * n + j], &c, &s);
                apply_givens(A, n, c, s, i - 1, i);
            }
        }
    }

    for (int i = 0; i < n-1; i++) {
        d[i] = A[i * n + i];
        if (i < n - 1)
            e[i] = A[i * n + i + 1];
    }
    d[n-1] = A[(n-1) * n + (n-1)];
    e[n-1] = 0.0;
}

//-------------
// Corrected QR step for symmetric tridiagonal eigenvalue computation.
// We use the Wilkinson shift and then perform an implicit QR sweep that
// explicitly “chases the bulge” introduced by the rotation.
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

    double c1, s1;
    double a1 = A[0];
    double b1 = A[m];
    double hypot = sqrt(a1 * a1 + b1 * b1);
    c1 = a1 / hypot;
    s1 = -b1 / hypot;
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
        printf("A gauche\n");
        for (int row = 0; row < m; row++) {
            for (int col = 0; col < m; col++) {
                printf("%f ", A[row * m + col]);
            }
            printf("\n");
        }
    }

    double R_ik = A[0];
    double R_jk = A[1];
    A[0] = c1 * R_ik - s1 * R_jk;
    A[1] = s1 * R_ik + c1 * R_jk;
    
    double rR_ik = A[m];
    double rR_jk = A[m + 1];
    A[m] = c1 * rR_ik - s1 * rR_jk;
    A[m + 1] = c1 * rR_jk;

    printf("A droite\n");
    for (int row = 0; row < m; row++) {
        for (int col = 0; col < m; col++) {
            printf("%f ", A[row * m + col]);
        }
        printf("\n");
    }


    for(int i = 0; i < m-2; i++){
        //droite
        double a = A[i * m + i+1];
        double b = A[i * m + i+2];
        double hypot = sqrt(a * a + b * b);
        double c,s;
        if(hypot < 1e-10){
            c = 1.0;
            s = 0.0;
        }
        else{
            c = a / hypot;
            s = -b / hypot;
        }

        for (int k = i; k < i+3; k++) {  
                double R_ik = A[k * m + i+1];
                double R_jk = A[k * m + i+2];
                A[k * m + i+1] = c * R_ik - s * R_jk;
                A[k * m + i+2] = s * R_ik + c * R_jk;
    }
    printf("A droite\n");
        for (int row = 0; row < m; row++) {
            for (int col = 0; col < m; col++) {
                printf("%f ", A[row * m + col]);
            }
            printf("\n");
        }

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

    if(fabs(A[(m-2) * m + (m-1)]) > eps * (fabs(A[(m-1) * m + (m-1)]) + fabs(A[(m-2) * m + (m-2)]))){
        return m;
    }
    free(A);
    return m - 1;
}

//-------------
// Driver routine: note that we now allocate e with n elements because
// tridiagonalize_full writes to e[n-1].
int qr_eigs_band(double *A_compact, int n, int m, double eps, int max_iter, double *eigvals) {
    double *d = (double *)malloc(n * sizeof(double));
    double *e = (double *)malloc(n * sizeof(double));  // allocate n instead of (n-1)

    tridiagonalize_full(A_compact, n, m, d, e);

    int iter = 0;
    int current_size = n;
    while (current_size > 1 && iter < max_iter) {
        int ret = step_qr_tridiag(d, e, current_size, eps);
        if (ret == -1)
            break;
        current_size = ret;
        iter++;
    }

    for (int i = 0; i < n; i++)
        eigvals[i] = d[i];

    free(d);
    free(e);
    return iter;
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

void oui() {
    int nombre = 1;
    double *temps_arr = (double *)calloc(nombre, sizeof(double));
    double *taille_arr = (double *)calloc(nombre, sizeof(double));
    printf("Number of iterations: \n");
    clock_t start, end;
    for (int i = 0; i < nombre; i++) {
        printf("caca");
        double lx = 10.0;
        double ly = 10.0;
        int nx = 4;
        int ny = 4;
        double *A = NULL;
        double *d = (double *)calloc(nx*ny, sizeof(double));
        A = create_matrix(nx, ny, lx, ly, 0);
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
        double eps = 1e-10;
        int max_iter = 100;
        printf("ooo");
        start = clock();
        printf("start = %f\n", start);
        printf("start = %f\n", (double)start);
        int iteration = qr_eigs_band(A, n, k, eps, max_iter,d);
        end = clock();
        temps_arr[i] = (double)(end - start) / CLOCKS_PER_SEC;
        printf("d = Eigenvalue : ");
        for (int i = 0; i < n; i++) {
            printf("%f\n", d[i]);
        }
        free(d);
        free(A);
        printf("fin \n");
        printf("temps : %f\n", temps_arr[i]);
    }
    FILE *file = fopen("results.txt", "w");
    if (file != NULL) {
        for (int i = 0; i < nombre; i++) {
            fprintf(file, "%f %f\n", taille_arr[i], temps_arr[i]);
        }
        fclose(file);
    } else {
        printf("Error opening file!\n");
    }
    free(temps_arr);
    free(taille_arr);
}

int main() {
    oui();
    return 0;
}