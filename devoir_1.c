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

int step_qr_tridiag(double *d, double *e, int m, double eps){
    double *A = (double *)calloc(m * m, sizeof(double));
    for (int i = 0; i < m; i++) {
        A[i * m + i] = d[i];
        if (i < m - 1) {
            A[i * m + (i + 1)] = e[i];
            A[(i + 1) * m + i] = e[i];
        }
    }
    while(fabs(A[(m-1) * m + (m-2)]) > eps * (fabs(A[(m-1) * m + (m-1)]) + fabs(A[(m-2) * m + (m-2)]))) {
        iter++;
        double a = A[(m-2) * m + (m-2)];
        double b = A[(m-2) * m + (m-1)];
        double d_last = A[(m-1) * m + (m-1)];
        double trace = a + d_last;
        double determinant = SQUARE(a - d_last) + 4 * SQUARE(b);
        double mu;

        double racine = sqrt(determinant);

        double lambda1 = (trace + racine) / 2.0;
        double lambda2 = (trace - racine) / 2.0;
        if (determinant < 0) {
            mu = A[(m-1)*m + (m-1)];
        }
        else{
            double racine = sqrt(determinant);
            double lambda1 = (trace + racine) / 2.0;
            double lambda2 = (trace - racine) / 2.0;
            mu = (fabs(lambda1 - d_last) < fabs(lambda2 - d_last)) ? lambda1 : lambda2;
        }
        

        for (int i = 0; i < m; i++) {
            A[i * m + i] -= mu;
        }

        double *R = (double *)malloc(m * m * sizeof(double));
        double *Q = (double *)calloc(m * m, sizeof(double));
        memcpy(R, A, m * m * sizeof(double));

        for (int i = 0; i < m; i++) {
            Q[i * m + i] = 1.0;
        }

        for (int i = 0; i < m - 1; i++) {
            double c, s;
            givens_rotation(R[i * m + i], R[(i+1) * m + i], &c, &s);
            double a = R[i * m + i];
            double b = R[(i+1) * m + i];
            double hypot = sqrt(a * a + b * b);
            c = a / hypot;
            s = -b / hypot;


            for (int k = i; k < i+3 && k < m; k++) {
                if(k == m || k == i+2) {
                double R_jk = R[(i+1) * m + k];
                R[(i+1) * m + k] = c * R_jk;
                }
                else {
                double R_ik = R[i * m + k];
                double R_jk = R[(i+1) * m + k];
                R[i * m + k] = c * R_ik - s * R_jk;
                R[(i+1) * m + k] = s * R_ik + c * R_jk;
                }
                
            }
            

            for (int k = 0; k < m && k < i+2; k++) {
                double Q_ki = Q[k * m + i];
                double Q_kj = Q[k * m + (i+1)];
                Q[k * m + i] = c * Q_ki - s * Q_kj;
                Q[k * m + (i+1)] = s * Q_ki + c * Q_kj;
            }
        }

        double *temp = (double *)calloc(m * m, sizeof(double));

        //première ligne de temp

        for(int j = 0; j < m; j++) {
            for (int k = 0; k < 2; k++) {
                temp[j] += Q[k * m] * A[k * m + j];
            }
        }

        //mid lignes de temp

        for (int i = 1; i < m-1; i++) {
            for (int j = 0; j < m && j<i+3; j++) {
                for (int k = i-1; k < i+2; k++) {
                    temp[i * m + j] += Q[k * m + i] * A[k * m + j];
                }
            }
        }

        // dernière ligne de temp

        for(int j = 0; j < m; j++) {
            for (int k = m-2; k < m; k++) {
                temp[(m-1) *m + j] += Q[k * m + (m-1)] * A[k * m + j];
            }
        }

        //première ligne de A
        
        A[0] = 0;
        for (int k = 0; k < 3; k++) {
            A[0] += temp[k] * Q[k * m];
        }
        A[1] = 0;
        for (int k = 0; k < 3; k++) {
            A[1] += temp[k] * Q[k * m + 1];
        }

        //mid lignes de A
        for (int i = 1; i < m-1; i++) {
            for (int j = i-1; j<i+2 ; j++) {
                A[i * m + j] = 0;
                for (int k = 0; k < m; k++) {
                    A[i * m + j] += temp[i * m + k] * Q[k * m + j];
                }
            }
        }

        //dernières lignes de A

        A[m*m - 2] = 0;
        for (int k = 0; k < m; k++) {
            A[m*m - 2] += temp[(m-1)*m + k] * Q[k * m + m-2];
        }
        A[m*m - 1] = 0;
        for (int k = 0; k < m; k++) {
            A[m*m - 1] += temp[(m-1)*m + k] * Q[k * m + m-1];
        }

        //deshiftage de A

        for (int i = 0; i < m; i++) {
            A[i * m + i] += mu;
        }

        free(R);
        free(Q);
        free(temp);
    }
    for (int i = 0; i < m; i++) {
        d[i] = A[i * m + i];
    }
    for (int i = 0; i < m - 1; i++) {
        e[i] = A[i * m + (i + 1)];
    }
    free(A);
    return m - 1;
}

int qr_eigs_band(double *A, int n, int k, double eps, int max_iter, double *d){
    double *dr = (double *)calloc(n, sizeof(double));
    double *e = (double *)calloc(n-1, sizeof(double));
    tridiagonalize_full(A, n, k, dr, e);
    int m = n;
    for (int i = 0; i < n; i++) {
        m = step_qr_tridiag(dr, e, m, eps);
        d[n-1-i] = dr[m];
    }
    if(iter > max_iter) {
        return -1;
    }   
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
    int nombre = 6;
    double *temps_arr = (double *)calloc(nombre, sizeof(double));
    double *taille_arr = (double *)calloc(nombre, sizeof(double));
    printf("Number of iterations: \n");
    clock_t start, end;
    for (int i = 0; i < nombre; i++) {
        printf("caca");
        double lx = 10.0;
        double ly = 10.0;
        int nx = i + 1;
        int ny = i + 1;
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
    /*double lx = 10.0;
    double ly = 10.0;
    int nx = 5;
    int ny = 5;
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
    int k = 2;
    double eps = 1e-10;
    int max_iter = 100;
    int iteration = qr_eigs_band(A, n, k, eps, max_iter,d);
    printf("Number of iterations: %d\n", iteration);
    printf("Eigenvalues:\n");
    for (int i = 0; i < n; i++) {
        printf("%f\n", d[i]);
    }
    free(d);
    return 0;*/
}