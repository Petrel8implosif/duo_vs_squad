#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cblas.h>
#include <complex.h>
#define N 5
#define EPSILON 1e-9
#define SQUARE(x) ((x) * (x))
#define idxQ(i, j, m) ((i) + (j) * (m))
#define idxR(i, j) ((i) + (j) * n)
int iter = 0;



void givens_rotation(double a, double b, double *c, double *s) {
    double hypot = sqrt(a * a + b * b);
    *c = a / hypot;
    *s = -b / hypot;
}
double qr_givens(double *d, double *e, int n, double eps) {
    while(fabs(e[n-2]) > eps * (fabs(d[n-1]) + fabs(d[n-2]))) {
        iter++;
        double a = d[n-2];
        double b = e[n-2];
        double d_last = d[n-1];
        double trace = a + d_last;
        double determinant = SQUARE(a - d_last) + 4 * SQUARE(b);
        double mu;
        if (determinant < 0) {
            mu = d[n-1];
        }
        else{
            double racine = sqrt(determinant);
            double lambda1 = (trace + racine) / 2.0;
            double lambda2 = (trace - racine) / 2.0;
            mu = (fabs(lambda1 - d_last) < fabs(lambda2 - d_last)) ? lambda1 : lambda2;
        }
        
        for (int i = 0; i < n; i++) {
           d[i] -= mu;
        }

        double *d_R = (double *)malloc(n * sizeof(double));
        double *e_R_h = (double *)malloc((n-1) * sizeof(double));
        double *e_R_b = (double *)malloc((n-1) * sizeof(double));
        memcpy(d_R, d, n * sizeof(double)); 
        memcpy(e_R_h, e, (n-1) * sizeof(double));
        memcpy(e_R_b, e, (n-1) * sizeof(double));
        if(iter < 10) {
            printf("d_R vector after shift:\n");
            for (int i = 0; i < n; i++) {
                printf("%f ", d_R[i]);
            }
            printf("\n");

            printf("e_R vector after shift:\n");
            for (int i = 0; i < n - 1; i++) {
                printf("%f ", e_R_b[i]);
            }
            printf("\n");
        }
        double *Q = (double *)calloc(n * n, sizeof(double));

        for (int i = 0; i < n; i++) {
            Q[i * n + i] = 1.0;
        }
        
        double R_ik, R_jk;
        
        double Q_ki, Q_kj;

        for (int i = 0; i < n - 2; i++) {
            double c, s;
            givens_rotation(d_R[i], e_R_b[i], &c, &s);

            R_ik = d_R[i];
            R_jk = e_R_b[i];
            d_R[i] = c * R_ik - s * R_jk;
            e_R_b[i] = s * R_ik + c * R_jk;

            R_ik = e_R_h[i];
            R_jk = d_R[i + 1];
            e_R_h[i] = c * R_ik - s * R_jk;
            d_R[i + 1]  = s * R_ik + c * R_jk;

            double a = e_R_h[i+1];
            e_R_h[i + 1] = c * a;

            //iteration Q

            for (int k = 0; k < n && k < i+2; k++) {
                Q_ki = Q[k * n + i];
                Q_kj = Q[k * n + i+1];
                Q[k * n + i] = c * Q_ki - s * Q_kj;
                Q[k * n + i+1] = s * Q_ki + c * Q_kj;
            }
        }
        // iteration pour i = n-1
        double c, s;
        givens_rotation(d_R[n-2], e_R_b[n-2], &c, &s);

        R_ik = d_R[n-2];
        R_jk = e_R_b[n-2];
        d_R[n-2] = c * R_ik - s * R_jk;
        e_R_b[n-2] = s * R_ik + c * R_jk;

        R_ik = e_R_h[n-2];
        R_jk = d_R[n-1];
        e_R_h[n-2] = c * R_ik - s * R_jk;
        d_R[n-1] = s * R_ik + c * R_jk;

        //iteration Q pour i = n-1
        for (int k = 0; k < n; k++) {
            Q_ki = Q[k * n + n-2];
            Q_kj = Q[k * n + n-1];
            Q[k * n + n-2] = c * Q_ki - s * Q_kj;
            Q[k * n + n-1] = s * Q_ki + c * Q_kj;
        }


        double *temp = (double *)calloc(n * n, sizeof(double));

        //première ligne de temp

        for(int j = 0; j < n; j++) {
            temp[j] = d[0] * Q[0 + j] + e[0] * Q[n + j];
        }

        //mid lignes de temp

        for (int i = 1; i < n-1; i++) {
            for (int j = 0; j < n; j++) {
                temp[i * n + j] = e[i-1] * Q[n * (i-1) + j] + d[i] * Q[n * i + j] + e[i] * Q[n * (i+1) + j];
            }
        }

        // dernière ligne de temp

        temp[n*n - 3] = e[n-2] * Q[n * (n-2) + n-3];
        temp[n*n - 2] = d[n-1] * Q[n * (n-1) + n-2] + e[n-2] * Q[n * (n-2) + n-2];
        temp[n*n - 1] = d[n-1] * Q[n * (n-1) + n-1] + e[n-2] * Q[n * (n-2) + n-1];

        //première ligne de A
        
        d[0] = 0;
        for (int k = 0; k < 3; k++) {
            d[0] += temp[k] * Q[k * n];
        }
        e[0] = 0;
        for (int k = 0; k < 3; k++) {
            e[0] += temp[k] * Q[k * n + 1];
        }

        //mid lignes de A
        for (int i = 1; i < n-1; i++) {
            for (int k = 0; k < n; k++) {
                d[i] += temp[i * n + k] * Q[k * n + i];
                e[i] += temp[i * n + k] * Q[k * n + i + 1];
            }
        }

        //dernières lignes de A

        d[n-1] = 0;
        for (int k = 0; k < n; k++) {
            d[n-1] += temp[(n-1)*n + k] * Q[k * n + n-1];
        }

        //deshiftage de A

        for (int i = 0; i < n; i++) {
            d[i] += mu;
        }
        free(Q);
        free(temp);
        free(d_R);
        free(e_R_b);
        free(e_R_h);

    }
    return n-1;
    
}

// version full avec A mais il faudrait la construire à partir de d et e.

/*double qr_givens(double *A, int n, double eps) {
    printf("%f\n", fabs(A[(n-2) * n + (n-1)]));
    printf("%f\n", eps * (fabs(A[(n-1) * n + (n-1)]) + fabs(A[(n-2) * n + (n-2)])));
    while(fabs(A[(n-1) * n + (n-2)]) > eps * (fabs(A[(n-1) * n + (n-1)]) + fabs(A[(n-2) * n + (n-2)]))) {
        iter++;
        double mu = trouver_mu(A, n);
        for (int i = 0; i < n; i++) {
            A[i * n + i] -= mu;
        }

        double *R = (double *)malloc(n * n * sizeof(double));
        double *Q = (double *)malloc(n * n * sizeof(double));
        memcpy(R, A, n * n * sizeof(double));

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                Q[i * n + j] = (i == j) ? 1.0 : 0.0;
            }
        }

        for (int i = 0; i < n - 1; i++) {
            for (int j = i + 1; j < i+2; j++) {
                double c, s;
                givens_rotation(R[i * n + i], R[j * n + i], &c, &s);

                for (int k = i; k < i+3 && k < n; k++) {
                    if(k == n || k == i+2) {
                        double R_jk = R[j * n + k];
                        R[j * n + k] = c * R_jk;
                    }
                    else {
                        double R_ik = R[i * n + k];
                        double R_jk = R[j * n + k];
                        R[i * n + k] = c * R_ik - s * R_jk;
                        R[j * n + k] = s * R_ik + c * R_jk;
                    }
                    
                }
                printf("Matrix R after Givens rotation at step (%d, %d):\n", i, j);
                for (int x = 0; x < n; x++) {
                    for (int y = 0; y < n; y++) {
                        printf("%f ", R[x * n + y]);
                    }
                    printf("\n");
                }
                

                for (int k = 0; k < n && k < i+2; k++) {
                    double Q_ki = Q[k * n + i];
                    double Q_kj = Q[k * n + j];
                    Q[k * n + i] = c * Q_ki - s * Q_kj;
                    Q[k * n + j] = s * Q_ki + c * Q_kj;
                }
                printf("Matrix Q after Givens rotation at step (%d, %d):\n", i, j);
                    for (int x = 0; x < n; x++) {
                        for (int y = 0; y < n; y++) {
                            printf("%f ", Q[x * n + y]);
                        }
                        printf("\n");
                    }
            }
        }

        double *temp = (double *)calloc(n * n, sizeof(double));

        //première ligne de temp

        for(int j = 0; j < n; j++) {
            for (int k = 0; k < 2; k++) {
                temp[j] += Q[k * n] * A[k * n + j];
            }
        }

        //mid lignes de temp

        for (int i = 1; i < n-1; i++) {
            for (int j = 0; j < n && j<i+3; j++) {
                for (int k = i-1; k < i+2; k++) {
                    temp[i * n + j] += Q[k * n + i] * A[k * n + j];
                }
            }
        }

        // dernière ligne de temp

        for(int j = 0; j < n; j++) {
            for (int k = n-2; k < n; k++) {
                temp[(n-1) *n + j] += Q[k * n + (n-1)] * A[k * n + j];
            }
        }

        //première ligne de A
        
        A[0] = 0;
        for (int k = 0; k < 3; k++) {
            A[0] += temp[k] * Q[k * n];
        }
        A[1] = 0;
        for (int k = 0; k < 3; k++) {
            A[1] += temp[k] * Q[k * n + 1];
        }

        //mid lignes de A
        for (int i = 1; i < n-1; i++) {
            for (int j = i-1; j<i+2 ; j++) {
                A[i * n + j] = 0;
                for (int k = 0; k < n; k++) {
                    A[i * n + j] += temp[i * n + k] * Q[k * n + j];
                }
            }
        }

        //dernières lignes de A

        A[n*n - 2] = 0;
        for (int k = 0; k < n; k++) {
            A[n*n - 2] += temp[(n-1)*n + k] * Q[k * n + n-2];
        }
        A[n*n - 1] = 0;
        for (int k = 0; k < n; k++) {
            A[n*n - 1] += temp[(n-1)*n + k] * Q[k * n + n-1];
        }

        //deshiftage de A

        for (int i = 0; i < n; i++) {
            A[i * n + i] += mu;
        }

        free(R);
        free(Q);
        free(temp);
    }

    printf("Matrix A after applying Givens rotation:\n");
    for (int x = 0; x < n; x++) {
        for (int y = 0; y < n; y++) {
            printf("%f ", A[x * n + y]);
        }
        printf("\n");
    }
    return n - 1;
}*/

// prob avec iter 
int spectral_decomposition(double *d, double* e, int n, double eps, int maxiter) {
    int m = n;
    for (int i = 0; i < n; i++) {
        m = qr_givens(d, e, m, eps);
    }
    if(iter > maxiter) {
        return -1;
    }
    else{
        return iter;
    }
}
/*int spectral_decomposition(double *A, int n, double eps, double* resultat) {
    int iteration = 0;
    int m = n;
    for (int i = 0; i < n; i++) {
        iteration++;
        printf("IIIIIIII = %d\n", i);
        m = qr_givens(A, m, eps);
        resultat[i] = A[(m) * (m+1) + (m)];
        double *new_A = (double *)malloc((m) * (m) * sizeof(double));
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
            new_A[i * (m) + j] = A[i * (m+1) + j];
            }
        }
        printf("Matrix new_A after reduction:\n");
        for (int x = 0; x < m; x++) {
            for (int y = 0; y < m; y++) {
            printf("%f ", new_A[x * m + y]);
            }
            printf("\n");
        }
        free(A);
        A = new_A;
        printf("Matrix A after reduction:\n");
        for (int x = 0; x < m; x++) {
            for (int y = 0; y < m; y++) {
            printf("%f ", A[x * m + y]);
            }
            printf("\n");
        }
        printf("%f\n", fabs(A[(m-2) * m + (m-1)]));
        printf("m = %d\n", m);
    }
    return m;
}*/

int main() {
    double *A = (double *)malloc(N * N * sizeof(double));
    double initial_A[N * N] = {
        1.0, 7.0, 0.0, 0.0, 0.0,
        7.0, 3.0, 2.0, 0.0, 0.0,
        0.0, 2.0, 4.0, 8.0, 0.0,
        0.0, 0.0, 8.0, 5.0, 11.0,
        0.0, 0.0, 0.0, 11.0, 6.0
    };
    memcpy(A, initial_A, N * N * sizeof(double));
    FILE *file_origin = fopen("A_origin.txt", "w");
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
    printf("d vector:\n");
    for (int i = 0; i < N; i++) {
        printf("%f ", d[i]);
    }
    printf("\n");

    printf("e vector:\n");
    for (int i = 0; i < N - 1; i++) {
        printf("%f ", e[i]);
    }
    printf("\n");

    printf("GOOOOO");
    spectral_decomposition(d, e, N, EPSILON, 1000);

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