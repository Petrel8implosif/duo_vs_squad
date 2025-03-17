#include "devoir_1.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
void tridiagonalize_band(double *A, int n, int k, double *d, double *e){
    
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
    printf("Matrix A before applying Givens rotation:\n");
    for (int x = 0; x < m; x++) {
        for (int y = 0; y < m; y++) {
            printf("%f ", A[x * m + y]);
        }
        printf("\n");
    }
    while(fabs(A[(m-1) * m + (m-2)]) > eps * (fabs(A[(m-1) * m + (m-1)]) + fabs(A[(m-2) * m + (m-2)]))) {
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
            printf("Matrix R after applying Givens rotation:\n");
            for (int x = 0; x < m; x++) {
                for (int y = 0; y < m; y++) {
                printf("%f ", R[x * m + y]);
                }
                printf("\n");
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

int qr_eigs_band(double *A, int n, int k, double eps, int max_iter){
    
}