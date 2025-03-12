#include "devoir_1.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cblas.h>
#include <complex.h>

#pragma GCC target("fma,avx,avx2,ssse3,sse,sse2")

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))
#define square(x) ((x) * (x))
#define SQUARE(x) ((x) * (x))
#define idxQ(i, j, m) ((i) + (j) * (m))
#define idxR(i, j) ((i) + (j) * n)  
#define EPSILON 1e-12

//gcc -Wall -Wextra step_trig_qr.c -lm -o main

int QR(double *A, double *R, int *P, int m, int n) {
    for (int i = 0; i < n; i++) {
        P[i] = i;
    }
    int rang = 0;
    for (int j = 0; j < n; j++) {
        double norm = 0.0;
        for (int i = 0; i < m; i++) {
            norm += SQUARE(A[idxQ(i, j, m)]);
        }
        R[idxR(j, j)] = sqrt(norm);
        

        if (R[idxR(j, j)] < EPSILON) {
            rang += 1;
            for (int i = 0; i < m; i++) {
                A[idxQ(i, j, m)] = 0.0;
            }
            
            continue;
            
        }
        for (int i = 0; i < m; i++) {
            A[idxQ(i, j, m)] /= R[idxR(j, j)];
        }

        for (int i = j + 1; i < n; i++) {
            R[idxR(j, i)] = 0.0;
            for (int k = 0; k < m; k++) {
                R[idxR(j, i)] += A[idxQ(k, j, m)] * A[idxQ(k, i, m)];
            }
            for (int k = 0; k < m; k++) {
                A[idxQ(k, i, m)] -= R[idxR(j, i)] * A[idxQ(k, j, m)];
            }
        }
    }
    
    return n - rang;
}

void construire_matrice_colonne_major(double *d, double *e, double *A, int n) {
    // Initialisation de la matrice à zéro
    for (int i = 0; i < n * n; i++) {
        A[i] = 0.0;
    }

    // Remplissage de la diagonale principale
    for (int i = 0; i < n; i++) {
        A[i + i * n] = d[i]; // A(i, i) = d[i]
    }

    // Remplissage des diagonales inférieure et supérieure
    for (int i = 0; i < n - 1; i++) {
        A[i + (i + 1) * n] = e[i];  // Supérieure A(i, i+1)
        A[(i + 1) + i * n] = e[i];  // Inférieure A(i+1, i) (symétrique)
    }
}

double trouver_mu(double *d, double *e, int m) {
    double a = d[m-2];
    double b = e[m-1];
    double d_last = d[m-1];
    double trace = a + d_last;
    double determinant = square(a - d_last) + 4 * square(b);
    double racine = sqrt(determinant);

    double lambda1 = (trace + racine) / 2.0;
    double lambda2 = (trace - racine) / 2.0;

    double mu = (fabs(lambda1 - d_last) < fabs(lambda2 - d_last)) ? lambda1 : lambda2;
    printf("mu = %f\n", mu);
    return mu;
}

void transpose_conjugue(int n, double *Q, double *Q_star) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Q_star[i * n + j] = Q[j * n + i]; // Conjugué et transposition
        }
    }
}

void multiplication_matrice(int n, double *A, double *B, double *C) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            C[i + j * n] = 0.0;
            for (int k = 0; k < n; k++) {
                C[i + j * n] += A[i + k * n] * B[k + j * n];
            }
        }
    }
}

void ajouter_muI(int n, double *M, double mu) {
    for (int i = 0; i < n; i++) {
        M[i  + i * n] += mu;
    }
}

int step_tridiag(double *d, double *e, double m, double eps, double* Q){
    
    int mi = m;
    double mu = trouver_mu(d, e, mi);
    double* shifted_d = (double*)malloc(mi * sizeof(double));
    
    //shift 
    for (int i = 0; i < mi; i++) {
        shifted_d[i] = d[i] - mu;
    }
    construire_matrice_colonne_major(shifted_d, e, Q, mi);

    printf("Matrice Q après construction:\n");
    for (int i = 0; i < mi; i++) {
        for (int j = 0; j < mi; j++) {
            printf("%f ", Q[j * mi + i]);
        }
        printf("\n");
    }

    double* R = (double*)malloc(mi * mi * sizeof(double));
    int *P = (int*)malloc(mi * sizeof(int));
    double * A_k_1 = Q;
    double * Q_conj = (double*)malloc(mi * mi * sizeof(double));
    
    QR(Q, R, P, mi, mi);
    transpose_conjugue(mi, Q, Q_conj);
    double *AQ = (double*)malloc(mi * mi * sizeof(double));
    multiplication_matrice(mi, A_k_1, Q, AQ);

    multiplication_matrice(mi, Q_conj, AQ, Q);
    ajouter_muI(mi, Q, mu);
    free(Q_conj);
    free(AQ);
    free(R);
    return 0;
}


int qr_eig_band(double *d, double *e,  int n, int k, double eps, int max_iter, double* Q){
    for (int i = 0; i < max_iter; i++) {
        step_tridiag(d, e, n, eps, Q);
    }
    return 0;
}


#define N 5

int main() {
    printf("working\n");
    fflush(stdout);  // S'assurer que le message est affiché immédiatement

    double *Q = (double *)calloc(N * N * sizeof(double), sizeof(double));
    if (Q == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        return 1;
    }

    double d[N] = {4.0, 3.0, 4.0, 5.0, 3.0};
    double e[N - 1] = {1.0, 1.0, 1.0, 1.0};
    FILE *file = fopen("matrix_Q.txt", "w");
    if (file == NULL) {
        fprintf(stderr, "Error opening file for writing.\n");
        free(Q);
        return 1;
    }

    setbuf(file, NULL);  // Disable buffering to ensure immediate writes

    fprintf(file, "Avant step_tridiag\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            fprintf(file, "%f ", Q[j * N + i]);
        }
        fprintf(file, "\n");
    }

    printf("Avant step_tridiag\n");
    fflush(stdout);

    qr_eig_band(d, e, N, 1, 1e-12, 1, Q);

    fprintf(file, "Après step_tridiag\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            fprintf(file, "%f ", Q[j * N + i]);
        }
        fprintf(file, "\n");
    }

    fclose(file);
    free(Q);
    return 0;
}
