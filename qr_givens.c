#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cblas.h>
#include <complex.h>
#define N 5
#define EPSILON 1e-12
#define SQUARE(x) ((x) * (x))
#define idxQ(i, j, m) ((i) + (j) * (m))
#define idxR(i, j) ((i) + (j) * n)

/*
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

// Function to apply Givens rotation to a matrix
void apply_givens(double A[N][N], double c, double s, int i, int j) {
    for (int k = 0; k < N; k++) {
        double temp = c * A[i][k] - s * A[j][k];
        A[j][k] = s * A[i][k] + c * A[j][k];
        A[i][k] = temp;
    }

    for (int k = 0; k < N; k++) {
        double temp = c * A[k][i] - s * A[k][j];
        A[k][j] = s * A[k][i] + c * A[k][j];
        A[k][i] = temp;
    }
}

void tridiagonalize(double A[N][N]) {
    int iter = 0;
    for (int j = 0; j < N - 1; j++) {  // Iterate through columns
        for (int i = N - 1; i > j; i--) {  // Work from the bottom row upward
            if (i == j+1) {
                printf(" i = j+1\n");
                printf("Matrix A before applying Givens rotation:\n");
                for (int x = 0; x < N; x++) {
                    for (int y = 0; y < N; y++) {
                        printf("%f ", A[x][y]);
                    }
                    printf("\n");
                }
            }
            if (fabs(A[i][j]) > 1e-10) {
                iter++;
                 // If significant, apply Givens rotation
                double c, s;
                givens_rotation(A[i - 1][j], A[i][j], &c, &s);

                // Apply Givens rotation to zero out A[i][j]
                apply_givens(A, c, s, i - 1, i);
                // Print the matrix A after applying Givens rotation
                printf("Matrix A after applying Givens rotation:\n");
                for (int x = 0; x < N; x++) {
                    for (int y = 0; y < N; y++) {
                        printf("%f ", A[x][y]);
                    }
                    printf("\n");
                }
            }
        }
    }
}
*/
/*
void fill_A(double *d, double *e, int n, double *A){
    for (int i = 0; i < n-1; i++) {
        A[i* 3 + 1] = d[i];
        A[i* 3 + 2] = e[i];
        A[i* 3 + 3] = e[i];
    }
    A[n*3 - 2] = d[n-1];
} */

double trouver_mu(double *d, double *e, int m) {
    double a = d[m-2];
    double b = e[m-1];
    double d_last = d[m-1];
    double trace = a + d_last;
    double determinant = SQUARE(a - d_last) + 4 * SQUARE(b);
    double racine = sqrt(determinant);

    double lambda1 = (trace + racine) / 2.0;
    double lambda2 = (trace - racine) / 2.0;

    double mu = (fabs(lambda1 - d_last) < fabs(lambda2 - d_last)) ? lambda1 : lambda2;
    printf("mu = %f\n", mu);
    return mu;
}
/*
void rotation_qr(double *d, double *e, int m, double c, double s) {
    double temp = d[0];
    d[0] = c * c * temp + 2 * c * s * e[0] + s * s * d[1];
    d[1] = s * s * temp - 2 * c * s * e[0] + c * c * d[1];
    e[0] = c * e[0] - s * temp;
} */


/*
int qr_step(double *d, double *e, double m, double eps, double* Q){
    double mu = trouver_mu(d, e, m);
    double* shifted_d = (double*)malloc(m * sizeof(double));
    for (int i = 0; i < m; i++) {
        shifted_d[i] = d[i] - mu;
    }

    double c, s;


    rotation_qr(shifted_d, e, m, c, s);

    for (int i = 0; i < m; i++) {
        shifted_d[i] = d[i] - mu;
    }

    return m-1;
} */


// Fonction pour calculer les coefficients c et s de la rotation de Givens
void givens_rotation(double a, double b, double *c, double *s) {
    double hypot = sqrt(a * a + b * b);
    *c = a / hypot;
    *s = -b / hypot;
}

// Applique la décomposition QR par Givens
void qr_givens(double A[N][N]) {
    // Initialisation de R avec A
    double R[N][N];
    double Q[N][N];
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            R[i][j] = A[i][j];
        }
    }

    // Initialisation de Q à l'identité
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            Q[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    // Application des rotations de Givens
    for (int i = 0; i < N - 1; i++) {
        for (int j = i + 1; j < N; j++) {
            double c, s;
            givens_rotation(R[i][i], R[j][i], &c, &s);

            // Appliquer la rotation à R (modifie lignes i et j)
            for (int k = 0; k < N; k++) {
                double R_ik = R[i][k];
                double R_jk = R[j][k];
                R[i][k] = c * R_ik - s * R_jk;
                R[j][k] = s * R_ik + c * R_jk;
            }

            // Appliquer la rotation à Q (modifie colonnes i et j)
            for (int k = 0; k < N; k++) {
                double Q_ki = Q[k][i];
                double Q_kj = Q[k][j];
                Q[k][i] = c * Q_ki - s * Q_kj;
                Q[k][j] = s * Q_ki + c * Q_kj;
            }
        }
    }
    // Calculate A = Q^T * A * Q
    double temp[N][N] = {0};

    // temp = Q^T * A
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                temp[i][j] += Q[k][i] * A[k][j];
            }
        }
    }

    // A = temp * Q
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A[i][j] = 0;
            for (int k = 0; k < N; k++) {
                A[i][j] += temp[i][k] * Q[k][j];
            }
        }
    }

    // Print the matrix A after applying Givens rotation
    printf("Matrix A after applying Givens rotation:\n");
    for (int x = 0; x < N; x++) {
        for (int y = 0; y < N; y++) {
            printf("%f ", A[x][y]);
        }
        printf("\n");
    }
}

int main(){
    /*
    double *d = (double *)malloc(N * sizeof(double));
    double *e = (double *)malloc((N - 1) * sizeof(double));

    d[0] = 5.0;
    d[1] = 6.0;
    d[2] = 7.0;
    d[3] = 8.0;
    d[4] = 9.0;

    e[0] = 1.0;
    e[1] = 2.0;
    e[2] = 3.0;
    e[3] = 4.0;
    double* A = (double *) calloc(N * 3 * sizeof(double), sizeof(double));
    fill_A(d, e, N, A);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%f ", A[i * 3 + j]);
        }
        printf("\n");
    }
    free(A);
    free(d);
    free(e);*/

    double A[N][N] = {
        {1.0, 7.0, 0.0, 9.0, 0.0},
        {7.0, 3.0, 2.0, 9.0, 0.0},
        {0.0, 2.0, 4.0, 8.0, 7.0},
        {9.0, 9.0, 8.0, 5.0, 11.0},
        {0.0, 0.0, 7.0, 11.0, 6.0}
    };

    FILE *file_origin = fopen("A_origin.txt", "w");
    if (file_origin != NULL) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                fprintf(file_origin, "%f ", A[i][j]);
            }
            fprintf(file_origin, "\n");
        }
        fclose(file_origin);
    } else {
        printf("Error opening file for writing.\n");
    }

    int iter = 20;

    for(int i = 0; i < iter; i++){
        qr_givens(A);
    }
    FILE *file = fopen("matrix_A.txt", "w");
    if (file != NULL) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                fprintf(file, "%f ", A[i][j]);
            }
            fprintf(file, "\n");
        }
        fclose(file);
    } else {
        printf("Error opening file for writing.\n");
    }
    return 0;
}