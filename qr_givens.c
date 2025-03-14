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
    int i;
    for (int j = 0; j < N - 1; j++) {  
        iter++;
        i = j+1;
        if (fabs(A[i][j]) > 1e-10) {
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
    printf("Number of iterations: %d\n", iter);
}


void fill_A(double *d, double *e, int n, double *A){
    for (int i = 0; i < n-1; i++) {
        A[i* 3 + 1] = d[i];
        A[i* 3 + 2] = e[i];
        A[i* 3 + 3] = e[i];
    }
    A[n*3 - 2] = d[n-1];
}

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

void rotation_qr(double *d, double *e, int m, double c, double s) {
    double temp = d[0];
    d[0] = c * c * temp + 2 * c * s * e[0] + s * s * d[1];
    d[1] = s * s * temp - 2 * c * s * e[0] + c * c * d[1];
    e[0] = c * e[0] - s * temp;
}


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
        {1.0, 7.0, 0.0, 0.0, 0.0},
        {7.0, 3.0, 2.0, 0.0, 0.0},
        {0.0, 2.0, 4.0, 8.0, 0.0},
        {0.0, 0.0, 8.0, 5.0, 11.0},
        {0.0, 0.0, 0.0, 11.0, 6.0}
    };

    tridiagonalize(A);
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