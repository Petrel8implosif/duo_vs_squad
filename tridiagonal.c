#include <stdio.h>
#include <math.h>
#include <stdlib.h>
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
    
    for (int i = 0; i < n; i++) {
        d[i] = A[i * n + i];
        if (i < n - 1) {
            e[i] = A[i * n + i + 1];
        }
    }
    e[4] = 0.0;

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

void print_matrix(double *A, int n) {
    double (*matrix)[n] = (double (*)[n])A; // Cast to 2D array for better readability
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%10.4f ", matrix[i][j]); // Print with better formatting
        }
        printf("\n");
    }
}

// Main function
int main() {
    int n = 5;  // Matrix size
    int k = 2;  // Bandwidth (assuming a tridiagonal matrix)
    double A[5 * 5] = {
        4.0, 1.0, 2.0, 0.0, 0.0,
        1.0, 3.0, 1.0, 2.0, 0.0,
        2.0, 1.0, 4.0, 1.0, 0.0,
        0.0, 2.0, 1.0, 5.0, 1.0,
        0.0, 0.0, 0.0, 1.0, 3.0
    };
    double d[5], e[5];

    tridiagonalize_full(A, n, k, d, e);

    printf("\nFinal Tridiagonal Matrix:\n");
    printf("Diagonal (d): ");
    for (int i = 0; i < n; i++) {
        printf("%lf ", d[i]);
    }
    printf("\nSubdiagonal (e): ");
    for (int i = 0; i < n ; i++) {
        printf("%lf ", e[i]);
    }
    printf("\n");

    double lx = 10.0;
    double ly = 10.0;
    int nx = 4;
    int ny = 4;
    double *matrixA = NULL;
    double *d_1 = (double *)calloc(nx*ny, sizeof(double));
    double *e_1 = (double *)calloc(nx*ny, sizeof(double));
    matrixA = create_matrix(nx, ny, lx, ly, 0);
    int taille = nx * ny;
    int bande = 0;
    for (int row = 0; row < n; row++) {
        for (int col = 0; col < n; col++) {
            if (A[row * n + col] != 0) {
                int band_width = abs(row - col);
                if (band_width > k) {
                    bande = band_width;
                }
            }
        }
    }
    tridiagonalize_full(matrixA, taille, bande, d_1, e_1);
    printf("Diagonal (d): ");
    for (int i = 0; i < taille; i++) {
        printf("%lf ", d_1[i]);
    }
    printf("\nSubdiagonal (e): ");
    for (int i = 0; i < taille; i++) {
        printf("%lf ", e_1[i]);
    }
    //print the final matrixA
    printf("\nFinal Matrix A:\n");
    print_matrix(matrixA, taille);
    free(d_1);
    free(e_1);
    free(matrixA);
    return 0;
}
