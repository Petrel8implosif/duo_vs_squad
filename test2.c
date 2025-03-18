#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define TOL 1e-12

// Retained function: adds a column of zeros at the beginning of each row.
// Input matrix A is in compact band form with n rows and k = m+1 columns.
double *add_column_of_zeros(double *A, int n, int k) {
    int ld = k + 1;  // new number of columns after adding the zero column
    double *B = malloc(n * ld * sizeof(double));
    if (B == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return NULL;
    }
    for (int i = 0; i < n; i++) {
        // Set the new (first) column to zero.
        B[i * ld] = 0.0;
        // Copy the original k entries into columns 1 .. k.
        for (int j = 0; j < k; j++) {
            B[i * ld + (j + 1)] = A[i * k + j];
        }
    }
    return B;
}


void givens_rotation(double a, double b, double *c, double *s) {
    double r = sqrt(a * a + b * b);
    if (r < TOL) {
        *c = 1.0;
        *s = 0.0;
    } else {
        *c = a / r;
        *s = -b / r;
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


void tridiagonalize_band(double *A, int n, int k, double *d, double *e) {
    int ld = k + 1;  
    double *B = add_column_of_zeros(A, n, k);
    if (B == NULL) {
        return;
    }
    
    for (int i = 1; i < n - 1; i++) {
        {
            double c, s;
            int col = ld - 2;  
            givens_rotation(B[i * ld + col], B[(i + 1) * ld + col], &c, &s);
            apply_givens(B, n, ld, c, s, i, i + 1);
        }
        {
            double c, s;
            int col = ld - 3;  
            givens_rotation(B[i * ld + col], B[(i + 1) * ld + col], &c, &s);
            apply_givens(B, n, ld, c, s, i, i + 1);
        }
    }
    
    for (int i = 0; i < n; i++) {
        d[i] = B[i * ld + (ld - 1)];
        if (i < n - 1) {
            e[i] = B[i * ld + (ld - 2)];
        }
    }
    
    free(B);
}


int main() {
    int n = 5;   
    int m = 2;   
    int k = m + 1;

    double A[] = {
        0.0, 0.0, 4.0,  // Row 0
        0.0, 1.0, 3.0,  // Row 1
        2.0, 1.0, 4.0,  // Row 2
        2.0, 1.0, 5.0,  // Row 3
        0.0, 1.0, 3.0   // Row 4
    };
    
    double d[n], e[n - 1];
    
    // Apply the tridiagonalization (QR reduction) sweep.
    tridiagonalize_band(A, n, k, d, e);
    
    // Output the tridiagonal matrix components.
    printf("Diagonal (d): ");
    for (int i = 0; i < n; i++) {
        printf("%lf ", d[i]);
    }
    printf("\nSubdiagonal (e): ");
    for (int i = 0; i < n - 1; i++) {
        printf("%lf ", e[i]);
    }
    printf("\n");
    
    return 0;
}
