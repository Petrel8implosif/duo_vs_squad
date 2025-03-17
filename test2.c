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

// Compute Givens rotation coefficients (c,s) so that applying the rotation
// will annihilate b (i.e. produce a zero in the second entry).
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

// Apply the Givens rotation to rows i and j of matrix A.
// A is stored with ld columns.
void apply_givens(double *A, int n, int ld, double c, double s, int i, int j) {
    for (int col = 0; col < ld; col++) {
        double temp = c * A[i * ld + col] - s * A[j * ld + col];
        A[j * ld + col] = s * A[i * ld + col] + c * A[j * ld + col];
        A[i * ld + col] = temp;
    }
}

/*
  Revised tridiagonalization algorithm for a symmetric band matrix.
  
  Input:
    - A: original n x k compact matrix (with k = m+1, m = number of superdiagonals).
         For example, each row of A holds: [A(i,i+2), A(i,i+1), A(i,i)]
    - n: order of the matrix.
    - k: number of columns in A.
  Output:
    - d: array (length n) that will hold the diagonal elements.
    - e: array (length n-1) that will hold the subdiagonal elements.
    
  Procedure:
    1. Expand A by adding a column of zeros to the left (new ld = k+1).
    2. For rows i = 1 to n–2, perform two Givens rotations (in the order indicated by the
       bold positions in the slides) to chase the bulge:
         • First, eliminate the entry in column (ld–2) (this column will become the subdiagonal).
         • Second, eliminate the entry in column (ld–3).
    3. Extract the tridiagonal matrix:
         - Main diagonal from column (ld–1)
         - Subdiagonal from column (ld–2)
*/
void tridiagonalize_band(double *A, int n, int k, double *d, double *e) {
    int ld = k + 1;  // extended leading dimension (for m = 2, ld = 4)
    double *B = add_column_of_zeros(A, n, k);
    if (B == NULL) {
        return;
    }
    
    // Bulge chasing loop.
    // For each row i (from 1 to n-2), perform two rotations.
    for (int i = 1; i < n - 1; i++) {
        // --- Rotation on column (ld-2) ---
        {
            double c, s;
            int col = ld - 2;  // for m=2, this is column 2.
            // The rotation is performed between row i and row i+1 at the chosen column.
            givens_rotation(B[i * ld + col], B[(i + 1) * ld + col], &c, &s);
            apply_givens(B, n, ld, c, s, i, i + 1);
        }
        // --- Rotation on column (ld-3) ---
        {
            double c, s;
            int col = ld - 3;  // for m=2, this is column 1.
            givens_rotation(B[i * ld + col], B[(i + 1) * ld + col], &c, &s);
            apply_givens(B, n, ld, c, s, i, i + 1);
        }
    }
    
    // Extract the tridiagonal matrix:
    // The main diagonal is stored in column (ld-1)
    // The subdiagonal is stored in column (ld-2)
    for (int i = 0; i < n; i++) {
        d[i] = B[i * ld + (ld - 1)];
        if (i < n - 1) {
            e[i] = B[i * ld + (ld - 2)];
        }
    }
    
    free(B);
}

///////////////////////
// Example usage:
///////////////////////
int main() {
    int n = 5;   // Matrix order.
    int m = 2;   // Number of original superdiagonals.
                 // (Thus the compact storage has k = m + 1 = 3 columns.)
    int k = m + 1;
    
    // Example symmetric band matrix in compact form:
    // Each row holds: [A(i,i+2), A(i,i+1), A(i,i)]
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
