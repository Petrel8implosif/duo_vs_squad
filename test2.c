#include <stdio.h>
#include <math.h>

#define TOL 1e-12

// Compute Givens rotation coefficients for (a, b)
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

/*
  In our compact band storage we assume that:
    - The matrix A is stored in an array of size n*(m+1)
    - Each row i (0-based) holds [A(i,i), A(i,i+1), ..., A(i,i+m)]
    - The main diagonal is in column 0 and the first superdiagonal in column 1.
  The goal is to eliminate the extra superdiagonals so that only the tridiagonal
  part (columns 0 and 1) remains.
*/

// Apply Givens rotation to eliminate the element in row i and column j
void apply_givens(double *A, int n, int k, double c, double s, int i, int j) {
    // Apply Givens rotation to the rows in the compact band storage
    for (int col = 0; col <= k; col++) {
        if (j + col < n) {
            double temp = c * A[i * (k + 1) + col] - s * A[j * (k + 1) + col];
            A[j * (k + 1) + col] = s * A[i * (k + 1) + col] + c * A[j * (k + 1) + col];
            A[i * (k + 1) + col] = temp;
        }
    }

    // Apply Givens rotation to the columns in the compact band storage
    for (int row = 0; row <= k; row++) {
        if (i + row < n) {
            double temp = c * A[(i + row) * (k + 1) + (k - row)] - s * A[(j + row) * (k + 1) + (k - row)];
            A[(j + row) * (k + 1) + (k - row)] = s * A[(i + row) * (k + 1) + (k - row)] + c * A[(j + row) * (k + 1) + (k - row)];
            A[(i + row) * (k + 1) + (k - row)] = temp;
        }
    }   

}



/*
  Chase the bulge downwards.
  Starting from a given row index, the saved bulge interacts with the next row's first
  superdiagonal (stored at column 1). A Givens rotation is computed to eliminate the bulge
  while updating that element; the new bulge (if any) is then carried to the following row.
  This continues until the bulge is negligible or we run out of rows.
*/
void chase_bulge(double *A, int n, int m, int start_row) {

    while(start_row < n - 1 && fabs(A[(start_row + 1) * (m + 1) + 1]) > TOL){
        int bulge_loc = 1;
    }
}

/*
  Tridiagonalize the symmetric band matrix A.
  Input:
    - A: an n x (m+1) compact array storing a band symmetric matrix.
    - n: order of A.
    - m: the band width (number of superdiagonals); the storage has m+1 columns.
  Output:
    - d: array of length n that will contain the diagonal elements.
    - e: array of length n-1 that will contain the subdiagonal (and superdiagonal) elements.
  The procedure applies a series of Givens rotations on each row to eliminate the extra
  off-diagonals. When a rotation would introduce a fill-in (bulge) outside the band,
  that bulge is saved and then chased down via further Givens rotations.
*/
void tridiagonalize_band(double *A, int n, int k, double *d, double *e){
    for(int i = 2; i< n; i++){
        for(int j = 0; j< -(k * k + 1); j-=k){
            if(k == 0){
                chase_bulge(A, n, k, i);
            }
            else{
                apply_givens(A,n,k,i,j);
            }
        }
    }
    // After the sweep, extract the diagonal and first superdiagonal.
    for (int i = 0; i < n; i++) {
        d[i] = A[i * (k+1)]; // main diagonal
        if (i < n - 1)
            e[i] = A[i * (k + 1) + 1]; // first superdiagonal (which equals subdiagonal)
    }
}

// Example usage
int main() {
    int n = 5;  // Matrix order
    int m = 2;  // Original number of superdiagonals (bandwidth m => storage width = m+1)
    
    // Example symmetric band matrix in compact form:
    // Each row holds: [2nd superdiagonal, 1st superdiagonal, diagonal]
    double A[] = {
        0.0, 0.0, 4.0,  // Row 0
        0.0, 1.0, 3.0,  // Row 1
        2.0, 1.0, 4.0,  // Row 2
        2.0, 1.0, 5.0,  // Row 3
        0.0, 1.0, 3.0   // Row 4
    };
    
    double d[n], e[n - 1];
    
    // Apply the tridiagonalization sweep.
    tridiagonalize_band(A, n, m, d, e);
    
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
