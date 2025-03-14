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

  When eliminating an element from column j (j >= 2) in row i, the corresponding
  Givens rotation acts on the pair (A[i, j-1], A[i, j]). The rotation is chosen so that
  the new combination in the j-th position becomes zero. However, the orthogonal transform
  produces a new element (the bulge) that “falls off” the band. We do not store that value in A;
  instead we save it in a temporary variable and then “chase” it downwards by applying
  successive Givens rotations.
*/

// Eliminates the element in row i, column j (j>=2) and returns the bulge.
double eliminate_element(double *A, int n, int m, int i, int j) {
    // Compute indices into A: row i has entries A[i*(m+1) + col], col = 0..m.
    int idx_prev = i * (m + 1) + (j - 1);  // element in column j-1 (remains in band)
    int idx_curr = i * (m + 1) + j;          // element to be eliminated

    double a_val = A[idx_prev];
    double b_val = A[idx_curr];
    double c, s;
    givens_rotation(a_val, b_val, &c, &s);
    
    // Update the element that stays in the band
    double new_val = c * a_val - s * b_val;
    A[idx_prev] = new_val;
    
    // The combination that would be created outside the band (the bulge)
    double bulge = s * a_val + c * b_val;
    
    // Eliminate the off–band element in A
    A[idx_curr] = 0.0;
    
    return bulge;
}

/*
  Chase the bulge downwards.
  Starting from a given row index, the saved bulge interacts with the next row's first
  superdiagonal (stored at column 1). A Givens rotation is computed to eliminate the bulge
  while updating that element; the new bulge (if any) is then carried to the following row.
  This continues until the bulge is negligible or we run out of rows.
*/
void chase_bulge(double *A, int n, int m, int start_row, double bulge) {
    int i = start_row;
    while (i < n - 1 && fabs(bulge) > TOL) {
        // For row i+1, the first superdiagonal is at column 1.
        int idx = (i + 1) * (m + 1) + 1;
        double a_val = A[idx];
        double c, s;
        givens_rotation(a_val, bulge, &c, &s);
        
        // Update the band element (which remains in A)
        double new_a = c * a_val - s * bulge;
        A[idx] = new_a;
        
        // Compute the new bulge (which is not stored in A)
        double new_bulge = s * a_val + c * bulge;
        bulge = new_bulge;
        
        i++;
    }
    // Any remaining bulge (if not chased out completely) is discarded.
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
void tridiagonalize_band(double *A, int n, int m, double *d, double *e) {
    // Loop over each row.
    for (int i = 0; i < n; i++) {
        // For each row, eliminate entries in columns 2..m (if any).
        // (Columns 0 and 1 will remain in the final tridiagonal form.)
        for (int j = m; j >= 2; j--) {
            int idx = i * (m + 1) + j;
            if (fabs(A[idx]) > TOL) {
                // Eliminate element at (i, j) to create a bulge.
                double bulge = eliminate_element(A, n, m, i, j);
                // Chase the bulge downwards starting from row i.
                chase_bulge(A, n, m, i, bulge);
            }
        }
    }
    // After the sweep, extract the diagonal and first superdiagonal.
    for (int i = 0; i < n; i++) {
        d[i] = A[i * (m + 1) + 0]; // main diagonal
        if (i < n - 1)
            e[i] = A[i * (m + 1) + 1]; // first superdiagonal (which equals subdiagonal)
    }
}

// Example usage
int main() {
    int n = 5;  // Matrix order
    int m = 2;  // Original number of superdiagonals (bandwidth m => storage width = m+1)
    
    // Example symmetric band matrix in compact form:
    // Each row holds: [diagonal, 1st superdiagonal, 2nd superdiagonal]
    double A[] = {
        4.0, 1.0, 2.0,  // Row 0
        3.0, 1.0, 2.0,  // Row 1
        4.0, 1.0, 0.0,  // Row 2
        5.0, 1.0, 0.0,  // Row 3
        3.0, 0.0, 0.0   // Row 4
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
