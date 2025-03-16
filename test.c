#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TOL 1e-12

/*
   In our compact row–major storage for symmetric band matrices,
   each row i (0 <= i < n) stores the entries:
      A(i,i), A(i,i+1), ..., A(i, i+bw-1)
   provided the column index does not exceed n-1.
   We allocate ld = bw + 1 columns so that any fill–in (bulge) that gets
   created by a rotation (near the “edge” of the band) can be temporarily stored.
*/

/* Get the (i,j) entry of the symmetric matrix stored compactly.
   If i <= j then the entry A(i,j) is stored in row i at index (j-i).
   If i > j, by symmetry: A(i,j)=A(j,i) stored in row j at index (i-j). */
double get_elem(double **A, int n, int ld, int i, int j) {
    if (i <= j) {
        int k = j - i;
        return (k < ld) ? A[i][k] : 0.0;
    } else {
        int k = i - j;
        return (k < ld) ? A[j][k] : 0.0;
    }
}

/* Set the (i,j) entry of the symmetric matrix stored compactly.
   When i <= j the element A(i,j) is stored in row i at index (j-i).
   When i > j, by symmetry A(i,j)=A(j,i) is stored in row j at index (i-j). */
void set_elem(double **A, int n, int ld, int i, int j, double val) {
    if (i <= j) {
        int k = j - i;
        if(k < ld) A[i][k] = val;
    } else {
        int k = i - j;
        if(k < ld) A[j][k] = val;
    }
}

/* Compute Givens rotation parameters (c,s) so that:
         [ c   s ] [ a ] = [ r ]
         [ -s  c ] [ b ]   [ 0 ]
   with r = sqrt(a^2+b^2) (unless both a and b are essentially 0).  */
void givens(double a, double b, double *c, double *s) {
    double r = hypot(a, b);
    if (r < TOL) {
        *c = 1.0;
        *s = 0.0;
    } else {
        *c = a / r;
        *s = b / r;
    }
}

/*
   Apply a symmetric Givens rotation – i.e. a congruence transformation
   A -> Q^T A Q – to rows and columns p and q of the full matrix.
   (We update only the stored (upper triangular) part.)
   Here we assume p < q.
*/
void apply_rotation_sym(double **A, int n, int ld, int p, int q, double c, double s) {
    /* Update the upper triangular part for rows p and q.
       For indices j from p to n-1 (for row p, these entries are stored),
       and similarly for row q. */
    for (int j = p; j < n; j++) {
        double a = get_elem(A, n, ld, p, j);
        double b = get_elem(A, n, ld, q, j);
        double new_a = c * a + s * b;
        double new_b = -s * a + c * b;
        set_elem(A, n, ld, p, j, new_a);
        set_elem(A, n, ld, q, j, new_b);
    }
    /* Also update the entries in the upper triangle that occur in rows i < p,
       i.e. the stored positions for A(i,p) and A(i,q). */
    for (int i = 0; i < p; i++) {
        double a = get_elem(A, n, ld, i, p);
        double b = get_elem(A, n, ld, i, q);
        double new_a = c * a + s * b;
        double new_b = -s * a + c * b;
        set_elem(A, n, ld, i, p, new_a);
        set_elem(A, n, ld, i, q, new_b);
    }
}

/*
   Compute the rotation that eliminates the element in row q and column "col"
   (which is (virtually) stored as A(q,col) via symmetry when q > col) using
   row p (which must hold the pivot entry). Then apply the symmetric rotation.
*/
void apply_givens_on_column(double **A, int n, int ld, int p, int q, int col) {
    double a = get_elem(A, n, ld, p, col);
    double b = get_elem(A, n, ld, q, col);
    double c, s;
    givens(a, b, &c, &s);
    apply_rotation_sym(A, n, ld, p, q, c, s);
}

/*
   Tridiagonalize a symmetric band matrix stored in compact form.
   The original matrix has nonzero entries for positions A(i,j) with j-i < orig_bw.
   We wish to zero out all “extra” off–diagonals (i.e. below the first subdiagonal).
   Bulge–chasing is performed following the procedure in the pdf.
*/
void tridiagonalize(double **A, int n, int ld, int orig_bw) {
    // Loop over columns 0,...,n-3
    for (int i = 0; i < n - 2; i++) {
        /* For column i, the original band extends up to row
           last = min(n-1, i+orig_bw-1). We wish to eliminate entries A(j,i) for j >= i+2.
           Note: For j > i, the element A(j, i) is stored as A(i, j-i). */
        int last = (i + orig_bw - 1 < n) ? i + orig_bw - 1 : n - 1;
        for (int j = i + 2; j <= last; j++) {
            if (fabs(get_elem(A, n, ld, j, i)) > TOL) {
                /* The pdf procedure calls for using the first subdiagonal entry as pivot.
                   Here we eliminate A(j,i) using a Givens rotation on rows (i+1) and j.
                   (Remember: A(j,i) = A(i,j) is stored in row i at index j-i.) */
                apply_givens_on_column(A, n, ld, i + 1, j, i);
                /* Bulge chasing: the rotation may have introduced a fill–in (bulge)
                   further down in column i. Chase it down until it is eliminated. */
                int m = j;
                while (m < n - 1 &&
                       fabs(get_elem(A, n, ld, m + 1, i)) > TOL) {
                    apply_givens_on_column(A, n, ld, m, m + 1, i);
                    m++;
                }
            }
        }
    }
}

/* Print the full symmetric matrix in standard (dense) format.
   This uses get_elem() to “virtually” expand the compact storage. */
void print_full_matrix(double **A, int n, int ld) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double val = get_elem(A, n, ld, i, j);
            printf("%8.3f ", val);
        }
        printf("\n");
    }
}

int main(void) {
    int n = 6;          // matrix dimension
    int orig_bw = 4;    // original band: nonzero entries appear for j-i < orig_bw
    int ld = orig_bw + 1; // allocate one extra column for possible bulge fill–in

    // Allocate the compact matrix: n rows each with ld columns.
    double **A = (double **)malloc(n * sizeof(double *));
    if (!A) { perror("malloc failed"); exit(EXIT_FAILURE); }
    for (int i = 0; i < n; i++) {
        A[i] = (double *)calloc(ld, sizeof(double));
        if (!A[i]) { perror("calloc failed"); exit(EXIT_FAILURE); }
    }

    // Initialize A as a symmetric band matrix stored in compact form.
    // For each row i, fill entries for j = i,..., min(n-1, i+orig_bw-1)
    // (Here we use A(i,j) = i + j + 1 as an example.)
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < orig_bw; k++) {
            int j = i + k;
            if (j < n) {
                A[i][k] = i + j + 1;
            }
        }
    }

    printf("Original matrix:\n");
    print_full_matrix(A, n, ld);
    printf("\n");

    // Tridiagonalize the matrix using Givens rotations with bulge chasing.
    tridiagonalize(A, n, ld, orig_bw);

    printf("Tridiagonalized matrix:\n");
    print_full_matrix(A, n, ld);

    // Free allocated memory.
    for (int i = 0; i < n; i++) {
        free(A[i]);
    }
    free(A);

    return 0;
}
