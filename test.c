#include <stdio.h>
#include <math.h>

#define TOL 1e-12

void givens_rotation(double a, double b, double *c, double *s) {
    double r = hypot(a, b);
    if (r < TOL) {
        *c = 1.0;
        *s = 0.0;
    } else {
        *c = a / r;
        *s = -b / r;
    }
}

double eliminate_element(double *A, int n, int m, int i, int j) {
    int idx_jm1 = i * (m + 1) + (j - 1);
    int idx_j = i * (m + 1) + j;

    double a = A[idx_jm1];
    double b = A[idx_j];
    double c, s;
    givens_rotation(a, b, &c, &s);

    // Apply rotation to row i
    A[idx_jm1] = c * a - s * b;
    A[idx_j] = 0.0;

    // Generate bulge in row i+1
    double bulge = 0.0;
    if (i + 1 < n) {
        int ip1 = i + 1;
        int col_in_ip1 = j - 1; // Corresponding column in next row
        if (col_in_ip1 >= 0) {
            int idx_ip1 = ip1 * (m + 1) + col_in_ip1;
            bulge = s * A[idx_ip1];
            A[idx_ip1] = c * A[idx_ip1];
        }
    }
    return bulge;
}

void chase_bulge(double *A, int n, int m, int start_row, double bulge) {
    int i = start_row;
    while (i < n - 1 && fabs(bulge) > TOL) {
        int row = i + 1;
        int col = 1; // Bulge always appears in first superdiagonal during chase
        
        // Get elements to eliminate bulge
        int idx = row * (m + 1) + 0; // Diagonal element
        double diag = A[idx];
        
        double c, s;
        givens_rotation(diag, bulge, &c, &s);
        
        // Update diagonal and superdiagonal
        A[idx] = c * diag - s * bulge;
        if (m >= 1) {
            A[row * (m + 1) + 1] = c * A[row * (m + 1) + 1];
        }
        
        // Propagate bulge to next row
        if (row + 1 < n) {
            int next_idx = (row + 1) * (m + 1) + 0;
            bulge = s * A[next_idx];
            A[next_idx] = c * A[next_idx];
        } else {
            bulge = 0.0;
        }
        
        i++;
    }
}

void tridiagonalize_band(double *A, int n, int m, double *d, double *e) {
    for (int i = 0; i < n; i++) {
        for (int j = m; j >= 2; j--) {
            if (i + j >= n) continue; // Stay within matrix bounds
            int idx = i * (m + 1) + j;
            if (fabs(A[idx]) > TOL) {
                double bulge = eliminate_element(A, n, m, i, j);
                chase_bulge(A, n, m, i, bulge);
            }
        }
    }

    // Extract results
    for (int i = 0; i < n; i++) {
        d[i] = A[i * (m + 1) + 0];
        if (i < n - 1) e[i] = A[i * (m + 1) + 1];
    }
}

int main() {
    int n = 5, m = 2;
    double A[] = {
        4.0, 1.0, 2.0,
        3.0, 1.0, 2.0,
        4.0, 1.0, 0.0,
        5.0, 1.0, 0.0,
        3.0, 0.0, 0.0
    };
    
    double d[n], e[n-1];
    tridiagonalize_band(A, n, m, d, e);

    printf("Diagonal (d): ");
    for (int i = 0; i < n; i++) printf("%.6f ", d[i]);
    printf("\nSubdiagonal (e): ");
    for (int i = 0; i < n-1; i++) printf("%.6f ", e[i]);
    printf("\n");
    
    return 0;
}