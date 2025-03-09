#include <stdio.h>

#define N 5   // Matrix size
#define K 2   // Bandwidth

void print_non_diag_non_subdiag(double *A) {
    printf("Elements not on the diagonal or first subdiagonal:\n");

    for (int i = 0; i < N; i++) {
        for (int j = 2; j <= K; j++) {  // Start from 2 to exclude diagonal and first subdiagonal
            int index = i * (K + 1) + j;  // Compute compact storage index
            if (index < N * (K + 1)) {  // Ensure index is within bounds
                printf("%lf ", A[index]);
            }
        }
    }
    printf("\n");
}

int main() {
    // Symmetric band matrix stored in row-major compact format (N × (K + 1))
    double A[] = {
        0, 0, 4,  // Row 0
        0, 1, 3,  // Row 1
        2, 1, 4,  // Row 2
        2, 1, 5,  // Row 3
        0, 1, 3   // Row 4
    };

    print_non_diag_non_subdiag(A);

    return 0;
}
