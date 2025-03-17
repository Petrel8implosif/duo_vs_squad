

#include <stdio.h>

int main() {
    double A[] = {
        0.0, 0.0, 4.0,  // Row 0
        0.0, 1.0, 3.0,  // Row 1
        2.0, 1.0, 4.0,  // Row 2
        2.0, 1.0, 5.0,  // Row 3
        0.0, 1.0, 3.0   // Row 4
    };

    int rows = 5;
    int cols = 3;

    // Create a new array with space for the additional column
    double B[rows * (cols + 1)];

    // Copy A into B and add a column of zeros
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            B[i * (cols + 1) + (j + 1)] = A[i * cols + j];
        }
        B[i * (cols + 1)] = 0.0; // Add zero at the beginning of each row
    }

    // Print the resulting matrix
    printf("Matrix after adding a column of zeros:\n");
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols + 1; j++) {
            printf("%.1f ", B[i * (cols + 1) + j]);
        }
        printf("\n");
    }

    return 0;
}
