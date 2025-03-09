#include "./code_1.h"

static int solve(double * const restrict A, double * const restrict b, int * const restrict p, int const n){

    // Initialize permutations 
    for (int i = 0; i < n; i++) p[i] = i;

    // Perform the LU in-place in A
    double absval;
    for (int i = 0; i < n; i++){
        double max = 0.0;
        int imax = i;

        for (int k = i; k < n; k++){
            if ((absval = fabs(A[k*n + i])) > max){ 
                max = absval;
                imax = k;
            }
        }
        if (max < 1e-9) return -1; // non-invertible

        if (imax != i){
            int tmp = p[i];
            p[i] = p[imax];
            p[imax] = tmp;

            // Pay the swapping costs upfront (and not during forward/backward solve)
            double temp;
            for (int j = 0; j < n; j++) {
                temp = A[i*n+j];
                A[i * n + j] = A[imax * n + j];
                A[imax * n + j] = temp;
            }

            temp = b[i];
            b[i] = b[imax];
            b[imax] = temp;
        }

        // Keep a read-only pointer to the row.
        // This is a constant pointer to constant data.
        // The restrict keyword ensure the compiler that 
        // the data can only be accessed through this pointer.
        double const * const restrict Ai = &A[i*n]; 
        double const Aii_inv = 1.0f/Ai[i]; // Invert the pivot only once
        for (int j = i+1; j < n; j++){
            // Do the same as Ai but allow to modify the data
            double * const restrict Aj = &A[j*n];

            Aj[i] *= Aii_inv;
            double const c = Aj[i];
            for (int k = i+1; k < n; k++) Aj[k] -= c * Ai[k];
        }
    }


    // forward substitution
    for (int i = 0; i < n; i++){
        double const * const restrict Ai = &A[i*n]; 
        double x = b[i];
        for (int k = 0; k < i; k++){
            x -= Ai[k] * b[k];
        }
        b[i] = x;
    }


    // Backward substitution

    for (int i = n-1; i >= 0; i--){
        double const * const restrict Ai = &A[i*n]; 
        double x = b[i];
        for (int k = i+1; k < n; k++){
            x -= Ai[k] * b[k];
        }
        b[i] = x / Ai[i];
    }

    return 0;
}


int solve_linear_system(double * A, double * b, int * p, int n){
    return solve(A, b, p, n);
}


