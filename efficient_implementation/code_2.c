#include "code_2.h"
#include <cblas.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// yes I cheated, but now you know this exists :)
// but beware, it can lead to incorrect computations when misused !
// (go read the documentation and understand why)
#pragma GCC optimize("-funroll-loops")
#pragma GCC optimize("-ffast-math")
#pragma GCC target("fma,avx2,avx,sse,sse2")

/**
 * @brief Solves Ax = b in-place (overwriting b with x)
 * @param A Symmetric matrix (Upper half)
 * @param b Right-hand-side
 * @param n Size of the linear system
 */
int solve_linear_system(double *A, double *b, int n) {
    // NOTE : the operator "a >> b" is a binary operator that will shift 
    // the binary representation of "a" by "b" bits to the right 
    // ex : 8 >> 1 = 4 (bin(8) = 0b1000 >> 1 = 0b0100 = 4)
    // Which is equivalent to dividing the integer by 2^b

    // Pre-check SDP so we don't have to perform an "if" in a loop
    for (int i = 0; i < n; i++) {
        if (A[((i*i+i) >> 1) + i] <= 0.0) return -1;
    }

    // Cholesky factorisation
    double sum;
    for (int i = 0; i < n; i++) {
        double *const restrict Ai = &A[(i*i+i) >> 1]; // Use constant arrays
        for (int j = 0; j < i; j++) {
            double const *const restrict Aj = &A[(j*j+j) >> 1];

            sum = 0.0; // Use an accumulator
            for (int k = 0; k < j; k++) sum += Ai[k] * Aj[k];

            Ai[j] = (1.0 / Aj[j]) * (Ai[j] - sum);
        }

        sum = 0.0;
        for (int k = 0; k < i; k++) sum += Ai[k] * Ai[k];

        Ai[i] = sqrt(Ai[i] - sum);
    }

    // forward substitution
    for (int i = 0; i < n; i++){
        double const * const restrict Ai = &A[(i*i+i) >> 1];

        double x = 0.0;
        for (int j = 0; j < i; j++) x += Ai[j] * b[j];

        b[i] = (b[i] - x) / Ai[i];
    }
    // Approx same performance
    /*cblas_dtpsv(CblasRowMajor, CblasLower, CblasNoTrans, CblasNonUnit, n, A, b, 1);*/

    for (int i = n-1; i >= 0; i--){
        double x = 0.0;
        for (int j = i+1; j < n; j++) x += A[((j*j + j) >> 1) + i] * b[j];
        b[i] = (b[i] - x)/A[((i*i+i)>>1) + i];
    }
    // Approx same performance
    /*cblas_dtpsv(CblasRowMajor, CblasLower, CblasTrans, CblasNonUnit, n, A, b, 1);*/

    return 0;
}

