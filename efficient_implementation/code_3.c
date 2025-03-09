#include <cblas.h>
#include <string.h>

#pragma GCC target("fma,avx,avx2,ssse3,sse,sse2")

/**
 * @brief Computes the QR decomposition of a matrix
 * @param A Matrix to be factorized, to be set with Q (!Column major!)
 * @param R Upper-triangular storage matrix (!Column major!)
 * @param P Permutation matrix of size n
 * @param m Number of rows of A
 * @param n Number of columns of A
 */
int QR(double *A, double *R, int *P, int m, int n) {
    int rank = n;
    for (int j = 0; j < n; j++) P[j] = j; // Setup the permutations

    for (int j = 0; j < rank; j++) {

        // ============== Gram-Schmidt : method 1 ====================
        // Compute all the scalar product at once !
        // but ... it produces non-negligible roundoff errors (bad for ill-conditioned matrices)

        /*
        // Rj = A^T*Aj ( R[k, j] = A[:, j] @ A[k, :] )
        cblas_dgemv(CblasColMajor, CblasTrans, m, j, 1.0, A, m, &A[j*m], 1, 0.0, &R[(j*(j+1))>>1], 1);

        // And subtract the linear combinations of columns at once too
        // Aj -= A*Rj
        cblas_dgemv(CblasColMajor, CblasNoTrans, m, j, -1.0, A, m, &R[(j*(j+1))>>1], 1, 1.0, &A[j*m], 1);
        */

        // ============== Gram-Schmidt : method 2 ====================
        // Slower, but does less roundoff errors
        
        for (int k = 0; k < j; k++) {
            R[((j*(j+1))>>1) + k] = cblas_ddot(m, A + k*m, 1, A + j*m, 1); 
            cblas_daxpy(m, -R[((j*(j+1))>>1) + k], A + k*m, 1, A + j*m, 1);
        }

        double const rjj = cblas_dnrm2(m, A + j*m, 1);
        if (rjj < 1e-12) {
            --rank;
            cblas_dswap(m, A + j*m, 1, A + rank*m, 1);  // swap A cols
            memset(&A[rank*m], 0, sizeof(*A)*m);

            cblas_dcopy(j, &R[(j*(j+1))>>1], 1, &R[(rank*(rank+1))>>1], 1);
            memset(&R[j + ((rank*(rank+1))>>1)], 0, sizeof(*R)*(rank-j));

            int tmp = P[j];
            P[j] = P[rank];
            P[rank] = tmp;
            j--;
            continue;
        }

        R[((j*(j+1))>>1) + j] = rjj;
        cblas_dscal(m, 1./rjj, A + j*m, 1);
    }

    return rank;
}
