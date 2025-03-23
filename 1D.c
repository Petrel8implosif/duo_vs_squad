#include "devoir_1.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#define PI 3.14159265358979323846
#define SQUARE(x) ((x) * (x))
double vector_norm(double *v, int n) {
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += v[i] * v[i];
    }
    return sqrt(sum);
}

// Function to perform Householder transformation
void tridiagonalize_full(double *A, int n, int k, double *d, double *e) {
    for (int k = 0; k < n - 2; k++) {
        double *v = (double *)calloc(n, sizeof(double));
        double norm_x = 0.0;

        // Compute the Householder vector
        for (int i = k + 1; i < n; i++) {
            v[i] = A[i *n + k];
        }
        norm_x = vector_norm(&v[k + 1], n - k - 1);
        if (v[k + 1] >= 0) {
            v[k + 1] += norm_x;
        } else {
            v[k + 1] -= norm_x;
        }

        double vtv = 0.0;
        for (int i = k + 1; i < n; i++) {
            vtv += v[i] * v[i];
        }
        if (vtv < 1e-10) { // Avoid division by zero
            free(v);
            continue;
        }
        double beta = 2.0 / vtv;

        // Compute p = beta * A * v
        double *p = (double *)calloc(n, sizeof(double));
        for (int i = 0; i < n; i++) {
            for (int j = k + 1; j < n; j++) {
                p[i] += A[i *n + j] * v[j];
            }
            p[i] *= beta;
        }

        // Compute t = beta * (v ⋅ p) / 2
        double t = 0.0;
        for (int j = k + 1; j < n; j++) {
            t += v[j] * p[j];
        }
        t *= beta;
        t /= 2.0;

        // Compute q = p - t * v
        double *q = (double *)malloc(n * sizeof(double));
        for (int i = 0; i < n; i++) {
            q[i] = p[i] - t * v[i];
        }

        // Update A: A = A - v*q^T - q*v^T
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                A[i *n + j] -= v[i] * q[j] + q[i] * v[j];
            }
        }

        

        free(v);
        free(p);
        free(q);
    }
    // update d and e with the diagonal and subdiagonal elements respectively
    for (int k = 0; k < n-1; k++) {
    d[k] = A[k *n + k];
    e[k] = A[(k + 1) *n + k];
    }
    d[n-1] = A[(n-1) *n + (n-1)];
}
int step_qr_tridiag(double *d, double *e, int m, double eps){
    double *A = (double *)calloc(m * m, sizeof(double));
    for (int i = 0; i < m; i++) {
        A[i * m + i] = d[i];
        if (i < m - 1) {
            A[i * m + (i + 1)] = e[i];
            A[(i + 1) * m + i] = e[i];
        }
    }
    // Print matrix A
    printf("Matrix A:\n");
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            printf("%f ", A[i * m + j]);
        }
        printf("\n");
    }
    double t_nn = A[m * m-1];  // calcul du shift de wilkinson
    double t_n1n1 = A[(m-2) * m + (m-2)];
    double t_n1n = A[(m-1) * m + (m-2)];

    double dd = (t_n1n1 - t_nn) / 2.0;
    if (dd == 0) {
        double mu = t_nn - fabs(t_n1n);
    }
    double sign_d = (dd > 0) - (dd < 0); // Equivalent à sign(d)
    double mu = t_nn + dd - sign_d * sqrt(dd * dd + t_n1n * t_n1n);

 
    

    for (int i = 0; i < m; i++) {
        A[i * m + i] -= mu;
    }

    double *cos = (double *)calloc(m-1, sizeof(double));
    double *sin = (double *)calloc(m-1, sizeof(double));

    for (int i = 0; i < m - 1; i++) {
        double c = 1.0;
        double s = 0.0;
        double a = A[i * m + i];
        double b = A[(i+1) * m + i];
        double hypot = sqrt(a * a + b * b);
        if(hypot < 1e-10){
            c = 1.0;
            s = 0.0;
        }
        else{
            c = a / hypot;
            s = -b / hypot;
        }
        cos[i] = c;
        sin[i] = s;


        for (int k = i; k < i+3 && k < m; k++) {  // gauche 
            if(k == m || k == i+2) {
            double R_ik = A[i * m + k];
            double R_jk = A[(i+1) * m + k];
            A[i * m + k] = c * R_ik - s * R_jk;
            A[(i+1) * m + k] = c * R_jk;
            }
            else {
            double R_ik = A[i * m + k];
            double R_jk = A[(i+1) * m + k];
            A[i * m + k] = c * R_ik - s * R_jk;
            A[(i+1) * m + k] = s * R_ik + c * R_jk;
            }
        }
    }

        double c1 = cos[0];
        double s1 = sin[0];
    
        double R_ik = A[0];
        double R_jk = A[1];
        A[0] = c1 * R_ik - s1 * R_jk;
        A[1] = s1 * R_ik + c1 * R_jk;
        
        double rR_ik = A[m];
        double rR_jk = A[m + 1];
        A[m] = c1 * rR_ik - s1 * rR_jk;
        A[m + 1] = c1 * rR_jk;
    



    for(int i = 0; i < m-2; i++){
        //droite
            double c = 1.0;
            double s = 0.0;
            c = cos[i+1];
            s = sin[i+1];
            for (int k = i; k < i+3; k++) {  
                    double R_ik = A[k * m + i+1];
                    double R_jk = A[k * m + i+2];
                    A[k * m + i+1] = c * R_ik - s * R_jk;
                    A[k * m + i+2] = s * R_ik + c * R_jk;
        }
    }
    
    free(cos);
    free(sin);
    // Print matrix A after QR step
    printf("Matrix A after QR step:\n");
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            printf("%f ", A[i * m + j]);
        }
        printf("\n");
    }
    for(int i = 0; i < m; i++){
        A[i * m + i] += mu;
    }
    for (int i = 0; i < m; i++) {
        d[i] = A[i * m + i];
    }
    for (int i = 0; i < m - 1; i++) {
        e[i] = A[i * m + (i+1)];
    }

    if(fabs(A[(m-1) * m + (m-2)]) > eps * (fabs(A[(m-1) * m + (m-1)]) + fabs(A[(m-2) * m + (m-2)]))){
        free(A);
        return m;
    }
    free(A);
    return m - 1;
}

int qr_eigs_full( double *e, int n, int k, double eps, int max_iter, double *d){
    int m = n;
    int i = 0;
    for (i = 0; i < max_iter && m > 1; i++) {
        m = step_qr_tridiag(d, e, m, eps);
    }
    if (i == max_iter) {
        return -1;
    }
    free(e);
    return i;
}

int main(void) {
    // Choose the number of interior nodes 👎 so that there are n+2 grid points
    int n = 100; // Feel free to change this value or loop over several n for convergence study
    double h = 2.0 / (n + 1);  // step size on [-1,1]

    // Allocate and build the discrete Laplacian matrix (n x n)
    double *dr = (double *)calloc(n, sizeof(double));
    double *er = (double *)calloc(n-1, sizeof(double));
    for (int i = 0; i < n-1; i++) {
        // Diagonal entries: -2/h^2
        dr[i] = -2.0 / (h * h);
        er[i] = 1.0 / (h * h);
    }
    dr[n-1] = -2.0 / (h * h);

    // Set parameters for the QR algorithm
    double eps = 1e-14;
    int max_iter = 1000;

    int iter = qr_eigs_full(er, n, 1, eps, max_iter, dr);
    if (iter < 0) {
        printf("QR algorithm did not converge.\n");
    } else {
        printf("QR algorithm converged in %d iterations.\n", iter);
    }

    // Print computed eigenvalues (these approximate the eigenvalues of the FD Laplacian)
    printf("\nComputed eigenvalues (discrete Laplacian):\n");
    for (int i = 0; i < n; i++) {
        printf("Eigenvalue %d: %e\n", i + 1, dr[i]);
    }

    // Compare with the analytical eigenvalues.
    // For the continuous problem u'' = λ u on [-1,1] with u(-1)=u(1)=0,
    // the eigenvalues are: λ_k = - (kπ/2)^2, k = 1, 2, ..., n.
    printf("\nAnalytical eigenvalues:\n");
    for (int k = 1; k <= n; k++) {
        double lambda_exact = -SQUARE(k * PI / 2.0);
        printf("Eigenvalue %d: %e\n", k, lambda_exact);
    }

    // Optionally, compute a total error norm for the first n eigenvalues.
    double error_norm = 0.0;
    for (int i = 0; i < n; i++) {
        // Note: In a more rigorous test, you might want to sort the eigenvalues
        // before comparing. Here we assume the computed eigenvalues are roughly in order.
        double lambda_exact = -SQUARE((i + 1) * PI / 2.0);
        error_norm += fabs(dr[i] - lambda_exact);
    }
    printf("\nTotal error norm (sum of absolute differences): %e\n", error_norm);

    free(er);
    free(dr);
    return 0;
}