#include "devoir_1.h"
#include <math.h>
#include <stdlib.h>

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
        e[k+1] = A[(k + 1) *n + k];
    }
    d[n-1] = A[(n-1) *n + (n-1)];
}
int step_qr_tridiag(double *d, double *e, int m, double eps){
    double *A = (double *)calloc(m * m, sizeof(double));
    for (int i = 0; i < m; i++) {
        A[i * m + i] = d[i];
        if (i < m - 1) {
            A[i * m + (i + 1)] = e[i+1];
            A[(i + 1) * m + i] = e[i+1];
        }
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
    for(int i = 0; i < m; i++){
        A[i * m + i] += mu;
    }
    for (int i = 0; i < m; i++) {
        d[i] = A[i * m + i];
    }
    for (int i = 0; i < m - 1; i++) {
        e[i+1] = A[i * m + (i+1)];
    }

    if(fabs(A[(m-1) * m + (m-2)]) > eps * (fabs(A[(m-1) * m + (m-1)]) + fabs(A[(m-2) * m + (m-2)]))){
        free(A);
        return m;
    }
    free(A);
    return m - 1;
}

int qr_eigs_full(double *A, int n, int k, double eps, int max_iter, double *d){
    double *dr = (double *)calloc(n, sizeof(double));
    double *e = (double *)calloc(n, sizeof(double));
    tridiagonalize_full(A, n, k, dr, e);
    int m = n;
    int i = 0;
    int m_temp = m;
    for (i = 0; i < max_iter && m > 1; i++) {
        m = step_qr_tridiag(dr, e, m, eps);
        if (m != m_temp) {
            d[m] = dr[m];
        }
        m_temp = m;
    }
    d[0] = dr[0];
    if (i == max_iter) {
        return -1;
    }
    free(dr);
    free(e);
    return i;
}