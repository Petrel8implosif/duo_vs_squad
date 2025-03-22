#ifndef CODE4_H
#define CODE4_H

void tridiagonalize_full(double *A, int n, int k, double *d, double *e);
int step_qr_tridiag(double *d, double *e, int m, double eps);
int qr_eigs_full(double *A, int n, int k, double eps, int max_iter, double *d);
#endif