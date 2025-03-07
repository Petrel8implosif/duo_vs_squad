#ifndef CODE4_H
#define CODE4_H

void tridiagonalize_band(double *A, int n, int k, double *d, double *e);
int step_qr_tridiag(double *d, double *e, double m, double eps);
int qr_eigs_band(double *A, int n, int k, double eps, int max_iter);
#endif