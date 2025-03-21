int step_qr_tridiag(double *d, double *e, int m, double eps){
    double *A = (double *)calloc(m * m, sizeof(double));
    for (int i = 0; i < m; i++) {
        A[i * m + i] = d[i];
        if (i < m - 1) {
            A[i * m + (i + 1)] = e[i];
            A[(i + 1) * m + i] = e[i];
        }
    }
    printf("Initial matrix A:\n");
    for (int row = 0; row < m; row++) {
        for (int col = 0; col < m; col++) {
            printf("%f ", A[row * m + col]);
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
        double c, s;
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
    printf("Gauche A:\n");
    for (int row = 0; row < m; row++) {
        for (int col = 0; col < m; col++) {
            printf("%f ", A[row * m + col]);
        }
        printf("\n");
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
            double c,s;
            c = cos[i+1];
            s = sin[i+1];
            for (int k = i; k < i+3; k++) {  
                    double R_ik = A[k * m + i+1];
                    double R_jk = A[k * m + i+2];
                    A[k * m + i+1] = c * R_ik - s * R_jk;
                    A[k * m + i+2] = s * R_ik + c * R_jk;
        }
        
        printf("Droite A:\n");
        for (int row = 0; row < m; row++) {
            for (int col = 0; col < m; col++) {
                printf("%f ", A[row * m + col]);
            }
            printf("\n");
        }
    }
    
    free(cos);
    free(sin);
    for(int i = 0; i < m; i++){
        A[i * m + i] += mu;
    }
    printf("Deshift A:\n");
    for (int row = 0; row < m; row++) {
        for (int col = 0; col < m; col++) {
            printf("%f ", A[row * m + col]);
        }
        printf("\n");
    }
    for (int i = 0; i < m; i++) {
        d[i] = A[i * m + i];
    }
    for (int i = 0; i < m - 1; i++) {
        e[i] = A[i * m + (i+1)];
    }
    printf("d: ");
    for (int i = 0; i < m; i++) {
        printf("%f ", d[i]);
    }
    printf("\ne: ");
    for (int i = 0; i < m - 1; i++) {
        printf("%f ", e[i]);
    }
    printf("\n");

    if(fabs(A[(m-1) * m + (m-2)]) > eps * (fabs(A[(m-1) * m + (m-1)]) + fabs(A[(m-2) * m + (m-2)]))){
        free(A);
        return m;
    }
    free(A);
    return m - 1;
}