#include"UO.h"

extern int T;
extern gsl_matrix *simRt;
extern gsl_vector *omega_t;
gsl_vector_view Rt;

double utility(double x,double A) {
   
    
    double uf = 0.0,tmp;
   
    vset(omega_t,0,x);
    vset(omega_t,1,1.0-x);
    for (int z = 0; z < T; z++)
    {
        Rt = gsl_matrix_row(simRt,z);
       // pvec(vv(Rt));
       // pvec(omega_t);
       // mypause();
        gsl_blas_ddot(omega_t,vv(Rt),&tmp);
        
        //printf("%g\n",(pow((1+tmp/100.0),1-A)/(1-A)));
        //mypause();
        //CRRA
         uf += 1.0/T *(pow((1+tmp/100.0),1-A)/(1-A));
        //CARA
        // uf += 1.0/T *((1.0-exp(-A*(1+tmp/100.0)))/(A));
    }
    
    return uf;
}


double utility_EU(double x,double A) {
   
    gsl_vector *cumEU = gsl_vector_calloc(T);
    double uf = 0.0,tmp;
    vset(omega_t,0,x);
    vset(omega_t,1,1.0-x);
    for (int z = 0; z < T; z++)
    {
        Rt = gsl_matrix_row(simRt,z);
       // pvec(vv(Rt));
       // pvec(omega_t);
       // mypause();
        gsl_blas_ddot(omega_t,vv(Rt),&tmp);
        //printf("%g\n",(pow((1+tmp/100.0),1-A)/(1-A)));
        //mypause();
        //CRRA
        uf += (pow((1+tmp/100.0),1-A)/(1-A));
        //CARA
        // uf += ((1.0-exp(-A*(1+tmp/100.0)))/(A));
        vset(cumEU,z,1.0/(z+1.0)*uf);
    }
    // char filename[50];  // Buffer to hold the filename
    // snprintf(filename, sizeof(filename), "cumEU(%d).out", n);  // Create filename

    // FILE *f = fopen(filename, "w");
    // if (f == NULL) {
    //     perror("Error opening file");
    //     exit(EXIT_FAILURE);
    // }
    // for (size_t l = 0; l < T; ++l) {
    //     fprintf(f, "%.10g\n", gsl_vector_get(cumEU, l));  // Write each element to the file
    // }

    // fclose(f);
    gsl_vector_free(cumEU);
    return 1.0/T *uf;
}




double golden_section_search(double (*func)(double,double),double A, double a, double b, double tol) {
    double gr = (sqrt(5) + 1) / 2; // golden ratio
    double c = b - (b - a) / gr;
    double d = a + (b - a) / gr;

    while (fabs(c - d) > tol) {
        if (func(c,A) < func(d,A)) {
            a = c;
        } else {
            b = d;
        }
        c = b - (b - a) / gr;
        d = a + (b - a) / gr;
    }

    return (b + a) / 2;
}





double R_utility(double x,double A,gsl_vector *row_y) {
   
    
    double uf = 0.0,tmp;
    vset(omega_t,0,x);
    vset(omega_t,1,1.0-x);
        gsl_blas_ddot(omega_t,row_y,&tmp);
        //printf("%g\n",(pow((1+tmp/100.0),1-A)/(1-A)));
        //mypause();
        //CRRA
        uf = (pow((1+tmp/100.0),1-A)/(1-A));
        //CARA
        //  uf = (1.0-exp(-A*(1+tmp/100.0)))/(A);
    
    
    return uf;
}