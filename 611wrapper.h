#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "matrix.h"
#include "mydefn.h"

#include "f2c.h"
#include "611new.h"


int uf(void);
int MAIN__(void);
int opt( gsl_vector *x0, int (*fun)(gsl_vector *x, double *f), int iprint, gsl_matrix *Lower_H );
int calcf(integer *p, double x[0], integer *nf, double *f, integer ui[0], double ur[0], int (*ufun)(void) );
