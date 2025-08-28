#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include "matrix.h"
#include "mydefn.h"

#include "ps.h"
#include "ran.h"
#include "distr.h"
#include "array.h"

#define GOLDEN_RATIO 1.61803398875
#define EPSILON 1e-6  // Precision of the solution

double utility(double x,double A);
double utility_EU(double x,double A);

double golden_section_search(double (*func)(double,double),double A, double a, double b, double tol);
double R_utility(double x,double A,gsl_vector *row_y);