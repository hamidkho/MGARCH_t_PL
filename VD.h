#include <stdio.h>
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

int log_prior(gsl_vector *theta, double *f);
void make_CCT_AAT_BBT (gsl_vector *theta,gsl_matrix *CCT,gsl_matrix *AAT,gsl_matrix *BBT);
int lgl(gsl_vector *theta, double *f);
double LFun(double rp,double VaR, double ES, double e);

