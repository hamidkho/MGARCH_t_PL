#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>   
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>

#include "mydefn.h"
#include "matrix.h"
#include "611wrapper.h"


int mode_hessian( int (*log_lgl)(gsl_vector *, double *), int (*log_prior)(gsl_vector *, double *),
		  int npar_block, gsl_vector_int *theta_ind, const gsl_vector *theta, 
		  gsl_vector *theta_new, gsl_matrix *chol_hessian );
/* int mode_hessian( int npar_block, gsl_vector_int *theta_ind, const gsl_vector *theta,  */
/* 		  gsl_vector *theta_new, gsl_matrix *chol_hessian); */
int post_sub(gsl_vector *theta_s, double *f);
void copy_block_to_full ( gsl_vector *vec_block, gsl_vector *vec_full, gsl_vector_int *vec_ind, int direction);
