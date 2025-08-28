/* code to optimiza a block of parameters in the posterior distribution */
#include "post-opt.h"



/* set vector views that are available for the whole file */
gsl_vector_view theta_fixed;
gsl_vector_int_view theta_ind_local;
gsl_vector *theta_full;

int npar_local;

/* function ptr to lgl and log_prior to be used locally  */
int (*log_lgl_local)(gsl_vector *,double *);
int (*log_prior_local)(gsl_vector *,double *);





/* Takes in int indicator vector (=1 for element to be optimated, =0 for fixed
   element) and returns estimate of mode and hessian */
/* Works on whole vector */
int mode_hessian( int (*log_lgl)(gsl_vector *, double *), int (*log_prior)(gsl_vector *, double *),
		  int npar_block, gsl_vector_int *theta_ind, const gsl_vector *theta, 
		  gsl_vector *theta_new, gsl_matrix *chol_hessian )
{
  int iflag;
  double f;
  gsl_vector_view theta_block;
  npar_local = theta -> size;
  gsl_vector *vec_work = gsl_vector_calloc(npar_local);
  gsl_matrix *Lower = gsl_matrix_calloc(npar_block,npar_block);
  
  log_lgl_local = log_lgl;
  log_prior_local = log_prior;

  theta_ind_local = gsl_vector_int_subvector( theta_ind ,0,npar_local );

  theta_fixed = gsl_vector_subvector( theta,0,npar_local);

  theta_full = gsl_vector_calloc( npar_local );

  gsl_vector_memcpy( theta_full , theta);

  /* copy theta to theta_new */
  gsl_vector_memcpy( theta_new,theta);

  /* copy theta parts to theta_block */
  theta_block = gsl_vector_subvector( vec_work, 0, npar_block);
    
  /* copy full to block based on theta_ind_local == 1 */
  copy_block_to_full ( vv(theta_block), theta_full , vv(theta_ind_local), 1);


  iflag = post_sub( vv(theta_block), &f); 

  if( iflag != 0){
    printf("******* Starting point gives iflag=%d in mode_hessian()******\n",iflag);
    exit(0);
  }

/*   printf("iflag = %d\n",iflag); */
/*   printf("f = %g\n",f); */

  /* call optimizor */
  //  printf("return before iflag \n");

  //  iflag = opt( vv(theta_block) , &post_sub , -1, Lower);

  iflag = opt( vv(theta_block) , &post_sub , 30, Lower);
  //  printf("return iflag = %d\n",iflag);

  /* get hessian if optimium achieved. */
/*   if( iflag == 0){ */
/*     //    printf("\nMode for the block\n"); */
/*     //    pvec( vv(theta_block) ); */
/*     //        num_hessian ( &post_sub , vv(theta_block), hessian ); */
/* /\*     inv_pd_sym( hessian ); *\/ */
/* /\*     printf("\nHessian for the block\n"); *\/ */
/* /\*     pmat( hessian ); *\/ */
/* /\*     inv_pd_sym( hessian ); *\/ */
/*     gsl_blas_dsyrk (CblasLower,CblasNoTrans,1.0,Lower,0.0,hessian); */
/*     copy_block_to_full ( vv(theta_block), theta_new , vv(theta_ind_local), -1); */
/*   }else{ */
/*     /\* FAILURE -- TRY TO GET SOME PROPOSAL *\/ */
/*     num_hessian ( &post_sub , vv(theta_block), hessian ); */
/*     copy_block_to_full ( vv(theta_block), theta_new , vv(theta_ind_local), -1); */
/*   } */
  
  if(iflag != 0){
    printf("Failure on OPT -- mode_hessian \n");
  }
  gsl_matrix_memcpy( chol_hessian, Lower);

  copy_block_to_full ( vv(theta_block), theta_new , vv(theta_ind_local), -1);

  gsl_vector_free( vec_work );
  gsl_vector_free ( theta_full );
  gsl_matrix_free (Lower);
  
  return iflag;
}




/* full posterior function to MIN */
int post_sub(gsl_vector *theta_s, double *f)
{
  int iflag,iflag2,i,j;
  gsl_vector *theta_full_local = gsl_vector_alloc(npar_local);
  double priorf;

  gsl_vector_memcpy(theta_full_local, vv(theta_fixed));
  j=0;
  for(i=0;i<npar_local;i++){
    if( vget_int( vv(theta_ind_local) ,i) == 1){
      vset(theta_full_local,i, vget(theta_s,j) );
      j++;
    }
  }  

  iflag = (*log_lgl_local)(theta_full_local,f);

  iflag2 = (*log_prior_local)( theta_full, &priorf );

  iflag = iflag + iflag2;
  *f =  - (*f + priorf);

  gsl_vector_free( theta_full_local );
  return iflag;
}


/* Copy vec_block elements to vec_full where indices of vec_ind = 1, if
   direction = -1 */
/* Copy vec_full elements using indices of vec_ind = 1 to vec_block, if
   direction = 1 */
void copy_block_to_full ( gsl_vector *vec_block, gsl_vector *vec_full, gsl_vector_int *vec_ind, int direction)
{
  size_t npar_block = vec_block->size;
  size_t npar_full = vec_full->size;
  int i,j,m;

  if( npar_block > npar_full ){
    printf("\n ********* ERROR in copy_block_to_full () *********\n");
    exit(-1);
  }


  if( direction == -1 ){		/* copy vec_block to vec_full */
    m=0;
    for(i=0;i<npar_full;i++){
      if( vget_int(vec_ind,i) == 1 ){
	vset( vec_full, i, vget(vec_block,m) );
	m++;
      }
    }
  }else if( direction == 1 ){
    m=0;
    for(i=0;i<npar_full;i++){
      if( vget_int(vec_ind,i) == 1 ){
	vset( vec_block, m, vget(vec_full,i) );
	m++;
      }
    }
  }else{
    printf("\n**** ERROR IN copy_block_to_full(), set direction =-1 or 1 ****\n");
    exit(-1);
  }

}






