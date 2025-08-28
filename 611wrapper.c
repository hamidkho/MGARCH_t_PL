/*  TOMS 611 MINIMIZOR WRAPPER USING f2c */
/* Make call to opt() function below */
/* For compiling 
-- in that order, at the end of the command line, as in
                cc *.o -lf2c -lm
*/

/* THIS USES f2c.h definition of integer which seems to version dependent */



#include "611wrapper.h"

int (*fun_local)(gsl_vector *x, double *f);

integer n_611,iv_print;
gsl_vector *x0;
int iflag;

/* Added */
gsl_matrix *Lower;


int opt( gsl_vector *theta, int (*fun)(gsl_vector *x, double *f), int iprint, gsl_matrix *Lower_H )
{
  n_611=theta->size;
  x0 = gsl_vector_calloc(n_611);
  gsl_vector_memcpy(x0,theta);

  if( iprint <= 0){
    iv_print = 0;
  }else{
    iv_print = iprint;
  }

  /* Added */
  Lower = gsl_matrix_calloc(n_611,n_611);

  fun_local = fun;
  MAIN__( );

  /* Added */
  gsl_matrix_memcpy( Lower_H, Lower);

  gsl_vector_memcpy(theta,x0);

  gsl_vector_free(x0);
  /* Added */
  gsl_matrix_free(Lower);

  return iflag;
}


int calcf(integer *p, double x[0], integer *nf, double *f, integer ui[0], double ur[0], int (*ufun)(void) )
{
  gsl_vector *theta = gsl_vector_alloc(*p);
  double funct_value;
  int ret_flag,i;

  for(i=0;i< *p;i++)
    vset(theta,i,x[i]);

  ret_flag = fun_local(theta, &funct_value);

  gsl_vector_free(theta);

  if( ret_flag == 0){
    *f = funct_value;
    return 0;
  }else{
    *nf = -1;
    return -19;
  }
}

/* Dummy routine to satisfy the optimizor call sequence */
int uf(void)
{
  return 0;
}

int MAIN__( void )
{
  integer i,j;
  integer p=n_611, dn=2;
  integer LIV=100 + 82+p,LV=100 + 77+p*(p+17)/2; 
  //  gsl_vector_long *IV = gsl_vector_long_calloc(LIV);
  integer *IV;

  IV = (integer *)malloc( sizeof(integer)*LIV );

  gsl_vector *V = gsl_vector_calloc(LV);
  gsl_vector *d = gsl_vector_calloc(p);

  integer UI[10];
  integer nf;
  double UR[10];
  double f;

  iflag=-19;

  for(i=0;i<p;i++){
    if( vget(x0,i) != 0.0 ){
      vset(d,i, fabs(1.0/vget(x0,i)) );
    }else{
      vset(d,i,1.0);
    }
  }
  //  gsl_vector_set_all(d,1.0); 

  deflt_(&dn,IV,&LIV,&LV,&(V->data[0]));

  if( iv_print > 0){
    IV[18]=iv_print;
  }else{
    IV[20]=0;
  }

  smsno_(&p,&(d->data[0]),&(x0->data[0]),calcf,IV,&LIV,&LV,&V->data[0],&UI[0],&UR[0],uf);

  int iv_code = IV[0];
  if( iv_print ){
    printf("IV(0)=%d\n",IV[0]);
  }

  if( iv_code >= 3 && iv_code <= 6){
    if( iv_print ){
      printf("\n **NORMAL CONVERGENCE** \n\n");
    }
    iflag=0;
  }else if( iv_code >= 7 && iv_code <= 10){
    /* Second Try */    
    printf("\n **FAIL** \n\n");
    IV[0]=0;
    deflt_(&dn,IV,&LIV,&LV,&(V->data[0]));

    if( iv_print > 0){
      IV[18]=iv_print;
    }else{
      IV[20]=0;
    }
    smsno_(&p,&(d->data[0]),&(x0->data[0]),calcf,IV,&LIV,&LV,&(V->data[0]),&UI[0],&UR[0],uf);
    iv_code = IV[0];
    if( iv_code >= 3 && iv_code <= 6){
      if( iv_print ){
	printf("\n **NORMAL CONVERGENCE** \n\n");
      }
      iflag=0;
    }
  }


  /* recover Lower triangle of Quasi-Newton Estimate of Hessian */
  int ij=0;
  for(i=0;i<p;i++){
    for(j=0;j<=i;j++){
      mset(Lower,i,j, V->data[ IV[41] + ij-1] );
      ij++;
    }
  }


  gsl_vector_free(d);
  gsl_vector_free(V);
  //  gsl_vector_long_free(IV);
  free( IV );

  return 0;
}
