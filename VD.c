#include "VD.h"

extern gsl_matrix *y;
extern gsl_matrix *e;
extern gsl_matrix *R_1;
extern gsl_matrix *sample_cov, *Htt1, *Ht1;
extern gsl_array *H;
extern int npar,st,ed,n,k,num_block;
extern gsl_rng *rng;

extern gsl_matrix *CCT;
extern gsl_matrix *AAT;
extern gsl_matrix *BBT;

void make_CCT_AAT_BBT (gsl_vector *theta,gsl_matrix *C,gsl_matrix *A,gsl_matrix *B)
{
    int i,j,m;
    gsl_vector_view block;
    gsl_matrix *C1=gsl_matrix_calloc(k,k);

    gsl_matrix_set_zero(C);
    gsl_matrix_set_zero(A);
    gsl_matrix_set_zero(B);
    gsl_matrix_set_zero(C1);

    block = gsl_vector_subvector ( theta, 0, k*(k+1)/2);
    vech_to_lower_matrix( vv(block) , C1 );
    gsl_matrix_memcpy(C,C1);
//    pmat(C);
    gsl_blas_dtrmm (CblasRight,CblasLower,CblasTrans,CblasNonUnit,1.0,C1,C);


    block = gsl_vector_subvector ( theta, k*(k+1)/2, k);
    gsl_blas_dsyr (CblasLower ,1.0, vv(block) , A);
    copy_lower_to_upper ( A );


    block = gsl_vector_subvector ( theta, k*(k+1)/2+k, k);
    gsl_blas_dsyr (CblasLower ,1.0, vv(block) , B);
    copy_lower_to_upper ( B );


    gsl_matrix_free( C1 );

}








int lgl(gsl_vector *theta, double *f)
{

    gsl_matrix *mat_work = gsl_matrix_calloc(k,k);
    gsl_matrix *eeT = gsl_matrix_calloc(k,k);
    gsl_vector *mu = gsl_vector_calloc(k);
    gsl_vector *eta = gsl_vector_calloc(k);
    gsl_vector *eps = gsl_vector_calloc(k);
    gsl_vector_view row,col_b,col;
    int i,j,m,err;
    double logl,logl_t,tmp,rp,psi;

    vset(mu,0,vget(theta,npar-2));
    vset(mu,1,vget(theta,npar-1));
    psi = vget(theta,npar-k-1);

    vset(eta,0,vget(theta,npar-2*k-1));
    vset(eta,1,vget(theta,npar-2*k));

    //  gsl_vector_set_zero( zero );

    make_CCT_AAT_BBT (theta,CCT,AAT,BBT);

    /* pmat(CCT); */
    /* pmat(AAT); */
    /* pmat(BBT); */
    /* pvec(mu); */

    logl=0.0;

    for(i=st;i<=ed;i++){

        if( i > st ){

            /* fill in V_t-I in the vector V_t */

            gsl_matrix_get_row ( eps, y, i-1 );
            gsl_vector_sub ( eps, eta );

            gsl_matrix_set_zero( eeT );
            gsl_blas_dsyr (CblasLower, 1.0, eps ,eeT);
            copy_lower_to_upper ( eeT );


            gsl_matrix_mul_elements( eeT, AAT);

            gsl_matrix_memcpy( mat_work, BBT);
            gsl_matrix_mul_elements( mat_work, H->mat[i-1]);

            gsl_matrix_add( mat_work, eeT);


            gsl_matrix_memcpy( H->mat[i] , CCT);
            gsl_matrix_add( H->mat[i], mat_work);
            if (i == ed) {

                gsl_matrix_get_row ( eps, y, i );
                gsl_vector_sub ( eps, eta );

//            gsl_vector_sub ( eps, vv(xi_v) );
//            pvec(vv(xi_v));
//            mypause();

                gsl_matrix_set_zero( eeT );
                gsl_blas_dsyr (CblasLower, 1.0, eps ,eeT);
                copy_lower_to_upper ( eeT );


                gsl_matrix_mul_elements( eeT, AAT);

                gsl_matrix_memcpy( mat_work, BBT);
                gsl_matrix_mul_elements( mat_work, H->mat[i]);

                gsl_matrix_add( mat_work, eeT);


                gsl_matrix_memcpy( Htt1 , CCT);
                gsl_matrix_add( Htt1, mat_work);



            }
        }else{
            /* *************** ARBITRARY startup *************** */
            gsl_matrix_memcpy( H->mat[i] , sample_cov );
            //      gsl_matrix_scale( H->mat[i], 1./(1.0 - mget(AAT,0,0)-mget(BBT,0,0)) );
        }

        row = gsl_matrix_row( y, i);

        //err = log_mnv_den_custom( vv(row) , mu , H->mat[i], &logl_t  );
        err = log_mvt_den_custom(vv(row) , mu ,psi,H->mat[i], &logl_t );
        //    printf("%d\n",i);
/*     logl_t = log_mvt_den( vv(row) , zero , nu , H->mat[i]); */
/*     err = GSL_SUCCESS; */

//    err = log_mvt_den_custom( vv(row) , zero , nu, H->mat[i], &logl_t  );

        //    pmat(H->mat[i]);
        /* logl_t = log_mnv_den( vv(row) , mu , H->mat[i] ); */
        /* err = GSL_SUCCESS; */

/*     //    printf("err=%d %d\n",err, GSL_SUCCESS); */
        if( err == GSL_SUCCESS ){
            logl += logl_t;
        }else{
            //      printf("err i = %d %d\n",err,i);
/*       pmat(H->mat[i]); */
/*       printf("CCT AAT BBT\n"); */
/*       pmat(CCT); */
/*       pmat(AAT); */
/*       pmat(BBT); */
            gsl_matrix_free(mat_work); 	/* FAILURE of FUNCTION EVALUATION -- cleanup */
            gsl_matrix_free(eeT);
            gsl_vector_free( mu );
            gsl_vector_free( eta );
            gsl_vector_free( eps );
            *f = -2.0e200;
            return 0;
            //      return -19;
        }

        //    logl += log_mnv_den( vv(row) , zero , H->mat[i] );


    }

    /* cleanup */
    gsl_matrix_free(mat_work);
    gsl_matrix_free(eeT);
    gsl_vector_free( mu );
    gsl_vector_free( eta );
    gsl_vector_free( eps );

    if( fabs(logl) < 2.0e300 && fabs(logl) > 2.0e-300 ){
        *f = logl;
        return 0;
    }else{
        //    printf("logl is NaN %g\n",logl);
        return -19;
        //    *f = 2.0e200;
    }
}


int log_prior(gsl_vector *theta, double *f)
{
    int i,j;
    double lprior;

    /* SET FOR k=2 ONLY************** */

    /* identification restriction */
    if(  vget(theta,0) < 0.0 || vget(theta,k*(k+1)/2) < 0.0 || vget(theta,k*(k+1)/2+k) < 0.0
         || vget(theta,2) < 0.0 ){
        return -19;
    }


    /* IGNORE stationarity */
    make_CCT_AAT_BBT (theta,CCT,AAT,BBT);
/*   gsl_matrix_scale( AAT, nu/(nu-2.0) ); */

/*   gsl_matrix_set_all( CCT, 1.0); */
/*   gsl_matrix_sub(CCT,AAT); */
/*   gsl_matrix_sub(CCT,BBT); */

/*   double m_min = gsl_matrix_min( CCT ); */

/*   if( m_min <= 0.0 ){ */
/*     //    pmat(CCT); */
/*     //    exit(0); */

/*     return -19; */
/*   } */




    lprior = 0.0;
    double psi = vget(theta,npar-k-1);
    for(i=0;i<npar;i++){
        lprior += log_nor( vget( theta,i), 0.0, 100.0 );
    }
    if (psi <= 100.0 && psi >= 2.0 )lprior += log(1.0/98.0)/log_nor( psi, 0.0, 100.0) ;
    else lprior = 2.0e-300;

    if( fabs(lprior) < 2.0e300 && fabs(lprior) > 2.0e-300 ){
        *f = lprior;
        //    printf(" log_prior = %g\n",lprior);
        return 0;
    }else{
      //  printf("lprior is NaN, %g\n",lprior);
       // pvec(theta);
        //exit(0);
        return -19;
        //    *f = 2.0e200;
    }
}

double LFun(double rp,double VaR, double ES, double e){

    double f = VaR/ES+log(-ES)-1.0;
    if (rp <= VaR) f += -1.0/(e*ES) * (VaR - rp);
    
    return f;
    
    }