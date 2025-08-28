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

#include "VD.h"
#include "ps.h"
#include "ran.h"
#include "distr.h"
#include "array.h"
#include "post-opt.h"
#include "mvtdstpack.h"
#include "UO.h"

#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "sys/time.h"

gsl_matrix *y,*Y;
gsl_matrix *sample_cov, *Htt1, *Ht1,*simRt;
gsl_array *H;
int T,npar,st,ed,n,k,num_block;
gsl_vector *omega_t;

gsl_matrix *CCT;
gsl_matrix *AAT;
gsl_matrix *BBT;


gsl_rng *rng;

int main(void) {
    /* getting the starting time */
    struct timeval start_t, stop;
    double secs = 0.0,sec_total = 0.0;

    int q,i,ii,j,flag,iter;
    gsl_vector_view col,theta_block;
    gsl_matrix_view y_block;
    char line[200];
    int err,n_b;
    double log_prior_0, lgl_0, log_post_0, log_prior_1, lgl_1, log_post_1;
    double MH,accept;
    double tau,plgl,lgi, plgl_p;
    FILE *f_af,*f_pd_port,*f_plgl, *f_prob_p, *f_prob_bv, *f_pd, *f_sim;


    
    n_b=10000;
    T=30000;

    n=2304; //number of observations


    st=0;
    ed=n-1;

    k=2;
    npar= k*(k+1)/2 + 2*k +k +1 +k  ; // CCT AAT BBT eta psi mu
    num_block = 1;
    gsl_vector *af = gsl_vector_calloc(num_block);

    printf("npar = %d\n",npar);
    sample_cov = gsl_matrix_calloc(k,k);
    gsl_vector_int *theta_ind = gsl_vector_int_calloc(npar);
    gsl_vector *prop = gsl_vector_calloc(npar);
    gsl_matrix *chol_hessian = gsl_matrix_calloc(npar,npar);
    gsl_matrix *chol_hessian_saved = gsl_matrix_calloc(npar,npar);
    gsl_vector *theta = gsl_vector_calloc(npar);
    gsl_vector *theta_new = gsl_vector_calloc(npar);
    double f;
    gsl_matrix_view mat_v;
    gsl_vector_view row_y, mu_t,row_Y;
    Ht1 = gsl_matrix_calloc(k,k);
    Htt1 = gsl_matrix_calloc(k,k);
    Y = gsl_matrix_calloc(n,k);
    CCT = gsl_matrix_calloc(k,k);
    AAT = gsl_matrix_calloc(k,k);
    BBT = gsl_matrix_calloc(k,k);

    FILE *fdata = fopen("sptlt.dat", "r");
    gsl_vector_int_set_all( theta_ind, 1);

    /* read in data */
    for(i=0;i<n;i++){

        fgets(line,200,fdata);
        /* fx.dat */
        //  sscanf(line, "%s %lf %lf %lf",gsl_vector_char_ptr(date,i),mptr(y,i,0),mptr(y,i,1),mptr(y,i,2));
        /* crsp.dat */
        sscanf(line, "%lf %lf",mptr(Y,i,0),mptr(Y,i,1));

     //   mset( Y,i,0, log( mget(Y,i,0)+1.0));
     //   mset( Y,i,1, log( mget(Y,i,1)+1.0));

    }

        // pmat(Y);

    gsl_matrix_scale(Y,100.0);


    /* initialize RNG */
    gsl_rng_env_setup();
    rng = gsl_rng_alloc (gsl_rng_taus);
    gsl_rng_set (rng,(unsigned long int)894529599);
    printf ("generator type: %s\n", gsl_rng_name (rng));


    int start = 1200 , end = 2304 ;
	gsl_vector *pl = gsl_vector_calloc(end-start);
    gsl_vector *pl_p = gsl_vector_calloc(end-start);
    gsl_vector *DR = gsl_vector_calloc(end-start);
    // probabilities
    gsl_vector *mut1 = gsl_vector_calloc(k);
     gsl_matrix *vt1 = gsl_matrix_calloc(k,k);
    gsl_vector *omega =gsl_vector_calloc(k);
    gsl_vector *prob_p = gsl_vector_calloc(7);
    gsl_vector *vec_work = gsl_vector_calloc(k);
    gsl_matrix *prob_p_mat = gsl_matrix_calloc(end-start,7);
    gsl_matrix *prob_bv_mat = gsl_matrix_calloc(end-start,7);
    vset(omega,0,0.6); vset(omega,1,0.4);
        

    double corr;

    gsl_vector *prob_bv = gsl_vector_calloc(7);
    gsl_vector_view diag;
    gsl_matrix *D = gsl_matrix_calloc(k,k);
    gsl_matrix *Cor = gsl_matrix_calloc(k,k);
    gsl_matrix *Cor2 = gsl_matrix_calloc(k,k);

    gsl_vector *ord = gsl_vector_calloc(200);
    gsl_matrix *pd = gsl_matrix_calloc(end-start,200);
    gsl_vector_view pd_v;

     //VaR ES
    gsl_vector *draw = gsl_vector_calloc(k);
    gsl_vector *sim = gsl_vector_calloc(T);
    gsl_matrix *VaR = gsl_matrix_calloc(end-start,3);
    gsl_matrix *ES = gsl_matrix_calloc(end-start,3);
    gsl_matrix *LF = gsl_matrix_calloc(end-start,3);
    gsl_rng *rng1 = gsl_rng_alloc (gsl_rng_taus);

     // GMVP
     gsl_matrix *Ecov = gsl_matrix_calloc(k,k);
    gsl_vector *Emu = gsl_vector_calloc(k);
    gsl_vector *post_mut1 = gsl_vector_calloc(k);
    gsl_matrix *tmp_mat = gsl_matrix_calloc(k,k);
    gsl_vector *omega_GMVP = gsl_vector_calloc(k);
    gsl_vector *id = gsl_vector_calloc(k);
    gsl_vector_set_all(id,1.0);
    gsl_matrix *covt1 = gsl_matrix_calloc(k,k);
    double sc, rp_GMVP;
    gsl_matrix *rp_GMV = gsl_matrix_calloc(end-start,2);
    gsl_matrix *omega_GMV = gsl_matrix_calloc(end-start,k);


     //UO
   
    simRt = gsl_matrix_calloc(T,k);
    omega_t = gsl_vector_calloc(k);
    gsl_matrix *opt_w = gsl_matrix_calloc(end-start,3);
    gsl_matrix *opt_u = gsl_matrix_calloc(end-start,3);
    gsl_matrix *realized_u = gsl_matrix_calloc(end-start,3);
    gsl_matrix *simRtbound = gsl_matrix_calloc(end-start,2);



    for ( i = 0; i < 200; i++)
    {
        vset(ord,i,-15.0+i*30.0/199);
        
    }


    for ( n = start ; n < end; ++n) {
        plgl = 0.0;
        plgl_p = 0.0;
        gsl_vector_set_zero(prob_p);
        gsl_vector_set_zero(prob_bv);
        gsl_matrix_set_zero(covt1);
        gsl_vector_set_zero(post_mut1);
        row_Y = gsl_matrix_row(Y,n);
        double rp1;
        gsl_blas_ddot(omega,vv(row_Y),&rp1);
        st=0;
        ed=n-1;
        gettimeofday(&start_t, NULL);
        y = gsl_matrix_calloc(n,k);
        H = gsl_array_calloc(n,k,k);

        mat_v = gsl_matrix_submatrix(Y,0,0,n,k);
        gsl_matrix_memcpy(y, mv(mat_v));
        cov_matrix_col(y, sample_cov);
       // gsl_rng_set (rng,(unsigned long int)894529599);
        //gsl_rng_set (rng1,(unsigned long int)894529599); 

        /* initial values for parameters */
        gsl_vector_set_all( theta,0.0 );

        /* CCT */
        theta_block = gsl_vector_subvector ( theta, 0,k*(k+1)/2);
        gsl_vector_set_all( vv(theta_block), 0.01 );
        /* AAT */
        theta_block = gsl_vector_subvector ( theta, k*(k+1)/2,k);
        gsl_vector_set_all( vv(theta_block), 0.20 );
        /* BBT */
        theta_block = gsl_vector_subvector( theta, k*(k+1)/2+k,k);
        gsl_vector_set_all( vv(theta_block), 0.95 );

        vset(theta,npar-k-1,5.0);
    //    pvec(theta);

        err = log_prior( theta, &log_prior_0);
        if(err != 0){
            printf("Starting Theta invalid\n");
            exit(0);
        }
        err = lgl( theta, &lgl_0);
        if(err != 0){
            printf("Starting Theta invalid\n");
            exit(0);
        }

        log_post_0 = log_prior_0 + lgl_0;
         gsl_vector_int_set_all( theta_ind, 1);
         
        /* Get initial Hessian Estimate for Block RW */

        // pvec(theta);
        flag = mode_hessian( &lgl, &log_prior, npar , theta_ind , theta, theta_new, chol_hessian_saved );

        theta_new->data[0] = abs(theta_new->data[0]);
        theta_new->data[2] = abs(theta_new->data[2]);

        printf("flag=%d\n",flag);
        //  pvec(theta);ƒ
        // pvec(theta_new);
        flag = lgl( theta_new, &f);
        // update_H();
//    printf("flag = %d\n",flag);
//    printf("f = %g\n",f);
//    sleep(5);

        gsl_vector_memcpy( theta, theta_new);


        gsl_blas_dsyrk( CblasLower, CblasNoTrans, 1.0,chol_hessian_saved, 0.0,chol_hessian);
        copy_lower_to_upper( chol_hessian);
        gsl_matrix_memcpy(chol_hessian_saved,chol_hessian);

        inv_pd_sym( chol_hessian_saved );
        gsl_linalg_cholesky_decomp (chol_hessian_saved );
        
       

        err = log_prior( theta, &log_prior_0);
        if(err != 0){
            pvec(theta);
            printf("Starting Theta invalid:=%d\n",ii);
            exit(0);
        }
        err = lgl( theta, &lgl_0);
        if(err != 0){
            printf("Starting Theta invalid:=%d\n",ii);
            exit(0);
        }

        log_post_0 = log_prior_0 + lgl_0;

        for(iter=0;iter<n_b+T;iter++){

	    gsl_matrix_memcpy( chol_hessian, chol_hessian_saved );

            tau=0.015; //-((n-start+1)/(end-start))*0.045;

            gsl_matrix_scale( chol_hessian, sqrt(tau) );
           

            if( gsl_rng_uniform(rng) > 0.1){
                ran_mvn_chol( theta , chol_hessian , prop);
            }else{
                gsl_matrix_scale( chol_hessian, 10.0 );
                ran_mvn_chol( theta , chol_hessian , prop);
            }

            //  ran_mvn_chol_precision(theta,chol_hessian,1,prop);


            accept=0.0;
            err = log_prior( prop, &log_prior_1);
            //      printf("err=%d\n",err);
            if( err == 0){
                err = lgl( prop, &lgl_1);
                if( err == 0){
                    /* get MH ratio */
                    log_post_1 = log_prior_1 + lgl_1;

                    MH = log_post_1 - log_post_0;

                    if( MH > 0.0 || gsl_rng_uniform(rng) < exp( MH ) ){
                        /* Accept */
                        gsl_vector_memcpy( theta, prop);
                        log_prior_0 = log_prior_1;
                        lgl_0 = lgl_1;
                        log_post_0 = log_post_1;
                        accept = 1.0;
                        gsl_matrix_memcpy(Ht1,Htt1);
                    }
                }
            }
                //  printf("accept=%g\n",accept);


                

            af->data[0] += accept/1000.0;



            if( iter % 1000 == 0){
                printf("iter:%d",iter);
                pvec(af);
                //fprintf(f_af,"%g ",vget(af,0));
                gsl_vector_set_zero( af );
            }


            

            if( iter >= n_b){
                j = iter - n_b;
                mu_t = gsl_vector_subvector(theta,npar-k,k);
                row_y = gsl_matrix_row(Y,n);
                double df = vget(theta,npar-k-1);
               // pvec(vv(mu_t));
                err = log_mvt_den_custom(vv(row_y),vv(mu_t),df,Ht1,&lgi);
                if (err == 0 ){
                    lgi = exp(lgi);
//                    printf("i:%d lgi:%g\n",iter,lgi);
                    plgl += 1.0/T *lgi;
//                    printf("i:%d lgi:%g\n",iter,plgl);
                }else{
                    printf("error!\n");
                }
                  // probabilities
                gsl_vector_memcpy(mut1,vv(mu_t));
                gsl_matrix_memcpy(vt1,Ht1);
                double mu_p; double sigma2_p;
                
                gsl_blas_ddot(omega,mut1,&mu_p);
                gsl_blas_dsymv(CblasLower,1.0,vt1,omega,0.0,vec_work);
                gsl_blas_ddot(vec_work,omega,&sigma2_p);
                double sigma_p = sqrt(sigma2_p);
                plgl_p += 1.0/T * den_st(rp1,mu_p,sigma2_p,df);
                double inv_stdev = 1.0/sqrt(sigma2_p);
                
                
                prob_p->data[0] += 1.0/T * gsl_cdf_tdist_P((-.0-mu_p)*inv_stdev,df);
                prob_p->data[1] += 1.0/T * gsl_cdf_tdist_P((-1.0-mu_p)*inv_stdev,df);
                prob_p->data[2] += 1.0/T * gsl_cdf_tdist_P((-2.0-mu_p)*inv_stdev,df);
                prob_p->data[3] += 1.0/T * gsl_cdf_tdist_P((-3.0-mu_p)*inv_stdev,df);
                prob_p->data[4] += 1.0/T * gsl_cdf_tdist_P((-4.0-mu_p)*inv_stdev,df);
                prob_p->data[5] += 1.0/T * gsl_cdf_tdist_P((-5.0-mu_p)*inv_stdev,df);
                prob_p->data[6] += 1.0/T * gsl_cdf_tdist_P((-10.0-mu_p)*inv_stdev,df);
                
                 // GMVP
            
                gsl_matrix_memcpy(Ecov,vt1);
                gsl_blas_dger(1.0,mut1,mut1,Ecov);
                gsl_matrix_scale(Ecov,1.0/T);
                gsl_matrix_add(covt1,Ecov);
                gsl_vector_memcpy(Emu, mut1);
                gsl_vector_scale(Emu,1.0/T);
                gsl_vector_add(post_mut1,Emu);

             
            
                // predictive density
                for ( i = 0; i < 200; i++)
                {
                    pd->data[(n-start) * pd->tda + i] += 1.0/T * den_st(vget(ord,i),mu_p,sigma2_p,df);
                }
                
                   // probabilities - BV dist
                   int nu = (int) df;
                mset(D,0,0,1/sqrt(mget(vt1,0,0)));
                mset(D,1,1,1/sqrt(mget(vt1,1,1)));
                gsl_blas_dsymm(CblasLeft,CblasLower,1.0,D,vt1,0.0,Cor);
                
                gsl_blas_dsymm(CblasLeft,CblasLower,1.0,Cor,D,0.0,Cor2);
                corr = mget(Cor2,1,0);
                double lower[] = {(-2.0-vget(mut1,0))*mget(D,0,0),(-2.0-vget(mut1,1))*mget(D,1,1)};
                double upper[] = {(2.0-vget(mut1,0))*mget(D,0,0),(2.0-vget(mut1,1))*mget(D,1,1)};
                int infi []= {2,2};
                //printf("l: %g, %g u: %g,%g corr %g func: %g\n",lower[0],lower[1],upper[0],upper[1],corr,bvnmvn_(lower,upper,infi,&corr));

                prob_bv ->data[0] += 1.0/T * (1-mvbvt_(&nu,lower,upper,infi,&corr));


                 lower[0] = (-3.0-vget(mut1,0))*mget(D,0,0); lower[1] = (-2.5-vget(mut1,1))*mget(D,1,1);
                 upper[0] = (3.0-vget(mut1,0))*mget(D,0,0); upper[1] = (2.5-vget(mut1,1))*mget(D,1,1);
                
                //printf("l: %g, %g u: %g,%g corr %g func: %g\n",lower[0],lower[1],upper[0],upper[1],corr,bvnmvn_(lower,upper,infi,&corr));

                prob_bv ->data[1] += 1.0/T * (1- mvbvt_(&nu,lower,upper,infi,&corr));

                lower[0] = (-4.0-vget(mut1,0))*mget(D,0,0); lower[1] = (-3.0-vget(mut1,1))*mget(D,1,1);
                upper[0] = (4.0-vget(mut1,0))*mget(D,0,0); upper[1] = (3.3-vget(mut1,1))*mget(D,1,1);
                
                //printf("l: %g, %g u: %g,%g corr %g func: %g\n",lower[0],lower[1],upper[0],upper[1],corr,bvnmvn_(lower,upper,infi,&corr));

                prob_bv ->data[2] += 1.0/T * (1-mvbvt_(&nu,lower,upper,infi,&corr));


                lower[0] = (-5.0-vget(mut1,0))*mget(D,0,0); lower[1] = (-3.5-vget(mut1,1))*mget(D,1,1);
                upper[0] = (5.0-vget(mut1,0))*mget(D,0,0); upper[1] = (3.5-vget(mut1,1))*mget(D,1,1);
               
                //printf("l: %g, %g u: %g,%g corr %g func: %g\n",lower[0],lower[1],upper[0],upper[1],corr,bvnmvn_(lower,upper,infi,&corr));

                prob_bv ->data[3] += 1.0/T * (1-mvbvt_(&nu,lower,upper,infi,&corr));


                lower[0] = (-10.0-vget(mut1,0))*mget(D,0,0); lower[1] = (-5.0-vget(mut1,1))*mget(D,1,1);
                upper[0] = (10.0-vget(mut1,0))*mget(D,0,0); upper[1] = (5.0-vget(mut1,1))*mget(D,1,1);

                //printf("l: %g, %g u: %g,%g corr %g func: %g\n",lower[0],lower[1],upper[0],upper[1],corr,bvnmvn_(lower,upper,infi,&corr));

                prob_bv ->data[4] += 1.0/T * (1-mvbvt_(&nu,lower,upper,infi,&corr));

                lower[0] = (-0.0-vget(mut1,0))*mget(D,0,0); lower[1] = (-0.0-vget(mut1,1))*mget(D,1,1);
                upper[0] = (-5.0-vget(mut1,0))*mget(D,0,0); upper[1] = (100.0-vget(mut1,1))*mget(D,1,1);
                infi[0] = 0 ; infi[1]= 0;
                //printf("l: %g, %g u: %g,%g corr %g func: %g\n",lower[0],lower[1],upper[0],upper[1],corr,bvnmvn_(lower,upper,infi,&corr));

                prob_bv ->data[5] += 1.0/T * (mvbvt_(&nu,lower,upper,infi,&corr));

                lower[0] = (-0.0-vget(mut1,0))*mget(D,0,0); lower[1] = (-0.0-vget(mut1,1))*mget(D,1,1);
                upper[0] = (-10.0-vget(mut1,0))*mget(D,0,0); upper[1] = (100.0-vget(mut1,1))*mget(D,1,1);
                infi[0] = 0 ; infi[1]= 0;
                //printf("l: %g, %g u: %g,%g corr %g func: %g\n",lower[0],lower[1],upper[0],upper[1],corr,bvnmvn_(lower,upper,infi,&corr));

                prob_bv ->data[6] += 1.0/T * (mvbvt_(&nu,lower,upper,infi,&corr));

                // Value at Risk - Expected Shortfall for portfolio
                //ran_mvn(mut1,vt1,draw); 
               // gsl_linalg_cholesky_decomp(vt1);
                ran_mvt(mut1,vt1,df,draw);
                gsl_matrix_set_row(simRt,j,draw);
                gsl_blas_ddot(omega,draw,&sim->data[j]);
            }

            
        }


        mset(opt_w,n-start,0,golden_section_search(utility,2.0, 0.0, 1.0, EPSILON));
        mset(opt_w,n-start,1,golden_section_search(utility,4.0, 0.0, 1.0, EPSILON));
        mset(opt_w,n-start,2,golden_section_search(utility,6.0, 0.0, 1.0, EPSILON));
        
        mset(opt_u,n-start,0,utility(mget(opt_w,n-start,0),2.0));
        mset(opt_u,n-start,1,utility(mget(opt_w,n-start,1),4.0));
        mset(opt_u,n-start,2,utility(mget(opt_w,n-start,2),6.0));


        mset(realized_u,n-start,0,R_utility(mget(opt_w,n-start,0),2.0,vv(row_y)));
        mset(realized_u,n-start,1,R_utility(mget(opt_w,n-start,1),4.0,vv(row_y)));
        mset(realized_u,n-start,2,R_utility(mget(opt_w,n-start,2),6.0,vv(row_y)));
        
        mset(simRtbound,n-start,0,gsl_matrix_min(simRt));
        mset(simRtbound,n-start,1,gsl_matrix_max(simRt));
        
        //utility_EU(0.75,2.0);
        
        f_sim = fopen("UO.out","w");
      
        for (int l = 0; l < n-start+1; ++l) {
            fprintf(f_sim,"%d  %g %g %g %g %g %g %g %g %g %g %g\n",l+start,mget(opt_w,l,0),mget(opt_w,l,1),mget(opt_w,l,2)
                                                                          ,mget(opt_u,l,0),mget(opt_u,l,1),mget(opt_u,l,2)
                                                                          ,mget(realized_u,l,0),mget(realized_u,l,1),mget(realized_u,l,2)
                                                                          ,mget(simRtbound,l,0),mget(simRtbound,l,1));
        }
        fclose(f_sim);

        // Downside Risk
        
        double count = 0.0;
        for (size_t zz = 0; zz < T; zz++)
        {
            if (mget(simRt,zz,0) > 0.0 && mget(simRt,zz,1)> 0.0) count ++;
        }
        vset (DR,n-start,count/T);
        f_sim = fopen("DR.out","w");
        for (int l = 0; l < n-start+1; l++)
        {
            fprintf(f_sim,"%g\n",vget(DR,l));
        }
        fclose(f_sim);



         char filename[50];
         snprintf(filename, sizeof(filename), "sim%d.out", n);

            // Open the file for writing
         f_sim = fopen(filename, "w");
         if (!f_sim) {
             fprintf(stderr, "Error: Unable to open file %s for writing\n", filename);
             return;
          }

            // Write the contents of simRt to the file
         
         for (size_t iii = 0; iii < T; iii++) {
                for (size_t jjj = 0; jjj < k; jjj++) {
                 fprintf(f_sim, "%.6f ", mget(simRt, iii, jjj));
             }
             fprintf(f_sim, "\n");
         }

         // Close the file
            fclose(f_sim);
            printf("Data saved to %s\n", filename);




        plgl = log(plgl);
        plgl_p = log(plgl_p);
        printf("n=%d, plgl = %g \n",n,plgl);
        vset(pl,n-start,plgl);
        f_plgl = fopen("Plgl.out","w");
        for (int l = 0; l < n-start+1; ++l) {
            fprintf(f_plgl,"%d  %g\n",l+start,vget(pl,l));
        }
        fclose(f_plgl);

         vset(pl_p,n-start,plgl_p);
        f_plgl = fopen("Plgl_p.out","w");
        for (int l = 0; l < n-start+1; ++l) {
            fprintf(f_plgl,"%d  %g\n",l+start,vget(pl_p,l));
        }
        fclose(f_plgl);

        printf("prob <0:%g,<1:%g,<2:%g,<3:%g,<4:%g\n",vget(prob_p,0),vget(prob_p,1),vget(prob_p,2),vget(prob_p,3),vget(prob_p,4));
        gsl_matrix_set_row(prob_p_mat,n-start,prob_p);
        f_prob_p = fopen("prob_p.out","w");
        for (int l = 0; l < n-start+1; ++l) {
            fprintf(f_prob_p,"%d  %g %g %g %g %g %g %g\n",l+start,mget(prob_p_mat,l,0),mget(prob_p_mat,l,1),mget(prob_p_mat,l,2),mget(prob_p_mat,l,3),mget(prob_p_mat,l,4),mget(prob_p_mat,l,5),mget(prob_p_mat,l,6));
        }
        fclose(f_prob_p);
         printf("prob <0:%g,<1:%g,<2:%g,<3:%g,<4:%g\n",vget(prob_bv,0),vget(prob_bv,1),vget(prob_bv,2),vget(prob_bv,3),vget(prob_bv,4));
        gsl_matrix_set_row(prob_bv_mat,n-start,prob_bv);
        f_prob_bv = fopen("prob_bv.out","w");
        for (int l = 0; l < n-start+1; ++l) {
            fprintf(f_prob_bv,"%d  %g %g %g %g %g %g %g\n",l+start,mget(prob_bv_mat,l,0),mget(prob_bv_mat,l,1),mget(prob_bv_mat,l,2),mget(prob_bv_mat,l,3),mget(prob_bv_mat,l,4),mget(prob_bv_mat,l,5),mget(prob_bv_mat,l,6));
        }
        fclose(f_prob_bv);
        f_pd = fopen("pd.out", "w");
        for (int l = 0; l < n-start+1; ++l) {
            for ( i = 0; i < 200; i++)
            {
                fprintf(f_pd, "%g %g\n",vget(ord,i),mget(pd,l,i) );
            }
            fprintf(f_pd,"\n\n");
            
        }
        fclose(f_pd);
        //VaR-ES-LF
        gsl_sort_vector(sim);
        // gsl_vector_scale(sim,0.01);
        int q1 = (int)(floor(T*1.0/100));
        int q5 = (int)(floor(T*5.0/100));
        int q10 = (int)(floor(T*10.0/100));
        mset(VaR,n-start,0,vget(sim,q10-1));
        mset(VaR,n-start,1,vget(sim,q5-1));
        mset(VaR,n-start,2,vget(sim,q1-1));
      
        double sum = 0;
        for ( i = 0; i< q10 ; i++)
        {
            sum += vget(sim,i);
            if (i == q1-1) mset(ES,n-start,2,sum/q1);
            if (i == q5-1) mset(ES,n-start,1,sum/q5);
        }
        
        mset(ES,n-start,0,sum/q10);

       
        // rp1 = rp1 * 0.01;
        printf("realized return for p : %g\n",rp1);
       
        mset(LF,n-start,0,LFun(rp1,mget(VaR,n-start,0),mget(ES,n-start,0),0.1));
        mset(LF,n-start,1,LFun(rp1,mget(VaR,n-start,1),mget(ES,n-start,1),0.05));
        mset(LF,n-start,2,LFun(rp1,mget(VaR,n-start,2),mget(ES,n-start,2),0.01));

        f_pd = fopen("VaR.out","w");
        for (int l = 0; l < n-start+1; ++l) {
            fprintf(f_pd,"%d  %g %g %g\n",l+start,mget(VaR,l,0),mget(VaR,l,1),mget(VaR,l,2));
        }
        fclose(f_pd);
        
        f_pd = fopen("ES.out","w");
        for (int l = 0; l < n-start+1; ++l) {
            fprintf(f_pd,"%d  %g %g %g\n",l+start,mget(ES,l,0),mget(ES,l,1),mget(ES,l,2));
        }
        fclose(f_pd);

        f_pd = fopen("LF.out","w");
        for (int l = 0; l < n-start+1; ++l) {
            fprintf(f_pd,"%d  %g %g %g\n",l+start,mget(LF,l,0),mget(LF,l,1),mget(LF,l,2));
        }
        fclose(f_pd);

        //GMVP
        gsl_matrix_set_zero(tmp_mat);
        gsl_blas_dger(1.0,post_mut1,post_mut1,tmp_mat);
        gsl_matrix_sub(covt1,tmp_mat);
        
        inv_pd_sym(covt1);
        gsl_blas_dsymv(CblasLower,1.0,covt1,id,0.0,vec_work);
        gsl_blas_ddot(id,vec_work,&sc);
      
       
        gsl_blas_dsymv(CblasLower,1.0,covt1,id,0.0,omega_GMVP);
        gsl_vector_scale(omega_GMVP,1.0/sc);
        
        gsl_blas_ddot(omega_GMVP,vv(row_Y),&rp_GMVP);
        pvec(omega_GMVP);
        printf("rp_GMVP:%g\n",rp_GMVP);
        mset(rp_GMV,n-start,0,rp1);
        mset(rp_GMV,n-start,1,rp_GMVP);
        
        f_pd = fopen("rp.out","w");
        for (int l = 0; l < n-start+1; ++l) {
            fprintf(f_pd,"%d  %g %g\n",l+start,mget(rp_GMV,l,0),mget(rp_GMV,l,1));
        }
        fclose(f_pd);
        gsl_matrix_set_row(omega_GMV,n-start,omega_GMVP);
        f_pd = fopen("omega.out","w");
        for (int l = 0; l < n-start+1; ++l) {
            fprintf(f_pd,"%d  %g %g\n",l+start,mget(omega_GMV,l,0),mget(omega_GMV,l,1));
        }
        fclose(f_pd);
        gsl_array_free(H);
        gsl_matrix_free(y);


        gettimeofday(&stop, NULL);
        secs = (double)(stop.tv_usec - start_t.tv_usec) / 1000000 + (double)(stop.tv_sec - start_t.tv_sec);
        sec_total += secs;
        printf("time taken for iteration %d %f secs\n",n, secs);
    }

    printf("Total time: %f mins\n",sec_total/60);

    return 0;
}
