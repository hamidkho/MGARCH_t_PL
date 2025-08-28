int mvndst_(integer *n, doublereal *lower, doublereal *upper,
	 integer *infin, doublereal *correl, integer *maxpts, doublereal *
	abseps, doublereal *releps, doublereal *error, doublereal *value, 
	integer *inform__);
	
int mvnlms_(doublereal *a, doublereal *b, integer *infin, 
	doublereal *lower, doublereal *upper);
	
 int covsrt_(integer *n, doublereal *lower, doublereal *upper,
	 doublereal *correl, integer *infin, doublereal *y, integer *infis, 
	doublereal *a, doublereal *b, doublereal *cov, integer *infi);
	
int dkswap_(doublereal *x, doublereal *y);

int rcswp_(integer *p, integer *q, doublereal *a, doublereal 
	*b, integer *infin, integer *n, doublereal *c__);

int dkbvrc_(integer *ndim, integer *minvls, integer *maxvls, 
	D_fp functn, doublereal *abseps, doublereal *releps, doublereal *
	abserr, doublereal *finest, integer *inform__);

int dksmrc_(integer *ndim, integer *klim, doublereal *sumkro,
	 integer *prime, doublereal *vk, D_fp functn, doublereal *x);

doublereal bvnmvn_(doublereal *lower, doublereal *upper, integer *infin, 
	doublereal *correl);
	 

	

