extern int mvtdst_(integer *n, integer *nu, doublereal *lower, 
	doublereal *upper, integer *infin, doublereal *correl, doublereal *
	delta, integer *maxpts, doublereal *abseps, doublereal *releps, 
	doublereal *error, doublereal *value, integer *inform__);
extern int mvsubr_0_(int n__, integer *n, doublereal *w, integer *
	nf, doublereal *f, integer *nuin, doublereal *correl, doublereal *
	lower, doublereal *upper, doublereal *delta, integer *infin, integer *
	nd, doublereal *vl, doublereal *er, integer *inform__);
extern int mvsubr_(integer *n, doublereal *w, integer *nf, 
	doublereal *f);
extern int mvsubr_(integer *n, doublereal *w, integer *nf, 
	doublereal *f);
extern int mvints_(integer *n, integer *nuin, doublereal *correl, 
	doublereal *lower, doublereal *upper, doublereal *delta, integer *
	infin, integer *nd, doublereal *vl, doublereal *er, integer *inform__);
extern  int mvspcl_(integer *nd, integer *nu, doublereal *a, 
	doublereal *b, doublereal *dl, doublereal *cov, integer *infi, 
	doublereal *snu, doublereal *vl, doublereal *er, integer *inform__);
extern int mvvlsb_(integer *n, doublereal *w, doublereal *r__, 
	doublereal *dl, integer *infi, doublereal *a, doublereal *b, 
	doublereal *cov, doublereal *y, doublereal *di, doublereal *ei, 
	integer *nd, doublereal *value);
extern int mvsort_(integer *n, doublereal *lower, doublereal *upper,
	 doublereal *delta, doublereal *correl, integer *infin, doublereal *y,
	 logical *pivot, integer *nd, doublereal *a, doublereal *b, 
	doublereal *dl, doublereal *cov, integer *infi, integer *inform__);
extern doublereal mvtdns_(integer *nu, doublereal *x);
extern int mvlims_(doublereal *a, doublereal *b, integer *infin, 
	doublereal *lower, doublereal *upper);
extern int mvsswp_(doublereal *x, doublereal *y);
extern int mvswap_(integer *p, integer *q, doublereal *a, 
	doublereal *b, doublereal *d__, integer *infin, integer *n, 
	doublereal *c__);
extern mvphi_(doublereal *z__);
extern mvphnv_(doublereal *p);
extern  mvbvn_(doublereal *lower, doublereal *upper, integer *infin, 
	doublereal *correl);
extern doublereal mvbvu_(doublereal *sh, doublereal *sk, doublereal *r__);
extern doublereal mvstdt_(integer *nu, doublereal *t);
extern doublereal mvbvt_(integer *nu, doublereal *lower, doublereal *upper, integer *
	infin, doublereal *correl);
extern mvbvtc_(integer *nu, doublereal *l, doublereal *u, integer *infin, 
	doublereal *rho);
extern doublereal mvbvtl_(integer *nu, doublereal *dh, doublereal *dk, doublereal *
	r__);
extern doublereal mvbvtl_(integer *nu, doublereal *dh, doublereal *dk, doublereal *
	r__);
extern doublereal mvchnv_(integer *n, doublereal *p);
extern doublereal mvchnc_(doublereal *lkn, integer *n, doublereal *p, doublereal *
	r__);
extern doublereal mvchnv_(integer *n, doublereal *p);
extern doublereal mvchnc_(doublereal *lkn, integer *n, doublereal *p, doublereal *
	r__);
extern int mvkbrv_(integer *ndim, integer *minvls, integer *maxvls, 
	integer *nf, U_fp funsub, doublereal *abseps, doublereal *releps, 
	doublereal *abserr, doublereal *finest, integer *inform__);
extern int mvkrsv_(integer *ndim, integer *kl, doublereal *values, 
	integer *prime, doublereal *vk, integer *nf, S_fp funsub, doublereal *
	x, doublereal *r__, integer *pr, doublereal *fs);
extern doublereal mvuni_(void);