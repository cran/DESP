/*
 *  DESP/src/srL_B.c by A. S. Dalalyan and S. Balmand  Copyright (C) 2015-
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License (version 3) as published by
 *  the Free Software Foundation.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  A copy of the GNU General Public License is available at
 *  http://www.r-project.org/Licenses/
 */

#include <R.h>
#include <Rinternals.h>
#include "scsSqR_Lasso_solve.h"
#include "cdSqR_Lasso_solve.h"
#include "R_ext/BLAS.h"


#ifdef _OPENMP
#include <omp.h>
#endif



SEXP srL_B(SEXP X, SEXP theta, SEXP lambda, SEXP sqRmethod, SEXP sto, SEXP nthreads){
	/* estimation of B */
	/* 
	X : data matrix
	theta : matrix of outliers, if NULL, the computation is similar to theta=0
	lambda : tuning parameter that promotes sparsity
	sqRmethod : choice of the method to solve the square-root Lasso problem ; 0 : coordinate descent ; 1 : solve SOCP using SCS
	sto : when using coordinate descent, whether a shuffled (2) or stochastic (1) coordinate descent must be performed or not (0)
	nthreads : number of parallel threads to be used
	*/

	SEXP B, dim, status, solution, solNames, t;
	double * beta, * z, * XX, * tt, * TTheta, * BB;
	int j, st, st_t, detect, nt, meth, stoc;
	int * stat;
	double lbda;

	double done = 1.0;
	int one = 1;
	double mone = -1.0;
	int nb;
  
	/* read the sample size and the number of variables */
	dim = getAttrib(X,R_DimSymbol);
	int n = INTEGER(dim)[0];               /* n is the sample size */
	int p = INTEGER(dim)[1];               /* p is the dimension */

	/* check whether the robust version has to be used or not */
	detect = 0;
	if(!isNull(theta)){ detect = 1;}

	PROTECT(solution = allocVector(VECSXP, 3));
	PROTECT(status = allocVector(INTSXP, 1));
	stat = INTEGER(status);

	/* the following lines are used to prevent false sharing. The columns of matrix BBB are aligned accordingly to cache line size (for a 64 bits processor). This matrix is used as a shared variable in replacement of BB in the parallel section. Doing this, we use more memory. In practice, it is just a little bit faster (not allways) when n becomes large. Note that it is interesting to create more threads than those available when n is not too large. */
	int colSize = p;
	if(p % (64/sizeof(double)) != 0){ colSize = (p/(64/sizeof(double)) + 1)*(64/sizeof(double)); }
	double * BBB;
	BBB = (double *) Calloc((size_t) colSize*p, double);
	memset(BBB, 0, sizeof(double)*colSize*p);

	PROTECT(t = allocVector(REALSXP,p));

	XX = REAL(X);
	if(detect){ 
	  TTheta = REAL(theta);
  }else{
    TTheta = NULL;
  }
	tt = REAL(t);
	lbda = (double) *(REAL(lambda));
	meth = (int) *(INTEGER(sqRmethod));
	if(meth == 0){ // use coordinate descent
		stoc = (int) *(INTEGER(sto));
    if (stoc!=0){
      GetRNGstate();
    }
	}
	nt = (int) *(INTEGER(nthreads));

	/* initialize status, 0 means success, 1 means failure */
	stat[0] = 0;
	st_t = 0;
 
	/* compute an estimator of each column by the square-root Lasso */
	#ifdef _OPENMP
	#pragma omp parallel if(p > 10) default(none) shared(XX, TTheta, BBB, colSize, tt, n, p, lbda, one, mone, done, detect, meth, stoc) private(j, beta, z, nb, st) reduction(+:st_t) num_threads(nt)
	#endif
	{ /* begin parallel */

	/* memory allocation for private threads'variables */
	if (meth == 1){// use SCS : conic
		beta = (double *) Calloc((size_t) p-1, double);
	}else{// use coordinate descent
		beta = (double *) Calloc((size_t) p, double);
	}
	z = (double *) Calloc((size_t) n, double);
	
	st = 0;
	nb = 0;

	#ifdef _OPENMP
	#pragma omp for nowait 
	#endif
	for(j=0; j<p; j++){
		if (meth == 1){// use SCS : conic
			st = scs_sqR_Lasso_solve_loop(XX, TTheta, lbda, n, p, beta, j);
			if (st != SCS_SOLVED){
				st_t += 1;
			}
			memcpy(BBB+(colSize*j), beta, j*sizeof(double));
			BBB[j+colSize*j] = 1;
			nb = p-j-1;

			memcpy(BBB+(colSize*j + j+1), beta+j, nb*sizeof(double));
		}else{// use coordinate descent
			st = cd_sqR_Lasso_solve_loop(XX, TTheta, lbda, n, p, beta, j, stoc);
			if (st != 0){
				st_t += 1;
			}
			memcpy(BBB+(colSize*j), beta, p*sizeof(double));
			BBB[j+colSize*j] = 1;
		}
		if(detect){
			memcpy(z, TTheta+(j*n), n*sizeof(double));

			/* we compute z = X*B[,j]-Theta[,j] */
			F77_CALL(dgemv)("N", &n, &p, &done, XX, &n, BBB+(j*colSize), &one, &mone, z, &one);
			tt[j] = F77_CALL(dnrm2)(&n, z, &one);
		}
	}

	/* free memory */
	if(beta) Free(beta);
	if(z) Free(z);

	} /* end parallel */

	if(st_t > 0){
		stat[0] = 1;
	}

	/* initialize the matrix B */
	PROTECT(B = allocMatrix(REALSXP,p,p));
	BB = REAL(B);
	for(j=0; j<p; j++){
		memcpy(BB+(p*j), BBB+(colSize*j), p*sizeof(double));
	}
	if(BBB) Free(BBB);

	SET_VECTOR_ELT(solution,0,B);
	SET_VECTOR_ELT(solution,1,status);
	SET_VECTOR_ELT(solution,2,t);
	
	PROTECT(solNames = allocVector(STRSXP,3));
	SET_STRING_ELT(solNames,0,mkChar("B"));
	SET_STRING_ELT(solNames,1,mkChar("status"));
	SET_STRING_ELT(solNames,2,mkChar("t"));
	setAttrib(solution,R_NamesSymbol,solNames);

	if(meth == 0){ // use coordinate descent
    if (stoc!=0){
      PutRNGstate();
    }
	}

	UNPROTECT(5);
	return solution;
}




