/*
 *  DESP/src/cdSqR_Lasso.c by A. S. Dalalyan and S. Balmand  Copyright (C) 2015-
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
#include "cdSqR_Lasso_solve.h"
#include "R_ext/BLAS.h"

#define DEFAULT_max_iters	5000	/* maximum iterations to take: 5000 */
#define DEFAULT_tol		1e-6	/* convergence tolerance: 1e-6 */
#define DEFAULT_tolGap		1e-8	/* duality gap tolerance: 1e-8 */

SEXP cd_sqR_Lasso(SEXP X, SEXP Y, SEXP lambda, SEXP sto){
	/* solves the square-root Lasso by coordinate descent */
  /* the original algorithm is due to Alexandre Belloni, Victor Chernozhukov and Lie Wang */
	/* 
	X : matrix of the explanatory variables
	Y : vector corresponding to the response variable
	lambda : tuning parameter that promotes sparsity
	sto : whether a stochastic coordinate descent must be performed or not
	*/

	double * XX, * YY, * bbeta, * bPrec, * e_j, * Xe_j, * e;
	int * stat;
	SEXP dim, status, beta, solution, solNames;
	int i, j, k, l, stop, stoc;
	double Q_j, xe_j, mb_j, ab, abp, srss, Ye_j, primObj, dualObj, lbda;

	/* read the sample size and the number of variables */
	dim = getAttrib(X,R_DimSymbol);
	int n = INTEGER(dim)[0];               /* n is the sample size */
	int p = INTEGER(dim)[1];               /* p is the dimension */

	double dzero = 0.0;
	double done = 1.0;
	int one = 1;

	int maxIter = (int) DEFAULT_max_iters;
	double tol = (double) DEFAULT_tol;
	double tolGap = (double) DEFAULT_tolGap;

	PROTECT(solution = allocVector(VECSXP, 2));
	PROTECT(beta = allocVector(REALSXP, p));
	PROTECT(status = allocVector(INTSXP, 1));
	bbeta = REAL(beta);
	stat = INTEGER(status);

	XX = REAL(X);
	YY = REAL(Y);
	lbda = (double) *(REAL(lambda))*sqrt(n);
	stoc = (int) *(INTEGER(sto));

	/* initialize beta to zero */
	memset(bbeta, 0, p*sizeof(double));

	if(stoc!=0){
    GetRNGstate();
	}

	/* for stochastic coordinate descent */
	if(stoc==1){
		maxIter *= p;
	}

	/* for shuffled coordinate descent */
	int * s = 0;
	if(stoc==2){
		s = (int *) Calloc((unsigned) p, int);
	}

	bPrec = (double *) Calloc((unsigned) p, double);
	Xe_j = (double *) Calloc((unsigned) p, double);

	e = (double *) Calloc((unsigned) p, double);
	e_j = (double *) Calloc((unsigned) n, double);
	memset(e, 0, p*sizeof(double));
	for(j=0; j<p; j++){
		for(i=0; i<n; i++){
			e[j] += pow(XX[i+j*n],2);
		}
		e[j] /= n;
	}
	
	memcpy(e_j, YY, n*sizeof(double));
	Q_j = F77_CALL(ddot)(&n,e_j,&one,e_j,&one);
	Q_j /= n;

	stop = 0;
	k = 1;
	while(stop != 1 && k <= maxIter){
		memcpy(bPrec, bbeta, p*sizeof(double));

		ab = 0; // sum of the absolute values of bbeta
		abp = 0; // sum of the absolute values of bbeta - bPrec
		srss = 0; // residual sum of squares

		/* for shuffled coordinate descent */
		if(stoc==2){
			sample(p, s);
		}
		l = 0;
		while(l<p){ 
			/* for stochastic coordinate descent */
			if(stoc==1 && k!=1){
				// stochastic
				j = (int) floor(unif_rand() * p);
				l = p;
			}else if(stoc==2){
				// shuffle
				j = s[l];
			}else{
				j = l;
			}
			if (bbeta[j] != 0){
				F77_CALL(daxpy)(&n, bbeta+j, XX+n*j, &one, e_j, &one);
				Q_j = F77_CALL(ddot)(&n,e_j,&one,e_j,&one);
				Q_j /= n;
			}
			xe_j = F77_CALL(ddot)(&n,XX+n*j,&one,e_j,&one);
			xe_j /= n;

			if (pow(n,2) < pow(lbda,2) / e[j]){
				bbeta[j] = 0;
        		}else if (xe_j > lbda/n * sqrt(Q_j)){
               		 	bbeta[j] = (xe_j - lbda*sqrt((Q_j-pow(xe_j,2)/e[j])/(pow(n,2)-(pow(lbda,2) / e[j]))))/e[j];
				mb_j = -bbeta[j];
				F77_CALL(daxpy)(&n, &mb_j, XX+n*j, &one, e_j, &one);
        		}else if(xe_j < - lbda/n * sqrt(Q_j)){
                		bbeta[j] = (xe_j + lbda*sqrt((Q_j-pow(xe_j,2)/e[j])/(pow(n,2)-(pow(lbda,2) / e[j]))))/e[j];
				mb_j = -bbeta[j];
				F77_CALL(daxpy)(&n, &mb_j, XX+n*j, &one, e_j, &one);
        		}else{
               	 		bbeta[j] = 0;
        		}
			l++;
		}
		for(j=0; j<p; j++){
			ab += fabs(bbeta[j]);
			abp += fabs(bbeta[j]-bPrec[j]);
		}

		/* primal objective */
		Q_j = F77_CALL(ddot)(&n,e_j,&one,e_j,&one);
		Q_j /= n;
		primObj = sqrt(Q_j) + lbda/n * ab;
		
		
		srss = F77_CALL(ddot)(&n,e_j,&one,e_j,&one);
		srss = sqrt(srss); // squared Residual Sum of Squares
		if (srss > 1e-10){
			Ye_j = F77_CALL(ddot)(&n,YY,&one,e_j,&one);
			F77_CALL(dgemv)("T", &n, &p, &done, XX, &n, e_j, &one, &dzero, Xe_j, &one);
			dualObj = 0;
			for(j=0; j<p; j++){
				dualObj -= fabs(lbda/n - fabs(Xe_j[j]/(sqrt(n)*srss))) * fabs(bbeta[j]);
			}
			dualObj += Ye_j/(sqrt(n)*srss);
		}else{
			dualObj = lbda/n * ab;
		}
       
		/* stopping criterion */
		if (abp < tol){
			if (primObj - dualObj < tolGap){
				stop = 1;
			}
		}
		k++;
	}

	stat[0] = 0;
	if (k == maxIter) stat[0] = -1;

	SET_VECTOR_ELT(solution,0,beta);
	SET_VECTOR_ELT(solution,1,status);
	
	PROTECT(solNames = allocVector(STRSXP,2));
	SET_STRING_ELT(solNames,0,mkChar("beta"));
	SET_STRING_ELT(solNames,1,mkChar("status"));
	setAttrib(solution,R_NamesSymbol,solNames);

	if(stoc!=0){
    PutRNGstate();
	}

	UNPROTECT(4);

	/* free memory */
	if(e) Free(e);
	if(e_j) Free(e_j);
	if(bPrec) Free(bPrec);
	if(Xe_j) Free(Xe_j);
	if(s) Free(s);

	return solution;
}

