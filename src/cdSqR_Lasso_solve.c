/*
 *  DESP/src/cdSqR_Lasso_solve_loop.c by A. S. Dalalyan and S. Balmand  Copyright (C) 2015-
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
#define DEFAULT_tol		1e-4	/* convergence tolerance: 1e-6 */
#define DEFAULT_tolGap		1e-6	/* duality gap tolerance: 1e-8 */

int cd_sqR_Lasso_solve_loop(double * X, double * theta, double lambda, int n, int p, double * beta, int col, int sto){
	/* 
	X : data matrix 
	theta : matrix of outliers
	lambda : tuning parameter that promotes sparsity
	n : sample size (number of rows of X)
	p : dimension (number of columns of X)
	beta : vector of the coefficients of regression
	col : column index in X of the response variable of the regression model
	sto : whether a stochastic coordinate descent must be performed or not
	*/

	double * bPrec, * e_j, * Xe_j, * e, * Y;
	int i, j, k, l, stop, st;
	double Q_j, xe_j, mb_j, ab, abp, srss, Ye_j, primObj, dualObj;


	double dzero = 0.0;
	double done = 1.0;
	double mone = -1.0;
	int one = 1;

	int maxIter = (int) DEFAULT_max_iters;
	double tol = (double) DEFAULT_tol;
	double tolGap = (double) DEFAULT_tolGap;

	/* initialize beta to zero */
	memset(beta, 0, p*sizeof(double));

	/* for stochastic coordinate descent */
	if(sto==1){
		maxIter *= p;
	}

	/* for shuffled coordinate descent */
	int * s = 0;
	if(sto==2){
		s = (int *) Calloc((unsigned) p, int);
	}

	bPrec = (double *) Calloc((unsigned) p, double);
	Xe_j = (double *) Calloc((unsigned) p, double);

	e = (double *) Calloc((unsigned) p, double);
	Y = (double *) Calloc((unsigned) n, double);
	e_j = (double *) Calloc((unsigned) n, double);
	memset(e, 0, p*sizeof(double));
	for(j=0; j<p; j++){
		if (j != col){
			for(i=0; i<n; i++){
			 	e[j] += pow(X[i+j*n],2);
			}
		}
		e[j] /= n;
	}
	
	if(theta != NULL){
		memcpy(Y, theta+col*n, n*sizeof(double));
	}else{
		memset(Y, 0, n*sizeof(double));
	}
	F77_CALL(daxpy)(&n, &mone, X+col*n, &one, Y, &one);

	memcpy(e_j, Y, n*sizeof(double));

	Q_j = F77_CALL(ddot)(&n,e_j,&one,e_j,&one);
	Q_j /= n;

	stop = 0;
	k = 1;
	while(stop != 1 && k <= maxIter){
		memcpy(bPrec, beta, p*sizeof(double));

		ab = 0; // sum of the absolute values of beta
		abp = 0; // sum of the absolute values of beta - bPrec
		srss = 0; // residual sum of squares
		
		/* for shuffled coordinate descent */
		if(sto==2){
			sample(p, s);
		}
		l = 0;
		while(l<p){
			/* for stochastic coordinate descent */
			if(sto==1 && k!=1){
				// stochastic
				j = (int) floor(unif_rand() * (p-1));
				if(j==col){
					j++;
				}
				l = p;
			}else if(sto==2){
				// shuffle
				j = s[l];
			}else{
				j = l;
			}
			if(j!=col){
				if (beta[j] != 0){
					F77_CALL(daxpy)(&n, beta+j, X+n*j, &one, e_j, &one);

					Q_j = F77_CALL(ddot)(&n,e_j,&one,e_j,&one);
					Q_j /= n;
				}
				xe_j = F77_CALL(ddot)(&n,X+n*j,&one,e_j,&one);
				xe_j /= n;

				if (pow(n,2) < pow(lambda,2) / e[j]){ 
					beta[j] = 0;
				}else if (xe_j > lambda/n * sqrt(Q_j)){
		       		 	beta[j] = (xe_j - lambda*sqrt((Q_j-pow(xe_j,2)/e[j])/(pow(n,2)-(pow(lambda,2) / e[j]))))/e[j];
					mb_j = -beta[j];
					F77_CALL(daxpy)(&n, &mb_j, X+n*j, &one, e_j, &one);
				}else if(xe_j < - lambda/n * sqrt(Q_j)){
		        		beta[j] = (xe_j + lambda*sqrt((Q_j-pow(xe_j,2)/e[j])/(pow(n,2)-(pow(lambda,2) / e[j]))))/e[j];
					mb_j = -beta[j];
					F77_CALL(daxpy)(&n, &mb_j, X+n*j, &one, e_j, &one);
				}else{
		       	 		beta[j] = 0;
				}
			}
			l++;
		}
		for(j=0; j<p; j++){
			if(j!=col){
				ab += fabs(beta[j]);
				abp += fabs(beta[j]-bPrec[j]);
			}
		}

		/* primal objective */
		Q_j = F77_CALL(ddot)(&n,e_j,&one,e_j,&one);
		Q_j /= n;
		primObj = sqrt(Q_j) + lambda/n * ab;
		
		
		srss = F77_CALL(ddot)(&n,e_j,&one,e_j,&one);
		srss = sqrt(srss); // squared Residual Sum of Squares
		if (srss > 1e-10){
			Ye_j = F77_CALL(ddot)(&n,Y,&one,e_j,&one);
			F77_CALL(dgemv)("T", &n, &p, &done, X, &n, e_j, &one, &dzero, Xe_j, &one);
			dualObj = 0;
			for(j=0; j<p; j++){
			if(j!=col){
				dualObj -= fabs(lambda/n - fabs(Xe_j[j]/(sqrt(n)*srss))) * fabs(beta[j]);
			}
			}
			dualObj += Ye_j/(sqrt(n)*srss);
		}else{
			dualObj = lambda/n * ab;
		}
       
		/* stopping criterion */
		if (abp < tol){
			if (primObj - dualObj < tolGap){
				stop = 1;
			}
		}
		k++;
	}

	st = 0;
	if (k == maxIter) st = -1;



	/* free memory */
	if(e) Free(e);
	if(e_j) Free(e_j);
	if(Y) Free(Y);
	if(bPrec) Free(bPrec);
	if(Xe_j) Free(Xe_j);
	if(s) Free(s);

	return st;
}


void sample(int n, int * s){
	/* Fisher-Yates shuffle */

	int i, j, tmp;
	
	for(i=0; i<n; i++){
		s[i] = i;
	}
	for(i=0; i<n; i++){
		j = (int) floor(unif_rand() * n);
		tmp = s[i];
		s[i] = s[j];
		s[j] = tmp;
	}
}

