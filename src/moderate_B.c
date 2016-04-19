/*
 *  DESP/src/moderate_B.c by A. S. Dalalyan and S. Balmand  Copyright (C) 2015-
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
#include "R_ext/BLAS.h"
#include "R_ext/Lapack.h"


SEXP moderate_B(SEXP X, SEXP theta){
	/* estimation of B */
	/* 
	X : data matrix
	theta : matrix of outliers, if NULL, the computation is similar to theta=0
	*/

	SEXP B, dim, solution, solNames, t;
	double * X_j, * xj, * work, * XX, * tt, * y, * y2, * z, * QQ , * BB, * TTheta;
	int j, info, lwork, nb, i, detect;
	double tmp;
	int * ipiv;

	/* read the sample size and the number of variables */
	dim = getAttrib(X,R_DimSymbol);
	int n = INTEGER(dim)[0];               /* n is the sample size */
	int p = INTEGER(dim)[1];               /* p is the dimension */

	/* check whether the robust version has to be used or not */
	detect = 0;
	if(!isNull(theta)){ detect = 1;}

	PROTECT(solution = allocVector(VECSXP, 2));

	/* initialize the matrix B */
	PROTECT(B = allocMatrix(REALSXP,p,p));
	BB = REAL(B);
	memset(BB, 0, sizeof(double)*p*p);

	PROTECT(t = allocVector(REALSXP,p));

	int one = 1;
	double done = 1.0;
	double mone = -1.0;
	double dzero = 0.0;

	X_j = (double *) Calloc((unsigned) n*(p-1), double);
	xj = (double *) Calloc((unsigned) n, double);
	QQ = (double *) Calloc((unsigned) (p-1)*(p-1), double);
	y = (double *) Calloc((unsigned) p-1, double);
	y2 = (double *) Calloc((unsigned) p-1, double);
	z = (double *) Calloc((unsigned) n, double);

	XX = REAL(X);
	if(detect){ 
	  TTheta = REAL(theta);
  }else{
    TTheta = NULL;
  }
	tt = REAL(t);

	ipiv = (int *) Calloc((unsigned) p-1, int);

	/* we perform a workspace query outside of the for loop, as the optimal size of the "work" array is obtained by a call to the Lapack ilaenv routine and because this routine only uses the size of the problem, not the actual data. */
	nb = p-1;
	lwork = -1;
	F77_CALL(dsytrf)("U", &nb, QQ, &nb, ipiv, &tmp, &lwork, &info);
	lwork = (int) tmp;

	work = (double *) Calloc((size_t) lwork, double);

	for(j=0; j<p; j++){
		nb = n*j;
		memcpy(X_j, XX, nb*sizeof(double));
		nb = n*(p-j-1);
		memcpy(X_j+(n*j), XX+(n*(j+1)), nb*sizeof(double));
		memcpy(xj, XX+(n*j), n*sizeof(double));

		if(detect == 1){
			for(i=0; i<n; i++){
				xj[i] -= TTheta[i+j*n];
			}
		}

		/* we compute X_j^T * X_j, the result is stored in QQ */
		nb = p-1;
		F77_CALL(dsyrk)("U", "T", &nb, &n, &done, X_j, &n, &dzero, QQ, &nb);
		/* we compute QQ^-1, the result is stored in QQ */
		F77_CALL(dsytrf)("U", &nb, QQ, &nb, ipiv, work, &lwork, &info);
		F77_CALL(dsytri)("U", &nb, QQ, &nb, ipiv, work, &info);

		/* we compute X_j^T * xj, the result is stored in y */
		F77_CALL(dgemv)("T", &n, &nb, &done, X_j, &n, xj, &one, &dzero, y, &one);
		/* we compute QQ * y, the result is stored in y2 */
		F77_CALL(dsymv)("U", &nb, &mone, QQ, &nb, y, &one, &dzero, y2, &one);

		memcpy(BB+(p*j), y2, j*sizeof(double));
		memcpy(BB+(p*j+j+1), y2+j, (p-j-1)*sizeof(double));
		BB[p*j+j] = 1.0;

		if(detect){
			memcpy(z, TTheta+(j*n), n*sizeof(double)); 

			/* compute z = X*B[,j]-Theta[,j] */
			F77_CALL(dgemv)("N", &n, &p, &done, XX, &n, BB+(j*p), &one, &mone, z, &one);
			tt[j] = F77_CALL(dnrm2)(&n, z, &one);
		}
	}

	SET_VECTOR_ELT(solution,0,B);
	SET_VECTOR_ELT(solution,1,t);
	
	PROTECT(solNames = allocVector(STRSXP,2));
	SET_STRING_ELT(solNames,0,mkChar("B"));
	SET_STRING_ELT(solNames,1,mkChar("t"));
	setAttrib(solution,R_NamesSymbol,solNames);

	/* free memory */
	if(xj) Free(xj);
	if(X_j) Free(X_j);
	if(work) Free(work);
	if(QQ) Free(QQ);
	if(ipiv) Free(ipiv);
	if(y) Free(y);
	if(y2) Free(y2);
	if(z) Free(z);

	UNPROTECT(4);
	return solution;
}


