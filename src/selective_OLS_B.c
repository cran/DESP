/*
 *  DESP/src/selective_OLS_B.c by A. S. Dalalyan and S. Balmand  Copyright (C) 2015-
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


SEXP selective_OLS_B(SEXP X, SEXP SPC){
	/* estimation of B */
	/* 
	X : data matrix
	SPC : matrix of squared partial correlations
	*/

	SEXP B, dim, solution, solNames, t;
	double * X_j, * xj, * work, * XX, * tt, * y, * y2, * QQ , * BB, * select;
	int j, info, lwork, i, k;
	double tmp;
	int * ipiv, * sel;

	/* read the sample size and the number of variables */
	dim = getAttrib(X,R_DimSymbol);
	int n = INTEGER(dim)[0];               /* n is the sample size */
	int p = INTEGER(dim)[1];               /* p is the dimension */


	/* initialize the matrix B */
	PROTECT(B = allocMatrix(REALSXP,p,p));
	BB = REAL(B);
	memset(BB, 0, sizeof(double)*p*p);



	int one = 1;
	double done = 1.0;
	double mone = -1.0;
	double dzero = 0.0;

	X_j = (double *) Calloc((unsigned) n*(p-1), double);
	xj = (double *) Calloc((unsigned) n, double);
	QQ = (double *) Calloc((unsigned) (p-1)*(p-1), double);
	y = (double *) Calloc((unsigned) p-1, double);
	y2 = (double *) Calloc((unsigned) p-1, double);
	sel = (int *) Calloc((unsigned) p-1, int);

	XX = REAL(X);
	select = REAL(SPC);

	ipiv = (int *) Calloc((unsigned) p-1, int);

	/* we perform a workspace query outside of the for loop, as the optimal size of the "work" array is obtained by a call to the Lapack ilaenv routine and because this routine only uses the size of the problem, not the actual data. */
	/* we allocate more space than if the computation of the size of the problem has been done in the loop, but only once */
	k = p-1; 
	lwork = -1;
	F77_CALL(dsytrf)("U", &k, QQ, &k, ipiv, &tmp, &lwork, &info);
	lwork = (int) tmp;

	work = (double *) Calloc((size_t) lwork, double);

	for(j=0; j<p; j++){
		k = 0;
		for(i=0; i<p; i++){
			if(i!=j && select[p*j+i]!=0){
				memcpy(X_j+(n*k), XX+(n*i), n*sizeof(double));
				sel[k] = i;
				k += 1;
			}
		}
		memcpy(xj, XX+(n*j), n*sizeof(double));
		
		if(k!=0){
			/* we compute X_j^T * X_j, the result is stored in QQ */
			F77_CALL(dsyrk)("U", "T", &k, &n, &done, X_j, &n, &dzero, QQ, &k);
			/* we compute QQ^-1, the result is stored in QQ */
			F77_CALL(dsytrf)("U", &k, QQ, &k, ipiv, work, &lwork, &info);
			F77_CALL(dsytri)("U", &k, QQ, &k, ipiv, work, &info);

			/* we compute X_j^T * xj, the result is stored in y */
			F77_CALL(dgemv)("T", &n, &k, &done, X_j, &n, xj, &one, &dzero, y, &one);
			/* we compute QQ * y, the result is stored in y2 */
			F77_CALL(dsymv)("U", &k, &mone, QQ, &k, y, &one, &dzero, y2, &one);

			for(i=0; i<k; i++){
				BB[p*j+sel[i]] = y2[i];
			}
		}
		BB[p*j+j] = 1.0;
	}


	/* free memory */
	if(xj) Free(xj);
	if(X_j) Free(X_j);
	if(work) Free(work);
	if(QQ) Free(QQ);
	if(ipiv) Free(ipiv);
	if(sel) Free(sel);
	if(y) Free(y);
	if(y2) Free(y2);

	UNPROTECT(1);
	return B;
}


