/*
 *  DESP/src/utils.c by A. S. Dalalyan and S. Balmand  Copyright (C) 2015-
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

#include "utils.h"


SEXP largestLambda(SEXP X) {
	// compute the largest value of the tuning parameter that encourages robustness, such that outliers can be detected
	/*
	X : data matrix
	*/

	//if (!isMatrix(X)) error(_("'X' must be a numeric matrix"));

	SEXP ans, dim;
	double * X_j, * xj, * work, * work2, * A, * A2, * vone, * d, * y, * y2, * QQ;
	int j, info, lwork, nb, i;
	double tmp;
	int * ipiv;

	/* read the sample size and the number of variables */
	dim = getAttrib(X,R_DimSymbol);
	int n = INTEGER(dim)[0];               /* n is the sample size */
	int p = INTEGER(dim)[1];               /* p is the dimension */

	int one = 1;
	double done = 1.0;
	double dzero = 0.0;

	X_j = (double *) Calloc((size_t) n*(p-1), double);
	xj = (double *) Calloc((size_t) n, double);
	vone = (double *) Calloc((size_t) n, double);
	d = (double *) Calloc((size_t) p, double);
	work2 = (double *) Calloc((size_t) n, double);
	QQ = (double *) Calloc((size_t) (p-1)*(p-1), double);
	y = (double *) Calloc((size_t) p-1, double);
	y2 = (double *) Calloc((size_t) p-1, double);

	A = REAL(X);

	A2 = (double *) Calloc((size_t) n*p, double);

	ans = allocVector(REALSXP, 1);

	ipiv = (int *) Calloc((size_t) p-1, int);
	
	/* we perform a workspace query outside of the for loop, as the optimal size of the "work" array is obtained by a call to the Lapack ilaenv routine and because this routine only uses the size of the problem, not the actual data. */
	nb = p-1;
	lwork = -1;
	F77_CALL(dsytrf)("U", &nb, QQ, &nb, ipiv, &tmp, &lwork, &info);
	lwork = (int) tmp;

	work = (double *) Calloc((size_t) lwork, double);

	for(j=0; j<p; j++){
		nb = n*j;
		memcpy(X_j, A, nb*sizeof(double));
		nb = n*(p-j-1);
		memcpy(X_j+(n*j), A+(n*(j+1)), nb*sizeof(double));
		memcpy(xj, A+(n*j), n*sizeof(double));

		/* we compute X_j^T * X_j, the result is stored in QQ */
		nb = p-1;
		F77_CALL(dsyrk)("U", "T", &nb, &n, &done, X_j, &n, &dzero, QQ, &nb);
		/* we compute QQ^-1, the result is stored in QQ */
		F77_CALL(dsytrf)("U", &nb, QQ, &nb, ipiv, work, &lwork, &info);
		F77_CALL(dsytri)("U", &nb, QQ, &nb, ipiv, work, &info);

		/* we compute X_j^T * xj, the result is stored in y */
		F77_CALL(dgemv)("T", &n, &nb, &done, X_j, &n, xj, &one, &dzero, y, &one);
		/* we compute QQ * y, the result is stored in y2 */
		F77_CALL(dsymv)("U", &nb, &done, QQ, &nb, y, &one, &dzero, y2, &one);
		/* we compute X_j * y2, the result is stored in A2+(n*j) */
		F77_CALL(dgemv)("N", &n, &nb, &done, X_j, &n, y2, &one, &dzero, A2+(n*j), &one);

	}
	for(i=0; i<n*p; i++){
		A2[i] = pow(A[i]-A2[i],2);
	}
	for (i = 0; i < n; i++) vone[i] = 1.0;
	F77_CALL(dgemv)("T", &n, &p, &done, A2, &n, vone, &one, &dzero, d, &one);

	for(j=0; j<p; j++){
		for(i=0; i<n; i++){
			A2[i+n*j] /= d[j];
		}
	}

	REAL(ans)[0] = sqrt(F77_CALL(dlange)("I", &n, &p, A2, &n, work2));
	
	/* free memory */
	if(work2) Free(work2);
	if(work) Free(work);
	if(ipiv) Free(ipiv);
	if(d) Free(d);
	if(vone) Free(vone);
	if(A2) Free(A2);
	if(y2) Free(y2);
	if(y) Free(y);
	if(QQ) Free(QQ);
	if(xj) Free(xj);
	if(X_j) Free(X_j);

	return ans;
}



SEXP update_theta(SEXP X, SEXP B, SEXP t, SEXP lambda){
	// compute the matrix Theta that corresponds to the outliers
	/*
	X : data matrix
	B : matrix of the coefficients of regression
	t : vector of the Euclidean norms of the columns of the matrix (X * B - Theta)
	lambda : tuning parameter that encourages robustness
	*/

	SEXP Theta, dim;
	double * XX, * BB, * TTheta, * y, * tt, * lbdat;
	int i, j;
	double normE2, lbda, lbda2, x;

	/* read the sample size and the number of variables */
	dim = getAttrib(X,R_DimSymbol);
	int n = INTEGER(dim)[0];               /* n is the sample size */
	int p = INTEGER(dim)[1];               /* p is the dimension */

	PROTECT(Theta = allocMatrix(REALSXP,n,p));
	TTheta = REAL(Theta);

	double done = 1.0;
	double dzero = 0.0;

	y = (double *) Calloc((unsigned) p*n, double);
	lbdat = (double *) Calloc((unsigned) p, double);

	XX = REAL(X);
	BB = REAL(B);
	tt = REAL(t);
	lbda = (double) *(REAL(lambda));

	memset(TTheta, 0, n*p*sizeof(double));

	/* we compute B^T * X^T, the result is stored in y */
	F77_CALL(dgemm)("T", "T", &p, &n, &p, &done, BB, &p, XX, &n, &dzero, y, &p);

	for(j=0; j<p; j++){
		lbdat[j] = lbda * tt[j];
	}
	lbda2 = pow(lbda, 2);
	for(i=0; i<n; i++){
		normE2 = 0.0;
		for(j=0; j<p; j++){
			normE2 += pow(y[i*p+j] / tt[j], 2);
		}
		if(normE2 > lbda2){
			x = thetaNorm2(y+(i*p), lbdat, p, 1e-6);
			for(j=0; j<p; j++){
				TTheta[i + j*n] = x * y[i*p+j] / (x + lbdat[j]);
			}
		}
	}
	if(y) Free(y);
	if(lbdat) Free(lbdat);

	UNPROTECT(1);
	return Theta;
}

double thetaNorm2(double * u, double * t, int p, double thres){
	// compute the Euclidean norm of the i-th row of the matrix Theta corresponding to outliers, using the method of Newton
	/*
	u : corresponds to the i-th row of the matrix X * B
	t : corresponds to lambda * t
	p : dimension (number of columns of Theta)
	thres : threshold depending on the desired accuracy
	*/

	double f_x, f_deriv, x, no;
	double * v;

	v = (double *) Calloc((unsigned) p, double);

	no = 0.0;
	for (int i = 0; i < p; i++) {
		no += pow(u[i]/t[i], 2);
	}
	x = 0.0;
	if (no > 1){
		memcpy(v, t, p*sizeof(double));
		f_x = no;
		while (f_x-1>=thres || 1-f_x>=thres){
			f_deriv = 0;
			for (int i = 0; i < p; i++) {
				f_deriv += pow(u[i], 2) / pow(v[i], 3);
			}
			f_deriv *= 2;

			x += (f_x-1)/f_deriv;
			for (int i = 0; i < p; i++) {
				v[i] = x+t[i];
			}
			no = 0.0;
			for (int i = 0; i < p; i++) {
				no += pow(u[i]/v[i], 2);
			}
			f_x = no;
		}
	}

	if(v) Free(v);

	return x;
}



