/*
 *  DESP/src/scsSqR_Lasso_solve.c by A. S. Dalalyan and S. Balmand  Copyright (C) 2015-
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

#include "scsSqR_Lasso_solve.h"
#include "R_ext/BLAS.h"

/* default settings values */ 
#define DEFAULT_normalize	1	/* boolean, heuristic data rescaling: 1 */
#define DEFAULT_scale		5	/* if normalized, rescales by this factor: 5 */
#define DEFAULT_rho_x		1e-3	/* x equality constraint scaling: 1e-3 */
#define DEFAULT_max_iters	2500	/* maximum iterations to take: 2500 */
#define DEFAULT_eps		1e-3	/* convergence tolerance: 1e-3 */
#define DEFAULT_alpha	 	1.8	/* relaxation parameter: 1.8 */
#define DEFAULT_cg_rate		2	/* for indirect, tolerance goes down like (1/iter)^cg_rate: 2 */
#define DEFAULT_verbose		1	/* boolean, write out progress: 1 */
#define DEFAULT_warm_start	0	/* boolean, warm start (put initial guess in Sol struct): 0 */


int scs_sqR_Lasso_solve(double * X, double * Y, double lambda, int n, int p, double * beta){

	Cone * k;
	Data * d;
	Sol * sol;
	Info info = { 0 };
	AMatrix * A;
	int i, j, st, idx, nnz;
	double * vn;
	scs_int * ii, * pp;
	scs_float * xx;

	int one = 1;
	int zero = 0;
	double done = 1.0;
	double dzero = 0.0;

	/* normalize the columns of X by dividing them by their norm */
	vn = (double *) calloc(p, sizeof(double));

	for(j=0; j<p; j++){
		vn[j] = sqrt(n) / F77_CALL(dnrm2)(&n, X+(j*n), &one);
	}

	/*
	 * set up the A matrix having 2p+2n+1 rows and 2p+n+1 columns
	 * ii integer vector of length nnz of row indices
	 * pp integer vector of row or column pointers
	 * xx vector of values
	 */	
	ii = (scs_int *) scs_calloc(4*p+2*n+n*p+1, sizeof(scs_int));
	xx = (scs_float *) scs_calloc(4*p+2*n+n*p+1, sizeof(scs_float));
	pp = (scs_int *) scs_calloc(2*p+n+2, sizeof(scs_int));

	idx = 0;
	pp[0] = 0;
	for(i=0; i<p; i++){
		xx[idx] = -1;
		ii[idx ++] = n+i;
		xx[idx] = -1;
		ii[idx ++] = n+p+i;
		pp[i+1] = pp[i] + 2;
	}
	for(j=0; j<p; j++){
		nnz = 0;
		for(i=0; i<n; i++){
			if (X[i+j*n] != 0){
				xx[idx] = X[i+j*n] * vn[j];
				ii[idx ++] = i;
				nnz += 1;
			}
		}
		xx[idx] = -1;
		ii[idx ++] = n+j;
		xx[idx] = 1;
		ii[idx ++] = n+p+j;
		pp[p+j+1] = pp[p+j] + nnz+2;
	}
	for(i=0; i<n; i++){
		xx[idx] = 1;
		ii[idx ++] = i;
		xx[idx] = 1;
		ii[idx ++] = n+2*p+1+i;
		pp[2*p+i+1] = pp[2*p+i] + 2;
	}
	xx[idx] = -1;
	ii[idx] = n+2*p;
	pp[2*p+n+1] = pp[2*p+n] + 1;

	k = (Cone *) scs_calloc(1, sizeof(Cone)); 
	k->q = (scs_int *) scs_calloc(1, sizeof(scs_int)); /* K vector of second-order cone constraints */
	d = (Data *) scs_calloc(1, sizeof(Data));
	d->stgs = (Settings *) scs_calloc(1, sizeof(Settings));
	d->b = (scs_float *) scs_calloc(2*p+2*n+1, sizeof(scs_float));
	d->c = (scs_float *) scs_calloc(2*p+n+1, sizeof(scs_float)); /* c vector : the cost function */
	A = (AMatrix *) scs_malloc(sizeof(AMatrix));
	sol = (Sol *) scs_calloc(1, sizeof(Sol));

  	/* set up SCS structures */
	d->m = (scs_int) 2*n+2*p+1; /* number of rows of A */
	d->n = (scs_int) 2*p+n+1; /* number of columns of A */
	F77_CALL(dcopy)(&n, Y, &one, d->b, &one);
	for(i=0; i<p; i++) (d->c)[i] = lambda; 
	(d->c)[2*p+n] = 1;

	A->n = d->n;
	A->m = d->m;
	A->x = xx;
	A->i = ii;
	A->p = pp;
	d->A = A;

	k->f = (scs_int) n;
	k->l = (scs_int) 2*p;
	k->qsize = (scs_int) 1;
	(k->q)[0] = 1+n;

	/* settings */
	d->stgs->alpha = (scs_float) DEFAULT_alpha;
	d->stgs->rho_x = (scs_float) DEFAULT_rho_x;
	d->stgs->max_iters = (scs_int) DEFAULT_max_iters;
	d->stgs->scale = (scs_float) DEFAULT_scale;
	d->stgs->eps = (scs_float) DEFAULT_eps;
	d->stgs->cg_rate = (scs_float) DEFAULT_cg_rate;
	d->stgs->verbose = 0; /* silent */
	d->stgs->normalize = 0;
	d->stgs->warm_start = (scs_int) DEFAULT_warm_start;

	/* solve! */
	st = scs(d, k, sol, &info);

	F77_CALL(dsbmv)("L", &p, &zero, &done, (sol->x)+p, &one, vn, &one, &dzero, beta, &one); 
	
	/* free memory */
	scs_free(ii);
	scs_free(xx);
	scs_free(pp);
	if (k->q) scs_free(k->q);
	scs_free(k);
	if (d->c) scs_free(d->c);
	if (d->b) scs_free(d->b);
	if (d->stgs) scs_free(d->stgs);
	scs_free(d);
	scs_free(A);
	if (sol->s) scs_free(sol->s);
	if (sol->y) scs_free(sol->y);
	if (sol->x) scs_free(sol->x);
	scs_free(sol);
	free(vn);

	return st;
}


int scs_sqR_Lasso_solve_loop(double * X, double * theta, double lambda, int n, int p, double * beta, int col){
	/* 
	X : data matrix 
	theta : matrix of outliers
	lambda : tuning parameter that promotes sparsity
	n : sample size (number of rows of X)
	p : dimension (number of columns of X)
	beta : vector of the coefficients of regression
	col : column index in X of the response variable of the regression model
	*/

	Cone * k;
	Data * d;
	Sol * sol;
	Info info = { 0 };
	AMatrix * A;
	int i, j, h, st, idx, nnz;
	double * vn;
	scs_int * ii, * pp;
	scs_float * xx;

	int one = 1;
	int zero = 0;
	double done = 1.0;
	double dzero = 0.0;

	/* normalize the columns of X by dividing them by their norm */
	vn = (double *) calloc(p-1, sizeof(double));


	for(j=0; j<p-1; j++){
		h = j;
		if(j >= col) h++;
		vn[j] = sqrt(n) / F77_CALL(dnrm2)(&n, X+(h*n), &one);
	}

	/*
	 * set up the A matrix having 2(p-1)+2n+1 rows and 2(p-1)+n+1 columns
	 * ii integer vector of length nnz of row indices
	 * pp integer vector of row or column pointers
	 * xx vector of values
	 */	
	ii = (scs_int *) scs_calloc(4*(p-1)+2*n+n*(p-1)+1, sizeof(scs_int));
	xx = (scs_float *) scs_calloc(4*(p-1)+2*n+n*(p-1)+1, sizeof(scs_float));
	pp = (scs_int *) scs_calloc(2*(p-1)+n+2, sizeof(scs_int));

	idx = 0;
	pp[0] = 0;
	for(i=0; i<p-1; i++){
		xx[idx] = -1;
		ii[idx ++] = n+i;
		xx[idx] = -1;
		ii[idx ++] = n+p-1+i;
		pp[i+1] = pp[i] + 2;
	}
	for(j=0; j<p-1; j++){
		h = j;
		if(j >= col) h++;
		nnz = 0;
		for(i=0; i<n; i++){
			if (X[i+h*n] != 0){
				xx[idx] = X[i+h*n] * vn[j];
				ii[idx ++] = i;
				nnz += 1;
			}
		}
		xx[idx] = -1;
		ii[idx ++] = n+j;
		xx[idx] = 1;
		ii[idx ++] = n+p-1+j;
		pp[p-1+j+1] = pp[p-1+j] + nnz+2;
	}
	for(i=0; i<n; i++){
		xx[idx] = 1;
		ii[idx ++] = i;
		xx[idx] = 1;
		ii[idx ++] = n+2*(p-1)+1+i;
		pp[2*(p-1)+i+1] = pp[2*(p-1)+i] + 2;
	}
	xx[idx] = -1;
	ii[idx] = n+2*(p-1);
	pp[2*(p-1)+n+1] = pp[2*(p-1)+n] + 1;

	k = (Cone *) scs_calloc(1, sizeof(Cone)); 
	k->q = (scs_int *) scs_calloc(1, sizeof(scs_int)); /* K vector of second-order cone constraints */
	d = (Data *) scs_calloc(1, sizeof(Data));
	d->stgs = (Settings *) scs_calloc(1, sizeof(Settings));
	d->b = (scs_float *) scs_calloc(2*(p-1)+2*n+1, sizeof(scs_float));
	d->c = (scs_float *) scs_calloc(2*(p-1)+n+1, sizeof(scs_float)); /* c vector : the cost function */
	A = (AMatrix *) scs_malloc(sizeof(AMatrix));
	sol = (Sol *) scs_calloc(1, sizeof(Sol));

  	/* set up SCS structures */
	d->m = (scs_int) 2*n+2*(p-1)+1; /* number of rows of A */
	d->n = (scs_int) 2*(p-1)+n+1; /* number of columns of A */
	if(theta != NULL){
		for(i=0; i<n; i++) (d->b)[i] = theta[i+col*n] - X[i+col*n]; 
	}else{
		for(i=0; i<n; i++) (d->b)[i] = - X[i+col*n]; 
	}
	for(i=0; i<p-1; i++) (d->c)[i] = lambda; 
	(d->c)[2*(p-1)+n] = 1;

	A->n = d->n;
	A->m = d->m;
	A->x = xx;
	A->i = ii;
	A->p = pp;
	d->A = A;

	k->f = (scs_int) n;
	k->l = (scs_int) 2*(p-1);
	k->qsize = (scs_int) 1;
	(k->q)[0] = 1+n;

	/* settings */
	d->stgs->alpha = (scs_float) DEFAULT_alpha;
	d->stgs->rho_x = (scs_float) DEFAULT_rho_x;
	d->stgs->max_iters = (scs_int) DEFAULT_max_iters;
	d->stgs->scale = (scs_float) DEFAULT_scale;
	d->stgs->eps = (scs_float) DEFAULT_eps;
	d->stgs->cg_rate = (scs_float) DEFAULT_cg_rate;
	d->stgs->verbose = 0; /* silent */
	d->stgs->normalize = 0;
	d->stgs->warm_start = (scs_int) DEFAULT_warm_start;

	/* solve! */
	st = scs(d, k, sol, &info);

	h = p-1;
	F77_CALL(dsbmv)("L", &h, &zero, &done, (sol->x)+(p-1), &one, vn, &one, &dzero, beta, &one); 
	
	/* free memory */
	scs_free(ii);
	scs_free(xx);
	scs_free(pp);
	if (k->q) scs_free(k->q);
	scs_free(k);
	if (d->c) scs_free(d->c);
	if (d->b) scs_free(d->b);
	if (d->stgs) scs_free(d->stgs);
	scs_free(d);
	scs_free(A);
	if (sol->s) scs_free(sol->s);
	if (sol->y) scs_free(sol->y);
	if (sol->x) scs_free(sol->x);
	scs_free(sol);
	free(vn);

	return st;
}

