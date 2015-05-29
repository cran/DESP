/*
 *  DESP/src/scsSolveSOCP.c by A. S. Dalalyan and S. Balmand  Copyright (C) 2015-
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
#include "scs.h"
#include "linsys/amatrix.h"

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

SEXP scs_SOCP_solve(SEXP Ax, SEXP Ai, SEXP Ap, SEXP Am, SEXP An, SEXP b, SEXP c, SEXP Kf, SEXP Kl, SEXP Kq, SEXP Kqsize, SEXP alpha, SEXP rho_x, SEXP max_iters, SEXP scale, SEXP eps, SEXP cg_rate, SEXP verbose, SEXP normalize, SEXP warm_start){
	Cone * k;
	Data * d;
	Sol * sol;
	Info info = { 0 };
	AMatrix * A;
	int i, st;
	SEXP x, status, solution, solNames;
	double *xx;
	int * stat;

	k = (Cone *) scs_calloc(1, sizeof(Cone));
	d = (Data *) scs_calloc(1, sizeof(Data));
	d->stgs = (Settings *) scs_calloc(1, sizeof(Settings));
	A = (AMatrix *) scs_malloc(sizeof(AMatrix));
	sol = (Sol *) scs_calloc(1, sizeof(Sol));

  	/* set up SCS structures */
	d->m = (scs_int) *(INTEGER(Am)); /* A has m rows */
	d->n = (scs_int) *(INTEGER(An)); /* A has n cols */
	d->b = (scs_float *) scs_malloc(d->m * sizeof(scs_float));
	d->b = (scs_float *) REAL(b);
	d->c = (scs_float *) scs_malloc(d->n * sizeof(scs_float));
	d->c = (scs_float *) REAL(c);

	A->n = d->n;
	A->m = d->m;
	A->x = (scs_float *) REAL(Ax);
	A->i = (scs_int *) INTEGER(Ai);
	A->p = (scs_int *) INTEGER(Ap);
	d->A = A;

	k->f = (scs_int) *(INTEGER(Kf));
	k->l = (scs_int) *(INTEGER(Kl));
	k->qsize = (scs_int) *(INTEGER(Kqsize));
	k->q = (scs_int *) scs_malloc(k->qsize * sizeof(scs_int));
	k->q = (scs_int *) INTEGER(Kq);

	/* settings */
	if (*(REAL(alpha)) != -1)
		d->stgs->alpha = (scs_float) *(REAL(alpha));
	else d->stgs->alpha = (scs_float) DEFAULT_alpha;
	if (*(REAL(rho_x)) != -1)
		d->stgs->rho_x = (scs_float) *(REAL(rho_x));
	else d->stgs->rho_x = (scs_float) DEFAULT_rho_x;
	if (*(INTEGER(max_iters)) != -1)
		d->stgs->max_iters = (scs_int) *(INTEGER(max_iters));
	else d->stgs->max_iters = (scs_int) DEFAULT_max_iters;
	if (*(REAL(scale)) != -1)
		d->stgs->scale = (scs_float) *(REAL(scale));
	else d->stgs->scale = (scs_float) DEFAULT_scale;
	if (*(REAL(eps)) != -1)
		d->stgs->eps = (scs_float) *(REAL(eps));
	else d->stgs->eps = (scs_float) DEFAULT_eps;
	if (*(REAL(cg_rate)) != -1)
		d->stgs->cg_rate = (scs_float) *(REAL(cg_rate));
	else d->stgs->cg_rate = (scs_float) DEFAULT_cg_rate;
	if (*(INTEGER(verbose)) != -1)
		d->stgs->verbose = (scs_int) *(INTEGER(verbose));
	else d->stgs->verbose = (scs_int) DEFAULT_verbose;
	if (*(INTEGER(normalize)) != -1)
		d->stgs->normalize = (scs_int) *(INTEGER(normalize));
	else d->stgs->normalize = (scs_int) DEFAULT_normalize;
	if (*(INTEGER(warm_start)) != -1)
		d->stgs->warm_start = (scs_int) *(INTEGER(warm_start));
	else d->stgs->warm_start = (scs_int) DEFAULT_warm_start;

	/* solve! */
	PROTECT(status = allocVector(INTSXP, 1));
	stat = INTEGER(status);
	st = scs(d, k, sol, &info);
	stat[0] = st;

	PROTECT(solution = allocVector(VECSXP, 2));
	PROTECT(x = allocVector(REALSXP, d->n));
	xx = REAL(x);

	for (i = 0; i < A->n; i++){
		xx[i] = (sol->x)[i];
	}
	
	SET_VECTOR_ELT(solution,0,x);
	SET_VECTOR_ELT(solution,1,status);
	
	PROTECT(solNames = allocVector(STRSXP,2));
	SET_STRING_ELT(solNames,0,mkChar("x"));
	SET_STRING_ELT(solNames,1,mkChar("status"));
	setAttrib(solution,R_NamesSymbol,solNames);

	UNPROTECT(4);
	return solution;

}


