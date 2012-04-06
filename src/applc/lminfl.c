/* ../appl/lminfl.f -- translated by f2c (version of 1 June 1993  23:00:00).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__10000 = 10000;
static integer c__1000 = 1000;
static integer c__1 = 1;

/* ----------------------------------------------------------------------- */

/*  R : A COMPUTER LANGAGE FOR STATISTICAL DATA ANALYSIS */
/*  COPYRIGHT (C) 1996  ROBERT GENTLEMAN AND ROSS IHAKA */

/*  THIS PROGRAM IS FREE SOFTWARE; YOU CAN REDISTRIBUTE IT AND/OR MODIFY */
/*  IT UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE AS PUBLISHED BY */
/*  THE FREE SOFTWARE FOUNDATION; EITHER VERSION 2 OF THE LICENSE, OR */
/*  (AT YOUR OPTION) ANY LATER VERSION. */

/*  THIS PROGRAM IS DISTRIBUTED IN THE HOPE THAT IT WILL BE USEFUL, */
/*  BUT WITHOUT ANY WARRANTY; WITHOUT EVEN THE IMPLIED WARRANTY OF */
/*  MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.  SEE THE */
/*  GNU GENERAL PUBLIC LICENSE FOR MORE DETAILS. */

/*  YOU SHOULD HAVE RECEIVED A COPY OF THE GNU GENERAL PUBLIC LICENSE */
/*  ALONG WITH THIS PROGRAM; IF NOT, WRITE TO THE FREE SOFTWARE */
/*  FOUNDATION, INC., 675 MASS AVE, CAMBRIDGE, MA 02139, USA. */

/* ----------------------------------------------------------------------- */

/*     LMINFL COMPUTES BASIC QUANTITIES USEFUL FOR COMPUTING */
/*     REGRESSION DIAGNOSTICS. */

/*     ON ENTRY */

/*         X         DOUBLE PRECISION(LDX,K) */
/*                   THE QR DECOMPOSITION AS COMPUTED BY DQRDC OR DQRDC2. */

/*         LDX       INTEGER */
/*                   THE LEADING DIMENSION OF THE ARRAY X. */

/*         N         INTEGER */
/*                   THE NUMBER OF ROWS OF THE MATRIX X. */

/*         K         INTEGER */
/*                   THE NUMBER OF COLUMNS IN THE MATRIX K. */

/*         QRAUX     DOUBLE PRECISION(K) */
/*                   AUXILIARY INFORMATION ABOUT THE QR DECOMPOSITION. */

/*         B         DOUBLE PRECISION(K) */
/*                   THE LEAST-SQUARES PARAMETER ESTIMATES. */

/*         RESID     DOUBLE PRECISION(K) */
/*                   THE RESIDUALS FROM THE REGRESSION. */

/*     ON RETURN */

/*         HAT       DOUBLE PRECISION(N) */
/*                   THE DIAGONAL OF THE HAT MATRIX. */

/*         COEF      DOUBLE PRECISION(N,P) */
/*                   A MATRIX WHICH HAS AS I-TH ROW CONTAINS THE ESTIMATED */
/*                   REGRESSION COEFFICIENTS WHEN THE I-TH CASE IS OMITTED */
/*                   FROM THE REGRESSION. */

/*         SIGMA     DOUBLE PRECISION(N) */
/*                   THE I-TH ELEMENT OF SIGMA CONTAINS AN ESTIMATE */
/*                   OF THE RESIDUAL STANDARD DEVIATION FOR THE MODEL WITH */
/*                   THE I-TH CASE OMITTED. */

/*     THIS VERSION DATED AUG 24, 1996. */
/*     ROSS IHAKA, UNIVERSITY OF AUCKLAND. */

/* Subroutine */ int lminfl_(doublereal *x, integer *ldx, integer *n, integer 
	*k, doublereal *qraux, doublereal *b, doublereal *resid, doublereal *
	hat, doublereal *coef, doublereal *sigma)
{
    /* System generated locals */
    integer x_dim1, x_offset, coef_dim1, coef_offset, i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer info, i, j;
    static doublereal denom;
    extern /* Subroutine */ int dqrsl_(doublereal *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, integer *, integer *), 
	    dtrsl_(doublereal *, integer *, integer *, doublereal *, integer *
	    , integer *);
    static doublereal dummy, sum;



/*     HAT MATRIX DIAGONAL */

    /* Parameter adjustments */
    --sigma;
    coef_dim1 = *n;
    coef_offset = coef_dim1 + 1;
    coef -= coef_offset;
    --hat;
    --resid;
    --b;
    --qraux;
    x_dim1 = *ldx;
    x_offset = x_dim1 + 1;
    x -= x_offset;

    /* Function Body */
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	hat[i] = 0.f;
/* L10: */
    }
    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i = 1; i <= i__2; ++i) {
	    sigma[i] = 0.f;
/* L20: */
	}
	sigma[j] = 1.f;
	dqrsl_(&x[x_offset], ldx, n, k, &qraux[1], &sigma[1], &sigma[1], &
		dummy, &dummy, &dummy, &dummy, &c__10000, &info);
	i__2 = *n;
	for (i = 1; i <= i__2; ++i) {
	    hat[i] += sigma[i] * sigma[i];
/* L30: */
	}
/* L40: */
    }

/*     CHANGES IN THE ESTIMATED COEFFICIENTS */

    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    sigma[j] = 0.f;
/* L50: */
	}
	sigma[i] = resid[i] / (1.f - hat[i]);
	dqrsl_(&x[x_offset], ldx, n, k, &qraux[1], &sigma[1], &dummy, &sigma[
		1], &dummy, &dummy, &dummy, &c__1000, &info);
	dtrsl_(&x[x_offset], ldx, k, &sigma[1], &c__1, &info);
	i__2 = *k;
	for (j = 1; j <= i__2; ++j) {
	    coef[i + j * coef_dim1] = sigma[j];
/* L60: */
	}
/* L70: */
    }

/*     ESTIMATED RESIDUAL STANDARD DEVIATION */

    denom = (doublereal) (*n - *k - 1);
    sum = 0.;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	sum += resid[i] * resid[i];
/* L80: */
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	sigma[i] = sqrt((sum - resid[i] * resid[i] / (1.f - hat[i])) / denom);
/* L90: */
    }
    return 0;
} /* lminfl_ */

