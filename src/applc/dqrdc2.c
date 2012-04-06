/* ../appl/dqrdc2.f -- translated by f2c (version of 1 June 1993  23:00:00).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int dqrdc2_(doublereal *x, integer *ldx, integer *n, integer 
	*p, doublereal *tol, integer *k, doublereal *qraux, integer *jpvt, 
	doublereal *work)
{
    /* System generated locals */
    integer x_dim1, x_offset, work_dim1, work_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *), sqrt(doublereal);

    /* Local variables */
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *), dnrm2_(integer *, doublereal *, integer *);
    static integer i, j, l;
    static doublereal t;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), daxpy_(integer *, doublereal *, doublereal *, integer 
	    *, doublereal *, integer *);
    static doublereal nrmxl, tt;
    static integer lp1, lup;
    static doublereal ttt;


/*     DQRDC2 USES HOUSEHOLDER TRANSFORMATIONS TO COMPUTE THE QR */
/*     FACTORIZATION OF AN N BY P MATRIX X.  A LIMITED COLUMN */
/*     PIVOTING STRATEGY BASED ON THE 2-NORMS OF THE REDUCED COLUMNS */
/*     MOVES COLUMNS WITH NEAR-ZERO NORM TO THE RIGHT-HAND EDGE OF */
/*     THE X MATRIX.  THIS STRATEGY MEANS THAT SEQUENTIAL ONE */
/*     DEGREE-OF-FREEDOM EFFECTS CAN BE COMPUTED IN A NATURAL WAY. */

/*     I AM VERY NERVOUS ABOUT MODIFYING LINPACK CODE IN THIS WAY. */
/*     IF YOU ARE A COMPUTATIONAL LINEAR ALGEBRA GURU AND YOU REALLY */
/*     UNDERSTAND HOW TO SOLVE THIS PROBLEM PLEASE FEEL FREE TO */
/*     SUGGEST IMPROVEMENTS TO THIS CODE. */

/*     ON ENTRY */

/*        X       DOUBLE PRECISION(LDX,P), WHERE LDX .GE. N. */
/*                X CONTAINS THE MATRIX WHOSE DECOMPOSITION IS TO BE */
/*                COMPUTED. */

/*        LDX     INTEGER. */
/*                LDX IS THE LEADING DIMENSION OF THE ARRAY X. */

/*        N       INTEGER. */
/*                N IS THE NUMBER OF ROWS OF THE MATRIX X. */

/*        P       INTEGER. */
/*                P IS THE NUMBER OF COLUMNS OF THE MATRIX X. */

/*        TOL     DOUBLE PRECISION */
/*                TOL IS THE NONNEGATIVE TOLERANCE USED TO */
/*                DETERMINE THE SUBSET OF THE COLUMNS OF X */
/*                INCLUDED IN THE SOLUTION. */

/*        JPVT    INTEGER(P). */
/*                INTEGERS WHICH ARE SWAPPED IN THE SAME WAY AS THE */
/*                THE COLUMNS OF X DURING PIVOTING.  ON ENTRY THESE */
/*                SHOULD BE SET EQUAL TO THE COLUMN INDICES OF THE */
/*                COLUMNS OF THE X MATRIX (TYPICALLY 1 to P). */

/*        WORK    DOUBLE PRECISION(P,2). */
/*                WORK IS A WORK ARRAY. */

/*     ON RETURN */

/*        X       X CONTAINS IN ITS UPPER TRIANGLE THE UPPER */
/*                TRIANGULAR MATRIX R OF THE QR FACTORIZATION. */
/*                BELOW ITS DIAGONAL X CONTAINS INFORMATION FROM */
/*                WHICH THE ORTHOGONAL PART OF THE DECOMPOSITION */
/*                CAN BE RECOVERED.  NOTE THAT IF PIVOTING HAS */
/*                BEEN REQUESTED, THE DECOMPOSITION IS NOT THAT */
/*                OF THE ORIGINAL MATRIX X BUT THAT OF X */
/*                WITH ITS COLUMNS PERMUTED AS DESCRIBED BY JPVT. */

/*        K       INTEGER. */
/*                K CONTAINS THE NUMBER OF COLUMNS OF X JUDGED */
/*                TO BE LINEARLY INDEPENDENT. */

/*        QRAUX   DOUBLE PRECISION(P). */
/*                QRAUX CONTAINS FURTHER INFORMATION REQUIRED TO RECOVER 
*/
/*                THE ORTHOGONAL PART OF THE DECOMPOSITION. */

/*        JPVT    JPVT(K) CONTAINS THE INDEX OF THE COLUMN OF THE */
/*                ORIGINAL MATRIX THAT HAS BEEN INTERCHANGED INTO */
/*                THE K-TH COLUMN, IF PIVOTING WAS REQUESTED. */

/*     THIS VERSION DATED 22 AUGUST 1995 */
/*     ROSS IHAKA */

/*     DQRDC USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS. */

/*     BLAS DAXPY,DDOT,DSCAL,DNRM2 */
/*     FORTRAN DABS,DMAX1,MIN0,DSQRT */

/*     INTERNAL VARIABLES */



/*     COMPUTE THE NORMS OF THE COLUMNS OF X. */

    /* Parameter adjustments */
    work_dim1 = *p;
    work_offset = work_dim1 + 1;
    work -= work_offset;
    --jpvt;
    --qraux;
    x_dim1 = *ldx;
    x_offset = x_dim1 + 1;
    x -= x_offset;

    /* Function Body */
    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
	qraux[j] = dnrm2_(n, &x[j * x_dim1 + 1], &c__1);
	work[j + work_dim1] = qraux[j];
	work[j + (work_dim1 << 1)] = qraux[j];
	if (work[j + (work_dim1 << 1)] == 0.) {
	    work[j + (work_dim1 << 1)] = 1.;
	}
/* L70: */
    }

/*     PERFORM THE HOUSEHOLDER REDUCTION OF X. */

    lup = min(*n,*p);
    *k = lup + 1;
    i__1 = lup;
    for (l = 1; l <= i__1; ++l) {

/*     CYCLE THE COLUMNS FROM L TO LUP LEFT-TO-RIGHT UNTIL ONE */
/*     WITH NON-NEGLIGIBLE NORM IS LOCATED.  A COLUMN IS CONSIDERED */
/*     TO HAVE BECOME NEGLIGIBLE IF ITS NORM HAS FALLEN BELOW */
/*     TOL TIMES ITS ORIGINAL NORM.  THE CHECK FOR L .LE. K */
/*     AVOIDS INFINITE CYCLING. */

L80:
	if (l >= *k || qraux[l] >= work[l + (work_dim1 << 1)] * *tol) {
	    goto L120;
	}
	lp1 = l + 1;
	i__2 = *n;
	for (i = 1; i <= i__2; ++i) {
	    t = x[i + l * x_dim1];
	    i__3 = lup;
	    for (j = lp1; j <= i__3; ++j) {
		x[i + (j - 1) * x_dim1] = x[i + j * x_dim1];
/* L90: */
	    }
	    x[i + lup * x_dim1] = t;
/* L100: */
	}
	i = jpvt[l];
	t = qraux[l];
	tt = work[l + work_dim1];
	ttt = work[l + (work_dim1 << 1)];
	i__2 = lup;
	for (j = lp1; j <= i__2; ++j) {
	    jpvt[j - 1] = jpvt[j];
	    qraux[j - 1] = qraux[j];
	    work[j - 1 + work_dim1] = work[j + work_dim1];
	    work[j - 1 + (work_dim1 << 1)] = work[j + (work_dim1 << 1)];
/* L110: */
	}
	jpvt[lup] = i;
	qraux[lup] = t;
	work[lup + work_dim1] = tt;
	work[lup + (work_dim1 << 1)] = ttt;
	--(*k);
	goto L80;
L120:
	if (l == *n) {
	    goto L190;
	}

/*           COMPUTE THE HOUSEHOLDER TRANSFORMATION FOR COLUMN L. */

	i__2 = *n - l + 1;
	nrmxl = dnrm2_(&i__2, &x[l + l * x_dim1], &c__1);
	if (nrmxl == 0.) {
	    goto L180;
	}
	if (x[l + l * x_dim1] != 0.) {
	    nrmxl = d_sign(&nrmxl, &x[l + l * x_dim1]);
	}
	i__2 = *n - l + 1;
	d__1 = 1. / nrmxl;
	dscal_(&i__2, &d__1, &x[l + l * x_dim1], &c__1);
	x[l + l * x_dim1] += 1.;

/*              APPLY THE TRANSFORMATION TO THE REMAINING COLUMNS, */
/*              UPDATING THE NORMS. */

	lp1 = l + 1;
	if (*p < lp1) {
	    goto L170;
	}
	i__2 = *p;
	for (j = lp1; j <= i__2; ++j) {
	    i__3 = *n - l + 1;
	    t = -ddot_(&i__3, &x[l + l * x_dim1], &c__1, &x[l + j * x_dim1], &
		    c__1) / x[l + l * x_dim1];
	    i__3 = *n - l + 1;
	    daxpy_(&i__3, &t, &x[l + l * x_dim1], &c__1, &x[l + j * x_dim1], &
		    c__1);
	    if (qraux[j] == 0.) {
		goto L150;
	    }
/* Computing 2nd power */
	    d__2 = (d__1 = x[l + j * x_dim1], abs(d__1)) / qraux[j];
	    tt = 1. - d__2 * d__2;
	    tt = max(tt,0.);
	    t = tt;
/* Computing 2nd power */
	    d__1 = qraux[j] / work[j + work_dim1];
	    tt = tt * .05 * (d__1 * d__1) + 1.;
	    if (tt == 1.) {
		goto L130;
	    }
	    qraux[j] *= sqrt(t);
	    goto L140;
L130:
	    i__3 = *n - l;
	    qraux[j] = dnrm2_(&i__3, &x[l + 1 + j * x_dim1], &c__1);
	    work[j + work_dim1] = qraux[j];
L140:
L150:
/* L160: */
	    ;
	}
L170:

/*              SAVE THE TRANSFORMATION. */

	qraux[l] = x[l + l * x_dim1];
	x[l + l * x_dim1] = -nrmxl;
L180:
L190:
/* L200: */
	;
    }
    --(*k);
    return 0;
} /* dqrdc2_ */

