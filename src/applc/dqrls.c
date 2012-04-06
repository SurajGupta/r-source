/* ../appl/dqrls.f -- translated by f2c (version of 1 June 1993  23:00:00).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1110 = 1110;

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

/* Subroutine */ int dqrls_(doublereal *x, integer *n, integer *p, doublereal 
	*y, integer *ny, doublereal *tol, doublereal *b, doublereal *rsd, 
	doublereal *qty, integer *k, integer *jpvt, doublereal *qraux, 
	doublereal *work)
{
    /* System generated locals */
    integer x_dim1, x_offset, y_dim1, y_offset, b_dim1, b_offset, rsd_dim1, 
	    rsd_offset, qty_dim1, qty_offset, i__1, i__2;

    /* Local variables */
    static integer info, j;
    extern /* Subroutine */ int dqrsl_(doublereal *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, integer *, integer *), 
	    dqrdc2_(doublereal *, integer *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *);
    static integer jj, kk;


/*     DQRFIT IS A SUBROUTINE TO COMPUTE LEAST SQUARES SOLUTIONS */
/*     TO THE SYSTEM */

/*     (1)               X * B = Y */

/*     WHICH MAY BE EITHER UNDER-DETERMINED OR OVER-DETERMINED. */
/*     THE USER MUST SUPPLY A TOLERANCE TO LIMIT THE COLUMNS OF */
/*     X USED IN COMPUTING THE SOLUTION.  IN EFFECT, A SET OF */
/*     COLUMNS WITH A CONDITION NUMBER APPROXIMATELY BOUNDED BY */
/*     1/TOL IS USED, THE OTHER COMPONENTS OF B BEING SET TO ZERO. */

/*     ON ENTRY */

/*        X      DOUBLE PRECISION(N,P). */
/*               X CONTAINS N-BY-P COEFFICIENT MATRIX OF */
/*               THE SYSTEM (1), X IS DESTROYED BY DQRFIT. */

/*        N      THE NUMBER OF ROWS OF THE MATRIX X. */

/*        P      THE NUMBER OF COLUMNS OF THE MATRIX X. */

/*        Y      DOUBLE PRECISION(N,NY) */
/*               Y CONTAINS THE RIGHT HAND SIDE(S) OF THE SYSTEM (1). */

/*        NY     THE NUMBER OF RIGHT HAND SIDES OF THE SYSTEM (1). */

/*        TOL    DOUBLE PRECISION */
/*               TOL IS THE NONNEGATIVE TOLERANCE USED TO */
/*               DETERMINE THE SUBSET OF COLUMNS OF X INCLUDED */
/*               IN THE SOLUTION.  COLUMNS ARE PIVOTED OUT OF */
/*               DECOMPOSITION IF */

/*        JPVT   INTEGER(P) */
/*               THE VALUES IN JPVT ARE PERMUTED IN THE SAME */
/*               WAY AS THE COLUMNS OF X.  THIS CAN BE USEFUL */
/*               IN UNSCRAMBLING COEFFICIENTS ETC. */

/*        WORK   DOUBLE PRECISION(2*P) */
/*               WORK IS AN ARRAY USED BY DQRDC2 AND DQRSL. */

/*     ON RETURN */

/*        X      CONTAINS THE OUTPUT ARRAY FROM DQRDC2. */
/*               NAMELY THE QR DECOMPOSITION OF X STORED IN */
/*               COMPACT FORM. */

/*        B      DOUBLE PRECISION(P,NY) */
/*               B CONTAINS THE SOLUTION VECTORS WITH ROWS PERMUTED */
/*               IN THE SAME WAY AS THE COLUMNS OF X.  COMPONENTS */
/*               CORRESPONDING TO COLUMNS NOT USED ARE SET TO ZERO. */

/*        RSD    DOUBLE PRECISION(N,NY) */
/*               RSD CONTAINS THE RESIDUAL VECTORS Y-X*B. */

/*        QTY    DOUBLE PRECISION(N,NY)     T */
/*               QTY CONTAINS THE VECTORS  Q Y.   NOTE THAT */
/*               THE INITIAL P ELEMENTS OF THIS VECTOR ARE */
/*               PERMUTED IN THE SAME WAY AS THE COLUMNS OF X. */

/*        K      INTEGER */
/*               K CONTAINS THE NUMBER OF COLUMNS USED IN THE */
/*               SOLUTION. */

/*        JPVT   HAS ITS CONTENTS PERMUTED AS DESCRIBED ABOVE. */

/*        QRAUX  DOUBLE PRECISION(P) */
/*               QRAUX CONTAINS AUXILIARY INFORMATION ON THE */
/*               QR DECOMPOSITION OF X. */


/*     ON RETURN THE ARRAYS X, JPVT AND QRAUX CONTAIN THE */
/*     USUAL OUTPUT FROM DQRDC, SO THAT THE QR DECOMPOSITION */
/*     OF X WITH PIVOTING IS FULLY AVAILABLE TO THE USER. */
/*     IN PARTICULAR, COLUMNS JPVT(1), JPVT(2),...,JPVT(K) */
/*     WERE USED IN THE SOLUTION, AND THE CONDITION NUMBER */
/*     ASSOCIATED WITH THOSE COLUMNS IS ESTIMATED BY */
/*     ABS(X(1,1)/X(K,K)). */

/*     DQRFIT USED THE LINPACK ROUTINES DQRDC AND DQRSL. */

/*     INTERNAL VARIABLES. */


/*     REDUCE X. */

    /* Parameter adjustments */
    --work;
    --qraux;
    --jpvt;
    qty_dim1 = *n;
    qty_offset = qty_dim1 + 1;
    qty -= qty_offset;
    rsd_dim1 = *n;
    rsd_offset = rsd_dim1 + 1;
    rsd -= rsd_offset;
    b_dim1 = *p;
    b_offset = b_dim1 + 1;
    b -= b_offset;
    y_dim1 = *n;
    y_offset = y_dim1 + 1;
    y -= y_offset;
    x_dim1 = *n;
    x_offset = x_dim1 + 1;
    x -= x_offset;

    /* Function Body */
    dqrdc2_(&x[x_offset], n, n, p, tol, k, &qraux[1], &jpvt[1], &work[1]);

/*     SOLVE THE TRUNCATED LEAST SQUARES PROBLEM FOR EACH RHS. */

    if (*k == 0) {
	goto L30;
    }
    i__1 = *ny;
    for (jj = 1; jj <= i__1; ++jj) {
/* L20: */
	dqrsl_(&x[x_offset], n, n, k, &qraux[1], &y[jj * y_dim1 + 1], &rsd[jj 
		* rsd_dim1 + 1], &qty[jj * qty_dim1 + 1], &b[jj * b_dim1 + 1],
		 &rsd[jj * rsd_dim1 + 1], &rsd[jj * rsd_dim1 + 1], &c__1110, &
		info);
    }
L30:

/*     SET THE UNUSED COMPONENTS OF B TO ZERO. */

    kk = *k + 1;
    i__1 = *p;
    for (j = kk; j <= i__1; ++j) {
	i__2 = *ny;
	for (jj = 1; jj <= i__2; ++jj) {
	    b[j + jj * b_dim1] = 0.f;
/* L40: */
	}
/* L50: */
    }
    return 0;
} /* dqrls_ */

