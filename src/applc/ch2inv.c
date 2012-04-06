/* ../appl/ch2inv.f -- translated by f2c (version of 1 June 1993  23:00:00).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

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

/*     CH2INV COMPUTES THE INVERSE OF A POSTIVE-DEFINITE SYMMETRIC */
/*     MATRIX FROM ITS CHOLESKI FACTORIZATION.  THIS CAN BE USED FOR */
/*     EXAMPLE TO COMPUTE THE DISPERSION MATRIX FOR THE ESTIMATED */
/*     PARAMETERS IN A REGRESSION ANALYSIS. */

/*     ON ENTRY */

/*         X         DOUBLE PRECISION(LDX,K) */
/*                   THE CHOLESKI DECOMPOSITION OR THE */
/*                   QR DECOMPOSITION AS COMPUTED BY DQRDC */
/*                   OR DQRDC2 */

/*         LDX       INTEGER */
/*                   THE LEADING DIMENSION OF THE ARRAY X */

/*         N         INTEGER */
/*                   THE NUMBER OF ROWS OF THE MATRIX X */

/*         K         INTEGER */
/*                   THE NUMBER OF COLUMNS IN THE MATRIX K */

/*     ON RETURN */

/*         V         DOUBLE PRECISION(K,K) */
/*                   THE VALUE OF INVERSE(X'X) */

/*     THIS VERSION DATED AUG 24, 1996. */
/*     ROSS IHAKA, UNIVERSITY OF AUCKLAND. */

/* Subroutine */ int ch2inv_(doublereal *x, integer *ldx, integer *n, 
	doublereal *v, integer *info)
{
    /* System generated locals */
    integer x_dim1, x_offset, v_dim1, v_offset, i__1, i__2;

    /* Local variables */
    static doublereal d;
    static integer i, j;
    extern /* Subroutine */ int dpodi_(doublereal *, integer *, integer *, 
	    doublereal *, integer *);
    static integer im1;



    /* Parameter adjustments */
    v_dim1 = *n;
    v_offset = v_dim1 + 1;
    v -= v_offset;
    x_dim1 = *ldx;
    x_offset = x_dim1 + 1;
    x -= x_offset;

    /* Function Body */
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	if (x[i + i * x_dim1] == 0.) {
	    *info = i;
	    return 0;
	}
	i__2 = *n;
	for (j = i; j <= i__2; ++j) {
	    v[i + j * v_dim1] = x[i + j * x_dim1];
/* L10: */
	}
/* L20: */
    }
    dpodi_(&v[v_offset], n, n, &d, &c__1);
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	im1 = i - 1;
	i__2 = im1;
	for (j = 1; j <= i__2; ++j) {
	    v[i + j * v_dim1] = v[j + i * v_dim1];
/* L30: */
	}
/* L40: */
    }
    return 0;
} /* ch2inv_ */

