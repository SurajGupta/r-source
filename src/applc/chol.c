/* ../appl/chol.f -- translated by f2c (version of 1 June 1993  23:00:00).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

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

/*     CHOL PERFORMS THE CHOLESKI DECOMPOSITION OF A SYMMETRIC */
/*     POSITIVE-DEFINITE MATRIX.  THIS IS JUST A WRAPPER FOR THE */
/*     LINPACK ROUTINE DPOFA. */

/*     ON ENTRY */

/*         A         DOUBLE PRECISION(LDA,N) */
/*                   THE UPPER TRIANGLE OF THE MATRIX TO BE FACTORIZED */
/*                   IS CONTAINED IN THE UPPER TRIANGLE OF A. */

/*         LDA       INTEGER */
/*                   THE LEADING DIMENSION OF A. */

/*         N         INTEGER */
/*                   THE NUMBER OR ROWS AND COLUMNS OF THE MATRIX */
/*                   TO BE FACTORIZED. */

/*     ON RETURN */

/*         V         DOUBLE PRECISION(N,N) */
/*                   THE SQUARE-ROOT (CHOLESKI) FACTOR. */

/*         INFO      INTEGER */
/*                   THE ERROR INDICATOR FROM DPOFA.  THIS WILL BE */
/*                   ZERO UNLESS THE MATRIX BEING FACTORIZED IS */
/*                   NOT POSITIVE DEFINITE. */

/*     THIS VERSION DATED AUG 25, 1996. */
/*     ROSS IHAKA, UNIVERSITY OF AUCKLAND. */

/* Subroutine */ int chol_(doublereal *a, integer *lda, integer *n, 
	doublereal *v, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, v_dim1, v_offset, i__1, i__2;

    /* Local variables */
    static integer i, j;
    extern /* Subroutine */ int dpofa_(doublereal *, integer *, integer *, 
	    integer *);


    /* Parameter adjustments */
    v_dim1 = *n;
    v_offset = v_dim1 + 1;
    v -= v_offset;
    a_dim1 = *lda;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	i__2 = *n;
	for (j = 1; j <= i__2; ++j) {
	    if (i > j) {
		v[i + j * v_dim1] = 0.;
	    } else {
		v[i + j * v_dim1] = a[i + j * a_dim1];
	    }
/* L10: */
	}
/* L20: */
    }
    dpofa_(&v[v_offset], n, n, info);
    return 0;
} /* chol_ */

