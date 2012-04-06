/* ../appl/dqrsl.f -- translated by f2c (version of 1 June 1993  23:00:00).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int dqrsl_(doublereal *x, integer *ldx, integer *n, integer *
	k, doublereal *qraux, doublereal *y, doublereal *qy, doublereal *qty, 
	doublereal *b, doublereal *rsd, doublereal *xb, integer *job, integer 
	*info)
{
    /* System generated locals */
    integer x_dim1, x_offset, i__1, i__2;

    /* Local variables */
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal temp;
    static logical cqty;
    static integer i, j;
    static doublereal t;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static logical cb;
    static integer jj;
    static logical cr;
    static integer ju, kp1;
    static logical cxb, cqy;


/*     DQRSL APPLIES THE OUTPUT OF DQRDC TO COMPUTE COORDINATE */
/*     TRANSFORMATIONS, PROJECTIONS, AND LEAST SQUARES SOLUTIONS. */
/*     FOR K .LE. MIN(N,P), LET XK BE THE MATRIX */

/*            XK = (X(JPVT(1)),X(JPVT(2)), ... ,X(JPVT(K))) */

/*     FORMED FROM COLUMNNS JPVT(1), ... ,JPVT(K) OF THE ORIGINAL */
/*     N X P MATRIX X THAT WAS INPUT TO DQRDC (IF NO PIVOTING WAS */
/*     DONE, XK CONSISTS OF THE FIRST K COLUMNS OF X IN THEIR */
/*     ORIGINAL ORDER).  DQRDC PRODUCES A FACTORED ORTHOGONAL MATRIX Q */
/*     AND AN UPPER TRIANGULAR MATRIX R SUCH THAT */

/*              XK = Q * (R) */
/*                       (0) */

/*     THIS INFORMATION IS CONTAINED IN CODED FORM IN THE ARRAYS */
/*     X AND QRAUX. */

/*     ON ENTRY */

/*        X      DOUBLE PRECISION(LDX,P). */
/*               X CONTAINS THE OUTPUT OF DQRDC. */

/*        LDX    INTEGER. */
/*               LDX IS THE LEADING DIMENSION OF THE ARRAY X. */

/*        N      INTEGER. */
/*               N IS THE NUMBER OF ROWS OF THE MATRIX XK.  IT MUST */
/*               HAVE THE SAME VALUE AS N IN DQRDC. */

/*        K      INTEGER. */
/*               K IS THE NUMBER OF COLUMNS OF THE MATRIX XK.  K */
/*               MUST NNOT BE GREATER THAN MIN(N,P), WHERE P IS THE */
/*               SAME AS IN THE CALLING SEQUENCE TO DQRDC. */

/*        QRAUX  DOUBLE PRECISION(P). */
/*               QRAUX CONTAINS THE AUXILIARY OUTPUT FROM DQRDC. */

/*        Y      DOUBLE PRECISION(N) */
/*               Y CONTAINS AN N-VECTOR THAT IS TO BE MANIPULATED */
/*               BY DQRSL. */

/*        JOB    INTEGER. */
/*               JOB SPECIFIES WHAT IS TO BE COMPUTED.  JOB HAS */
/*               THE DECIMAL EXPANSION ABCDE, WITH THE FOLLOWING */
/*               MEANING. */

/*                    IF A.NE.0, COMPUTE QY. */
/*                    IF B,C,D, OR E .NE. 0, COMPUTE QTY. */
/*                    IF C.NE.0, COMPUTE B. */
/*                    IF D.NE.0, COMPUTE RSD. */
/*                    IF E.NE.0, COMPUTE XB. */

/*               NOTE THAT A REQUEST TO COMPUTE B, RSD, OR XB */
/*               AUTOMATICALLY TRIGGERS THE COMPUTATION OF QTY, FOR */
/*               WHICH AN ARRAY MUST BE PROVIDED IN THE CALLING */
/*               SEQUENCE. */

/*     ON RETURN */

/*        QY     DOUBLE PRECISION(N). */
/*               QY CONNTAINS Q*Y, IF ITS COMPUTATION HAS BEEN */
/*               REQUESTED. */

/*        QTY    DOUBLE PRECISION(N). */
/*               QTY CONTAINS TRANS(Q)*Y, IF ITS COMPUTATION HAS */
/*               BEEN REQUESTED.  HERE TRANS(Q) IS THE */
/*               TRANSPOSE OF THE MATRIX Q. */

/*        B      DOUBLE PRECISION(K) */
/*               B CONTAINS THE SOLUTION OF THE LEAST SQUARES PROBLEM */

/*                    MINIMIZE NORM2(Y - XK*B), */

/*               IF ITS COMPUTATION HAS BEEN REQUESTED.  (NOTE THAT */
/*               IF PIVOTING WAS REQUESTED IN DQRDC, THE J-TH */
/*               COMPONENT OF B WILL BE ASSOCIATED WITH COLUMN JPVT(J) */
/*               OF THE ORIGINAL MATRIX X THAT WAS INPUT INTO DQRDC.) */

/*        RSD    DOUBLE PRECISION(N). */
/*               RSD CONTAINS THE LEAST SQUARES RESIDUAL Y - XK*B, */
/*               IF ITS COMPUTATION HAS BEEN REQUESTED.  RSD IS */
/*               ALSO THE ORTHOGONAL PROJECTION OF Y ONTO THE */
/*               ORTHOGONAL COMPLEMENT OF THE COLUMN SPACE OF XK. */

/*        XB     DOUBLE PRECISION(N). */
/*               XB CONTAINS THE LEAST SQUARES APPROXIMATION XK*B, */
/*               IF ITS COMPUTATION HAS BEEN REQUESTED.  XB IS ALSO */
/*               THE ORTHOGONAL PROJECTION OF Y ONTO THE COLUMN SPACE */
/*               OF X. */

/*        INFO   INTEGER. */
/*               INFO IS ZERO UNLESS THE COMPUTATION OF B HAS */
/*               BEEN REQUESTED AND R IS EXACTLY SINGULAR.  IN */
/*               THIS CASE, INFO IS THE INDEX OF THE FIRST ZERO */
/*               DIAGONAL ELEMENT OF R AND B IS LEFT UNALTERED. */

/*     THE PARAMETERS QY, QTY, B, RSD, AND XB ARE NOT REFERENCED */
/*     IF THEIR COMPUTATION IS NOT REQUESTED AND IN THIS CASE */
/*     CAN BE REPLACED BY DUMMY VARIABLES IN THE CALLING PROGRAM. */
/*     TO SAVE STORAGE, THE USER MAY IN SOME CASES USE THE SAME */
/*     ARRAY FOR DIFFERENT PARAMETERS IN THE CALLING SEQUENCE.  A */
/*     FREQUENTLY OCCURING EXAMPLE IS WHEN ONE WISHES TO COMPUTE */
/*     ANY OF B, RSD, OR XB AND DOES NOT NEED Y OR QTY.  IN THIS */
/*     CASE ONE MAY IDENTIFY Y, QTY, AND ONE OF B, RSD, OR XB, WHILE */
/*     PROVIDING SEPARATE ARRAYS FOR ANYTHING ELSE THAT IS TO BE */
/*     COMPUTED.  THUS THE CALLING SEQUENCE */

/*          CALL DQRSL(X,LDX,N,K,QRAUX,Y,DUM,Y,B,Y,DUM,110,INFO) */

/*     WILL RESULT IN THE COMPUTATION OF B AND RSD, WITH RSD */
/*     OVERWRITING Y.  MORE GENERALLY, EACH ITEM IN THE FOLLOWING */
/*     LIST CONTAINS GROUPS OF PERMISSIBLE IDENTIFICATIONS FOR */
/*     A SINGLE CALLINNG SEQUENCE. */

/*          1. (Y,QTY,B) (RSD) (XB) (QY) */

/*          2. (Y,QTY,RSD) (B) (XB) (QY) */

/*          3. (Y,QTY,XB) (B) (RSD) (QY) */

/*          4. (Y,QY) (QTY,B) (RSD) (XB) */

/*          5. (Y,QY) (QTY,RSD) (B) (XB) */

/*          6. (Y,QY) (QTY,XB) (B) (RSD) */

/*     IN ANY GROUP THE VALUE RETURNED IN THE ARRAY ALLOCATED TO */
/*     THE GROUP CORRESPONDS TO THE LAST MEMBER OF THE GROUP. */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB. */

/*     DQRSL USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS. */

/*     BLAS DAXPY,DCOPY,DDOT */
/*     FORTRAN DABS,MIN0,MOD */

/*     INTERNAL VARIABLES */



/*     SET INFO FLAG. */

    /* Parameter adjustments */
    --xb;
    --rsd;
    --b;
    --qty;
    --qy;
    --y;
    --qraux;
    x_dim1 = *ldx;
    x_offset = x_dim1 + 1;
    x -= x_offset;

    /* Function Body */
    *info = 0;

/*     DETERMINE WHAT IS TO BE COMPUTED. */

    cqy = *job / 10000 != 0;
    cqty = *job % 10000 != 0;
    cb = *job % 1000 / 100 != 0;
    cr = *job % 100 / 10 != 0;
    cxb = *job % 10 != 0;
/* Computing MIN */
    i__1 = *k, i__2 = *n - 1;
    ju = min(i__1,i__2);

/*     SPECIAL ACTION WHEN N=1. */

    if (ju != 0) {
	goto L40;
    }
    if (cqy) {
	qy[1] = y[1];
    }
    if (cqty) {
	qty[1] = y[1];
    }
    if (cxb) {
	xb[1] = y[1];
    }
    if (! cb) {
	goto L30;
    }
    if (x[x_dim1 + 1] != 0.) {
	goto L10;
    }
    *info = 1;
    goto L20;
L10:
    b[1] = y[1] / x[x_dim1 + 1];
L20:
L30:
    if (cr) {
	rsd[1] = 0.;
    }
    goto L250;
L40:

/*        SET UP TO COMPUTE QY OR QTY. */

    if (cqy) {
	dcopy_(n, &y[1], &c__1, &qy[1], &c__1);
    }
    if (cqty) {
	dcopy_(n, &y[1], &c__1, &qty[1], &c__1);
    }
    if (! cqy) {
	goto L70;
    }

/*           COMPUTE QY. */

    i__1 = ju;
    for (jj = 1; jj <= i__1; ++jj) {
	j = ju - jj + 1;
	if (qraux[j] == 0.) {
	    goto L50;
	}
	temp = x[j + j * x_dim1];
	x[j + j * x_dim1] = qraux[j];
	i__2 = *n - j + 1;
	t = -ddot_(&i__2, &x[j + j * x_dim1], &c__1, &qy[j], &c__1) / x[j + j 
		* x_dim1];
	i__2 = *n - j + 1;
	daxpy_(&i__2, &t, &x[j + j * x_dim1], &c__1, &qy[j], &c__1);
	x[j + j * x_dim1] = temp;
L50:
/* L60: */
	;
    }
L70:
    if (! cqty) {
	goto L100;
    }

/*           COMPUTE TRANS(Q)*Y. */

    i__1 = ju;
    for (j = 1; j <= i__1; ++j) {
	if (qraux[j] == 0.) {
	    goto L80;
	}
	temp = x[j + j * x_dim1];
	x[j + j * x_dim1] = qraux[j];
	i__2 = *n - j + 1;
	t = -ddot_(&i__2, &x[j + j * x_dim1], &c__1, &qty[j], &c__1) / x[j + 
		j * x_dim1];
	i__2 = *n - j + 1;
	daxpy_(&i__2, &t, &x[j + j * x_dim1], &c__1, &qty[j], &c__1);
	x[j + j * x_dim1] = temp;
L80:
/* L90: */
	;
    }
L100:

/*        SET UP TO COMPUTE B, RSD, OR XB. */

    if (cb) {
	dcopy_(k, &qty[1], &c__1, &b[1], &c__1);
    }
    kp1 = *k + 1;
    if (cxb) {
	dcopy_(k, &qty[1], &c__1, &xb[1], &c__1);
    }
    if (cr && *k < *n) {
	i__1 = *n - *k;
	dcopy_(&i__1, &qty[kp1], &c__1, &rsd[kp1], &c__1);
    }
    if (! cxb || kp1 > *n) {
	goto L120;
    }
    i__1 = *n;
    for (i = kp1; i <= i__1; ++i) {
	xb[i] = 0.;
/* L110: */
    }
L120:
    if (! cr) {
	goto L140;
    }
    i__1 = *k;
    for (i = 1; i <= i__1; ++i) {
	rsd[i] = 0.;
/* L130: */
    }
L140:
    if (! cb) {
	goto L190;
    }

/*           COMPUTE B. */

    i__1 = *k;
    for (jj = 1; jj <= i__1; ++jj) {
	j = *k - jj + 1;
	if (x[j + j * x_dim1] != 0.) {
	    goto L150;
	}
	*info = j;
/*           ......EXIT */
	goto L180;
L150:
	b[j] /= x[j + j * x_dim1];
	if (j == 1) {
	    goto L160;
	}
	t = -b[j];
	i__2 = j - 1;
	daxpy_(&i__2, &t, &x[j * x_dim1 + 1], &c__1, &b[1], &c__1);
L160:
/* L170: */
	;
    }
L180:
L190:
    if (! cr && ! cxb) {
	goto L240;
    }

/*           COMPUTE RSD OR XB AS REQUIRED. */

    i__1 = ju;
    for (jj = 1; jj <= i__1; ++jj) {
	j = ju - jj + 1;
	if (qraux[j] == 0.) {
	    goto L220;
	}
	temp = x[j + j * x_dim1];
	x[j + j * x_dim1] = qraux[j];
	if (! cr) {
	    goto L200;
	}
	i__2 = *n - j + 1;
	t = -ddot_(&i__2, &x[j + j * x_dim1], &c__1, &rsd[j], &c__1) / x[j + 
		j * x_dim1];
	i__2 = *n - j + 1;
	daxpy_(&i__2, &t, &x[j + j * x_dim1], &c__1, &rsd[j], &c__1);
L200:
	if (! cxb) {
	    goto L210;
	}
	i__2 = *n - j + 1;
	t = -ddot_(&i__2, &x[j + j * x_dim1], &c__1, &xb[j], &c__1) / x[j + j 
		* x_dim1];
	i__2 = *n - j + 1;
	daxpy_(&i__2, &t, &x[j + j * x_dim1], &c__1, &xb[j], &c__1);
L210:
	x[j + j * x_dim1] = temp;
L220:
/* L230: */
	;
    }
L240:
L250:
    return 0;
} /* dqrsl_ */

