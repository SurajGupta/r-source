/* ../appl/uncmin.f -- translated by f2c (version of 1 June 1993  23:00:00).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static doublereal c_b2 = .33333333333333331;
static doublereal c_b3 = 10.;
static doublereal c_b33 = 0.;
static integer c__4 = 4;
static integer c__2 = 2;
static integer c__1 = 1;
static integer c__3 = 3;
static doublereal c_b146 = 1.;
static integer c__0 = 0;

/* Subroutine */ int fdhess_(integer *n, doublereal *x, doublereal *fval, 
	S_fp fun, doublereal *h, integer *nfd, doublereal *step, doublereal *
	f, integer *ndigit, doublereal *typx)
{
    /* System generated locals */
    integer h_dim1, h_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_di(doublereal *, integer *), pow_dd(doublereal *, doublereal *)
	    ;

    /* Local variables */
    static integer i, j;
    static doublereal tempi, tempj, fii, eta, fij;


/* This subroutine calculates a numerical approximation to the upper */
/* triangular portion of the second derivative matrix (the Hessian). */
/* Algorithm A5.6.2 from Dennis and Schnable (1983), Numerical Methods */
/* for Unconstrained Optimization and Nonlinear Equations, */
/* Prentice-Hall, 321-322. */
/* Programmed by Richard H. Jones, January 11, 1989 */

/* Input to subroutine */
/*     n.....The number of parameters */
/*     x.....Vector of parameter values */
/*     fval..Double precision value of function at X */
/*     fun...A function provided by the user which must be declared as */
/*           EXTERNAL in the calling program.  Its call must be of the */
/*           CALL FUN(N,X,FVAL) where FVAL is the computed value of the */
/*           function */
/*     nfd...First dimension of H in the calling program */
/* Output from subroutine */
/*     h.....An n by n matrix of the approximate Hessian */
/* Work space */
/*     step....A real array of length n */
/*     f.......A double precision array of length n */

    /* Parameter adjustments */
    --typx;
    --f;
    --step;
    h_dim1 = *nfd;
    h_offset = h_dim1 + 1;
    h -= h_offset;
    --x;

    /* Function Body */
    i__1 = -(*ndigit);
    d__1 = pow_di(&c_b3, &i__1);
    eta = pow_dd(&d__1, &c_b2);
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
/* Computing MAX */
	d__1 = x[i], d__2 = typx[i];
	step[i] = eta * max(d__1,d__2);
	if (typx[i] < 0.) {
	    step[i] = -step[i];
	}
	tempi = x[i];
	x[i] += step[i];
	step[i] = x[i] - tempi;
	(*fun)(n, &x[1], &f[i]);
	x[i] = tempi;
/* L10: */
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	tempi = x[i];
	x[i] += step[i] * 2.;
	(*fun)(n, &x[1], &fii);
	h[i + i * h_dim1] = (*fval - f[i] + (fii - f[i])) / (step[i] * step[i]
		);
	x[i] = tempi + step[i];
	if (i < *n) {
	    i__2 = *n;
	    for (j = i + 1; j <= i__2; ++j) {
		tempj = x[j];
		x[j] += step[j];
		(*fun)(n, &x[1], &fij);
		h[i + j * h_dim1] = (*fval - f[i] + (fij - f[j])) / (step[i] *
			 step[j]);
		x[j] = tempj;
/* L20: */
	    }
	}
	x[i] = tempi;
/* L30: */
    }
    return 0;
} /* fdhess_ */


/* Subroutine */ int bakslv_(integer *nr, integer *n, doublereal *a, 
	doublereal *x, doublereal *b)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    static integer i, j, ip1;
    static doublereal sum;


/* PURPOSE */
/* ------- */
/* SOLVE  AX=B  WHERE A IS UPPER TRIANGULAR MATRIX. */
/* NOTE THAT A IS INPUT AS A LOWER TRIANGULAR MATRIX AND */
/* THAT THIS ROUTINE TAKES ITS TRANSPOSE IMPLICITLY. */

/* PARAMETERS */
/* ---------- */
/* NR           --> ROW DIMENSION OF MATRIX */
/* N            --> DIMENSION OF PROBLEM */
/* A(N,N)       --> LOWER TRIANGULAR MATRIX (PRESERVED) */
/* X(N)        <--  SOLUTION VECTOR */
/* B(N)         --> RIGHT-HAND SIDE VECTOR */

/* NOTE */
/* ---- */
/* IF B IS NO LONGER REQUIRED BY CALLING ROUTINE, */
/* THEN VECTORS B AND X MAY SHARE THE SAME STORAGE. */


/* SOLVE (L-TRANSPOSE)X=B. (BACK SOLVE) */

    /* Parameter adjustments */
    --b;
    --x;
    a_dim1 = *nr;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    i = *n;
    x[i] = b[i] / a[i + i * a_dim1];
    if (*n == 1) {
	return 0;
    }
L30:
    ip1 = i;
    --i;
    sum = 0.;
    i__1 = *n;
    for (j = ip1; j <= i__1; ++j) {
	sum += a[j + i * a_dim1] * x[j];
/* L40: */
    }
    x[i] = (b[i] - sum) / a[i + i * a_dim1];
    if (i > 1) {
	goto L30;
    }
    return 0;
} /* bakslv_ */

/* Subroutine */ int chlhsn_(integer *nr, integer *n, doublereal *a, 
	doublereal *epsm, doublereal *sx, doublereal *udiag)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i, j;
    static doublereal evmin, evmax;
    extern /* Subroutine */ int choldc_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static doublereal addmax, diagmn, diagmx, offmax, offrow, posmax;
    static integer im1, jm1, ip1, jp1;
    static doublereal sdd, amu, tol;


/* PURPOSE */
/* ------- */
/* FIND THE L(L-TRANSPOSE) [WRITTEN LL+] DECOMPOSITION OF THE PERTURBED */
/* MODEL HESSIAN MATRIX A+MU*I(WHERE MU\0 AND I IS THE IDENTITY MATRIX) */
/* WHICH IS SAFELY POSITIVE DEFINITE.  IF A IS SAFELY POSITIVE DEFINITE */
/* UPON ENTRY, THEN MU=0. */

/* PARAMETERS */
/* ---------- */
/* NR           --> ROW DIMENSION OF MATRIX */
/* N            --> DIMENSION OF PROBLEM */
/* A(N,N)      <--> ON ENTRY; "A" IS MODEL HESSIAN (ONLY LOWER */
/*                  TRIANGULAR PART AND DIAGONAL STORED) */
/*                  ON EXIT:  A CONTAINS L OF LL+ DECOMPOSITION OF */
/*                  PERTURBED MODEL HESSIAN IN LOWER TRIANGULAR */
/*                  PART AND DIAGONAL AND CONTAINS HESSIAN IN UPPER */
/*                  TRIANGULAR PART AND UDIAG */
/* EPSM         --> MACHINE EPSILON */
/* SX(N)        --> DIAGONAL SCALING MATRIX FOR X */
/* UDIAG(N)    <--  ON EXIT: CONTAINS DIAGONAL OF HESSIAN */

/* INTERNAL VARIABLES */
/* ------------------ */
/* TOL              TOLERANCE */
/* DIAGMN           MINIMUM ELEMENT ON DIAGONAL OF A */
/* DIAGMX           MAXIMUM ELEMENT ON DIAGONAL OF A */
/* OFFMAX           MAXIMUM OFF-DIAGONAL ELEMENT OF A */
/* OFFROW           SUM OF OFF-DIAGONAL ELEMENTS IN A ROW OF A */
/* EVMIN            MINIMUM EIGENVALUE OF A */
/* EVMAX            MAXIMUM EIGENVALUE OF A */

/* DESCRIPTION */
/* ----------- */
/* 1. IF "A" HAS ANY NEGATIVE DIAGONAL ELEMENTS, THEN CHOOSE MU>0 */
/* SUCH THAT THE DIAGONAL OF A:=A+MU*I IS ALL POSITIVE */
/* WITH THE RATIO OF ITS SMALLEST TO LARGEST ELEMENT ON THE */
/* ORDER OF SQRT(EPSM). */

/* 2. "A" UNDERGOES A PERTURBED CHOLESKY DECOMPOSITION WHICH */
/* RESULTS IN AN LL+ DECOMPOSITION OF A+D, WHERE D IS A */
/* NON-NEGATIVE DIAGONAL MATRIX WHICH IS IMPLICITLY ADDED TO */
/* "A" DURING THE DECOMPOSITION IF "A" IS NOT POSITIVE DEFINITE. */
/* "A" IS RETAINED AND NOT CHANGED DURING THIS PROCESS BY */
/* COPYING L INTO THE UPPER TRIANGULAR PART OF "A" AND THE */
/* DIAGONAL INTO UDIAG.  THEN THE CHOLESKY DECOMPOSITION ROUTINE */
/* IS CALLED.  ON RETURN, ADDMAX CONTAINS MAXIMUM ELEMENT OF D. */

/* 3. IF ADDMAX=0, "A" WAS POSITIVE DEFINITE GOING INTO STEP 2 */
/* AND RETURN IS MADE TO CALLING PROGRAM.  OTHERWISE, */
/* THE MINIMUM NUMBER SDD WHICH MUST BE ADDED TO THE */
/* DIAGONAL OF A TO MAKE IT SAFELY STRICTLY DIAGONALLY DOMINANT */
/* IS CALCULATED.  SINCE A+ADDMAX*I AND A+SDD*I ARE SAFELY */
/* POSITIVE DEFINITE, CHOOSE MU=MIN(ADDMAX,SDD) AND DECOMPOSE */
/* A+MU*I TO OBTAIN L. */


/* SCALE HESSIAN */
/* PRE- AND POST- MULTIPLY "A" BY INV(SX) */

    /* Parameter adjustments */
    --udiag;
    --sx;
    a_dim1 = *nr;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i = j; i <= i__2; ++i) {
	    a[i + j * a_dim1] /= sx[i] * sx[j];
/* L10: */
	}
/* L20: */
    }

/* STEP1 */
/* ----- */
/* NOTE:  IF A DIFFERENT TOLERANCE IS DESIRED THROUGHOUT THIS */
/* ALGORITHM, CHANGE TOLERANCE HERE: */
    tol = sqrt(*epsm);

    diagmx = a[a_dim1 + 1];
    diagmn = a[a_dim1 + 1];
    if (*n == 1) {
	goto L35;
    }
    i__1 = *n;
    for (i = 2; i <= i__1; ++i) {
	if (a[i + i * a_dim1] < diagmn) {
	    diagmn = a[i + i * a_dim1];
	}
	if (a[i + i * a_dim1] > diagmx) {
	    diagmx = a[i + i * a_dim1];
	}
/* L30: */
    }
L35:
    posmax = max(diagmx,0.);

/* DIAGMN .LE. 0 */

    if (diagmn > posmax * tol) {
	goto L100;
    }
/*     IF(DIAGMN.LE.POSMAX*TOL) */
/*     THEN */
    amu = tol * (posmax - diagmn) - diagmn;
    if (amu != 0.) {
	goto L60;
    }
/*       IF(AMU.EQ.0.0D0) */
/*       THEN */

/* FIND LARGEST OFF-DIAGONAL ELEMENT OF A */
    offmax = 0.;
    if (*n == 1) {
	goto L50;
    }
    i__1 = *n;
    for (i = 2; i <= i__1; ++i) {
	im1 = i - 1;
	i__2 = im1;
	for (j = 1; j <= i__2; ++j) {
	    if ((d__1 = a[i + j * a_dim1], abs(d__1)) > offmax) {
		offmax = (d__2 = a[i + j * a_dim1], abs(d__2));
	    }
/* L40: */
	}
/* L45: */
    }
L50:
    amu = offmax;
    if (amu != 0.) {
	goto L55;
    }
/*         IF(AMU.EQ.0.0D0) */
/*         THEN */
    amu = 1.;
    goto L60;
/*         ELSE */
L55:
    amu *= tol + 1.;
/*         ENDIF */
/*       ENDIF */

/* A=A + MU*I */

L60:
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	a[i + i * a_dim1] += amu;
/* L65: */
    }
    diagmx += amu;
/*     ENDIF */

/* STEP2 */
/* ----- */
/* COPY LOWER TRIANGULAR PART OF "A" TO UPPER TRIANGULAR PART */
/* AND DIAGONAL OF "A" TO UDIAG */

L100:
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	udiag[j] = a[j + j * a_dim1];
	if (j == *n) {
	    goto L110;
	}
	jp1 = j + 1;
	i__2 = *n;
	for (i = jp1; i <= i__2; ++i) {
	    a[j + i * a_dim1] = a[i + j * a_dim1];
/* L105: */
	}
L110:
	;
    }

    choldc_(nr, n, &a[a_offset], &diagmx, &tol, &addmax);


/* STEP3 */
/* ----- */
/* IF ADDMAX=0, "A" WAS POSITIVE DEFINITE GOING INTO STEP 2, */
/* THE LL+ DECOMPOSITION HAS BEEN DONE, AND WE RETURN. */
/* OTHERWISE, ADDMAX>0.  PERTURB "A" SO THAT IT IS SAFELY */
/* DIAGONALLY DOMINANT AND FIND LL+ DECOMPOSITION */

    if (addmax <= 0.) {
	goto L170;
    }
/*     IF(ADDMAX.GT.0.0D0) */
/*     THEN */

/* RESTORE ORIGINAL "A" (LOWER TRIANGULAR PART AND DIAGONAL) */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	a[j + j * a_dim1] = udiag[j];
	if (j == *n) {
	    goto L120;
	}
	jp1 = j + 1;
	i__2 = *n;
	for (i = jp1; i <= i__2; ++i) {
	    a[i + j * a_dim1] = a[j + i * a_dim1];
/* L115: */
	}
L120:
	;
    }

/* FIND SDD SUCH THAT A+SDD*I IS SAFELY POSITIVE DEFINITE */
/* NOTE:  EVMIN<0 SINCE A IS NOT POSITIVE DEFINITE; */

    evmin = 0.;
    evmax = a[a_dim1 + 1];
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	offrow = 0.;
	if (i == 1) {
	    goto L135;
	}
	im1 = i - 1;
	i__2 = im1;
	for (j = 1; j <= i__2; ++j) {
	    offrow += (d__1 = a[i + j * a_dim1], abs(d__1));
/* L130: */
	}
L135:
	if (i == *n) {
	    goto L145;
	}
	ip1 = i + 1;
	i__2 = *n;
	for (j = ip1; j <= i__2; ++j) {
	    offrow += (d__1 = a[j + i * a_dim1], abs(d__1));
/* L140: */
	}
L145:
/* Computing MIN */
	d__1 = evmin, d__2 = a[i + i * a_dim1] - offrow;
	evmin = min(d__1,d__2);
/* Computing MAX */
	d__1 = evmax, d__2 = a[i + i * a_dim1] + offrow;
	evmax = max(d__1,d__2);
/* L150: */
    }
    sdd = tol * (evmax - evmin) - evmin;

/* PERTURB "A" AND DECOMPOSE AGAIN */

    amu = min(sdd,addmax);
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	a[i + i * a_dim1] += amu;
	udiag[i] = a[i + i * a_dim1];
/* L160: */
    }

/* "A" NOW GUARANTEED SAFELY POSITIVE DEFINITE */

    choldc_(nr, n, &a[a_offset], &c_b33, &tol, &addmax);
/*     ENDIF */

/* UNSCALE HESSIAN AND CHOLESKY DECOMPOSITION MATRIX */

L170:
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i = j; i <= i__2; ++i) {
	    a[i + j * a_dim1] = sx[i] * a[i + j * a_dim1];
/* L175: */
	}
	if (j == 1) {
	    goto L185;
	}
	jm1 = j - 1;
	i__2 = jm1;
	for (i = 1; i <= i__2; ++i) {
	    a[i + j * a_dim1] = sx[i] * sx[j] * a[i + j * a_dim1];
/* L180: */
	}
L185:
	udiag[j] = udiag[j] * sx[j] * sx[j];
/* L190: */
    }
    return 0;
} /* chlhsn_ */

/* Subroutine */ int choldc_(integer *nr, integer *n, doublereal *a, 
	doublereal *diagmx, doublereal *tol, doublereal *addmax)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal temp;
    static integer i, j, k;
    static doublereal aminl, offmax, amnlsq;
    static integer jm1, jp1;
    static doublereal sum;


/* PURPOSE */
/* ------- */
/* FIND THE PERTURBED L(L-TRANSPOSE) [WRITTEN LL+] DECOMPOSITION */
/* OF A+D, WHERE D IS A NON-NEGATIVE DIAGONAL MATRIX ADDED TO A IF */
/* NECESSARY TO ALLOW THE CHOLESKY DECOMPOSITION TO CONTINUE. */

/* PARAMETERS */
/* ---------- */
/* NR           --> ROW DIMENSION OF MATRIX */
/* N            --> DIMENSION OF PROBLEM */
/* A(N,N)      <--> ON ENTRY: MATRIX FOR WHICH TO FIND PERTURBED */
/*                       CHOLESKY DECOMPOSITION */
/*                  ON EXIT:  CONTAINS L OF LL+ DECOMPOSITION */
/*                  IN LOWER TRIANGULAR PART AND DIAGONAL OF "A" */
/* DIAGMX       --> MAXIMUM DIAGONAL ELEMENT OF "A" */
/* TOL          --> TOLERANCE */
/* ADDMAX      <--  MAXIMUM AMOUNT IMPLICITLY ADDED TO DIAGONAL OF "A" */
/*                  IN FORMING THE CHOLESKY DECOMPOSITION OF A+D */
/* INTERNAL VARIABLES */
/* ------------------ */
/* AMINL    SMALLEST ELEMENT ALLOWED ON DIAGONAL OF L */
/* AMNLSQ   =AMINL**2 */
/* OFFMAX   MAXIMUM OFF-DIAGONAL ELEMENT IN COLUMN OF A */


/* DESCRIPTION */
/* ----------- */
/* THE NORMAL CHOLESKY DECOMPOSITION IS PERFORMED.  HOWEVER, IF AT ANY */
/* POINT THE ALGORITHM WOULD ATTEMPT TO SET L(I,I)=SQRT(TEMP) */
/* WITH TEMP < TOL*DIAGMX, THEN L(I,I) IS SET TO SQRT(TOL*DIAGMX) */
/* INSTEAD.  THIS IS EQUIVALENT TO ADDING TOL*DIAGMX-TEMP TO A(I,I) */



    /* Parameter adjustments */
    a_dim1 = *nr;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    *addmax = 0.;
    aminl = sqrt(*diagmx * *tol);
    amnlsq = aminl * aminl;

/* FORM COLUMN J OF L */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* FIND DIAGONAL ELEMENTS OF L */
	sum = 0.;
	if (j == 1) {
	    goto L20;
	}
	jm1 = j - 1;
	i__2 = jm1;
	for (k = 1; k <= i__2; ++k) {
	    sum += a[j + k * a_dim1] * a[j + k * a_dim1];
/* L10: */
	}
L20:
	temp = a[j + j * a_dim1] - sum;
	if (temp < amnlsq) {
	    goto L30;
	}
/*       IF(TEMP.GE.AMINL**2) */
/*       THEN */
	a[j + j * a_dim1] = sqrt(temp);
	goto L40;
/*       ELSE */

/* FIND MAXIMUM OFF-DIAGONAL ELEMENT IN COLUMN */
L30:
	offmax = 0.;
	if (j == *n) {
	    goto L37;
	}
	jp1 = j + 1;
	i__2 = *n;
	for (i = jp1; i <= i__2; ++i) {
	    if ((d__1 = a[i + j * a_dim1], abs(d__1)) > offmax) {
		offmax = (d__2 = a[i + j * a_dim1], abs(d__2));
	    }
/* L35: */
	}
L37:
	if (offmax <= amnlsq) {
	    offmax = amnlsq;
	}

/* ADD TO DIAGONAL ELEMENT  TO ALLOW CHOLESKY DECOMPOSITION TO CONTINU
E */
	a[j + j * a_dim1] = sqrt(offmax);
/* Computing MAX */
	d__1 = *addmax, d__2 = offmax - temp;
	*addmax = max(d__1,d__2);
/*       ENDIF */

/* FIND I,J ELEMENT OF LOWER TRIANGULAR MATRIX */
L40:
	if (j == *n) {
	    goto L100;
	}
	jp1 = j + 1;
	i__2 = *n;
	for (i = jp1; i <= i__2; ++i) {
	    sum = 0.;
	    if (j == 1) {
		goto L60;
	    }
	    jm1 = j - 1;
	    i__3 = jm1;
	    for (k = 1; k <= i__3; ++k) {
		sum += a[i + k * a_dim1] * a[j + k * a_dim1];
/* L50: */
	    }
L60:
	    a[i + j * a_dim1] = (a[i + j * a_dim1] - sum) / a[j + j * a_dim1];
/* L70: */
	}
L100:
	;
    }
    return 0;
} /* choldc_ */

/* Subroutine */ int d1fcn_(integer *n, doublereal *x, doublereal *g)
{

/* PURPOSE */
/* ------- */
/* DUMMY ROUTINE TO PREVENT UNSATISFIED EXTERNAL DIAGNOSTIC */
/* WHEN SPECIFIC ANALYTIC GRADIENT FUNCTION NOT SUPPLIED. */

    /* Parameter adjustments */
    --g;
    --x;

    /* Function Body */
    g[*n] = g[*n];
    x[*n] = x[*n];
    return 0;
} /* d1fcn_ */

/* Subroutine */ int d2fcn_(integer *nr, integer *n, doublereal *x, 
	doublereal *h)
{
    /* System generated locals */
    integer h_dim1, h_offset;


/* PURPOSE */
/* ------- */
/* DUMMY ROUTINE TO PREVENT UNSATISFIED EXTERNAL DIAGNOSTIC */
/* WHEN SPECIFIC ANALYTIC HESSIAN FUNCTION NOT SUPPLIED. */

    /* Parameter adjustments */
    h_dim1 = *nr;
    h_offset = h_dim1 + 1;
    h -= h_offset;
    --x;

    /* Function Body */
    h[*nr + h_dim1] = h[*nr + h_dim1];
    x[*n] = x[*n];
    return 0;
} /* d2fcn_ */

/* Subroutine */ int dfault_(integer *n, doublereal *x, doublereal *typsiz, 
	doublereal *fscale, integer *method, integer *iexp, integer *msg, 
	integer *ndigit, integer *itnlim, integer *iagflg, integer *iahflg, 
	integer *ipr, doublereal *dlt, doublereal *gradtl, doublereal *stepmx,
	 doublereal *steptl)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal);

    /* Local variables */
    static doublereal epsm;
    static integer i;
    extern doublereal d1mach_(integer *);
    extern integer i1mach_(integer *);


/* PURPOSE */
/* ------- */
/* SET DEFAULT VALUES FOR EACH INPUT VARIABLE TO */
/* MINIMIZATION ALGORITHM. */

/* PARAMETERS */
/* ---------- */
/* N            --> DIMENSION OF PROBLEM */
/* X(N)         --> INITIAL GUESS TO SOLUTION (TO COMPUTE MAX STEP SIZE) 
*/
/* TYPSIZ(N)   <--  TYPICAL SIZE FOR EACH COMPONENT OF X */
/* FSCALE      <--  ESTIMATE OF SCALE OF MINIMIZATION FUNCTION */
/* METHOD      <--  ALGORITHM TO USE TO SOLVE MINIMIZATION PROBLEM */
/* IEXP        <--  =0 IF MINIMIZATION FUNCTION NOT EXPENSIVE TO EVALUATE 
*/
/* MSG         <--  MESSAGE TO INHIBIT CERTAIN AUTOMATIC CHECKS + OUTPUT 
*/
/* NDIGIT      <--  NUMBER OF GOOD DIGITS IN MINIMIZATION FUNCTION */
/* ITNLIM      <--  MAXIMUM NUMBER OF ALLOWABLE ITERATIONS */
/* IAGFLG      <--  =0 IF ANALYTIC GRADIENT NOT SUPPLIED */
/* IAHFLG      <--  =0 IF ANALYTIC HESSIAN NOT SUPPLIED */
/* IPR         <--  DEVICE TO WHICH TO SEND OUTPUT */
/* DLT         <--  TRUST REGION RADIUS */
/* GRADTL      <--  TOLERANCE AT WHICH GRADIENT CONSIDERED CLOSE ENOUGH */
/*                  TO ZERO TO TERMINATE ALGORITHM */
/* STEPMX      <--  VALUE OF ZERO TO TRIP DEFAULT MAXIMUM IN OPTCHK */
/* STEPTL      <--  TOLERANCE AT WHICH SUCCESSIVE ITERATES CONSIDERED */
/*                  CLOSE ENOUGH TO TERMINATE ALGORITHM */

    /* Parameter adjustments */
    --typsiz;
    --x;

    /* Function Body */
    x[*n] = x[*n];

/* SET TYPICAL SIZE OF X AND MINIMIZATION FUNCTION */
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	typsiz[i] = 1.;
/* L10: */
    }
    *fscale = 1.;

/* SET TOLERANCES */
    *dlt = -1.;
    epsm = d1mach_(&c__4);
    *gradtl = pow_dd(&epsm, &c_b2);
    *stepmx = 0.;
    *steptl = sqrt(epsm);

/* SET FLAGS */
    *method = 1;
    *iexp = 1;
    *msg = 0;
    *ndigit = -1;
    *itnlim = 150;
    *iagflg = 0;
    *iahflg = 0;
    *ipr = i1mach_(&c__2);

    return 0;
} /* dfault_ */

/* Subroutine */ int dogdrv_(integer *nr, integer *n, doublereal *x, 
	doublereal *f, doublereal *g, doublereal *a, doublereal *p, 
	doublereal *xpls, doublereal *fpls, U_fp fcn, doublereal *sx, 
	doublereal *stepmx, doublereal *steptl, doublereal *dlt, integer *
	iretcd, logical *mxtake, doublereal *sc, doublereal *wrk1, doublereal 
	*wrk2, doublereal *wrk3, integer *ipr)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i;
    static doublereal fplsp;
    static logical fstdog, nwtake;
    extern /* Subroutine */ int dogstp_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, logical *, logical *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *)
	    , tregup_(integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, U_fp, doublereal *, doublereal *, 
	    logical *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, logical *,
	     integer *, integer *, doublereal *);
    static doublereal rnwtln, eta, cln, tmp;


/* PURPOSE */
/* ------- */
/* FIND A NEXT NEWTON ITERATE (XPLS) BY THE DOUBLE DOGLEG METHOD */

/* PARAMETERS */
/* ---------- */
/* NR           --> ROW DIMENSION OF MATRIX */
/* N            --> DIMENSION OF PROBLEM */
/* X(N)         --> OLD ITERATE X[K-1] */
/* F            --> FUNCTION VALUE AT OLD ITERATE, F(X) */
/* G(N)         --> GRADIENT  AT OLD ITERATE, G(X), OR APPROXIMATE */
/* A(N,N)       --> CHOLESKY DECOMPOSITION OF HESSIAN */
/*                  IN LOWER TRIANGULAR PART AND DIAGONAL */
/* P(N)         --> NEWTON STEP */
/* XPLS(N)     <--  NEW ITERATE X[K] */
/* FPLS        <--  FUNCTION VALUE AT NEW ITERATE, F(XPLS) */
/* FCN          --> NAME OF SUBROUTINE TO EVALUATE FUNCTION */
/* SX(N)        --> DIAGONAL SCALING MATRIX FOR X */
/* STEPMX       --> MAXIMUM ALLOWABLE STEP SIZE */
/* STEPTL       --> RELATIVE STEP SIZE AT WHICH SUCCESSIVE ITERATES */
/*                  CONSIDERED CLOSE ENOUGH TO TERMINATE ALGORITHM */
/* DLT         <--> TRUST REGION RADIUS */
/*                  [RETAIN VALUE BETWEEN SUCCESSIVE CALLS] */
/* IRETCD      <--  RETURN CODE */
/*                    =0 SATISFACTORY XPLS FOUND */
/*                    =1 FAILED TO FIND SATISFACTORY XPLS SUFFICIENTLY */
/*                       DISTINCT FROM X */
/* MXTAKE      <--  BOOLEAN FLAG INDICATING STEP OF MAXIMUM LENGTH USED */
/* SC(N)        --> WORKSPACE [CURRENT STEP] */
/* WRK1(N)      --> WORKSPACE (AND PLACE HOLDING ARGUMENT TO TREGUP) */
/* WRK2(N)      --> WORKSPACE */
/* WRK3(N)      --> WORKSPACE */
/* IPR          --> DEVICE TO WHICH TO SEND OUTPUT */


    /* Parameter adjustments */
    --wrk3;
    --wrk2;
    --wrk1;
    --sc;
    --sx;
    --xpls;
    --p;
    a_dim1 = *nr;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --g;
    --x;

    /* Function Body */
    *iretcd = 4;
    fstdog = TRUE_;
    tmp = 0.;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	tmp += sx[i] * sx[i] * p[i] * p[i];
/* L5: */
    }
    rnwtln = sqrt(tmp);
/* $    WRITE(IPR,954) RNWTLN */

L100:

/* FIND NEW STEP BY DOUBLE DOGLEG ALGORITHM */
    dogstp_(nr, n, &g[1], &a[a_offset], &p[1], &sx[1], &rnwtln, dlt, &nwtake, 
	    &fstdog, &wrk1[1], &wrk2[1], &cln, &eta, &sc[1], ipr, stepmx);

/* CHECK NEW POINT AND UPDATE TRUST REGION */
    tregup_(nr, n, &x[1], f, &g[1], &a[a_offset], (U_fp)fcn, &sc[1], &sx[1], &
	    nwtake, stepmx, steptl, dlt, iretcd, &wrk3[1], &fplsp, &xpls[1], 
	    fpls, mxtake, ipr, &c__2, &wrk1[1]);
    if (*iretcd <= 1) {
	return 0;
    }
    goto L100;
/* %950 FORMAT(42H DOGDRV    INITIAL TRUST REGION NOT GIVEN., */
/* %   +       22H  COMPUTE CAUCHY STEP.) */
/* %951 FORMAT(18H DOGDRV    ALPHA =,E20.13/ */
/* %   +       18H DOGDRV    BETA  =,E20.13/ */
/* %   +       18H DOGDRV    DLT   =,E20.13/ */
/* %   +       18H DOGDRV    NWTAKE=,L1    ) */
/* %952 FORMAT(28H DOGDRV    CURRENT STEP (SC)) */
/* %954 FORMAT(18H0DOGDRV    RNWTLN=,E20.13) */
/* %955 FORMAT(14H DOGDRV       ,5(E20.13,3X)) */
} /* dogdrv_ */

/* Subroutine */ int dogstp_(integer *nr, integer *n, doublereal *g, 
	doublereal *a, doublereal *p, doublereal *sx, doublereal *rnwtln, 
	doublereal *dlt, logical *nwtake, logical *fstdog, doublereal *ssd, 
	doublereal *v, doublereal *cln, doublereal *eta, doublereal *sc, 
	integer *ipr, doublereal *stepmx)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal alam, beta;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer i, j;
    static doublereal alpha, tmp, dot1, dot2;


/* PURPOSE */
/* ------- */
/* FIND NEW STEP BY DOUBLE DOGLEG ALGORITHM */

/* PARAMETERS */
/* ---------- */
/* NR           --> ROW DIMENSION OF MATRIX */
/* N            --> DIMENSION OF PROBLEM */
/* G(N)         --> GRADIENT AT CURRENT ITERATE, G(X) */
/* A(N,N)       --> CHOLESKY DECOMPOSITION OF HESSIAN IN */
/*                  LOWER PART AND DIAGONAL */
/* P(N)         --> NEWTON STEP */
/* SX(N)        --> DIAGONAL SCALING MATRIX FOR X */
/* RNWTLN       --> NEWTON STEP LENGTH */
/* DLT         <--> TRUST REGION RADIUS */
/* NWTAKE      <--> BOOLEAN, =.TRUE. IF NEWTON STEP TAKEN */
/* FSTDOG      <--> BOOLEAN, =.TRUE. IF ON FIRST LEG OF DOGLEG */
/* SSD(N)      <--> WORKSPACE [CAUCHY STEP TO THE MINIMUM OF THE */
/*                  QUADRATIC MODEL IN THE SCALED STEEPEST DESCENT */
/*                  DIRECTION] [RETAIN VALUE BETWEEN SUCCESSIVE CALLS] */
/* V(N)        <--> WORKSPACE  [RETAIN VALUE BETWEEN SUCCESSIVE CALLS] */
/* CLN         <--> CAUCHY LENGTH */
/*                  [RETAIN VALUE BETWEEN SUCCESSIVE CALLS] */
/* ETA              [RETAIN VALUE BETWEEN SUCCESSIVE CALLS] */
/* SC(N)       <--  CURRENT STEP */
/* IPR          --> DEVICE TO WHICH TO SEND OUTPUT */
/* STEPMX       --> MAXIMUM ALLOWABLE STEP SIZE */

/* INTERNAL VARIABLES */
/* ------------------ */
/* CLN              LENGTH OF CAUCHY STEP */

    /* Parameter adjustments */
    --sc;
    --v;
    --ssd;
    --sx;
    --p;
    a_dim1 = *nr;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --g;

    /* Function Body */
    *ipr = *ipr;

/* CAN WE TAKE NEWTON STEP */

    if (*rnwtln > *dlt) {
	goto L100;
    }
/*     IF(RNWTLN.LE.DLT) */
/*     THEN */
    *nwtake = TRUE_;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	sc[i] = p[i];
/* L10: */
    }
    *dlt = *rnwtln;
/* $      WRITE(IPR,951) */
    goto L700;
/*     ELSE */

/* NEWTON STEP TOO LONG */
/* CAUCHY STEP IS ON DOUBLE DOGLEG CURVE */

L100:
    *nwtake = FALSE_;
    if (! (*fstdog)) {
	goto L200;
    }
/*       IF(FSTDOG) */
/*       THEN */

/*         CALCULATE DOUBLE DOGLEG CURVE (SSD) */
    *fstdog = FALSE_;
    alpha = 0.;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	alpha += g[i] * g[i] / (sx[i] * sx[i]);
/* L110: */
    }
    beta = 0.;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	tmp = 0.;
	i__2 = *n;
	for (j = i; j <= i__2; ++j) {
	    tmp += a[j + i * a_dim1] * g[j] / (sx[j] * sx[j]);
/* L120: */
	}
	beta += tmp * tmp;
/* L130: */
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	ssd[i] = -(alpha / beta) * g[i] / sx[i];
/* L140: */
    }
    *cln = alpha * sqrt(alpha) / beta;
    *eta = alpha * .8f * alpha / (-beta * ddot_(n, &g[1], &c__1, &p[1], &c__1)
	    ) + .2f;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	v[i] = *eta * sx[i] * p[i] - ssd[i];
/* L150: */
    }
    if (*dlt == -1.) {
	*dlt = min(*cln,*stepmx);
    }
/* $        WRITE(IPR,954) ALPHA,BETA,CLN,ETA */
/* $        WRITE(IPR,955) */
/* $        WRITE(IPR,960) (SSD(I),I=1,N) */
/* $        WRITE(IPR,956) */
/* $        WRITE(IPR,960) (V(I),I=1,N) */
/*       ENDIF */
L200:
    if (*eta * *rnwtln > *dlt) {
	goto L220;
    }
/*       IF(ETA*RNWTLN .LE. DLT) */
/*       THEN */

/*         TAKE PARTIAL STEP IN NEWTON DIRECTION */

    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	sc[i] = *dlt / *rnwtln * p[i];
/* L210: */
    }
/* $        WRITE(IPR,957) */
    goto L700;
/*       ELSE */
L220:
    if (*cln < *dlt) {
	goto L240;
    }
/*         IF(CLN.GE.DLT) */
/*         THEN */
/*           TAKE STEP IN STEEPEST DESCENT DIRECTION */

    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	sc[i] = *dlt / *cln * ssd[i] / sx[i];
/* L230: */
    }
/* $          WRITE(IPR,958) */
    goto L700;
/*         ELSE */

/*           CALCULATE CONVEX COMBINATION OF SSD AND ETA*P */
/*           WHICH HAS SCALED LENGTH DLT */

L240:
    dot1 = ddot_(n, &v[1], &c__1, &ssd[1], &c__1);
    dot2 = ddot_(n, &v[1], &c__1, &v[1], &c__1);
    alam = (-dot1 + sqrt(dot1 * dot1 - dot2 * (*cln * *cln - *dlt * *dlt))) / 
	    dot2;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	sc[i] = (ssd[i] + alam * v[i]) / sx[i];
/* L250: */
    }
/* $          WRITE(IPR,959) */
/*         ENDIF */
/*       ENDIF */
/*     ENDIF */
L700:
/* $    WRITE(IPR,952) FSTDOG,NWTAKE,RNWTLN,DLT */
/* $    WRITE(IPR,953) */
/* $    WRITE(IPR,960) (SC(I),I=1,N) */
    return 0;

/* %951 FORMAT(27H0DOGSTP    TAKE NEWTON STEP) */
/* %952 FORMAT(18H DOGSTP    FSTDOG=,L1/ */
/* %   +       18H DOGSTP    NWTAKE=,L1/ */
/* %   +       18H DOGSTP    RNWTLN=,E20.13/ */
/* %   +       18H DOGSTP    DLT   =,E20.13) */
/* %953 FORMAT(28H DOGSTP    CURRENT STEP (SC)) */
/* %954 FORMAT(18H DOGSTP    ALPHA =,E20.13/ */
/* %   +       18H DOGSTP    BETA  =,E20.13/ */
/* %   +       18H DOGSTP    CLN   =,E20.13/ */
/* %   +       18H DOGSTP    ETA   =,E20.13) */
/* %955 FORMAT(28H DOGSTP    CAUCHY STEP (SSD)) */
/* %956 FORMAT(12H DOGSTP    V) */
/* %957 FORMAT(48H0DOGSTP    TAKE PARTIAL STEP IN NEWTON DIRECTION) */
/* %958 FORMAT(50H0DOGSTP    TAKE STEP IN STEEPEST DESCENT DIRECTION) */
/* %959 FORMAT(39H0DOGSTP    TAKE CONVEX COMBINATION STEP) */
/* %960 FORMAT(14H DOGSTP       ,5(E20.13,3X)) */
} /* dogstp_ */

/* Subroutine */ int forslv_(integer *nr, integer *n, doublereal *a, 
	doublereal *x, doublereal *b)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i, j, im1;
    static doublereal sum;


/* PURPOSE */
/* ------- */
/* SOLVE  AX=B  WHERE A IS LOWER TRIANGULAR MATRIX */

/* PARAMETERS */
/* ---------- */
/* NR           --> ROW DIMENSION OF MATRIX */
/* N            --> DIMENSION OF PROBLEM */
/* A(N,N)       --> LOWER TRIANGULAR MATRIX (PRESERVED) */
/* X(N)        <--  SOLUTION VECTOR */
/* B(N)         --> RIGHT-HAND SIDE VECTOR */

/* NOTE */
/* ---- */
/* IF B IS NO LONGER REQUIRED BY CALLING ROUTINE, */
/* THEN VECTORS B AND X MAY SHARE THE SAME STORAGE. */


/* SOLVE LX=B. (FOREWARD SOLVE) */

    /* Parameter adjustments */
    --b;
    --x;
    a_dim1 = *nr;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    x[1] = b[1] / a[a_dim1 + 1];
    if (*n == 1) {
	return 0;
    }
    i__1 = *n;
    for (i = 2; i <= i__1; ++i) {
	sum = 0.;
	im1 = i - 1;
	i__2 = im1;
	for (j = 1; j <= i__2; ++j) {
	    sum += a[i + j * a_dim1] * x[j];
/* L10: */
	}
	x[i] = (b[i] - sum) / a[i + i * a_dim1];
/* L20: */
    }
    return 0;
} /* forslv_ */

/* Subroutine */ int fstocd_(integer *n, doublereal *x, S_fp fcn, doublereal *
	sx, doublereal *rnoise, doublereal *g)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static integer i;
    static doublereal third, stepi, fplus, fminus, xtempi;

/* PURPOSE */
/* ------- */
/* FIND CENTRAL DIFFERENCE APPROXIMATION G TO THE FIRST DERIVATIVE */
/* (GRADIENT) OF THE FUNCTION DEFINED BY FCN AT THE POINT X. */

/* PARAMETERS */
/* ---------- */
/* N            --> DIMENSION OF PROBLEM */
/* X            --> POINT AT WHICH GRADIENT IS TO BE APPROXIMATED. */
/* FCN          --> NAME OF SUBROUTINE TO EVALUATE FUNCTION. */
/* SX           --> DIAGONAL SCALING MATRIX FOR X. */
/* RNOISE       --> RELATIVE NOISE IN FCN [F(X)]. */
/* G           <--  CENTRAL DIFFERENCE APPROXIMATION TO GRADIENT. */



/* FIND I TH  STEPSIZE, EVALUATE TWO NEIGHBORS IN DIRECTION OF I TH */
/* UNIT VECTOR, AND EVALUATE I TH  COMPONENT OF GRADIENT. */

    /* Parameter adjustments */
    --g;
    --sx;
    --x;

    /* Function Body */
    third = .33333333333333331;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
/* Computing MAX */
	d__2 = (d__1 = x[i], abs(d__1)), d__3 = 1. / sx[i];
	stepi = pow_dd(rnoise, &third) * max(d__2,d__3);
	xtempi = x[i];
	x[i] = xtempi + stepi;
	(*fcn)(n, &x[1], &fplus);
	x[i] = xtempi - stepi;
	(*fcn)(n, &x[1], &fminus);
	x[i] = xtempi;
	g[i] = (fplus - fminus) / (stepi * 2.);
/* L10: */
    }
    return 0;
} /* fstocd_ */

/* Subroutine */ int fstofd_(integer *nr, integer *m, integer *n, doublereal *
	xpls, S_fp fcn, doublereal *fpls, doublereal *a, doublereal *sx, 
	doublereal *rnoise, doublereal *fhat, integer *icase)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i, j;
    static doublereal xtmpj;
    static integer jp1, nm1;
    static doublereal stepsz;

/* PURPOSE */
/* ------- */
/* FIND FIRST ORDER FORWARD FINITE DIFFERENCE APPROXIMATION "A" TO THE */
/* FIRST DERIVATIVE OF THE FUNCTION DEFINED BY THE SUBPROGRAM "FNAME" */
/* EVALUATED AT THE NEW ITERATE "XPLS". */


/* FOR OPTIMIZATION USE THIS ROUTINE TO ESTIMATE: */
/* 1) THE FIRST DERIVATIVE (GRADIENT) OF THE OPTIMIZATION FUNCTION "FCN */
/*    ANALYTIC USER ROUTINE HAS BEEN SUPPLIED; */
/* 2) THE SECOND DERIVATIVE (HESSIAN) OF THE OPTIMIZATION FUNCTION */
/*    IF NO ANALYTIC USER ROUTINE HAS BEEN SUPPLIED FOR THE HESSIAN BUT */
/*    ONE HAS BEEN SUPPLIED FOR THE GRADIENT ("FCN") AND IF THE */
/*    OPTIMIZATION FUNCTION IS INEXPENSIVE TO EVALUATE */

/* NOTE */
/* ---- */
/* _M=1 (OPTIMIZATION) ALGORITHM ESTIMATES THE GRADIENT OF THE FUNCTION */
/*      (FCN).   FCN(X) # F: R(N)-->R(1) */
/* _M=N (SYSTEMS) ALGORITHM ESTIMATES THE JACOBIAN OF THE FUNCTION */
/*      FCN(X) # F: R(N)-->R(N). */
/* _M=N (OPTIMIZATION) ALGORITHM ESTIMATES THE HESSIAN OF THE OPTIMIZATIO 
*/
/*      FUNCTION, WHERE THE HESSIAN IS THE FIRST DERIVATIVE OF "FCN" */

/* PARAMETERS */
/* ---------- */
/* NR           --> ROW DIMENSION OF MATRIX */
/* M            --> NUMBER OF ROWS IN A */
/* N            --> NUMBER OF COLUMNS IN A; DIMENSION OF PROBLEM */
/* XPLS(N)      --> NEW ITERATE:  X[K] */
/* FCN          --> NAME OF SUBROUTINE TO EVALUATE FUNCTION */
/* FPLS(M)      --> _M=1 (OPTIMIZATION) FUNCTION VALUE AT NEW ITERATE: */
/*                       FCN(XPLS) */
/*                  _M=N (OPTIMIZATION) VALUE OF FIRST DERIVATIVE */
/*                       (GRADIENT) GIVEN BY USER FUNCTION FCN */
/*                  _M=N (SYSTEMS)  FUNCTION VALUE OF ASSOCIATED */
/*                       MINIMIZATION FUNCTION */
/* A(NR,N)     <--  FINITE DIFFERENCE APPROXIMATION (SEE NOTE).  ONLY */
/*                  LOWER TRIANGULAR MATRIX AND DIAGONAL ARE RETURNED */
/* SX(N)        --> DIAGONAL SCALING MATRIX FOR X */
/* RNOISE       --> RELATIVE NOISE IN FCN [F(X)] */
/* FHAT(M)      --> WORKSPACE */
/* ICASE        --> =1 OPTIMIZATION (GRADIENT) */
/*                  =2 SYSTEMS */
/*                  =3 OPTIMIZATION (HESSIAN) */

/* INTERNAL VARIABLES */
/* ------------------ */
/* STEPSZ - STEPSIZE IN THE J-TH VARIABLE DIRECTION */


/* FIND J-TH COLUMN OF A */
/* EACH COLUMN IS DERIVATIVE OF F(FCN) WITH RESPECT TO XPLS(J) */

    /* Parameter adjustments */
    --fhat;
    --sx;
    a_dim1 = *nr;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --fpls;
    --xpls;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	d__2 = (d__1 = xpls[j], abs(d__1)), d__3 = 1. / sx[j];
	stepsz = sqrt(*rnoise) * max(d__2,d__3);
	xtmpj = xpls[j];
	xpls[j] = xtmpj + stepsz;
	(*fcn)(n, &xpls[1], &fhat[1]);
	xpls[j] = xtmpj;
	i__2 = *m;
	for (i = 1; i <= i__2; ++i) {
	    a[i + j * a_dim1] = (fhat[i] - fpls[i]) / stepsz;
/* L20: */
	}
/* L30: */
    }
    if (*icase != 3) {
	return 0;
    }

/* IF COMPUTING HESSIAN, A MUST BE SYMMETRIC */

    if (*n == 1) {
	return 0;
    }
    nm1 = *n - 1;
    i__1 = nm1;
    for (j = 1; j <= i__1; ++j) {
	jp1 = j + 1;
	i__2 = *m;
	for (i = jp1; i <= i__2; ++i) {
	    a[i + j * a_dim1] = (a[i + j * a_dim1] + a[j + i * a_dim1]) / 2.;
/* L40: */
	}
/* L50: */
    }
    return 0;
} /* fstofd_ */

/* Subroutine */ int grdchk_(integer *n, doublereal *x, S_fp fcn, doublereal *
	f, doublereal *g, doublereal *typsiz, doublereal *sx, doublereal *
	fscale, doublereal *rnf, doublereal *analtl, doublereal *wrk1, 
	integer *msg, integer *ipr)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3, d__4;

    /* Local variables */
    static integer i;
    static doublereal gs;
    extern /* Subroutine */ int fstofd_(integer *, integer *, integer *, 
	    doublereal *, S_fp, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *);
    static integer ker;
    static doublereal wrk;


/* PURPOSE */
/* ------- */
/* CHECK ANALYTIC GRADIENT AGAINST ESTIMATED GRADIENT */

/* PARAMETERS */
/* ---------- */
/* N            --> DIMENSION OF PROBLEM */
/* X(N)         --> ESTIMATE TO A ROOT OF FCN */
/* FCN          --> NAME OF SUBROUTINE TO EVALUATE OPTIMIZATION FUNCTION 
*/
/*                  MUST BE DECLARED EXTERNAL IN CALLING ROUTINE */
/*                       FCN:  R(N) --> R(1) */
/* F            --> FUNCTION VALUE:  FCN(X) */
/* G(N)         --> GRADIENT:  G(X) */
/* TYPSIZ(N)    --> TYPICAL SIZE FOR EACH COMPONENT OF X */
/* SX(N)        --> DIAGONAL SCALING MATRIX:  SX(I)=1./TYPSIZ(I) */
/* FSCALE       --> ESTIMATE OF SCALE OF OBJECTIVE FUNCTION FCN */
/* RNF          --> RELATIVE NOISE IN OPTIMIZATION FUNCTION FCN */
/* ANALTL       --> TOLERANCE FOR COMPARISON OF ESTIMATED AND */
/*                  ANALYTICAL GRADIENTS */
/* WRK1(N)      --> WORKSPACE */
/* MSG         <--  MESSAGE OR ERROR CODE */
/*                    ON OUTPUT: =-21, PROBABLE CODING ERROR OF GRADIENT 
*/
/* IPR          --> DEVICE TO WHICH TO SEND OUTPUT */


/* COMPUTE FIRST ORDER FINITE DIFFERENCE GRADIENT AND COMPARE TO */
/* ANALYTIC GRADIENT. */

    /* Parameter adjustments */
    --wrk1;
    --sx;
    --typsiz;
    --g;
    --x;

    /* Function Body */
    fstofd_(&c__1, &c__1, n, &x[1], (S_fp)fcn, f, &wrk1[1], &sx[1], rnf, &wrk,
	     &c__1);
    ker = 0;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
/* Computing MAX */
	d__2 = abs(*f);
/* Computing MAX */
	d__3 = (d__1 = x[i], abs(d__1)), d__4 = typsiz[i];
	gs = max(d__2,*fscale) / max(d__3,d__4);
/* Computing MAX */
	d__3 = (d__2 = g[i], abs(d__2));
	if ((d__1 = g[i] - wrk1[i], abs(d__1)) > max(d__3,gs) * *analtl) {
	    ker = 1;
	}
/* L5: */
    }
    if (ker == 0) {
	goto L20;
    }
/* %      WRITE(IPR,901) */
/* %      WRITE(IPR,902) (I,G(I),WRK1(I),I=1,N) */
    *msg = -21;
L20:
    return 0;
/* %901 FORMAT(47H0GRDCHK    PROBABLE ERROR IN CODING OF ANALYTIC, */
/* %   +       19H GRADIENT FUNCTION./ */
/* %   +       16H GRDCHK     COMP,12X,8HANALYTIC,12X,8HESTIMATE) */
/* %902 FORMAT(11H GRDCHK    ,I5,3X,E20.13,3X,E20.13) */
} /* grdchk_ */

/* Subroutine */ int heschk_(integer *nr, integer *n, doublereal *x, U_fp fcn,
	 S_fp d1fcn, S_fp d2fcn, doublereal *f, doublereal *g, doublereal *a, 
	doublereal *typsiz, doublereal *sx, doublereal *rnf, doublereal *
	analtl, integer *iagflg, doublereal *udiag, doublereal *wrk1, 
	doublereal *wrk2, integer *msg, integer *ipr)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2, d__3, d__4, d__5;

    /* Local variables */
    static integer i, j;
    static doublereal hs;
    extern /* Subroutine */ int sndofd_(integer *, integer *, doublereal *, 
	    U_fp, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), fstofd_(integer *, integer *, 
	    integer *, doublereal *, S_fp, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *);
    static integer jp1, ker;


/* PURPOSE */
/* ------- */
/* CHECK ANALYTIC HESSIAN AGAINST ESTIMATED HESSIAN */
/*  (THIS MAY BE DONE ONLY IF THE USER SUPPLIED ANALYTIC HESSIAN */
/*   D2FCN FILLS ONLY THE LOWER TRIANGULAR PART AND DIAGONAL OF A) */

/* PARAMETERS */
/* ---------- */
/* NR           --> ROW DIMENSION OF MATRIX */
/* N            --> DIMENSION OF PROBLEM */
/* X(N)         --> ESTIMATE TO A ROOT OF FCN */
/* FCN          --> NAME OF SUBROUTINE TO EVALUATE OPTIMIZATION FUNCTION 
*/
/*                  MUST BE DECLARED EXTERNAL IN CALLING ROUTINE */
/*                       FCN:  R(N) --> R(1) */
/* D1FCN        --> NAME OF SUBROUTINE TO EVALUATE GRADIENT OF FCN. */
/*                  MUST BE DECLARED EXTERNAL IN CALLING ROUTINE */
/* D2FCN        --> NAME OF SUBROUTINE TO EVALUATE HESSIAN OF FCN. */
/*                  MUST BE DECLARED EXTERNAL IN CALLING ROUTINE */
/* F            --> FUNCTION VALUE:  FCN(X) */
/* G(N)        <--  GRADIENT:  G(X) */
/* A(N,N)      <--  ON EXIT:  HESSIAN IN LOWER TRIANGULAR PART AND DIAG */
/* TYPSIZ(N)    --> TYPICAL SIZE FOR EACH COMPONENT OF X */
/* SX(N)        --> DIAGONAL SCALING MATRIX:  SX(I)=1./TYPSIZ(I) */
/* RNF          --> RELATIVE NOISE IN OPTIMIZATION FUNCTION FCN */
/* ANALTL       --> TOLERANCE FOR COMPARISON OF ESTIMATED AND */
/*                  ANALYTICAL GRADIENTS */
/* IAGFLG       --> =1 IF ANALYTIC GRADIENT SUPPLIED */
/* UDIAG(N)     --> WORKSPACE */
/* WRK1(N)      --> WORKSPACE */
/* WRK2(N)      --> WORKSPACE */
/* MSG         <--> MESSAGE OR ERROR CODE */
/*                    ON INPUT : IF =1XX DO NOT COMPARE ANAL + EST HESS */
/*                    ON OUTPUT: =-22, PROBABLE CODING ERROR OF HESSIAN */
/* IPR          --> DEVICE TO WHICH TO SEND OUTPUT */


/* COMPUTE FINITE DIFFERENCE APPROXIMATION A TO THE HESSIAN. */

    /* Parameter adjustments */
    --wrk2;
    --wrk1;
    --udiag;
    --sx;
    --typsiz;
    a_dim1 = *nr;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --g;
    --x;

    /* Function Body */
    if (*iagflg == 1) {
	fstofd_(nr, n, n, &x[1], (S_fp)d1fcn, &g[1], &a[a_offset], &sx[1], 
		rnf, &wrk1[1], &c__3);
    }
    if (*iagflg != 1) {
	sndofd_(nr, n, &x[1], (U_fp)fcn, f, &a[a_offset], &sx[1], rnf, &wrk1[
		1], &wrk2[1]);
    }

    ker = 0;

/* COPY LOWER TRIANGULAR PART OF "A" TO UPPER TRIANGULAR PART */
/* AND DIAGONAL OF "A" TO UDIAG */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	udiag[j] = a[j + j * a_dim1];
	if (j == *n) {
	    goto L30;
	}
	jp1 = j + 1;
	i__2 = *n;
	for (i = jp1; i <= i__2; ++i) {
	    a[j + i * a_dim1] = a[i + j * a_dim1];
/* L25: */
	}
L30:
	;
    }

/* COMPUTE ANALYTIC HESSIAN AND COMPARE TO FINITE DIFFERENCE */
/* APPROXIMATION. */

    (*d2fcn)(nr, n, &x[1], &a[a_offset]);
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	d__3 = (d__1 = g[j], abs(d__1));
/* Computing MAX */
	d__4 = (d__2 = x[j], abs(d__2)), d__5 = typsiz[j];
	hs = max(d__3,1.) / max(d__4,d__5);
/* Computing MAX */
	d__3 = (d__2 = udiag[j], abs(d__2));
	if ((d__1 = a[j + j * a_dim1] - udiag[j], abs(d__1)) > max(d__3,hs) * 
		*analtl) {
	    ker = 1;
	}
	if (j == *n) {
	    goto L40;
	}
	jp1 = j + 1;
	i__2 = *n;
	for (i = jp1; i <= i__2; ++i) {
/* Computing MAX */
	    d__3 = (d__2 = a[i + j * a_dim1], abs(d__2));
	    if ((d__1 = a[i + j * a_dim1] - a[j + i * a_dim1], abs(d__1)) > 
		    max(d__3,hs) * *analtl) {
		ker = 1;
	    }
/* L35: */
	}
L40:
	;
    }

    if (ker == 0) {
	goto L90;
    }
/* %      WRITE(IPR,901) */
/* %      DO 50 I=1,N */
/* %        IF(I.EQ.1) GO TO 45 */
/* %        IM1=I-1 */
/* %        DO 43 J=1,IM1 */
/* %          WRITE(IPR,902) I,J,A(I,J),A(J,I) */
/* % 43     CONTINUE */
/* % 45     WRITE(IPR,902) I,I,A(I,I),UDIAG(I) */
/* % 50   CONTINUE */
    *msg = -22;
/*     ENDIF */
L90:
    return 0;
/* %901 FORMAT(47H HESCHK    PROBABLE ERROR IN CODING OF ANALYTIC, */
/* %   +       18H HESSIAN FUNCTION./ */
/* %   +       21H HESCHK      ROW  COL,14X,8HANALYTIC,14X,10H(ESTIMATE)) 
*/
/* %902 FORMAT(11H HESCHK    ,2I5,2X,E20.13,2X,1H(,E20.13,1H)) */
} /* heschk_ */

/* Subroutine */ int hookdr_(integer *nr, integer *n, doublereal *x, 
	doublereal *f, doublereal *g, doublereal *a, doublereal *udiag, 
	doublereal *p, doublereal *xpls, doublereal *fpls, U_fp fcn, 
	doublereal *sx, doublereal *stepmx, doublereal *steptl, doublereal *
	dlt, integer *iretcd, logical *mxtake, doublereal *amu, doublereal *
	dltp, doublereal *phi, doublereal *phip0, doublereal *sc, doublereal *
	xplsp, doublereal *wrk0, doublereal *epsm, integer *itncnt, integer *
	ipr)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal beta;
    static integer i, j;
    static doublereal alpha, fplsp;
    static logical fstime, nwtake;
    extern /* Subroutine */ int tregup_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, U_fp, doublereal *, 
	    doublereal *, logical *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, doublereal *, doublereal *, doublereal *
	    , logical *, integer *, integer *, doublereal *), hookst_(integer 
	    *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, logical *,
	     doublereal *, logical *, doublereal *, doublereal *, integer *);
    static doublereal rnwtln, tmp;


/* PURPOSE */
/* ------- */
/* FIND A NEXT NEWTON ITERATE (XPLS) BY THE MORE-HEBDON METHOD */

/* PARAMETERS */
/* ---------- */
/* NR           --> ROW DIMENSION OF MATRIX */
/* N            --> DIMENSION OF PROBLEM */
/* X(N)         --> OLD ITERATE X[K-1] */
/* F            --> FUNCTION VALUE AT OLD ITERATE, F(X) */
/* G(N)         --> GRADIENT AT OLD ITERATE, G(X), OR APPROXIMATE */
/* A(N,N)       --> CHOLESKY DECOMPOSITION OF HESSIAN IN LOWER */
/*                  TRIANGULAR PART AND DIAGONAL. */
/*                  HESSIAN IN UPPER TRIANGULAR PART AND UDIAG. */
/* UDIAG(N)     --> DIAGONAL OF HESSIAN IN A(.,.) */
/* P(N)         --> NEWTON STEP */
/* XPLS(N)     <--  NEW ITERATE X[K] */
/* FPLS        <--  FUNCTION VALUE AT NEW ITERATE, F(XPLS) */
/* FCN          --> NAME OF SUBROUTINE TO EVALUATE FUNCTION */
/* SX(N)        --> DIAGONAL SCALING MATRIX FOR X */
/* STEPMX       --> MAXIMUM ALLOWABLE STEP SIZE */
/* STEPTL       --> RELATIVE STEP SIZE AT WHICH SUCCESSIVE ITERATES */
/*                  CONSIDERED CLOSE ENOUGH TO TERMINATE ALGORITHM */
/* DLT         <--> TRUST REGION RADIUS */
/* IRETCD      <--  RETURN CODE */
/*                    =0 SATISFACTORY XPLS FOUND */
/*                    =1 FAILED TO FIND SATISFACTORY XPLS SUFFICIENTLY */
/*                       DISTINCT FROM X */
/* MXTAKE      <--  BOOLEAN FLAG INDICATING STEP OF MAXIMUM LENGTH USED */
/* AMU         <--> [RETAIN VALUE BETWEEN SUCCESSIVE CALLS] */
/* DLTP        <--> [RETAIN VALUE BETWEEN SUCCESSIVE CALLS] */
/* PHI         <--> [RETAIN VALUE BETWEEN SUCCESSIVE CALLS] */
/* PHIP0       <--> [RETAIN VALUE BETWEEN SUCCESSIVE CALLS] */
/* SC(N)        --> WORKSPACE */
/* XPLSP(N)     --> WORKSPACE */
/* WRK0(N)      --> WORKSPACE */
/* EPSM         --> MACHINE EPSILON */
/* ITNCNT       --> ITERATION COUNT */
/* IPR          --> DEVICE TO WHICH TO SEND OUTPUT */


    /* Parameter adjustments */
    --wrk0;
    --xplsp;
    --sc;
    --sx;
    --xpls;
    --p;
    --udiag;
    a_dim1 = *nr;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --g;
    --x;

    /* Function Body */
    *iretcd = 4;
    fstime = TRUE_;
    tmp = 0.;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	tmp += sx[i] * sx[i] * p[i] * p[i];
/* L5: */
    }
    rnwtln = sqrt(tmp);
/* $    WRITE(IPR,954) RNWTLN */

    if (*itncnt > 1) {
	goto L100;
    }
/*     IF(ITNCNT.EQ.1) */
/*     THEN */
    *amu = 0.;

/*       IF FIRST ITERATION AND TRUST REGION NOT PROVIDED BY USER, */
/*       COMPUTE INITIAL TRUST REGION. */

    if (*dlt != -1.) {
	goto L100;
    }
/*       IF(DLT.EQ. (-1.0D0)) */
/*       THEN */
    alpha = 0.;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	alpha += g[i] * g[i] / (sx[i] * sx[i]);
/* L10: */
    }
    beta = 0.;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	tmp = 0.;
	i__2 = *n;
	for (j = i; j <= i__2; ++j) {
	    tmp += a[j + i * a_dim1] * g[j] / (sx[j] * sx[j]);
/* L20: */
	}
	beta += tmp * tmp;
/* L30: */
    }
    *dlt = alpha * sqrt(alpha) / beta;
    *dlt = min(*dlt,*stepmx);
/* $        WRITE(IPR,950) */
/* $        WRITE(IPR,951) ALPHA,BETA,DLT */
/*       ENDIF */
/*     ENDIF */

L100:

/* FIND NEW STEP BY MORE-HEBDON ALGORITHM */
    hookst_(nr, n, &g[1], &a[a_offset], &udiag[1], &p[1], &sx[1], &rnwtln, 
	    dlt, amu, dltp, phi, phip0, &fstime, &sc[1], &nwtake, &wrk0[1], 
	    epsm, ipr);
    *dltp = *dlt;

/* CHECK NEW POINT AND UPDATE TRUST REGION */
    tregup_(nr, n, &x[1], f, &g[1], &a[a_offset], (U_fp)fcn, &sc[1], &sx[1], &
	    nwtake, stepmx, steptl, dlt, iretcd, &xplsp[1], &fplsp, &xpls[1], 
	    fpls, mxtake, ipr, &c__3, &udiag[1]);
    if (*iretcd <= 1) {
	return 0;
    }
    goto L100;

/* %950 FORMAT(43H HOOKDR    INITIAL TRUST REGION NOT GIVEN. , */
/* %   +       21H COMPUTE CAUCHY STEP.) */
/* %951 FORMAT(18H HOOKDR    ALPHA =,E20.13/ */
/* %   +       18H HOOKDR    BETA  =,E20.13/ */
/* %   +       18H HOOKDR    DLT   =,E20.13) */
/* %952 FORMAT(28H HOOKDR    CURRENT STEP (SC)) */
/* %954 FORMAT(18H0HOOKDR    RNWTLN=,E20.13) */
/* %955 FORMAT(14H HOOKDR       ,5(E20.13,3X)) */
} /* hookdr_ */

/* Subroutine */ int hookst_(integer *nr, integer *n, doublereal *g, 
	doublereal *a, doublereal *udiag, doublereal *p, doublereal *sx, 
	doublereal *rnwtln, doublereal *dlt, doublereal *amu, doublereal *
	dltp, doublereal *phi, doublereal *phip0, logical *fstime, doublereal 
	*sc, logical *nwtake, doublereal *wrk0, doublereal *epsm, integer *
	ipr)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static logical done;
    static doublereal phip;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static integer i, j;
    static doublereal amulo, amuup, hi;
    extern /* Subroutine */ int choldc_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static doublereal addmax, stepln;
    extern /* Subroutine */ int forslv_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *);
    static integer jp1;
    extern /* Subroutine */ int lltslv_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal alo;


/* PURPOSE */
/* ------- */
/* FIND NEW STEP BY MORE-HEBDON ALGORITHM */

/* PARAMETERS */
/* ---------- */
/* NR           --> ROW DIMENSION OF MATRIX */
/* N            --> DIMENSION OF PROBLEM */
/* G(N)         --> GRADIENT AT CURRENT ITERATE, G(X) */
/* A(N,N)       --> CHOLESKY DECOMPOSITION OF HESSIAN IN */
/*                  LOWER TRIANGULAR PART AND DIAGONAL. */
/*                  HESSIAN OR APPROX IN UPPER TRIANGULAR PART */
/* UDIAG(N)     --> DIAGONAL OF HESSIAN IN A(.,.) */
/* P(N)         --> NEWTON STEP */
/* SX(N)        --> DIAGONAL SCALING MATRIX FOR N */
/* RNWTLN       --> NEWTON STEP LENGTH */
/* DLT         <--> TRUST REGION RADIUS */
/* AMU         <--> [RETAIN VALUE BETWEEN SUCCESSIVE CALLS] */
/* DLTP         --> TRUST REGION RADIUS AT LAST EXIT FROM THIS ROUTINE */
/* PHI         <--> [RETAIN VALUE BETWEEN SUCCESSIVE CALLS] */
/* PHIP0       <--> [RETAIN VALUE BETWEEN SUCCESSIVE CALLS] */
/* FSTIME      <--> BOOLEAN. =.TRUE. IF FIRST ENTRY TO THIS ROUTINE */
/*                  DURING K-TH ITERATION */
/* SC(N)       <--  CURRENT STEP */
/* NWTAKE      <--  BOOLEAN, =.TRUE. IF NEWTON STEP TAKEN */
/* WRK0         --> WORKSPACE */
/* EPSM         --> MACHINE EPSILON */
/* IPR          --> DEVICE TO WHICH TO SEND OUTPUT */


/* HI AND ALO ARE CONSTANTS USED IN THIS ROUTINE. */
/* CHANGE HERE IF OTHER VALUES ARE TO BE SUBSTITUTED. */
    /* Parameter adjustments */
    --wrk0;
    --sc;
    --sx;
    --p;
    --udiag;
    a_dim1 = *nr;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --g;

    /* Function Body */
    *ipr = *ipr;
    hi = 1.5;
    alo = .75;
/* ----- */
    if (*rnwtln > hi * *dlt) {
	goto L15;
    }
/*     IF(RNWTLN.LE.HI*DLT) */
/*     THEN */

/*       TAKE NEWTON STEP */

    *nwtake = TRUE_;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	sc[i] = p[i];
/* L10: */
    }
    *dlt = min(*dlt,*rnwtln);
    *amu = 0.;
/* $      WRITE(IPR,951) */
    return 0;
/*     ELSE */

/*       NEWTON STEP NOT TAKEN */

L15:
/* $      WRITE(IPR,952) */
    *nwtake = FALSE_;
    if (*amu <= 0.) {
	goto L20;
    }
/*       IF(AMU.GT.0.0D0) */
/*       THEN */
    *amu -= (*phi + *dltp) * (*dltp - *dlt + *phi) / (*dlt * phip);
/* $        WRITE(IPR,956) AMU */
/*       ENDIF */
L20:
    *phi = *rnwtln - *dlt;
    if (! (*fstime)) {
	goto L28;
    }
/*       IF(FSTIME) */
/*       THEN */
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	wrk0[i] = sx[i] * sx[i] * p[i];
/* L25: */
    }

/*         SOLVE L*Y = (SX**2)*P */

    forslv_(nr, n, &a[a_offset], &wrk0[1], &wrk0[1]);
/* Computing 2nd power */
    d__1 = dnrm2_(n, &wrk0[1], &c__1);
    *phip0 = -(d__1 * d__1) / *rnwtln;
    *fstime = FALSE_;
/*       ENDIF */
L28:
    phip = *phip0;
    amulo = -(*phi) / phip;
    amuup = 0.;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	amuup += g[i] * g[i] / (sx[i] * sx[i]);
/* L30: */
    }
    amuup = sqrt(amuup) / *dlt;
    done = FALSE_;
/* $      WRITE(IPR,956) AMU */
/* $      WRITE(IPR,959) PHI */
/* $      WRITE(IPR,960) PHIP */
/* $      WRITE(IPR,957) AMULO */
/* $      WRITE(IPR,958) AMUUP */

/*       TEST VALUE OF AMU; GENERATE NEXT AMU IF NECESSARY */

L100:
    if (done) {
	return 0;
    }
/* $      WRITE(IPR,962) */
    if (*amu >= amulo && *amu <= amuup) {
	goto L110;
    }
/*       IF(AMU.LT.AMULO .OR.  AMU.GT.AMUUP) */
/*       THEN */
/* Computing MAX */
    d__1 = sqrt(amulo * amuup), d__2 = amuup * .001;
    *amu = max(d__1,d__2);
/* $        WRITE(IPR,956) AMU */
/*       ENDIF */
L110:

/*       COPY (H,UDIAG) TO L */
/*       WHERE H <-- H+AMU*(SX**2) [DO NOT ACTUALLY CHANGE (H,UDIAG)] */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	a[j + j * a_dim1] = udiag[j] + *amu * sx[j] * sx[j];
	if (j == *n) {
	    goto L130;
	}
	jp1 = j + 1;
	i__2 = *n;
	for (i = jp1; i <= i__2; ++i) {
	    a[i + j * a_dim1] = a[j + i * a_dim1];
/* L120: */
	}
L130:
	;
    }

/*       FACTOR H=L(L+) */

    d__1 = sqrt(*epsm);
    choldc_(nr, n, &a[a_offset], &c_b33, &d__1, &addmax);

/*       SOLVE H*P = L(L+)*SC = -G */

    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	wrk0[i] = -g[i];
/* L140: */
    }
    lltslv_(nr, n, &a[a_offset], &sc[1], &wrk0[1]);
/* $      WRITE(IPR,955) */
/* $      WRITE(IPR,963) (SC(I),I=1,N) */

/*       RESET H.  NOTE SINCE UDIAG HAS NOT BEEN DESTROYED WE NEED DO */
/*       NOTHING HERE.  H IS IN THE UPPER PART AND IN UDIAG, STILL INTACT 
*/

    stepln = 0.;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	stepln += sx[i] * sx[i] * sc[i] * sc[i];
/* L150: */
    }
    stepln = sqrt(stepln);
    *phi = stepln - *dlt;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	wrk0[i] = sx[i] * sx[i] * sc[i];
/* L160: */
    }
    forslv_(nr, n, &a[a_offset], &wrk0[1], &wrk0[1]);
/* Computing 2nd power */
    d__1 = dnrm2_(n, &wrk0[1], &c__1);
    phip = -(d__1 * d__1) / stepln;
/* $      WRITE(IPR,961) DLT,STEPLN */
/* $      WRITE(IPR,959) PHI */
/* $      WRITE(IPR,960) PHIP */
    if ((alo * *dlt > stepln || stepln > hi * *dlt) && amuup - amulo > 0.) {
	goto L170;
    }
/*       IF((ALO*DLT.LE.STEPLN .AND. STEPLN.LE.HI*DLT) .OR. */
/*            (AMUUP-AMULO.LE.0.0D0)) */
/*       THEN */

/*         SC IS ACCEPTABLE HOOKSTEP */

/* $        WRITE(IPR,954) */
    done = TRUE_;
    goto L100;
/*       ELSE */

/*         SC NOT ACCEPTABLE HOOKSTEP.  SELECT NEW AMU */

L170:
/* $        WRITE(IPR,953) */
/* Computing MAX */
    d__1 = amulo, d__2 = *amu - *phi / phip;
    amulo = max(d__1,d__2);
    if (*phi < 0.) {
	amuup = min(amuup,*amu);
    }
    *amu -= stepln * *phi / (*dlt * phip);
/* $        WRITE(IPR,956) AMU */
/* $        WRITE(IPR,957) AMULO */
/* $        WRITE(IPR,958) AMUUP */
    goto L100;
/*       ENDIF */
/*     ENDIF */

/* %951 FORMAT(27H0HOOKST    TAKE NEWTON STEP) */
/* %952 FORMAT(32H0HOOKST    NEWTON STEP NOT TAKEN) */
/* %953 FORMAT(31H HOOKST    SC IS NOT ACCEPTABLE) */
/* %954 FORMAT(27H HOOKST    SC IS ACCEPTABLE) */
/* %955 FORMAT(28H HOOKST    CURRENT STEP (SC)) */
/* %956 FORMAT(18H HOOKST    AMU   =,E20.13) */
/* %957 FORMAT(18H HOOKST    AMULO =,E20.13) */
/* %958 FORMAT(18H HOOKST    AMUUP =,E20.13) */
/* %959 FORMAT(18H HOOKST    PHI   =,E20.13) */
/* %960 FORMAT(18H HOOKST    PHIP  =,E20.13) */
/* %961 FORMAT(18H HOOKST    DLT   =,E20.13/ */
/* %   +       18H HOOKST    STEPLN=,E20.13) */
/* %962 FORMAT(23H0HOOKST    FIND NEW AMU) */
/* %963 FORMAT(14H HOOKST       ,5(E20.13,3X)) */
} /* hookst_ */

/* Subroutine */ int hsnint_(integer *nr, integer *n, doublereal *a, 
	doublereal *sx, integer *method)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i, j, jp1;


/* PURPOSE */
/* ------- */
/* PROVIDE INITIAL HESSIAN WHEN USING SECANT UPDATES */

/* PARAMETERS */
/* ---------- */
/* NR           --> ROW DIMENSION OF MATRIX */
/* N            --> DIMENSION OF PROBLEM */
/* A(N,N)      <--  INITIAL HESSIAN (LOWER TRIANGULAR MATRIX) */
/* SX(N)        --> DIAGONAL SCALING MATRIX FOR X */
/* METHOD       --> ALGORITHM TO USE TO SOLVE MINIMIZATION PROBLEM */
/*                    =1,2 FACTORED SECANT METHOD USED */
/*                    =3   UNFACTORED SECANT METHOD USED */


    /* Parameter adjustments */
    --sx;
    a_dim1 = *nr;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	if (*method == 3) {
	    a[j + j * a_dim1] = sx[j] * sx[j];
	}
	if (*method != 3) {
	    a[j + j * a_dim1] = sx[j];
	}
	if (j == *n) {
	    goto L100;
	}
	jp1 = j + 1;
	i__2 = *n;
	for (i = jp1; i <= i__2; ++i) {
	    a[i + j * a_dim1] = 0.;
/* L90: */
	}
L100:
	;
    }
    return 0;
} /* hsnint_ */

/* Subroutine */ int lltslv_(integer *nr, integer *n, doublereal *a, 
	doublereal *x, doublereal *b)
{
    /* System generated locals */
    integer a_dim1, a_offset;

    /* Local variables */
    extern /* Subroutine */ int bakslv_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *), forslv_(integer *, integer *, 
	    doublereal *, doublereal *, doublereal *);


/* PURPOSE */
/* ------- */
/* SOLVE AX=B WHERE A HAS THE FORM L(L-TRANSPOSE) */
/* BUT ONLY THE LOWER TRIANGULAR PART, L, IS STORED. */

/* PARAMETERS */
/* ---------- */
/* NR           --> ROW DIMENSION OF MATRIX */
/* N            --> DIMENSION OF PROBLEM */
/* A(N,N)       --> MATRIX OF FORM L(L-TRANSPOSE). */
/*                  ON RETURN A IS UNCHANGED. */
/* X(N)        <--  SOLUTION VECTOR */
/* B(N)         --> RIGHT-HAND SIDE VECTOR */

/* NOTE */
/* ---- */
/* IF B IS NOT REQUIRED BY CALLING PROGRAM, THEN */
/* B AND X MAY SHARE THE SAME STORAGE. */


/* FORWARD SOLVE, RESULT IN X */

    /* Parameter adjustments */
    --b;
    --x;
    a_dim1 = *nr;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    forslv_(nr, n, &a[a_offset], &x[1], &b[1]);

/* BACK SOLVE, RESULT IN X */

    bakslv_(nr, n, &a[a_offset], &x[1], &x[1]);
    return 0;
} /* lltslv_ */

/* Subroutine */ int lnsrch_(integer *n, doublereal *x, doublereal *f, 
	doublereal *g, doublereal *p, doublereal *xpls, doublereal *fpls, 
	S_fp fcn, logical *mxtake, integer *iretcd, doublereal *stepmx, 
	doublereal *steptl, doublereal *sx, integer *ipr)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;

    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    static doublereal disc;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal a, b;
    static integer i;
    static doublereal pfpls, t1, t2, t3, almbda, plmbda, tlmbda, rmnlmb;
    extern /* Subroutine */ int sclmul_(integer *, doublereal *, doublereal *,
	     doublereal *);
    static doublereal scl, rln, sln, slp, tmp;

/* PURPOSE */
/* ------- */
/* FIND A NEXT NEWTON ITERATE BY LINE SEARCH. */

/* PARAMETERS */
/* ---------- */
/* N            --> DIMENSION OF PROBLEM */
/* X(N)         --> OLD ITERATE:   X[K-1] */
/* F            --> FUNCTION VALUE AT OLD ITERATE, F(X) */
/* G(N)         --> GRADIENT AT OLD ITERATE, G(X), OR APPROXIMATE */
/* P(N)         --> NON-ZERO NEWTON STEP */
/* XPLS(N)     <--  NEW ITERATE X[K] */
/* FPLS        <--  FUNCTION VALUE AT NEW ITERATE, F(XPLS) */
/* FCN          --> NAME OF SUBROUTINE TO EVALUATE FUNCTION */
/* IRETCD      <--  RETURN CODE */
/* MXTAKE      <--  BOOLEAN FLAG INDICATING STEP OF MAXIMUM LENGTH USED */
/* STEPMX       --> MAXIMUM ALLOWABLE STEP SIZE */
/* STEPTL       --> RELATIVE STEP SIZE AT WHICH SUCCESSIVE ITERATES */
/*                  CONSIDERED CLOSE ENOUGH TO TERMINATE ALGORITHM */
/* SX(N)        --> DIAGONAL SCALING MATRIX FOR X */
/* IPR          --> DEVICE TO WHICH TO SEND OUTPUT */

/* INTERNAL VARIABLES */
/* ------------------ */
/* SLN              NEWTON LENGTH */
/* RLN              RELATIVE LENGTH OF NEWTON STEP */


    /* Parameter adjustments */
    --sx;
    --xpls;
    --p;
    --g;
    --x;

    /* Function Body */
    *ipr = *ipr;
    *mxtake = FALSE_;
    *iretcd = 2;
/* $    WRITE(IPR,954) */
/* $    WRITE(IPR,955) (P(I),I=1,N) */
    tmp = 0.;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	tmp += sx[i] * sx[i] * p[i] * p[i];
/* L5: */
    }
    sln = sqrt(tmp);
    if (sln <= *stepmx) {
	goto L10;
    }

/* NEWTON STEP LONGER THAN MAXIMUM ALLOWED */
    scl = *stepmx / sln;
    sclmul_(n, &scl, &p[1], &p[1]);
    sln = *stepmx;
/* $      WRITE(IPR,954) */
/* $      WRITE(IPR,955) (P(I),I=1,N) */
L10:
    slp = ddot_(n, &g[1], &c__1, &p[1], &c__1);
    rln = 0.;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
/* Computing MAX */
/* Computing MAX */
	d__5 = (d__2 = x[i], abs(d__2)), d__6 = 1. / sx[i];
	d__3 = rln, d__4 = (d__1 = p[i], abs(d__1)) / max(d__5,d__6);
	rln = max(d__3,d__4);
/* L15: */
    }
    rmnlmb = *steptl / rln;
    almbda = 1.;
/* $    WRITE(IPR,952) SLN,SLP,RMNLMB,STEPMX,STEPTL */

/* LOOP */
/* CHECK IF NEW ITERATE SATISFACTORY.  GENERATE NEW LAMBDA IF NECESSARY. 
*/

L100:
    if (*iretcd < 2) {
	return 0;
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	xpls[i] = x[i] + almbda * p[i];
/* L105: */
    }
    (*fcn)(n, &xpls[1], fpls);
/* $    WRITE(IPR,950) ALMBDA */
/* $    WRITE(IPR,951) */
/* $    WRITE(IPR,955) (XPLS(I),I=1,N) */
/* $    WRITE(IPR,953) FPLS */
    if (*fpls > *f + slp * 1e-4 * almbda) {
	goto L130;
    }
/*     IF(FPLS.LE. F+SLP*1.D-4*ALMBDA) */
/*     THEN */

/* SOLUTION FOUND */

    *iretcd = 0;
    if (almbda == 1. && sln > *stepmx * .99) {
	*mxtake = TRUE_;
    }
    goto L100;

/* SOLUTION NOT (YET) FOUND */

/*     ELSE */
L130:
    if (almbda >= rmnlmb) {
	goto L140;
    }
/*       IF(ALMBDA .LT. RMNLMB) */
/*       THEN */

/* NO SATISFACTORY XPLS FOUND SUFFICIENTLY DISTINCT FROM X */

    *iretcd = 1;
    goto L100;
/*       ELSE */

/* CALCULATE NEW LAMBDA */

L140:
    if (almbda != 1.) {
	goto L150;
    }
/*         IF(ALMBDA.EQ.1.0D0) */
/*         THEN */

/* FIRST BACKTRACK: QUADRATIC FIT */

    tlmbda = -slp / ((*fpls - *f - slp) * 2.);
    goto L170;
/*         ELSE */

/* ALL SUBSEQUENT BACKTRACKS: CUBIC FIT */

L150:
    t1 = *fpls - *f - almbda * slp;
    t2 = pfpls - *f - plmbda * slp;
    t3 = 1. / (almbda - plmbda);
    a = t3 * (t1 / (almbda * almbda) - t2 / (plmbda * plmbda));
    b = t3 * (t2 * almbda / (plmbda * plmbda) - t1 * plmbda / (almbda * 
	    almbda));
    disc = b * b - a * 3. * slp;
    if (disc <= b * b) {
	goto L160;
    }
/*           IF(DISC.GT. B*B) */
/*           THEN */

/* ONLY ONE POSITIVE CRITICAL POINT, MUST BE MINIMUM */

    tlmbda = (-b + d_sign(&c_b146, &a) * sqrt(disc)) / (a * 3.);
    goto L165;
/*           ELSE */

/* BOTH CRITICAL POINTS POSITIVE, FIRST IS MINIMUM */

L160:
    tlmbda = (-b - d_sign(&c_b146, &a) * sqrt(disc)) / (a * 3.);
/*           ENDIF */
L165:
    if (tlmbda > almbda * .5f) {
	tlmbda = almbda * .5f;
    }
/*         ENDIF */
L170:
    plmbda = almbda;
    pfpls = *fpls;
    if (tlmbda >= almbda * .1f) {
	goto L180;
    }
/*         IF(TLMBDA.LT.ALMBDA/10.0D0) */
/*         THEN */
    almbda *= .1f;
    goto L190;
/*         ELSE */
L180:
    almbda = tlmbda;
/*         ENDIF */
/*       ENDIF */
/*     ENDIF */
L190:
    goto L100;
/* %950 FORMAT(18H LNSRCH    ALMBDA=,E20.13) */
/* %951 FORMAT(29H LNSRCH    NEW ITERATE (XPLS)) */
/* %952 FORMAT(18H LNSRCH    SLN   =,E20.13/ */
/* %   +       18H LNSRCH    SLP   =,E20.13/ */
/* %   +       18H LNSRCH    RMNLMB=,E20.13/ */
/* %   +       18H LNSRCH    STEPMX=,E20.13/ */
/* %   +       18H LNSRCH    STEPTL=,E20.13) */
/* %953 FORMAT(19H LNSRCH    F(XPLS)=,E20.13) */
/* %954 FORMAT(26H0LNSRCH    NEWTON STEP (P)) */
/* %955 FORMAT(14H LNSRCH       ,5(E20.13,3X)) */
} /* lnsrch_ */

/* Subroutine */ int mvmltl_(integer *nr, integer *n, doublereal *a, 
	doublereal *x, doublereal *y)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i, j;
    static doublereal sum;


/* PURPOSE */
/* ------- */
/* COMPUTE Y=LX */
/* WHERE L IS A LOWER TRIANGULAR MATRIX STORED IN A */

/* PARAMETERS */
/* ---------- */
/* NR           --> ROW DIMENSION OF MATRIX */
/* N            --> DIMENSION OF PROBLEM */
/* A(N,N)       --> LOWER TRIANGULAR (N*N) MATRIX */
/* X(N)         --> OPERAND VECTOR */
/* Y(N)        <--  RESULT VECTOR */

/* NOTE */
/* ---- */
/* X AND Y CANNOT SHARE STORAGE */

    /* Parameter adjustments */
    --y;
    --x;
    a_dim1 = *nr;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	sum = 0.;
	i__2 = i;
	for (j = 1; j <= i__2; ++j) {
	    sum += a[i + j * a_dim1] * x[j];
/* L10: */
	}
	y[i] = sum;
/* L30: */
    }
    return 0;
} /* mvmltl_ */

/* Subroutine */ int mvmlts_(integer *nr, integer *n, doublereal *a, 
	doublereal *x, doublereal *y)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i, j, ip1;
    static doublereal sum;


/* PURPOSE */
/* ------- */
/* COMPUTE Y=AX */
/* WHERE "A" IS A SYMMETRIC (N*N) MATRIX STORED IN ITS LOWER */
/* TRIANGULAR PART AND X,Y ARE N-VECTORS */

/* PARAMETERS */
/* ---------- */
/* NR           --> ROW DIMENSION OF MATRIX */
/* N            --> DIMENSION OF PROBLEM */
/* A(N,N)       --> SYMMETRIC (N*N) MATRIX STORED IN */
/*                  LOWER TRIANGULAR PART AND DIAGONAL */
/* X(N)         --> OPERAND VECTOR */
/* Y(N)        <--  RESULT VECTOR */

/* NOTE */
/* ---- */
/* X AND Y CANNOT SHARE STORAGE. */

    /* Parameter adjustments */
    --y;
    --x;
    a_dim1 = *nr;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	sum = 0.;
	i__2 = i;
	for (j = 1; j <= i__2; ++j) {
	    sum += a[i + j * a_dim1] * x[j];
/* L10: */
	}
	if (i == *n) {
	    goto L25;
	}
	ip1 = i + 1;
	i__2 = *n;
	for (j = ip1; j <= i__2; ++j) {
	    sum += a[j + i * a_dim1] * x[j];
/* L20: */
	}
L25:
	y[i] = sum;
/* L30: */
    }
    return 0;
} /* mvmlts_ */

/* Subroutine */ int mvmltu_(integer *nr, integer *n, doublereal *a, 
	doublereal *x, doublereal *y)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i, j;
    static doublereal sum;


/* PURPOSE */
/* ------- */
/* COMPUTE Y=(L+)X */
/* WHERE L IS A LOWER TRIANGULAR MATRIX STORED IN A */
/* (L-TRANSPOSE (L+) IS TAKEN IMPLICITLY) */

/* PARAMETERS */
/* ---------- */
/* NR           --> ROW DIMENSION OF MATRIX */
/* N            --> DIMENSION OF PROBLEM */
/* A(NR,1)       --> LOWER TRIANGULAR (N*N) MATRIX */
/* X(N)         --> OPERAND VECTOR */
/* Y(N)        <--  RESULT VECTOR */

/* NOTE */
/* ---- */
/* X AND Y CANNOT SHARE STORAGE */

    /* Parameter adjustments */
    --y;
    --x;
    a_dim1 = *nr;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	sum = 0.;
	i__2 = *n;
	for (j = i; j <= i__2; ++j) {
	    sum += a[j + i * a_dim1] * x[j];
/* L10: */
	}
	y[i] = sum;
/* L30: */
    }
    return 0;
} /* mvmltu_ */

/* Subroutine */ int optchk_(integer *n, doublereal *x, doublereal *typsiz, 
	doublereal *sx, doublereal *fscale, doublereal *gradtl, integer *
	itnlim, integer *ndigit, doublereal *epsm, doublereal *dlt, integer *
	method, integer *iexp, integer *iagflg, integer *iahflg, doublereal *
	stepmx, integer *msg, integer *ipr)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal), d_lg10(doublereal *);

    /* Local variables */
    static integer i;
    static doublereal stpsiz;


/* PURPOSE */
/* ------- */
/* CHECK INPUT FOR REASONABLENESS */

/* PARAMETERS */
/* ---------- */
/* N            --> DIMENSION OF PROBLEM */
/* X(N)         --> ON ENTRY, ESTIMATE TO ROOT OF FCN */
/* TYPSIZ(N)   <--> TYPICAL SIZE OF EACH COMPONENT OF X */
/* SX(N)       <--  DIAGONAL SCALING MATRIX FOR X */
/* FSCALE      <--> ESTIMATE OF SCALE OF OBJECTIVE FUNCTION FCN */
/* GRADTL       --> TOLERANCE AT WHICH GRADIENT CONSIDERED CLOSE */
/*                  ENOUGH TO ZERO TO TERMINATE ALGORITHM */
/* ITNLIM      <--> MAXIMUM NUMBER OF ALLOWABLE ITERATIONS */
/* NDIGIT      <--> NUMBER OF GOOD DIGITS IN OPTIMIZATION FUNCTION FCN */
/* EPSM         --> MACHINE EPSILON */
/* DLT         <--> TRUST REGION RADIUS */
/* METHOD      <--> ALGORITHM INDICATOR */
/* IEXP        <--> EXPENSE FLAG */
/* IAGFLG      <--> =1 IF ANALYTIC GRADIENT SUPPLIED */
/* IAHFLG      <--> =1 IF ANALYTIC HESSIAN SUPPLIED */
/* STEPMX      <--> MAXIMUM STEP SIZE */
/* MSG         <--> MESSAGE AND ERROR CODE */
/* IPR          --> DEVICE TO WHICH TO SEND OUTPUT */


/* CHECK THAT PARAMETERS ONLY TAKE ON ACCEPTABLE VALUES. */
/* IF NOT, SET THEM TO DEFAULT VALUES. */
    /* Parameter adjustments */
    --sx;
    --typsiz;
    --x;

    /* Function Body */
    if (*method < 1 || *method > 3) {
	*method = 1;
    }
    if (*iagflg != 1) {
	*iagflg = 0;
    }
    if (*iahflg != 1) {
	*iahflg = 0;
    }
    if (*iexp != 0) {
	*iexp = 1;
    }
    if (*msg / 2 % 2 == 1 && *iagflg == 0) {
	goto L830;
    }
    if (*msg / 4 % 2 == 1 && *iahflg == 0) {
	goto L835;
    }

/* CHECK DIMENSION OF PROBLEM */

    if (*n <= 0) {
	goto L805;
    }
    if (*n == 1 && *msg % 2 == 0) {
	goto L810;
    }

/* COMPUTE SCALE MATRIX */

    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	if (typsiz[i] == 0.) {
	    typsiz[i] = 1.f;
	}
	if (typsiz[i] < 0.) {
	    typsiz[i] = -typsiz[i];
	}
	sx[i] = 1. / typsiz[i];
/* L10: */
    }

/* CHECK MAXIMUM STEP SIZE */

    if (*stepmx > 0.) {
	goto L20;
    }
    stpsiz = 0.;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	stpsiz += x[i] * x[i] * sx[i] * sx[i];
/* L15: */
    }
    stpsiz = sqrt(stpsiz);
/* Computing MAX */
    d__1 = stpsiz * 1e3;
    *stepmx = max(d__1,1e3);
L20:
/* CHECK FUNCTION SCALE */
    if (*fscale == 0.) {
	*fscale = 1.;
    }
    if (*fscale < 0.) {
	*fscale = -(*fscale);
    }

/* CHECK GRADIENT TOLERANCE */
    if (*gradtl < 0.) {
	goto L815;
    }

/* CHECK ITERATION LIMIT */
    if (*itnlim <= 0) {
	goto L820;
    }

/* CHECK NUMBER OF DIGITS OF ACCURACY IN FUNCTION FCN */
    if (*ndigit == 0) {
	goto L825;
    }
    if (*ndigit < 0) {
	*ndigit = (integer) (-d_lg10(epsm));
    }

/* CHECK TRUST REGION RADIUS */
    if (*dlt <= 0.) {
	*dlt = -1.;
    }
    if (*dlt > *stepmx) {
	*dlt = *stepmx;
    }
    return 0;

/* ERROR EXITS */

/* %805 WRITE(IPR,901) N */
/* %    MSG=-1 */
L805:
    *msg = -1;
    goto L895;
/* %810 WRITE(IPR,902) */
/* %    MSG=-2 */
L810:
    *msg = -2;
    goto L895;
/* %815 WRITE(IPR,903) GRADTL */
/* %    MSG=-3 */
L815:
    *msg = -3;
    goto L895;
/* %820 WRITE(IPR,904) ITNLIM */
/* %    MSG=-4 */
L820:
    *msg = -4;
    goto L895;
/* %825 WRITE(IPR,905) NDIGIT */
/* %    MSG=-5 */
L825:
    *msg = -5;
    goto L895;
/* %830 WRITE(IPR,906) MSG,IAGFLG */
/* %    MSG=-6 */
L830:
    *msg = -6;
    goto L895;
/* %835 WRITE(IPR,907) MSG,IAHFLG */
/* %    MSG=-7 */
L835:
    *msg = -7;
L895:
    return 0;
/* %901 FORMAT(32H0OPTCHK    ILLEGAL DIMENSION, N=,I5) */
/* %902 FORMAT(55H0OPTCHK    +++ WARNING +++  THIS PACKAGE IS INEFFICIENT,
 */
/* %   +       26H FOR PROBLEMS OF SIZE N=1./ */
/* %   +       48H OPTCHK    CHECK INSTALLATION LIBRARIES FOR MORE, */
/* %   +       22H APPROPRIATE ROUTINES./ */
/* %   +       41H OPTCHK    IF NONE, SET MSG AND RESUBMIT.) */
/* %903 FORMAT(38H0OPTCHK    ILLEGAL TOLERANCE.  GRADTL=,E20.13) */
/* %904 FORMAT(44H0OPTCHK    ILLEGAL ITERATION LIMIT.  ITNLIM=,I5) */
/* %905 FORMAT(52H0OPTCHK    MINIMIZATION FUNCTION HAS NO GOOD DIGITS., */
/* %   +        9H  NDIGIT=,I5) */
/* %906 FORMAT(50H0OPTCHK    USER REQUESTS THAT ANALYTIC GRADIENT BE, */
/* %   +       33H ACCEPTED AS PROPERLY CODED (MSG=,I5, 2H),/ */
/* %   +       45H OPTCHK    BUT ANALYTIC GRADIENT NOT SUPPLIED, */
/* %   +        9H (IAGFLG=,I5, 2H).) */
/* %907 FORMAT(49H0OPTCHK    USER REQUESTS THAT ANALYTIC HESSIAN BE, */
/* %   +       33H ACCEPTED AS PROPERLY CODED (MSG=,I5, 2H),/ */
/* %   +       44H OPTCHK    BUT ANALYTIC HESSIAN NOT SUPPLIED, */
/* %   +        9H (IAHFLG=,I5, 2H).) */
} /* optchk_ */

/* Subroutine */ int optdrv_(integer *nr, integer *n, doublereal *x, S_fp fcn,
	 S_fp d1fcn, S_fp d2fcn, doublereal *typsiz, doublereal *fscale, 
	integer *method, integer *iexp, integer *msg, integer *ndigit, 
	integer *itnlim, integer *iagflg, integer *iahflg, integer *ipr, 
	doublereal *dlt, doublereal *gradtl, doublereal *stepmx, doublereal *
	steptl, doublereal *xpls, doublereal *fpls, doublereal *gpls, integer 
	*itrmcd, doublereal *a, doublereal *udiag, doublereal *g, doublereal *
	p, doublereal *sx, doublereal *wrk0, doublereal *wrk1, doublereal *
	wrk2, doublereal *wrk3)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_di(doublereal *, integer *), sqrt(doublereal);

    /* Local variables */
    static doublereal dltp, epsm, phip0, f;
    static integer i;
    extern doublereal d1mach_(integer *);
    extern /* Subroutine */ int secfac_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, logical *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), grdchk_(
	    integer *, doublereal *, S_fp, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *), heschk_(
	    integer *, integer *, doublereal *, S_fp, S_fp, S_fp, doublereal *
	    , doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *);
    static integer iretcd;
    static doublereal analtl;
    extern /* Subroutine */ int sndofd_(integer *, integer *, doublereal *, 
	    S_fp, doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), chlhsn_(integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), fstocd_(
	    integer *, doublereal *, S_fp, doublereal *, doublereal *, 
	    doublereal *), secunf_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    logical *, doublereal *, doublereal *, doublereal *), fstofd_(
	    integer *, integer *, integer *, doublereal *, S_fp, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, integer *
	    );
    static integer icscmx;
    extern /* Subroutine */ int dogdrv_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, S_fp, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, logical *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *), optchk_(
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, integer *);
    static logical mxtake;
    static doublereal dlpsav, phisav, dltsav, amusav;
    static integer itncnt;
    extern /* Subroutine */ int lnsrch_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, S_fp, 
	    logical *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *);
    static doublereal phpsav;
    extern /* Subroutine */ int hookdr_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, S_fp, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, logical *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *), hsnint_(integer *, integer *, doublereal *, 
	    doublereal *, integer *);
    static logical noupdt;
    extern /* Subroutine */ int result_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, integer *), lltslv_(integer *, integer *, doublereal *
	    , doublereal *, doublereal *), optstp_(integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, integer *, logical *, integer *, integer *);
    static doublereal phi, amu, rnf, wrk;


/* PURPOSE */
/* ------- */
/* DRIVER FOR NON-LINEAR OPTIMIZATION PROBLEM */

/* PARAMETERS */
/* ---------- */
/* NR           --> ROW DIMENSION OF MATRIX */
/* N            --> DIMENSION OF PROBLEM */
/* X(N)         --> ON ENTRY: ESTIMATE TO A ROOT OF FCN */
/* FCN          --> NAME OF SUBROUTINE TO EVALUATE OPTIMIZATION FUNCTION 
*/
/*                  MUST BE DECLARED EXTERNAL IN CALLING ROUTINE */
/*                            FCN: R(N) --> R(1) */
/* D1FCN        --> (OPTIONAL) NAME OF SUBROUTINE TO EVALUATE GRADIENT */
/*                  OF FCN.  MUST BE DECLARED EXTERNAL IN CALLING ROUTINE 
*/
/* D2FCN        --> (OPTIONAL) NAME OF SUBROUTINE TO EVALUATE HESSIAN OF 
*/
/*                  OF FCN.  MUST BE DECLARED EXTERNAL IN CALLING ROUTINE 
*/
/* TYPSIZ(N)    --> TYPICAL SIZE FOR EACH COMPONENT OF X */
/* FSCALE       --> ESTIMATE OF SCALE OF OBJECTIVE FUNCTION */
/* METHOD       --> ALGORITHM TO USE TO SOLVE MINIMIZATION PROBLEM */
/*                    =1 LINE SEARCH */
/*                    =2 DOUBLE DOGLEG */
/*                    =3 MORE-HEBDON */
/* IEXP         --> =1 IF OPTIMIZATION FUNCTION FCN IS EXPENSIVE TO */
/*                  EVALUATE, =0 OTHERWISE.  IF SET THEN HESSIAN WILL */
/*                  BE EVALUATED BY SECANT UPDATE INSTEAD OF */
/*                  ANALYTICALLY OR BY FINITE DIFFERENCES */
/* MSG         <--> ON INPUT:  (.GT.0) MESSAGE TO INHIBIT CERTAIN */
/*                    AUTOMATIC CHECKS */
/*                  ON OUTPUT: (.LT.0) ERROR CODE; =0 NO ERROR */
/* NDIGIT       --> NUMBER OF GOOD DIGITS IN OPTIMIZATION FUNCTION FCN */
/* ITNLIM       --> MAXIMUM NUMBER OF ALLOWABLE ITERATIONS */
/* IAGFLG       --> =1 IF ANALYTIC GRADIENT SUPPLIED */
/* IAHFLG       --> =1 IF ANALYTIC HESSIAN SUPPLIED */
/* IPR          --> DEVICE TO WHICH TO SEND OUTPUT */
/* DLT          --> TRUST REGION RADIUS */
/* GRADTL       --> TOLERANCE AT WHICH GRADIENT CONSIDERED CLOSE */
/*                  ENOUGH TO ZERO TO TERMINATE ALGORITHM */
/* STEPMX       --> MAXIMUM ALLOWABLE STEP SIZE */
/* STEPTL       --> RELATIVE STEP SIZE AT WHICH SUCCESSIVE ITERATES */
/*                  CONSIDERED CLOSE ENOUGH TO TERMINATE ALGORITHM */
/* XPLS(N)     <--> ON EXIT:  XPLS IS LOCAL MINIMUM */
/* FPLS        <--> ON EXIT:  FUNCTION VALUE AT SOLUTION, XPLS */
/* GPLS(N)     <--> ON EXIT:  GRADIENT AT SOLUTION XPLS */
/* ITRMCD      <--  TERMINATION CODE */
/* A(N,N)       --> WORKSPACE FOR HESSIAN (OR ESTIMATE) */
/*                  AND ITS CHOLESKY DECOMPOSITION */
/* UDIAG(N)     --> WORKSPACE [FOR DIAGONAL OF HESSIAN] */
/* G(N)         --> WORKSPACE (FOR GRADIENT AT CURRENT ITERATE) */
/* P(N)         --> WORKSPACE FOR STEP */
/* SX(N)        --> WORKSPACE (FOR DIAGONAL SCALING MATRIX) */
/* WRK0(N)      --> WORKSPACE */
/* WRK1(N)      --> WORKSPACE */
/* WRK2(N)      --> WORKSPACE */
/* WRK3(N)      --> WORKSPACE */


/* INTERNAL VARIABLES */
/* ------------------ */
/* ANALTL           TOLERANCE FOR COMPARISON OF ESTIMATED AND */
/*                  ANALYTICAL GRADIENTS AND HESSIANS */
/* EPSM             MACHINE EPSILON */
/* F                FUNCTION VALUE: FCN(X) */
/* ITNCNT           CURRENT ITERATION, K */
/* RNF              RELATIVE NOISE IN OPTIMIZATION FUNCTION FCN. */
/*                       NOISE=10.**(-NDIGIT) */


/* INITIALIZATION */
/* -------------- */
    /* Parameter adjustments */
    --wrk3;
    --wrk2;
    --wrk1;
    --wrk0;
    --sx;
    --p;
    --g;
    --udiag;
    a_dim1 = *nr;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --gpls;
    --xpls;
    --typsiz;
    --x;

    /* Function Body */
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	p[i] = 0.;
/* L10: */
    }
    itncnt = 0;
    iretcd = -1;
    epsm = d1mach_(&c__4);
    optchk_(n, &x[1], &typsiz[1], &sx[1], fscale, gradtl, itnlim, ndigit, &
	    epsm, dlt, method, iexp, iagflg, iahflg, stepmx, msg, ipr);
    if (*msg < 0) {
	return 0;
    }
/* Computing MAX */
    i__1 = -(*ndigit);
    d__1 = pow_di(&c_b3, &i__1);
    rnf = max(d__1,epsm);
/* Computing MAX */
    d__1 = .01, d__2 = sqrt(rnf);
    analtl = max(d__1,d__2);

/* %    IF(MOD(MSG/8,2).EQ.1) GO TO 15 */
/* %    WRITE(IPR,901) */
/* %    WRITE(IPR,900) (TYPSIZ(I),I=1,N) */
/* %    WRITE(IPR,902) */
/* %    WRITE(IPR,900) (SX(I),I=1,N) */
/* %    WRITE(IPR,903) FSCALE */
/* %    WRITE(IPR,904) NDIGIT,IAGFLG,IAHFLG,IEXP,METHOD,ITNLIM,EPSM */
/* %    WRITE(IPR,905) STEPMX,STEPTL,GRADTL,DLT,RNF,ANALTL */
/* % 15 CONTINUE */

/* EVALUATE FCN(X) */

    (*fcn)(n, &x[1], &f);

/* EVALUATE ANALYTIC OR FINITE DIFFERENCE GRADIENT AND CHECK ANALYTIC */
/* GRADIENT, IF REQUESTED. */

    if (*iagflg == 1) {
	goto L20;
    }
/*     IF (IAGFLG .EQ. 0) */
/*     THEN */
    fstofd_(&c__1, &c__1, n, &x[1], (S_fp)fcn, &f, &g[1], &sx[1], &rnf, &wrk, 
	    &c__1);
    goto L25;

L20:
    (*d1fcn)(n, &x[1], &g[1]);
    if (*msg / 2 % 2 == 1) {
	goto L25;
    }
/*     IF (MOD(MSG/2,2).EQ.0) */
/*     THEN */
    grdchk_(n, &x[1], (S_fp)fcn, &f, &g[1], &typsiz[1], &sx[1], fscale, &rnf, 
	    &analtl, &wrk1[1], msg, ipr);
    if (*msg < 0) {
	return 0;
    }
L25:

    optstp_(n, &x[1], &f, &g[1], &wrk1[1], &itncnt, &icscmx, itrmcd, gradtl, 
	    steptl, &sx[1], fscale, itnlim, &iretcd, &mxtake, ipr, msg);
    if (*itrmcd != 0) {
	goto L700;
    }

    if (*iexp != 1) {
	goto L80;
    }

/* IF OPTIMIZATION FUNCTION EXPENSIVE TO EVALUATE (IEXP=1), THEN */
/* HESSIAN WILL BE OBTAINED BY SECANT UPDATES.  GET INITIAL HESSIAN. */

    hsnint_(nr, n, &a[a_offset], &sx[1], method);
    goto L90;
L80:

/* EVALUATE ANALYTIC OR FINITE DIFFERENCE HESSIAN AND CHECK ANALYTIC */
/* HESSIAN IF REQUESTED (ONLY IF USER-SUPPLIED ANALYTIC HESSIAN */
/* ROUTINE D2FCN FILLS ONLY LOWER TRIANGULAR PART AND DIAGONAL OF A). */

    if (*iahflg == 1) {
	goto L82;
    }
/*     IF (IAHFLG .EQ. 0) */
/*     THEN */
    if (*iagflg == 1) {
	fstofd_(nr, n, n, &x[1], (S_fp)d1fcn, &g[1], &a[a_offset], &sx[1], &
		rnf, &wrk1[1], &c__3);
    }
    if (*iagflg != 1) {
	sndofd_(nr, n, &x[1], (S_fp)fcn, &f, &a[a_offset], &sx[1], &rnf, &
		wrk1[1], &wrk2[1]);
    }
    goto L88;

/*     ELSE */
L82:
    if (*msg / 4 % 2 == 0) {
	goto L85;
    }
/*        IF (MOD(MSG/4, 2) .EQ. 1) */
/*        THEN */
    (*d2fcn)(nr, n, &x[1], &a[a_offset]);
    goto L88;

/*        ELSE */
L85:
    heschk_(nr, n, &x[1], (S_fp)fcn, (S_fp)d1fcn, (S_fp)d2fcn, &f, &g[1], &a[
	    a_offset], &typsiz[1], &sx[1], &rnf, &analtl, iagflg, &udiag[1], &
	    wrk1[1], &wrk2[1], msg, ipr);

/*           HESCHK EVALUATES D2FCN AND CHECKS IT AGAINST THE FINITE */
/*           DIFFERENCE HESSIAN WHICH IT CALCULATES BY CALLING FSTOFD */
/*           (IF IAGFLG .EQ. 1) OR SNDOFD (OTHERWISE). */

    if (*msg < 0) {
	return 0;
    }
L88:

L90:
    if (*msg / 8 % 2 == 0) {
	result_(nr, n, &x[1], &f, &g[1], &a[a_offset], &p[1], &itncnt, &c__1, 
		ipr);
    }


/* ITERATION */
/* --------- */
L100:
    ++itncnt;

/* FIND PERTURBED LOCAL MODEL HESSIAN AND ITS LL+ DECOMPOSITION */
/* (SKIP THIS STEP IF LINE SEARCH OR DOGSTEP TECHNIQUES BEING USED WITH */
/* SECANT UPDATES.  CHOLESKY DECOMPOSITION L ALREADY OBTAINED FROM */
/* SECFAC.) */

    if (*iexp == 1 && *method != 3) {
	goto L105;
    }
L103:
    chlhsn_(nr, n, &a[a_offset], &epsm, &sx[1], &udiag[1]);
L105:

/* SOLVE FOR NEWTON STEP:  AP=-G */

    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	wrk1[i] = -g[i];
/* L110: */
    }
    lltslv_(nr, n, &a[a_offset], &p[1], &wrk1[1]);

/* DECIDE WHETHER TO ACCEPT NEWTON STEP  XPLS=X + P */
/* OR TO CHOOSE XPLS BY A GLOBAL STRATEGY. */

    if (*iagflg != 0 || *method == 1) {
	goto L111;
    }
    dltsav = *dlt;
    if (*method == 2) {
	goto L111;
    }
    amusav = amu;
    dlpsav = dltp;
    phisav = phi;
    phpsav = phip0;
L111:
    if (*method == 1) {
	lnsrch_(n, &x[1], &f, &g[1], &p[1], &xpls[1], fpls, (S_fp)fcn, &
		mxtake, &iretcd, stepmx, steptl, &sx[1], ipr);
    }
    if (*method == 2) {
	dogdrv_(nr, n, &x[1], &f, &g[1], &a[a_offset], &p[1], &xpls[1], fpls, 
		(S_fp)fcn, &sx[1], stepmx, steptl, dlt, &iretcd, &mxtake, &
		wrk0[1], &wrk1[1], &wrk2[1], &wrk3[1], ipr);
    }
    if (*method == 3) {
	hookdr_(nr, n, &x[1], &f, &g[1], &a[a_offset], &udiag[1], &p[1], &
		xpls[1], fpls, (S_fp)fcn, &sx[1], stepmx, steptl, dlt, &
		iretcd, &mxtake, &amu, &dltp, &phi, &phip0, &wrk0[1], &wrk1[1]
		, &wrk2[1], &epsm, &itncnt, ipr);
    }

/* IF COULD NOT FIND SATISFACTORY STEP AND FORWARD DIFFERENCE */
/* GRADIENT WAS USED, RETRY USING CENTRAL DIFFERENCE GRADIENT. */

    if (iretcd != 1 || *iagflg != 0) {
	goto L112;
    }
/*     IF (IRETCD .EQ. 1 .AND. IAGFLG .EQ. 0) */
/*     THEN */

/*        SET IAGFLG FOR CENTRAL DIFFERENCES */

    *iagflg = -1;
/* %       WRITE(IPR,906) ITNCNT */

    fstocd_(n, &x[1], (S_fp)fcn, &sx[1], &rnf, &g[1]);
    if (*method == 1) {
	goto L105;
    }
    *dlt = dltsav;
    if (*method == 2) {
	goto L105;
    }
    amu = amusav;
    dltp = dlpsav;
    phi = phisav;
    phip0 = phpsav;
    goto L103;
/*     ENDIF */

/* CALCULATE STEP FOR OUTPUT */

L112:
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	p[i] = xpls[i] - x[i];
/* L114: */
    }

/* CALCULATE GRADIENT AT XPLS */

    if (*iagflg == -1) {
	goto L116;
    }
    if (*iagflg == 0) {
	goto L118;
    }

/* ANALYTIC GRADIENT */
    (*d1fcn)(n, &xpls[1], &gpls[1]);
    goto L120;

/* CENTRAL DIFFERENCE GRADIENT */
L116:
    fstocd_(n, &xpls[1], (S_fp)fcn, &sx[1], &rnf, &gpls[1]);
    goto L120;

/* FORWARD DIFFERENCE GRADIENT */
L118:
    fstofd_(&c__1, &c__1, n, &xpls[1], (S_fp)fcn, fpls, &gpls[1], &sx[1], &
	    rnf, &wrk, &c__1);
L120:

/* CHECK WHETHER STOPPING CRITERIA SATISFIED */

    optstp_(n, &xpls[1], fpls, &gpls[1], &x[1], &itncnt, &icscmx, itrmcd, 
	    gradtl, steptl, &sx[1], fscale, itnlim, &iretcd, &mxtake, ipr, 
	    msg);
    if (*itrmcd != 0) {
	goto L690;
    }

/* EVALUATE HESSIAN AT XPLS */

    if (*iexp == 0) {
	goto L130;
    }
    if (*method == 3) {
	secunf_(nr, n, &x[1], &g[1], &a[a_offset], &udiag[1], &xpls[1], &gpls[
		1], &epsm, &itncnt, &rnf, iagflg, &noupdt, &wrk1[1], &wrk2[1],
		 &wrk3[1]);
    }
    if (*method != 3) {
	secfac_(nr, n, &x[1], &g[1], &a[a_offset], &xpls[1], &gpls[1], &epsm, 
		&itncnt, &rnf, iagflg, &noupdt, &wrk0[1], &wrk1[1], &wrk2[1], 
		&wrk3[1]);
    }
    goto L150;
L130:
    if (*iahflg == 1) {
	goto L140;
    }
    if (*iagflg == 1) {
	fstofd_(nr, n, n, &xpls[1], (S_fp)d1fcn, &gpls[1], &a[a_offset], &sx[
		1], &rnf, &wrk1[1], &c__3);
    }
    if (*iagflg != 1) {
	sndofd_(nr, n, &xpls[1], (S_fp)fcn, fpls, &a[a_offset], &sx[1], &rnf, 
		&wrk1[1], &wrk2[1]);
    }
    goto L150;
L140:
    (*d2fcn)(nr, n, &xpls[1], &a[a_offset]);
L150:
    if (*msg / 16 % 2 == 1) {
	result_(nr, n, &xpls[1], fpls, &gpls[1], &a[a_offset], &p[1], &itncnt,
		 &c__1, ipr);
    }

/* X <-- XPLS  AND  G <-- GPLS  AND  F <-- FPLS */

    f = *fpls;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	x[i] = xpls[i];
	g[i] = gpls[i];
/* L160: */
    }
    goto L100;

/* TERMINATION */
/* ----------- */
/* RESET XPLS,FPLS,GPLS,  IF PREVIOUS ITERATE SOLUTION */

L690:
    if (*itrmcd != 3) {
	goto L710;
    }
L700:
    *fpls = f;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	xpls[i] = x[i];
	gpls[i] = g[i];
/* L705: */
    }

/* PRINT RESULTS */

L710:
    if (*msg / 8 % 2 == 0) {
	result_(nr, n, &xpls[1], fpls, &gpls[1], &a[a_offset], &p[1], &itncnt,
		 &c__0, ipr);
    }
    *msg = 0;
    return 0;

/* %900 FORMAT(14H OPTDRV       ,5(E20.13,3X)) */
/* %901 FORMAT(20H0OPTDRV    TYPICAL X) */
/* %902 FORMAT(40H OPTDRV    DIAGONAL SCALING MATRIX FOR X) */
/* %903 FORMAT(22H OPTDRV    TYPICAL F =,E20.13) */
/* %904 FORMAT(40H0OPTDRV    NUMBER OF GOOD DIGITS IN FCN=,I5/ */
/* %   +       27H OPTDRV    GRADIENT FLAG  =,I5,18H   (=1 IF ANALYTIC, */
/* %   +       19H GRADIENT SUPPLIED)/ */
/* %   +       27H OPTDRV    HESSIAN FLAG   =,I5,18H   (=1 IF ANALYTIC, */
/* %   +       18H HESSIAN SUPPLIED)/ */
/* %   +       27H OPTDRV    EXPENSE FLAG   =,I5, 9H   (=1 IF, */
/* %   +       45H MINIMIZATION FUNCTION EXPENSIVE TO EVALUATE)/ */
/* %   +       27H OPTDRV    METHOD TO USE  =,I5,19H   (=1,2,3 FOR LINE, 
*/
/* %   +       49H SEARCH, DOUBLE DOGLEG, MORE-HEBDON RESPECTIVELY)/ */
/* %   +       27H OPTDRV    ITERATION LIMIT=,I5/ */
/* %   +       27H OPTDRV    MACHINE EPSILON=,E20.13) */
/* %905 FORMAT(30H0OPTDRV    MAXIMUM STEP SIZE =,E20.13/ */
/* %   +       30H OPTDRV    STEP TOLERANCE    =,E20.13/ */
/* %   +       30H OPTDRV    GRADIENT TOLERANCE=,E20.13/ */
/* %   +       30H OPTDRV    TRUST REG RADIUS  =,E20.13/ */
/* %   +       30H OPTDRV    REL NOISE IN FCN  =,E20.13/ */
/* %   +       30H OPTDRV    ANAL-FD TOLERANCE =,E20.13) */
/* %906 FORMAT(52H OPTDRV    SHIFT FROM FORWARD TO CENTRAL DIFFERENCES, */
/* %   1   14H IN ITERATION , I5) */
} /* optdrv_ */

/* Subroutine */ int optif0_(integer *nr, integer *n, doublereal *x, S_fp fcn,
	 doublereal *xpls, doublereal *fpls, doublereal *gpls, integer *
	itrmcd, doublereal *a, doublereal *wrk)
{
    /* System generated locals */
    integer a_dim1, a_offset, wrk_dim1, wrk_offset;

    /* Local variables */
    static integer iexp;
    extern /* Subroutine */ int d1fcn_(integer *, doublereal *, doublereal *),
	     d2fcn_(integer *, integer *, doublereal *, doublereal *);
    static integer iagflg, iahflg;
    static doublereal fscale, gradtl;
    static integer ndigit;
    extern /* Subroutine */ int dfault_(integer *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static integer method, itnlim;
    static doublereal steptl;
    extern /* Subroutine */ int optdrv_(integer *, integer *, doublereal *, 
	    S_fp, S_fp, S_fp, doublereal *, doublereal *, integer *, integer *
	    , integer *, integer *, integer *, integer *, integer *, integer *
	    , doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static doublereal stepmx, dlt;
    static integer msg, ipr;


/* PURPOSE */
/* ------- */
/* PROVIDE SIMPLEST INTERFACE TO MINIMIZATION PACKAGE. */
/* USER HAS NO CONTROL OVER OPTIONS. */

/* PARAMETERS */
/* ---------- */
/* NR           --> ROW DIMENSION OF MATRIX */
/* N            --> DIMENSION OF PROBLEM */
/* X(N)         --> INITIAL ESTIMATE OF MINIMUM */
/* FCN          --> NAME OF ROUTINE TO EVALUATE MINIMIZATION FUNCTION. */
/*                  MUST BE DECLARED EXTERNAL IN CALLING ROUTINE. */
/* XPLS(N)     <--  LOCAL MINIMUM */
/* FPLS        <--  FUNCTION VALUE AT LOCAL MINIMUM XPLS */
/* GPLS(N)     <--  GRADIENT AT LOCAL MINIMUM XPLS */
/* ITRMCD      <--  TERMINATION CODE */
/* A(N,N)       --> WORKSPACE */
/* WRK(N,9)     --> WORKSPACE */


/* EQUIVALENCE WRK(N,1) = UDIAG(N) */
/*             WRK(N,2) = G(N) */
/*             WRK(N,3) = P(N) */
/*             WRK(N,4) = TYPSIZ(N) */
/*             WRK(N,5) = SX(N) */
/*             WRK(N,6) = WRK0(N) */
/*             WRK(N,7) = WRK1(N) */
/*             WRK(N,8) = WRK2(N) */
/*             WRK(N,9) = WRK3(N) */

    /* Parameter adjustments */
    wrk_dim1 = *nr;
    wrk_offset = wrk_dim1 + 1;
    wrk -= wrk_offset;
    a_dim1 = *nr;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --gpls;
    --xpls;
    --x;

    /* Function Body */
    dfault_(n, &x[1], &wrk[(wrk_dim1 << 2) + 1], &fscale, &method, &iexp, &
	    msg, &ndigit, &itnlim, &iagflg, &iahflg, &ipr, &dlt, &gradtl, &
	    stepmx, &steptl);
    optdrv_(nr, n, &x[1], (S_fp)fcn, (S_fp)d1fcn_, (S_fp)d2fcn_, &wrk[(
	    wrk_dim1 << 2) + 1], &fscale, &method, &iexp, &msg, &ndigit, &
	    itnlim, &iagflg, &iahflg, &ipr, &dlt, &gradtl, &stepmx, &steptl, &
	    xpls[1], fpls, &gpls[1], itrmcd, &a[a_offset], &wrk[wrk_dim1 + 1],
	     &wrk[(wrk_dim1 << 1) + 1], &wrk[wrk_dim1 * 3 + 1], &wrk[wrk_dim1 
	    * 5 + 1], &wrk[wrk_dim1 * 6 + 1], &wrk[wrk_dim1 * 7 + 1], &wrk[(
	    wrk_dim1 << 3) + 1], &wrk[wrk_dim1 * 9 + 1]);
    return 0;
} /* optif0_ */

/* Subroutine */ int optif9_(integer *nr, integer *n, doublereal *x, S_fp fcn,
	 S_fp d1fcn, S_fp d2fcn, doublereal *typsiz, doublereal *fscale, 
	integer *method, integer *iexp, integer *msg, integer *ndigit, 
	integer *itnlim, integer *iagflg, integer *iahflg, integer *ipr, 
	doublereal *dlt, doublereal *gradtl, doublereal *stepmx, doublereal *
	steptl, doublereal *xpls, doublereal *fpls, doublereal *gpls, integer 
	*itrmcd, doublereal *a, doublereal *wrk)
{
    /* System generated locals */
    integer a_dim1, a_offset, wrk_dim1, wrk_offset;

    /* Local variables */
    extern /* Subroutine */ int optdrv_(integer *, integer *, doublereal *, 
	    S_fp, S_fp, S_fp, doublereal *, doublereal *, integer *, integer *
	    , integer *, integer *, integer *, integer *, integer *, integer *
	    , doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);


/* PURPOSE */
/* ------- */
/* PROVIDE COMPLETE INTERFACE TO MINIMIZATION PACKAGE. */
/* USER HAS FULL CONTROL OVER OPTIONS. */

/* PARAMETERS */
/* ---------- */
/* NR           --> ROW DIMENSION OF MATRIX */
/* N            --> DIMENSION OF PROBLEM */
/* X(N)         --> ON ENTRY: ESTIMATE TO A ROOT OF FCN */
/* FCN          --> NAME OF SUBROUTINE TO EVALUATE OPTIMIZATION FUNCTION 
*/
/*                  MUST BE DECLARED EXTERNAL IN CALLING ROUTINE */
/*                            FCN: R(N) --> R(1) */
/* D1FCN        --> (OPTIONAL) NAME OF SUBROUTINE TO EVALUATE GRADIENT */
/*                  OF FCN.  MUST BE DECLARED EXTERNAL IN CALLING ROUTINE 
*/
/* D2FCN        --> (OPTIONAL) NAME OF SUBROUTINE TO EVALUATE HESSIAN OF 
*/
/*                  OF FCN.  MUST BE DECLARED EXTERNAL IN CALLING ROUTINE 
*/
/* TYPSIZ(N)    --> TYPICAL SIZE FOR EACH COMPONENT OF X */
/* FSCALE       --> ESTIMATE OF SCALE OF OBJECTIVE FUNCTION */
/* METHOD       --> ALGORITHM TO USE TO SOLVE MINIMIZATION PROBLEM */
/*                    =1 LINE SEARCH */
/*                    =2 DOUBLE DOGLEG */
/*                    =3 MORE-HEBDON */
/* IEXP         --> =1 IF OPTIMIZATION FUNCTION FCN IS EXPENSIVE TO */
/*                  EVALUATE, =0 OTHERWISE.  IF SET THEN HESSIAN WILL */
/*                  BE EVALUATED BY SECANT UPDATE INSTEAD OF */
/*                  ANALYTICALLY OR BY FINITE DIFFERENCES */
/* MSG         <--> ON INPUT:  (.GT.0) MESSAGE TO INHIBIT CERTAIN */
/*                    AUTOMATIC CHECKS */
/*                  ON OUTPUT: (.LT.0) ERROR CODE; =0 NO ERROR */
/* NDIGIT       --> NUMBER OF GOOD DIGITS IN OPTIMIZATION FUNCTION FCN */
/* ITNLIM       --> MAXIMUM NUMBER OF ALLOWABLE ITERATIONS */
/* IAGFLG       --> =1 IF ANALYTIC GRADIENT SUPPLIED */
/* IAHFLG       --> =1 IF ANALYTIC HESSIAN SUPPLIED */
/* IPR          --> DEVICE TO WHICH TO SEND OUTPUT */
/* DLT          --> TRUST REGION RADIUS */
/* GRADTL       --> TOLERANCE AT WHICH GRADIENT CONSIDERED CLOSE */
/*                  ENOUGH TO ZERO TO TERMINATE ALGORITHM */
/* STEPMX       --> MAXIMUM ALLOWABLE STEP SIZE */
/* STEPTL       --> RELATIVE STEP SIZE AT WHICH SUCCESSIVE ITERATES */
/*                  CONSIDERED CLOSE ENOUGH TO TERMINATE ALGORITHM */
/* XPLS(N)     <--> ON EXIT:  XPLS IS LOCAL MINIMUM */
/* FPLS        <--> ON EXIT:  FUNCTION VALUE AT SOLUTION, XPLS */
/* GPLS(N)     <--> ON EXIT:  GRADIENT AT SOLUTION XPLS */
/* ITRMCD      <--  TERMINATION CODE */
/* A(N,N)       --> WORKSPACE FOR HESSIAN (OR ESTIMATE) */
/*                  AND ITS CHOLESKY DECOMPOSITION */
/* WRK(N,8)     --> WORKSPACE */


/* EQUIVALENCE WRK(N,1) = UDIAG(N) */
/*             WRK(N,2) = G(N) */
/*             WRK(N,3) = P(N) */
/*             WRK(N,4) = SX(N) */
/*             WRK(N,5) = WRK0(N) */
/*             WRK(N,6) = WRK1(N) */
/*             WRK(N,7) = WRK2(N) */
/*             WRK(N,8) = WRK3(N) */

    /* Parameter adjustments */
    wrk_dim1 = *nr;
    wrk_offset = wrk_dim1 + 1;
    wrk -= wrk_offset;
    a_dim1 = *nr;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --gpls;
    --xpls;
    --typsiz;
    --x;

    /* Function Body */
    optdrv_(nr, n, &x[1], (S_fp)fcn, (S_fp)d1fcn, (S_fp)d2fcn, &typsiz[1], 
	    fscale, method, iexp, msg, ndigit, itnlim, iagflg, iahflg, ipr, 
	    dlt, gradtl, stepmx, steptl, &xpls[1], fpls, &gpls[1], itrmcd, &a[
	    a_offset], &wrk[wrk_dim1 + 1], &wrk[(wrk_dim1 << 1) + 1], &wrk[
	    wrk_dim1 * 3 + 1], &wrk[(wrk_dim1 << 2) + 1], &wrk[wrk_dim1 * 5 + 
	    1], &wrk[wrk_dim1 * 6 + 1], &wrk[wrk_dim1 * 7 + 1], &wrk[(
	    wrk_dim1 << 3) + 1]);
    return 0;
} /* optif9_ */

/* Subroutine */ int optstp_(integer *n, doublereal *xpls, doublereal *fpls, 
	doublereal *gpls, doublereal *x, integer *itncnt, integer *icscmx, 
	integer *itrmcd, doublereal *gradtl, doublereal *steptl, doublereal *
	sx, doublereal *fscale, integer *itnlim, integer *iretcd, logical *
	mxtake, integer *ipr, integer *msg)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3, d__4;

    /* Local variables */
    static doublereal d;
    static integer i;
    static doublereal relgrd;
    static integer jtrmcd;
    static doublereal relstp, rgx, rsx;


/* UNCONSTRAINED MINIMIZATION STOPPING CRITERIA */
/* -------------------------------------------- */
/* FIND WHETHER THE ALGORITHM SHOULD TERMINATE, DUE TO ANY */
/* OF THE FOLLOWING: */
/* 1) PROBLEM SOLVED WITHIN USER TOLERANCE */
/* 2) CONVERGENCE WITHIN USER TOLERANCE */
/* 3) ITERATION LIMIT REACHED */
/* 4) DIVERGENCE OR TOO RESTRICTIVE MAXIMUM STEP (STEPMX) SUSPECTED */


/* PARAMETERS */
/* ---------- */
/* N            --> DIMENSION OF PROBLEM */
/* XPLS(N)      --> NEW ITERATE X[K] */
/* FPLS         --> FUNCTION VALUE AT NEW ITERATE F(XPLS) */
/* GPLS(N)      --> GRADIENT AT NEW ITERATE, G(XPLS), OR APPROXIMATE */
/* X(N)         --> OLD ITERATE X[K-1] */
/* ITNCNT       --> CURRENT ITERATION K */
/* ICSCMX      <--> NUMBER CONSECUTIVE STEPS .GE. STEPMX */
/*                  [RETAIN VALUE BETWEEN SUCCESSIVE CALLS] */
/* ITRMCD      <--  TERMINATION CODE */
/* GRADTL       --> TOLERANCE AT WHICH RELATIVE GRADIENT CONSIDERED CLOSE 
*/
/*                  ENOUGH TO ZERO TO TERMINATE ALGORITHM */
/* STEPTL       --> RELATIVE STEP SIZE AT WHICH SUCCESSIVE ITERATES */
/*                  CONSIDERED CLOSE ENOUGH TO TERMINATE ALGORITHM */
/* SX(N)        --> DIAGONAL SCALING MATRIX FOR X */
/* FSCALE       --> ESTIMATE OF SCALE OF OBJECTIVE FUNCTION */
/* ITNLIM       --> MAXIMUM NUMBER OF ALLOWABLE ITERATIONS */
/* IRETCD       --> RETURN CODE */
/* MXTAKE       --> BOOLEAN FLAG INDICATING STEP OF MAXIMUM LENGTH USED */
/* IPR          --> DEVICE TO WHICH TO SEND OUTPUT */
/* MSG          --> IF MSG INCLUDES A TERM 8, SUPPRESS OUTPUT */



    /* Parameter adjustments */
    --sx;
    --x;
    --gpls;
    --xpls;

    /* Function Body */
    *itrmcd = 0;

/* LAST GLOBAL STEP FAILED TO LOCATE A POINT LOWER THAN X */
    if (*iretcd != 1) {
	goto L50;
    }
/*     IF(IRETCD.EQ.1) */
/*     THEN */
    jtrmcd = 3;
    goto L600;
/*     ENDIF */
L50:

/* FIND DIRECTION IN WHICH RELATIVE GRADIENT MAXIMUM. */
/* CHECK WHETHER WITHIN TOLERANCE */

/* Computing MAX */
    d__1 = abs(*fpls);
    d = max(d__1,*fscale);
    rgx = 0.;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
/* Computing MAX */
	d__3 = (d__2 = xpls[i], abs(d__2)), d__4 = 1. / sx[i];
	relgrd = (d__1 = gpls[i], abs(d__1)) * max(d__3,d__4) / d;
	rgx = max(rgx,relgrd);
/* L100: */
    }
    jtrmcd = 1;
    if (rgx <= *gradtl) {
	goto L600;
    }

    if (*itncnt == 0) {
	return 0;
    }

/* FIND DIRECTION IN WHICH RELATIVE STEPSIZE MAXIMUM */
/* CHECK WHETHER WITHIN TOLERANCE. */

    rsx = 0.;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
/* Computing MAX */
	d__3 = (d__2 = xpls[i], abs(d__2)), d__4 = 1. / sx[i];
	relstp = (d__1 = xpls[i] - x[i], abs(d__1)) / max(d__3,d__4);
	rsx = max(rsx,relstp);
/* L120: */
    }
    jtrmcd = 2;
    if (rsx <= *steptl) {
	goto L600;
    }

/* CHECK ITERATION LIMIT */

    jtrmcd = 4;
    if (*itncnt >= *itnlim) {
	goto L600;
    }

/* CHECK NUMBER OF CONSECUTIVE STEPS \ STEPMX */

    if (*mxtake) {
	goto L140;
    }
/*     IF(.NOT.MXTAKE) */
/*     THEN */
    *icscmx = 0;
    return 0;
/*     ELSE */
L140:
/* %      IF (MOD(MSG/8,2) .EQ. 0) WRITE(IPR,900) */
    ++(*icscmx);
    if (*icscmx < 5) {
	return 0;
    }
    jtrmcd = 5;
/*     ENDIF */


/* PRINT TERMINATION CODE */

L600:
    *itrmcd = jtrmcd;
/* %    IF (MOD(MSG/8,2) .EQ. 0) GO TO(601,602,603,604,605), ITRMCD */
/* %    GO TO 700 */
/* %601 WRITE(IPR,901) */
/* %    GO TO 700 */
/* %602 WRITE(IPR,902) */
/* %    GO TO 700 */
/* %603 WRITE(IPR,903) */
/* %    GO TO 700 */
/* %604 WRITE(IPR,904) */
/* %    GO TO 700 */
/* %605 WRITE(IPR,905) */

/* L700: */
    return 0;

/* %900 FORMAT(48H0OPTSTP    STEP OF MAXIMUM LENGTH (STEPMX) TAKEN) */
/* %901 FORMAT(43H0OPTSTP    RELATIVE GRADIENT CLOSE TO ZERO./ */
/* %   +       48H OPTSTP    CURRENT ITERATE IS PROBABLY SOLUTION.) */
/* %902 FORMAT(48H0OPTSTP    SUCCESSIVE ITERATES WITHIN TOLERANCE./ */
/* %   +       48H OPTSTP    CURRENT ITERATE IS PROBABLY SOLUTION.) */
/* %903 FORMAT(52H0OPTSTP    LAST GLOBAL STEP FAILED TO LOCATE A POINT, */
/* %   +       14H LOWER THAN X./ */
/* %   +       51H OPTSTP    EITHER X IS AN APPROXIMATE LOCAL MINIMUM, */
/* %   +       17H OF THE FUNCTION,/ */
/* %   +       50H OPTSTP    THE FUNCTION IS TOO NON-LINEAR FOR THIS, */
/* %   +       11H ALGORITHM,/ */
/* %   +       34H OPTSTP    OR STEPTL IS TOO LARGE.) */
/* %904 FORMAT(36H0OPTSTP    ITERATION LIMIT EXCEEDED./ */
/* %   +       28H OPTSTP    ALGORITHM FAILED.) */
/* %905 FORMAT(39H0OPTSTP    MAXIMUM STEP SIZE EXCEEDED 5, */
/* %   +       19H CONSECUTIVE TIMES./ */
/* %   +       50H OPTSTP    EITHER THE FUNCTION IS UNBOUNDED BELOW,/ */
/* %   +       47H OPTSTP    BECOMES ASYMPTOTIC TO A FINITE VALUE, */
/* %   +       30H FROM ABOVE IN SOME DIRECTION,/ */
/* %   +       33H OPTSTP    OR STEPMX IS TOO SMALL) */
} /* optstp_ */

/* Subroutine */ int qraux1_(integer *nr, integer *n, doublereal *r, integer *
	i)
{
    /* System generated locals */
    integer r_dim1, r_offset, i__1;

    /* Local variables */
    static integer j;
    static doublereal tmp;


/* PURPOSE */
/* ------- */
/* INTERCHANGE ROWS I,I+1 OF THE UPPER HESSENBERG MATRIX R, */
/* COLUMNS I TO N */

/* PARAMETERS */
/* ---------- */
/* NR           --> ROW DIMENSION OF MATRIX */
/* N            --> DIMENSION OF MATRIX */
/* R(N,N)      <--> UPPER HESSENBERG MATRIX */
/* I            --> INDEX OF ROW TO INTERCHANGE (I.LT.N) */

    /* Parameter adjustments */
    r_dim1 = *nr;
    r_offset = r_dim1 + 1;
    r -= r_offset;

    /* Function Body */
    i__1 = *n;
    for (j = *i; j <= i__1; ++j) {
	tmp = r[*i + j * r_dim1];
	r[*i + j * r_dim1] = r[*i + 1 + j * r_dim1];
	r[*i + 1 + j * r_dim1] = tmp;
/* L10: */
    }
    return 0;
} /* qraux1_ */

/* Subroutine */ int qraux2_(integer *nr, integer *n, doublereal *r, integer *
	i, doublereal *a, doublereal *b)
{
    /* System generated locals */
    integer r_dim1, r_offset, i__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal c;
    static integer j;
    static doublereal s, y, z, den;


/* PURPOSE */
/* ------- */
/* PRE-MULTIPLY R BY THE JACOBI ROTATION J(I,I+1,A,B) */

/* PARAMETERS */
/* ---------- */
/* NR           --> ROW DIMENSION OF MATRIX */
/* N            --> DIMENSION OF MATRIX */
/* R(N,N)      <--> UPPER HESSENBERG MATRIX */
/* I            --> INDEX OF ROW */
/* A            --> SCALAR */
/* B            --> SCALAR */

    /* Parameter adjustments */
    r_dim1 = *nr;
    r_offset = r_dim1 + 1;
    r -= r_offset;

    /* Function Body */
    den = sqrt(*a * *a + *b * *b);
    c = *a / den;
    s = *b / den;
    i__1 = *n;
    for (j = *i; j <= i__1; ++j) {
	y = r[*i + j * r_dim1];
	z = r[*i + 1 + j * r_dim1];
	r[*i + j * r_dim1] = c * y - s * z;
	r[*i + 1 + j * r_dim1] = s * y + c * z;
/* L10: */
    }
    return 0;
} /* qraux2_ */

/* Subroutine */ int qrupdt_(integer *nr, integer *n, doublereal *a, 
	doublereal *u, doublereal *v)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i, j, k;
    static doublereal t1, t2;
    extern /* Subroutine */ int qraux1_(integer *, integer *, doublereal *, 
	    integer *);
    static integer ii;
    extern /* Subroutine */ int qraux2_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *);
    static integer km1;


/* PURPOSE */
/* ------- */
/* FIND AN ORTHOGONAL (N*N) MATRIX (Q*) AND AN UPPER TRIANGULAR (N*N) */
/* MATRIX (R*) SUCH THAT (Q*)(R*)=R+U(V+) */

/* PARAMETERS */
/* ---------- */
/* NR           --> ROW DIMENSION OF MATRIX */
/* N            --> DIMENSION OF PROBLEM */
/* A(N,N)      <--> ON INPUT:  CONTAINS R */
/*                  ON OUTPUT: CONTAINS (R*) */
/* U(N)         --> VECTOR */
/* V(N)         --> VECTOR */


/* DETERMINE LAST NON-ZERO IN U(.) */

    /* Parameter adjustments */
    --v;
    --u;
    a_dim1 = *nr;
    a_offset = a_dim1 + 1;
    a -= a_offset;

    /* Function Body */
    k = *n;
L10:
    if (u[k] != 0. || k == 1) {
	goto L20;
    }
/*     IF(U(K).EQ.0.0D0 .AND. K.GT.1) */
/*     THEN */
    --k;
    goto L10;
/*     ENDIF */

/* (K-1) JACOBI ROTATIONS TRANSFORM */
/*     R + U(V+) --> (R*) + (U(1)*E1)(V+) */
/* WHICH IS UPPER HESSENBERG */

L20:
    if (k <= 1) {
	goto L40;
    }
    km1 = k - 1;
    i__1 = km1;
    for (ii = 1; ii <= i__1; ++ii) {
	i = km1 - ii + 1;
	if (u[i] != 0.) {
	    goto L25;
	}
/*         IF(U(I).EQ.0.0D0) */
/*         THEN */
	qraux1_(nr, n, &a[a_offset], &i);
	u[i] = u[i + 1];
	goto L30;
/*         ELSE */
L25:
	d__1 = -u[i + 1];
	qraux2_(nr, n, &a[a_offset], &i, &u[i], &d__1);
	u[i] = sqrt(u[i] * u[i] + u[i + 1] * u[i + 1]);
/*         ENDIF */
L30:
	;
    }
/*     ENDIF */

/* R <-- R + (U(1)*E1)(V+) */

L40:
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	a[j * a_dim1 + 1] += u[1] * v[j];
/* L50: */
    }

/* (K-1) JACOBI ROTATIONS TRANSFORM UPPER HESSENBERG R */
/* TO UPPER TRIANGULAR (R*) */

    if (k <= 1) {
	goto L100;
    }
    km1 = k - 1;
    i__1 = km1;
    for (i = 1; i <= i__1; ++i) {
	if (a[i + i * a_dim1] != 0.) {
	    goto L70;
	}
/*         IF(A(I,I).EQ.0.0D0) */
/*         THEN */
	qraux1_(nr, n, &a[a_offset], &i);
	goto L80;
/*         ELSE */
L70:
	t1 = a[i + i * a_dim1];
	t2 = -a[i + 1 + i * a_dim1];
	qraux2_(nr, n, &a[a_offset], &i, &t1, &t2);
/*         ENDIF */
L80:
	;
    }
/*     ENDIF */
L100:
    return 0;
} /* qrupdt_ */

/* Subroutine */ int sclmul_(integer *n, doublereal *s, doublereal *v, 
	doublereal *z)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i;


/* PURPOSE */
/* ------- */
/* MULTIPLY VECTOR BY SCALAR */
/* RESULT VECTOR MAY BE OPERAND VECTOR */

/* PARAMETERS */
/* ---------- */
/* N            --> DIMENSION OF VECTORS */
/* S            --> SCALAR */
/* V(N)         --> OPERAND VECTOR */
/* Z(N)        <--  RESULT VECTOR */
    /* Parameter adjustments */
    --z;
    --v;

    /* Function Body */
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	z[i] = *s * v[i];
/* L100: */
    }
    return 0;
} /* sclmul_ */

/* Subroutine */ int secfac_(integer *nr, integer *n, doublereal *x, 
	doublereal *g, doublereal *a, doublereal *xpls, doublereal *gpls, 
	doublereal *epsm, integer *itncnt, doublereal *rnf, integer *iagflg, 
	logical *noupdt, doublereal *s, doublereal *y, doublereal *u, 
	doublereal *w)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2, d__3, d__4, d__5;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *), dnrm2_(integer *, doublereal *, integer *);
    static doublereal ynrm2;
    static integer i, j;
    static doublereal snorm2, reltol;
    static logical skpupd;
    static integer im1;
    extern /* Subroutine */ int mvmltl_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *), qrupdt_(integer *, integer *, 
	    doublereal *, doublereal *, doublereal *), mvmltu_(integer *, 
	    integer *, doublereal *, doublereal *, doublereal *);
    static doublereal alp, den1, den2;


/* PURPOSE */
/* ------- */
/* UPDATE HESSIAN BY THE BFGS FACTORED METHOD */

/* PARAMETERS */
/* ---------- */
/* NR           --> ROW DIMENSION OF MATRIX */
/* N            --> DIMENSION OF PROBLEM */
/* X(N)         --> OLD ITERATE, X[K-1] */
/* G(N)         --> GRADIENT OR APPROXIMATE AT OLD ITERATE */
/* A(N,N)      <--> ON ENTRY: CHOLESKY DECOMPOSITION OF HESSIAN IN */
/*                    LOWER PART AND DIAGONAL. */
/*                  ON EXIT:  UPDATED CHOLESKY DECOMPOSITION OF HESSIAN */
/*                    IN LOWER TRIANGULAR PART AND DIAGONAL */
/* XPLS(N)      --> NEW ITERATE, X[K] */
/* GPLS(N)      --> GRADIENT OR APPROXIMATE AT NEW ITERATE */
/* EPSM         --> MACHINE EPSILON */
/* ITNCNT       --> ITERATION COUNT */
/* RNF          --> RELATIVE NOISE IN OPTIMIZATION FUNCTION FCN */
/* IAGFLG       --> =1 IF ANALYTIC GRADIENT SUPPLIED, =0 ITHERWISE */
/* NOUPDT      <--> BOOLEAN: NO UPDATE YET */
/*                  [RETAIN VALUE BETWEEN SUCCESSIVE CALLS] */
/* S(N)         --> WORKSPACE */
/* Y(N)         --> WORKSPACE */
/* U(N)         --> WORKSPACE */
/* W(N)         --> WORKSPACE */


    /* Parameter adjustments */
    --w;
    --u;
    --y;
    --s;
    --gpls;
    --xpls;
    a_dim1 = *nr;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --g;
    --x;

    /* Function Body */
    if (*itncnt == 1) {
	*noupdt = TRUE_;
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	s[i] = xpls[i] - x[i];
	y[i] = gpls[i] - g[i];
/* L10: */
    }
    den1 = ddot_(n, &s[1], &c__1, &y[1], &c__1);
    snorm2 = dnrm2_(n, &s[1], &c__1);
    ynrm2 = dnrm2_(n, &y[1], &c__1);
    if (den1 < sqrt(*epsm) * snorm2 * ynrm2) {
	goto L110;
    }
/*     IF(DEN1.GE.SQRT(EPSM)*SNORM2*YNRM2) */
/*     THEN */
    mvmltu_(nr, n, &a[a_offset], &s[1], &u[1]);
    den2 = ddot_(n, &u[1], &c__1, &u[1], &c__1);

/*       L <-- SQRT(DEN1/DEN2)*L */

    alp = sqrt(den1 / den2);
    if (! (*noupdt)) {
	goto L50;
    }
/*       IF(NOUPDT) */
/*       THEN */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	u[j] = alp * u[j];
	i__2 = *n;
	for (i = j; i <= i__2; ++i) {
	    a[i + j * a_dim1] = alp * a[i + j * a_dim1];
/* L20: */
	}
/* L30: */
    }
    *noupdt = FALSE_;
    den2 = den1;
    alp = 1.;
/*       ENDIF */
L50:
    skpupd = TRUE_;

/*       W = L(L+)S = HS */

    mvmltl_(nr, n, &a[a_offset], &u[1], &w[1]);
    i = 1;
    if (*iagflg != 0) {
	goto L55;
    }
/*       IF(IAGFLG.EQ.0) */
/*       THEN */
    reltol = sqrt(*rnf);
    goto L60;
/*       ELSE */
L55:
    reltol = *rnf;
/*       ENDIF */
L60:
    if (i > *n || ! skpupd) {
	goto L70;
    }
/*       IF(I.LE.N .AND. SKPUPD) */
/*       THEN */
/* Computing MAX */
    d__4 = (d__2 = g[i], abs(d__2)), d__5 = (d__3 = gpls[i], abs(d__3));
    if ((d__1 = y[i] - w[i], abs(d__1)) < reltol * max(d__4,d__5)) {
	goto L65;
    }
/*         IF(ABS(Y(I)-W(I)) .GE. RELTOL*DMAX1(ABS(G(I)),ABS(GPLS(I)))) */
/*         THEN */
    skpupd = FALSE_;
    goto L60;
/*         ELSE */
L65:
    ++i;
    goto L60;
/*         ENDIF */
/*       ENDIF */
L70:
    if (skpupd) {
	goto L110;
    }
/*       IF(.NOT.SKPUPD) */
/*       THEN */

/*         W=Y-ALP*L(L+)S */

    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	w[i] = y[i] - alp * w[i];
/* L75: */
    }

/*         ALP=1/SQRT(DEN1*DEN2) */

    alp /= den1;

/*         U=(L+)/SQRT(DEN1*DEN2) = (L+)S/SQRT((Y+)S * (S+)L(L+)S) */

    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	u[i] = alp * u[i];
/* L80: */
    }

/*         COPY L INTO UPPER TRIANGULAR PART.  ZERO L. */

    if (*n == 1) {
	goto L93;
    }
    i__1 = *n;
    for (i = 2; i <= i__1; ++i) {
	im1 = i - 1;
	i__2 = im1;
	for (j = 1; j <= i__2; ++j) {
	    a[j + i * a_dim1] = a[i + j * a_dim1];
	    a[i + j * a_dim1] = 0.;
/* L85: */
	}
/* L90: */
    }

/*         FIND Q, (L+) SUCH THAT  Q(L+) = (L+) + U(W+) */

L93:
    qrupdt_(nr, n, &a[a_offset], &u[1], &w[1]);

/*         UPPER TRIANGULAR PART AND DIAGONAL OF A NOW CONTAIN UPDATED */
/*         CHOLESKY DECOMPOSITION OF HESSIAN.  COPY BACK TO LOWER */
/*         TRIANGULAR PART. */

    if (*n == 1) {
	goto L110;
    }
    i__1 = *n;
    for (i = 2; i <= i__1; ++i) {
	im1 = i - 1;
	i__2 = im1;
	for (j = 1; j <= i__2; ++j) {
	    a[i + j * a_dim1] = a[j + i * a_dim1];
/* L95: */
	}
/* L100: */
    }
/*       ENDIF */
/*     ENDIF */
L110:
    return 0;
} /* secfac_ */

/* Subroutine */ int secunf_(integer *nr, integer *n, doublereal *x, 
	doublereal *g, doublereal *a, doublereal *udiag, doublereal *xpls, 
	doublereal *gpls, doublereal *epsm, integer *itncnt, doublereal *rnf, 
	integer *iagflg, logical *noupdt, doublereal *s, doublereal *y, 
	doublereal *t)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *), dnrm2_(integer *, doublereal *, integer *);
    static doublereal ynrm2;
    static integer i, j;
    static doublereal snorm2;
    static logical skpupd;
    static integer jp1;
    extern /* Subroutine */ int mvmlts_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal gam, tol, den1, den2;


/* PURPOSE */
/* ------- */
/* UPDATE HESSIAN BY THE BFGS UNFACTORED METHOD */

/* PARAMETERS */
/* ---------- */
/* NR           --> ROW DIMENSION OF MATRIX */
/* N            --> DIMENSION OF PROBLEM */
/* X(N)         --> OLD ITERATE, X[K-1] */
/* G(N)         --> GRADIENT OR APPROXIMATE AT OLD ITERATE */
/* A(N,N)      <--> ON ENTRY: APPROXIMATE HESSIAN AT OLD ITERATE */
/*                    IN UPPER TRIANGULAR PART (AND UDIAG) */
/*                  ON EXIT:  UPDATED APPROX HESSIAN AT NEW ITERATE */
/*                    IN LOWER TRIANGULAR PART AND DIAGONAL */
/*                  [LOWER TRIANGULAR PART OF SYMMETRIC MATRIX] */
/* UDIAG        --> ON ENTRY: DIAGONAL OF HESSIAN */
/* XPLS(N)      --> NEW ITERATE, X[K] */
/* GPLS(N)      --> GRADIENT OR APPROXIMATE AT NEW ITERATE */
/* EPSM         --> MACHINE EPSILON */
/* ITNCNT       --> ITERATION COUNT */
/* RNF          --> RELATIVE NOISE IN OPTIMIZATION FUNCTION FCN */
/* IAGFLG       --> =1 IF ANALYTIC GRADIENT SUPPLIED, =0 OTHERWISE */
/* NOUPDT      <--> BOOLEAN: NO UPDATE YET */
/*                  [RETAIN VALUE BETWEEN SUCCESSIVE CALLS] */
/* S(N)         --> WORKSPACE */
/* Y(N)         --> WORKSPACE */
/* T(N)         --> WORKSPACE */


/* COPY HESSIAN IN UPPER TRIANGULAR PART AND UDIAG TO */
/* LOWER TRIANGULAR PART AND DIAGONAL */

    /* Parameter adjustments */
    --t;
    --y;
    --s;
    --gpls;
    --xpls;
    --udiag;
    a_dim1 = *nr;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --g;
    --x;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	a[j + j * a_dim1] = udiag[j];
	if (j == *n) {
	    goto L5;
	}
	jp1 = j + 1;
	i__2 = *n;
	for (i = jp1; i <= i__2; ++i) {
	    a[i + j * a_dim1] = a[j + i * a_dim1];
/* L4: */
	}
L5:
	;
    }

    if (*itncnt == 1) {
	*noupdt = TRUE_;
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	s[i] = xpls[i] - x[i];
	y[i] = gpls[i] - g[i];
/* L10: */
    }
    den1 = ddot_(n, &s[1], &c__1, &y[1], &c__1);
    snorm2 = dnrm2_(n, &s[1], &c__1);
    ynrm2 = dnrm2_(n, &y[1], &c__1);
    if (den1 < sqrt(*epsm) * snorm2 * ynrm2) {
	goto L100;
    }
/*     IF(DEN1.GE.SQRT(EPSM)*SNORM2*YNRM2) */
/*     THEN */
    mvmlts_(nr, n, &a[a_offset], &s[1], &t[1]);
    den2 = ddot_(n, &s[1], &c__1, &t[1], &c__1);
    if (! (*noupdt)) {
	goto L50;
    }
/*       IF(NOUPDT) */
/*       THEN */

/*         H <-- [(S+)Y/(S+)HS]H */

    gam = den1 / den2;
    den2 = gam * den2;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	t[j] = gam * t[j];
	i__2 = *n;
	for (i = j; i <= i__2; ++i) {
	    a[i + j * a_dim1] = gam * a[i + j * a_dim1];
/* L20: */
	}
/* L30: */
    }
    *noupdt = FALSE_;
/*       ENDIF */
L50:
    skpupd = TRUE_;

/*       CHECK UPDATE CONDITION ON ROW I */

    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
/* Computing MAX */
	d__3 = (d__1 = g[i], abs(d__1)), d__4 = (d__2 = gpls[i], abs(d__2));
	tol = *rnf * max(d__3,d__4);
	if (*iagflg == 0) {
	    tol /= sqrt(*rnf);
	}
	if ((d__1 = y[i] - t[i], abs(d__1)) < tol) {
	    goto L60;
	}
/*         IF(ABS(Y(I)-T(I)).GE.TOL) */
/*         THEN */
	skpupd = FALSE_;
	goto L70;
/*         ENDIF */
L60:
	;
    }
L70:
    if (skpupd) {
	goto L100;
    }
/*       IF(.NOT.SKPUPD) */
/*       THEN */

/*         BFGS UPDATE */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *n;
	for (i = j; i <= i__2; ++i) {
	    a[i + j * a_dim1] = a[i + j * a_dim1] + y[i] * y[j] / den1 - t[i] 
		    * t[j] / den2;
/* L80: */
	}
/* L90: */
    }
/*       ENDIF */
/*     ENDIF */
L100:
    return 0;
} /* secunf_ */

/* Subroutine */ int sndofd_(integer *nr, integer *n, doublereal *xpls, S_fp 
	fcn, doublereal *fpls, doublereal *a, doublereal *sx, doublereal *
	rnoise, doublereal *stepsz, doublereal *anbr)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static doublereal fhat;
    static integer i, j;
    static doublereal xtmpi, xtmpj;
    static integer ip1;
    static doublereal ov3;

/* PURPOSE */
/* ------- */
/* FIND SECOND ORDER FORWARD FINITE DIFFERENCE APPROXIMATION "A" */
/* TO THE SECOND DERIVATIVE (HESSIAN) OF THE FUNCTION DEFINED BY THE SUBP 
*/
/* "FCN" EVALUATED AT THE NEW ITERATE "XPLS" */

/* FOR OPTIMIZATION USE THIS ROUTINE TO ESTIMATE */
/* 1) THE SECOND DERIVATIVE (HESSIAN) OF THE OPTIMIZATION FUNCTION */
/*    IF NO ANALYTICAL USER FUNCTION HAS BEEN SUPPLIED FOR EITHER */
/*    THE GRADIENT OR THE HESSIAN AND IF THE OPTIMIZATION FUNCTION */
/*    "FCN" IS INEXPENSIVE TO EVALUATE. */

/* PARAMETERS */
/* ---------- */
/* NR           --> ROW DIMENSION OF MATRIX */
/* N            --> DIMENSION OF PROBLEM */
/* XPLS(N)      --> NEW ITERATE:   X[K] */
/* FCN          --> NAME OF SUBROUTINE TO EVALUATE FUNCTION */
/* FPLS         --> FUNCTION VALUE AT NEW ITERATE, F(XPLS) */
/* A(N,N)      <--  FINITE DIFFERENCE APPROXIMATION TO HESSIAN */
/*                  ONLY LOWER TRIANGULAR MATRIX AND DIAGONAL */
/*                  ARE RETURNED */
/* SX(N)        --> DIAGONAL SCALING MATRIX FOR X */
/* RNOISE       --> RELATIVE NOISE IN FNAME [F(X)] */
/* STEPSZ(N)    --> WORKSPACE (STEPSIZE IN I-TH COMPONENT DIRECTION) */
/* ANBR(N)      --> WORKSPACE (NEIGHBOR IN I-TH DIRECTION) */



/* FIND I-TH STEPSIZE AND EVALUATE NEIGHBOR IN DIRECTION */
/* OF I-TH UNIT VECTOR. */

    /* Parameter adjustments */
    --anbr;
    --stepsz;
    --sx;
    a_dim1 = *nr;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --xpls;

    /* Function Body */
    ov3 = .33333333333333331;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
/* Computing MAX */
	d__2 = (d__1 = xpls[i], abs(d__1)), d__3 = 1. / sx[i];
	stepsz[i] = pow_dd(rnoise, &ov3) * max(d__2,d__3);
	xtmpi = xpls[i];
	xpls[i] = xtmpi + stepsz[i];
	(*fcn)(n, &xpls[1], &anbr[i]);
	xpls[i] = xtmpi;
/* L10: */
    }

/* CALCULATE COLUMN I OF A */

    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	xtmpi = xpls[i];
	xpls[i] = xtmpi + stepsz[i] * 2.;
	(*fcn)(n, &xpls[1], &fhat);
	a[i + i * a_dim1] = (*fpls - anbr[i] + (fhat - anbr[i])) / (stepsz[i] 
		* stepsz[i]);

/* CALCULATE SUB-DIAGONAL ELEMENTS OF COLUMN */
	if (i == *n) {
	    goto L25;
	}
	xpls[i] = xtmpi + stepsz[i];
	ip1 = i + 1;
	i__2 = *n;
	for (j = ip1; j <= i__2; ++j) {
	    xtmpj = xpls[j];
	    xpls[j] = xtmpj + stepsz[j];
	    (*fcn)(n, &xpls[1], &fhat);
	    a[j + i * a_dim1] = (*fpls - anbr[i] + (fhat - anbr[j])) / (
		    stepsz[i] * stepsz[j]);
	    xpls[j] = xtmpj;
/* L20: */
	}
L25:
	xpls[i] = xtmpi;
/* L30: */
    }
    return 0;
} /* sndofd_ */

/* Subroutine */ int tregup_(integer *nr, integer *n, doublereal *x, 
	doublereal *f, doublereal *g, doublereal *a, S_fp fcn, doublereal *sc,
	 doublereal *sx, logical *nwtake, doublereal *stepmx, doublereal *
	steptl, doublereal *dlt, integer *iretcd, doublereal *xplsp, 
	doublereal *fplsp, doublereal *xpls, doublereal *fpls, logical *
	mxtake, integer *ipr, integer *method, doublereal *udiag)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;

    /* Local variables */
    static doublereal dltf;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal temp;
    static integer i, j;
    static doublereal dltfp, dltmp;
    static integer ip1;
    static doublereal rln, slp;


/* PURPOSE */
/* ------- */
/* DECIDE WHETHER TO ACCEPT XPLS=X+SC AS THE NEXT ITERATE AND UPDATE THE 
*/
/* TRUST REGION DLT. */

/* PARAMETERS */
/* ---------- */
/* NR           --> ROW DIMENSION OF MATRIX */
/* N            --> DIMENSION OF PROBLEM */
/* X(N)         --> OLD ITERATE X[K-1] */
/* F            --> FUNCTION VALUE AT OLD ITERATE, F(X) */
/* G(N)         --> GRADIENT AT OLD ITERATE, G(X), OR APPROXIMATE */
/* A(N,N)       --> CHOLESKY DECOMPOSITION OF HESSIAN IN */
/*                  LOWER TRIANGULAR PART AND DIAGONAL. */
/*                  HESSIAN OR APPROX IN UPPER TRIANGULAR PART */
/* FCN          --> NAME OF SUBROUTINE TO EVALUATE FUNCTION */
/* SC(N)        --> CURRENT STEP */
/* SX(N)        --> DIAGONAL SCALING MATRIX FOR X */
/* NWTAKE       --> BOOLEAN, =.TRUE. IF NEWTON STEP TAKEN */
/* STEPMX       --> MAXIMUM ALLOWABLE STEP SIZE */
/* STEPTL       --> RELATIVE STEP SIZE AT WHICH SUCCESSIVE ITERATES */
/*                  CONSIDERED CLOSE ENOUGH TO TERMINATE ALGORITHM */
/* DLT         <--> TRUST REGION RADIUS */
/* IRETCD      <--> RETURN CODE */
/*                    =0 XPLS ACCEPTED AS NEXT ITERATE; */
/*                       DLT TRUST REGION FOR NEXT ITERATION. */
/*                    =1 XPLS UNSATISFACTORY BUT ACCEPTED AS NEXT ITERATE 
*/
/*                       BECAUSE XPLS-X .LT. SMALLEST ALLOWABLE */
/*                       STEP LENGTH. */
/*                    =2 F(XPLS) TOO LARGE.  CONTINUE CURRENT ITERATION */
/*                       WITH NEW REDUCED DLT. */
/*                    =3 F(XPLS) SUFFICIENTLY SMALL, BUT QUADRATIC MODEL 
*/
/*                       PREDICTS F(XPLS) SUFFICIENTLY WELL TO CONTINUE */
/*                       CURRENT ITERATION WITH NEW DOUBLED DLT. */
/* XPLSP(N)    <--> WORKSPACE [VALUE NEEDS TO BE RETAINED BETWEEN */
/*                  SUCCESIVE CALLS OF K-TH GLOBAL STEP] */
/* FPLSP       <--> [RETAIN VALUE BETWEEN SUCCESSIVE CALLS] */
/* XPLS(N)     <--  NEW ITERATE X[K] */
/* FPLS        <--  FUNCTION VALUE AT NEW ITERATE, F(XPLS) */
/* MXTAKE      <--  BOOLEAN FLAG INDICATING STEP OF MAXIMUM LENGTH USED */
/* IPR          --> DEVICE TO WHICH TO SEND OUTPUT */
/* METHOD       --> ALGORITHM TO USE TO SOLVE MINIMIZATION PROBLEM */
/*                    =1 LINE SEARCH */
/*                    =2 DOUBLE DOGLEG */
/*                    =3 MORE-HEBDON */
/* UDIAG(N)     --> DIAGONAL OF HESSIAN IN A(.,.) */


    /* Parameter adjustments */
    --udiag;
    --xpls;
    --xplsp;
    --sx;
    --sc;
    a_dim1 = *nr;
    a_offset = a_dim1 + 1;
    a -= a_offset;
    --g;
    --x;

    /* Function Body */
    *ipr = *ipr;
    *mxtake = FALSE_;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	xpls[i] = x[i] + sc[i];
/* L100: */
    }
    (*fcn)(n, &xpls[1], fpls);
    dltf = *fpls - *f;
    slp = ddot_(n, &g[1], &c__1, &sc[1], &c__1);

/* NEXT STATEMENT ADDED FOR CASE OF COMPILERS WHICH DO NOT OPTIMIZE */
/* EVALUATION OF NEXT "IF" STATEMENT (IN WHICH CASE FPLSP COULD BE */
/* UNDEFINED). */
    if (*iretcd == 4) {
	*fplsp = 0.;
    }
/* $    WRITE(IPR,961) IRETCD,FPLS,FPLSP,DLTF,SLP */
    if (*iretcd != 3 || *fpls < *fplsp && dltf <= slp * 1e-4) {
	goto L130;
    }
/*     IF(IRETCD.EQ.3 .AND. (FPLS.GE.FPLSP .OR. DLTF.GT. 1.D-4*SLP)) */
/*     THEN */

/*       RESET XPLS TO XPLSP AND TERMINATE GLOBAL STEP */

    *iretcd = 0;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	xpls[i] = xplsp[i];
/* L110: */
    }
    *fpls = *fplsp;
    *dlt *= .5f;
/* $      WRITE(IPR,951) */
    goto L230;
/*     ELSE */

/*       FPLS TOO LARGE */

L130:
    if (dltf <= slp * 1e-4) {
	goto L170;
    }
/*       IF(DLTF.GT. 1.D-4*SLP) */
/*       THEN */
/* $        WRITE(IPR,952) */
    rln = 0.;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
/* Computing MAX */
/* Computing MAX */
	d__5 = (d__2 = xpls[i], abs(d__2)), d__6 = 1. / sx[i];
	d__3 = rln, d__4 = (d__1 = sc[i], abs(d__1)) / max(d__5,d__6);
	rln = max(d__3,d__4);
/* L140: */
    }
/* $        WRITE(IPR,962) RLN */
    if (rln >= *steptl) {
	goto L150;
    }
/*         IF(RLN.LT.STEPTL) */
/*         THEN */

/*           CANNOT FIND SATISFACTORY XPLS SUFFICIENTLY DISTINCT FROM X */

    *iretcd = 1;
/* $          WRITE(IPR,954) */
    goto L230;
/*         ELSE */

/*           REDUCE TRUST REGION AND CONTINUE GLOBAL STEP */

L150:
    *iretcd = 2;
    dltmp = -slp * *dlt / ((dltf - slp) * 2.);
/* $          WRITE(IPR,963) DLTMP */
    if (dltmp >= *dlt * .1f) {
	goto L155;
    }
/*           IF(DLTMP.LT. .1*DLT) */
/*           THEN */
    *dlt *= .1f;
    goto L160;
/*           ELSE */
L155:
    *dlt = dltmp;
/*           ENDIF */
L160:
/* $          WRITE(IPR,955) */
    goto L230;
/*         ENDIF */
/*       ELSE */

/*         FPLS SUFFICIENTLY SMALL */

L170:
/* $        WRITE(IPR,958) */
    dltfp = 0.;
    if (*method == 3) {
	goto L180;
    }

/*         IF (METHOD .EQ. 2) */
/*         THEN */

    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	temp = 0.;
	i__2 = *n;
	for (j = i; j <= i__2; ++j) {
	    temp += a[j + i * a_dim1] * sc[j];
/* L173: */
	}
	dltfp += temp * temp;
/* L177: */
    }
    goto L190;

/*         ELSE */

L180:
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	dltfp += udiag[i] * sc[i] * sc[i];
	if (i == *n) {
	    goto L187;
	}
	temp = 0.;
	ip1 = i + 1;
	i__2 = *n;
	for (j = ip1; j <= i__2; ++j) {
	    temp += a[i + j * a_dim1] * sc[i] * sc[j];
/* L183: */
	}
	dltfp += temp * 2.;
L187:
	;
    }

/*         END IF */

L190:
    dltfp = slp + dltfp / 2.;
/* $        WRITE(IPR,964) DLTFP,NWTAKE */
    if (*iretcd == 2 || (d__1 = dltfp - dltf, abs(d__1)) > abs(dltf) * .1f || 
	    *nwtake || *dlt > *stepmx * .99f) {
	goto L210;
    }
/*         IF(IRETCD.NE.2 .AND. (ABS(DLTFP-DLTF) .LE. .1*ABS(DLTF)) */
/*    +         .AND. (.NOT.NWTAKE) .AND. (DLT.LE. .99*STEPMX)) */
/*         THEN */

/*           DOUBLE TRUST REGION AND CONTINUE GLOBAL STEP */

    *iretcd = 3;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	xplsp[i] = xpls[i];
/* L200: */
    }
    *fplsp = *fpls;
/* Computing MIN */
    d__1 = *dlt * 2.;
    *dlt = min(d__1,*stepmx);
/* $          WRITE(IPR,959) */
    goto L230;
/*         ELSE */

/*           ACCEPT XPLS AS NEXT ITERATE.  CHOOSE NEW TRUST REGION. */

L210:
/* $          WRITE(IPR,960) */
    *iretcd = 0;
    if (*dlt > *stepmx * .99f) {
	*mxtake = TRUE_;
    }
    if (dltf < dltfp * .1f) {
	goto L220;
    }
/*           IF(DLTF.GE. .1*DLTFP) */
/*           THEN */

/*             DECREASE TRUST REGION FOR NEXT ITERATION */

    *dlt *= .5f;
    goto L230;
/*           ELSE */

/*             CHECK WHETHER TO INCREASE TRUST REGION FOR NEXT ITERATION 
*/

L220:
    if (dltf <= dltfp * .75f) {
/* Computing MIN */
	d__1 = *dlt * 2.;
	*dlt = min(d__1,*stepmx);
    }
/*           ENDIF */
/*         ENDIF */
/*       ENDIF */
/*     ENDIF */
L230:
/* $    WRITE(IPR,953) */
/* $    WRITE(IPR,956) IRETCD,MXTAKE,DLT,FPLS */
/* $    WRITE(IPR,957) */
/* $    WRITE(IPR,965) (XPLS(I),I=1,N) */
    return 0;

/* %951 FORMAT(55H TREGUP    RESET XPLS TO XPLSP. TERMINATION GLOBAL STEP)
 */
/* %952 FORMAT(26H TREGUP    FPLS TOO LARGE.) */
/* %953 FORMAT(38H0TREGUP    VALUES AFTER CALL TO TREGUP) */
/* %954 FORMAT(54H TREGUP    CANNOT FIND SATISFACTORY XPLS DISTINCT FROM, 
*/
/* %   +       27H X.  TERMINATE GLOBAL STEP.) */
/* %955 FORMAT(53H TREGUP    REDUCE TRUST REGION. CONTINUE GLOBAL STEP.) 
*/
/* %956 FORMAT(21H TREGUP       IRETCD=,I3/ */
/* %   +       21H TREGUP       MXTAKE=,L1/ */
/* %   +       21H TREGUP       DLT   =,E20.13/ */
/* %   +       21H TREGUP       FPLS  =,E20.13) */
/* %957 FORMAT(32H TREGUP       NEW ITERATE (XPLS)) */
/* %958 FORMAT(35H TREGUP    FPLS SUFFICIENTLY SMALL.) */
/* %959 FORMAT(54H TREGUP    DOUBLE TRUST REGION.  CONTINUE GLOBAL STEP.) 
*/
/* %960 FORMAT(50H TREGUP    ACCEPT XPLS AS NEW ITERATE.  CHOOSE NEW, */
/* %   +       38H TRUST REGION.  TERMINATE GLOBAL STEP.) */
/* %961 FORMAT(18H TREGUP    IRETCD=,I5/ */
/* %   +       18H TREGUP    FPLS  =,E20.13/ */
/* %   +       18H TREGUP    FPLSP =,E20.13/ */
/* %   +       18H TREGUP    DLTF  =,E20.13/ */
/* %   +       18H TREGUP    SLP   =,E20.13) */
/* %962 FORMAT(18H TREGUP    RLN   =,E20.13) */
/* %963 FORMAT(18H TREGUP    DLTMP =,E20.13) */
/* %964 FORMAT(18H TREGUP    DLTFP =,E20.13/ */
/* %   +       18H TREGUP    NWTAKE=,L1) */
/* %965 FORMAT(14H TREGUP       ,5(E20.13,3X)) */
} /* tregup_ */

