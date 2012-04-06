/*
 *  R : A Computer Langage for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include "Defn.h"
#include "Mathlib.h"
#include "Print.h"

void formatString(SEXP *x, int n, int *fieldwidth, int quote)
{
	int xmax = 0, naflag = 0;
	int i, l;

	for (i = 0; i < n; i++) {
		if (CHAR(x[i]) == NULL)
			naflag = 1;
		else {
			l = Rstrlen(CHAR(x[i]));
			if (l > xmax)
				xmax = l;
		}
	}
	*fieldwidth = xmax;
	if (quote)
		*fieldwidth += 2;
	if (naflag && *fieldwidth < print_na_width)
		*fieldwidth = print_na_width;
}

void formatLogical(int *x, int n, int *fieldwidth)
{
	int i;

	*fieldwidth = 1;
	for(i=0 ; i<n ; i++) {
		if (x[i] == 1 && *fieldwidth < 4) 
			*fieldwidth = 4;
		if (x[i] == 0 && *fieldwidth < 5 ) {
			*fieldwidth = 5;
			break;
			/* this is the widest it can be so stop */
		}
		if (x[i] == NA_LOGICAL && *fieldwidth <  print_na_width) 
			*fieldwidth =  print_na_width;
	}
}

void formatFactor(int *x, int n, int *fieldwidth, SEXP levels, int nlevs)
{
	int xmax = INT_MIN, naflag = 0;
	int i, l = 0;

	if(isNull(levels)) {
		for(i=0 ; i<n ; i++) {
			if (x[i] == NA_INTEGER || x[i] < 1 || x[i] > nlevs)
				naflag = 1;
			else if (x[i] > xmax)
				xmax = x[i];
		}
		if (xmax > 0)
			l = IndexWidth(xmax);
	}
	else {
		l = 0;
		for(i=0 ; i<n ; i++) {
			if (x[i] == NA_INTEGER || x[i] < 1 || x[i] > nlevs)
				naflag = 1;
			else {
				xmax = strlen(CHAR(STRING(levels)[x[i]-1]));
				if (xmax > l) l = xmax;
			}
		}
	}
	if (naflag) *fieldwidth = print_na_width;
	else *fieldwidth = 1;
	if (l > *fieldwidth) *fieldwidth = l;
}

void formatInteger(int *x, int n, int *fieldwidth)
{
	int xmin = INT_MAX, xmax = INT_MIN, naflag = 0;
	int i, l;

	for (i = 0; i < n; i++) {
		if (x[i] == NA_INTEGER)
			naflag = 1;
		else {
			if (x[i] < xmin) xmin = x[i];
			if (x[i] > xmax) xmax = x[i];
		}
	}

	if (naflag) *fieldwidth = print_na_width;
	else *fieldwidth = 1;

	if (xmin < 0) {
		l = IndexWidth(-xmin) + 1;	/* +1 for sign */
		if (l > *fieldwidth) *fieldwidth = l;
	}
	if (xmax > 0) {
		l = IndexWidth(xmax);
		if (l > *fieldwidth) *fieldwidth = l;
	}
}


	/*  scientific format determination for real numbers.
	 *  This is time-critical code.  It is worth optimizing.
	 *    nsig digits altogether
	 *    kpower+1 digits to the left of .
	 *    kpower+1+sgn including sign  */

	/* F Format is used if all three following conditions hold:
	 * 1) The maximum exponent is less than 7
	 * 2) There are no more than 5 leading zeros after .
	 * 3) There are at most MAXDIG digits in all  */

	/* E Format has the form
	 * [S]X[.XXX]E+XX[X]
	 * This is indicated by setting *e to non-zero (usually 1)
	 * If the additional exponent digit is required *e is set to 2 */


#define MAXDIG		11

static double tbl[] =
{
	0.e0, 1.e0, 1.e1, 1.e2, 1.e3, 1.e4, 1.e5, 1.e6, 1.e7, 1.e8, 1.e9
};

static double eps;


static void scientific(double *x, int *sgn, int *kpower, int *nsig)
{
	register double alpha;
	register double r;
	int j;

	if (*x == 0.0) {
		*kpower = 0;
		*nsig = 1;
		*sgn = 0;
	}
	else {
		if(*x < 0.0) {
			*sgn = 1;
			r = -*x;
		}
		else {
			*sgn = 0;
			r = *x;
		}
		*kpower = floor(log10(r));
		if (abs(*kpower) < 10) {
			if (*kpower >= 0)
				alpha = r / tbl[*kpower + 1];	/* division slow ? */
			else
				alpha = r * tbl[-*kpower + 1];
		}
		else alpha = r / pow(10.0, (double)*kpower);

			/* make sure that alpha is in [1,10) */

		if (10.0 - alpha < eps) {
			alpha = alpha / 10.0;
			*kpower = *kpower + 1;
		}

			/* compute number of digits */

		*nsig = print_digits;
		for (j=1; j <= *nsig; j++) {
			if (fabs(alpha - floor(alpha+0.5)) / alpha < eps) {
				*nsig = j;
				break;
			}
			alpha *= 10.0;
		}
	}
}

#ifdef UNUSED
static void FormatDbl(double *x, int l, int incr, int *m, int *n, int *e)
{
	int left, right, sleft;
	int mnl, mxl, rt, mxs, mxe;
	int neg, sgn;
	int i, kpower, nsig;
	int naflag;

	eps = pow(10.0, -(double)print_digits);

	neg = 0;
	mxl = INT_MIN;
	rt = INT_MAX;
	mxs = INT_MIN;
	mnl = INT_MAX;
	mxe = INT_MIN;
	naflag = 0;

	for (i = 0; i < l; i += incr) {
		if (!FINITE(x[i])) {
			naflag = 1;
		}
		else {
			scientific(&x[i], &sgn, &kpower, &nsig);

			left = kpower + 1;
			sleft = sgn + ((left <= 0) ? 1 : left);
			right = left - nsig;
			if (sgn) neg = 1;

			if (right < rt) rt = right;	/* digits to right of . */
			if (left > mxl) mxl = left;	/* max digits to left of . */
			if (left < mnl) mnl = left;	/* min digits to left of . */
			if (sleft > mxs) mxs = sleft;	/* max left including sign */
			if (nsig >= mxe) mxe = nsig;	/* max sig digits */

		}
	}


	/* F Format is used if all three following conditions hold: */
	/* 1) The maximum exponent is less than 7 */
	/* 2) There are no more than 5 leading zeros after . */
	/* 3) There are at most MAXDIG digits in all */

	/* E Format has the form */
	/* [S]X[.XXX]E+XX[X] */
	/* This is indicated by setting *e to non-zero (usually 1) */
	/* If the additional exponent digit is required *e is set to 2 */

	if (mxl < 8 && mnl > -5 && mxl - rt <= MAXDIG) {
		*e = 0;
		if (mxl != mnl && mxl - rt > MAXDIG)
			rt = mxl - MAXDIG;
		if (mxl < 0)
			mxs = 1 + neg;
		if (rt > 0)
			rt = 0;
		if (rt < -MAXDIG)
			rt = -MAXDIG;
		*n = -rt;
		*m = mxs - rt + (rt != 0);
	}
	else {
		*e = 1;
		*n = mxe - 1;
		*m = neg + (*n > 0) + *n + 5;
		if (mxl > 100 || mnl < -99) {
			*m += 1;
			*e = 2;
		}
	}
	if (naflag && *m < print_na_width)
		*m = print_na_width;
}
#endif

void formatReal(double *x, int l, int *m, int *n, int *e)
{
	int left, right, sleft;
	int mnl, mxl, rt, mxs, mxe;
	int neg, sgn;
	int i, kpower, nsig;
	int naflag;

	eps = pow(10.0, -(double)print_digits);

	neg = 0;
	mxl = INT_MIN;
	rt = INT_MAX;
	mxs = INT_MIN;
	mnl = INT_MAX;
	mxe = INT_MIN;
	naflag = 0;

	for (i=0; i<l; i++) {
		if (!FINITE(x[i])) {
			naflag = 1;
		}
		else {
			scientific(&x[i], &sgn, &kpower, &nsig);

			left = kpower + 1;
			sleft = sgn + ((left <= 0) ? 1 : left);
			right = left - nsig;
			if (sgn) neg = 1;

			if (right < rt) rt = right;	/* digits to right of . */
			if (left > mxl) mxl = left;	/* max digits to left of . */
			if (left < mnl) mnl = left;	/* min digits to left of . */
			if (sleft > mxs) mxs = sleft;	/* max left including sign */
			if (nsig >= mxe) mxe = nsig;	/* max sig digits */
		}
	}
	if (mxl < 8 && mnl > -5 && mxl - rt <= MAXDIG) {
		*e = 0;
		if (mxl != mnl && mxl - rt > MAXDIG)
			rt = mxl - MAXDIG;
		if (mxl < 0)
			mxs = 1 + neg;
		if (rt > 0)
			rt = 0;
		if (rt < -MAXDIG)
			rt = -MAXDIG;
		*n = -rt;
		*m = mxs - rt + (rt != 0);
	}
	else {
		*e = 1;
		*n = mxe - 1;
		*m = neg + (*n > 0) + *n + 5;
		if (mxl > 100 || mnl < -99) {
			*m += 1;
			*e = 2;
		}
	}
	if (naflag && *m < print_na_width)
		*m = print_na_width;
}

#ifdef COMPLEX_DATA
void formatComplex(complex *x, int l,
	int *mr, int *nr, int *er,
	int *mi, int *ni, int *ei)
{
	int left, right, sleft;
	int rt, mnl, mxl, mxs, mxe;
	int i_rt, i_mnl, i_mxl, i_mxs, i_mxe;
	int neg, sgn;
	int i, kpower, nsig;
	int naflag;

	eps = pow(10.0, -(double)print_digits);

	naflag = 0;
	neg = 0;

	rt = INT_MAX;
	mxl = INT_MIN;
	mnl = INT_MAX;
	mxs = INT_MIN;
	mxe = INT_MIN;

	i_rt = INT_MAX;
	i_mxl = INT_MIN;
	i_mnl = INT_MAX;
	i_mxs = INT_MIN;
	i_mxe = INT_MIN;

	for (i=0; i<l; i++) {
		if (!FINITE(x[i].r) || !FINITE(x[i].i)) {
			naflag = 1;
		}
		else {
				/* real part */

			scientific(&(x[i].r), &sgn, &kpower, &nsig);

			left = kpower + 1;
			sleft = sgn + ((left <= 0) ? 1 : left);
			right = left - nsig;
			if (sgn) neg = 1;

			if (right < rt) rt = right;		/* digits to right of . */
			if (left > mxl) mxl = left;		/* max digits to left of . */
			if (left < mnl) mnl = left;		/* min digits to left of . */
			if (sleft > mxs) mxs = sleft;		/* max left including sign */
			if (nsig >= mxe) mxe = nsig;		/* max sig digits */

				/* imaginary part */
				/* this is always unsigned */
				/* we explicitly put the sign in */
				/* when we print */

			scientific(&(x[i].i), &sgn, &kpower, &nsig);

			left = kpower + 1;
			sleft = ((left <= 0) ? 1 : left);
			right = left - nsig;

			if (right < i_rt) i_rt = right;		/* digits to right of . */
			if (left > i_mxl) i_mxl = left;		/* max digits to left of . */
			if (left < i_mnl) i_mnl = left;		/* min digits to left of . */
			if (sleft > i_mxs) i_mxs = sleft;	/* max left including sign */
			if (nsig >= i_mxe) i_mxe = nsig;	/* max sig digits */
		}
	}

		/* overall format for real part */

	if (mnl == INT_MAX) {
		er = 0;
		ei = 0;
		*mr = print_na_width - 2;
		*mi = 0;
		*nr = 0;
		*ni = 0;
		return;
	}

	if (mxl < 8 && mnl > -5 && mxl - rt <= MAXDIG) {
		*er = 0;
		if (mxl != mnl && mxl - rt > MAXDIG) rt = mxl - MAXDIG;
		if (mxl < 0) mxs = 1 + neg;
		if (rt > 0) rt = 0;
		if (rt < -MAXDIG) rt = -MAXDIG;
		*nr = -rt;
		*mr = mxs - rt + (rt != 0);
	}
	else {
		*er = 1;
		*nr = mxe - 1;
		*mr = neg + (*nr > 0) + *nr + 5;
		if (mxl > 100 || mnl < -99) {
			*mr += 1;
			*er = 2;
		}
	}

		/* overall format for imaginary part */

	if (i_mxl < 8 && i_mnl > -5 && i_mxl - i_rt <= MAXDIG) {
		*ei = 0;
		if (i_mxl != i_mnl && i_mxl - i_rt > MAXDIG) i_rt = i_mxl - MAXDIG;
		if (i_mxl < 0) i_mxs = 1;
		if (i_rt > 0) i_rt = 0;
		if (i_rt < -MAXDIG) i_rt = -MAXDIG;
		*ni = -i_rt;
		*mi = i_mxs - i_rt + (i_rt != 0);
	}
	else {
		*ei = 1;
		*ni = i_mxe - 1;
		*mi = (*ni > 0) + *ni + 5;
		if (i_mxl > 100 || i_mnl < -99) {
			*mi += 1;
			*ei = 2;
		}
	}
}
#endif
