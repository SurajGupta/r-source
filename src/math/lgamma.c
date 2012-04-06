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

#include "Mathlib.h"

int signgam = 0;

/* log(2*pi)/2 and pi */

static double hl2pi = 0.9189385332046727417803297;
static double xpi = M_PI;

 /* Coefficients from Cheney and Hart */

#define M 6
#define N 8
static double p1[] =
{
	0.83333333333333101837e-1,
	-.277777777735865004e-2,
	0.793650576493454e-3,
	-.5951896861197e-3,
	0.83645878922e-3,
	-.1633436431e-2,
};
static double p2[] =
{
	-.42353689509744089647e5,
	-.20886861789269887364e5,
	-.87627102978521489560e4,
	-.20085274013072791214e4,
	-.43933044406002567613e3,
	-.50108693752970953015e2,
	-.67449507245925289918e1,
	0.0,
};
static double q2[] =
{
	-.42353689509744090010e5,
	-.29803853309256649932e4,
	0.99403074150827709015e4,
	-.15286072737795220248e4,
	-.49902852662143904834e3,
	0.18949823415702801641e3,
	-.23081551524580124562e2,
	0.10000000000000000000e1,
};

static double posarg(double);
static double negarg(double);
static double asform(double);

double lgamma(double arg)
{
	signgam = 1.0;
	if (arg <= 0.0)
		return (negarg(arg));
	if (arg > 8.0)
		return (asform(arg));
	return (log(posarg(arg)));
}

/* Equation 6.1.41 Abramowitz and Stegun */
/* See also ACM algorithm 291 */

static double asform(arg)
double arg;
{
	double log();
	double n, argsq;
	int i;

	argsq = 1. / (arg * arg);
	for (n = 0, i = M - 1; i >= 0; i--) {
		n = n * argsq + p1[i];
	}
	return ((arg - .5) * log(arg) - arg + hl2pi + n / arg);
}

static double negarg(arg)
double arg;
{
	double temp;
	double log(), sin(), posarg();

	arg = -arg;
	temp = sin(xpi * arg);
	if (temp == 0.0)
		DOMAIN_ERROR;
	if (temp < 0.0)
		temp = -temp;
	else
		signgam = -1;
	return (-log(arg * posarg(arg) * temp / xpi));
}

static double posarg(arg)
double arg;
{
	double n, d, s;
	register i;

	if (arg < 2.0)
		return (posarg(arg + 1.0) / arg);
	if (arg > 3.0)
		return ((arg - 1.0) * posarg(arg - 1.0));

	s = arg - 2.;
	for (n = 0, d = 0, i = N - 1; i >= 0; i--) {
		n = n * s + p2[i];
		d = d * s + q2[i];
	}
	return (n / d);
}
