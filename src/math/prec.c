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

#define MAXPLACES 14

double prec(double x, double digits)
{
	double pow10, sgn;

	digits = floor(digits+0.5);

	if (digits > MAXPLACES)
		digits = MAXPLACES;
	else if (digits < 1)
		digits = 1;

	if(x == 0.0) return x;
	sgn = 1.0;
	if(x < 0.0) {
		sgn = -sgn;
		x = -x;
	}
	pow10 = pow(10.0, digits-floor(log10(x))-1.0);
	return sgn*floor(x*pow10+0.5)/pow10;
}
