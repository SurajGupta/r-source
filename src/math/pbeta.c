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

/* Reference:
 * Bosten, N. E. and E. L. Battiste (1974).
 * Remark on Algorithm 179,
 * Communications of the ACM, 17, 153.
 *
 * x = value to which function is to be integrated. x must be in (0,1).
 * pin = input (1st) parameter (must be greater than 0)
 * qin = input (2nd) parameter (must be greater than 0)
 */

#include "Mathlib.h"

double pbeta(double x, double pin, double qin)
{
	int i, ib, n;
	double c, finsum, p, p1, ps, q, term, xb, y;
	double betai;
	static double eps = 0.0;
	static double alneps = 0.0;
	static double sml = 0.0;
	static double alnsml = 0.0;

	if (eps == 0.0) {
		eps = DBL_EPSILON;
		sml = DBL_MIN;
		alneps = log(eps);
		alnsml = log(sml);
	}
	if (pin <= 0.0 || qin <= 0.0)
		DOMAIN_ERROR;
	if (x <= 0.0)
		return 0.0;
	if (x >= 1.0)
		return 1.0;

	y = x;
	p = pin;
	q = qin;
	if (q > p || x >= 0.8)
		if (x >= 0.2) {
			y = 1.0 - y;
			p = qin;
			q = pin;
		}
	if ((p + q) * y / (p + 1.) < eps) {
		betai = 0.0;
		xb = p * log(fmax2(y, sml)) - log(p) - lbeta(p, q);
		if (xb > alnsml && y != 0.0)
			betai = exp(xb);
		if (y != x || p != pin)
			betai = 1.0 - betai;
	} else {
		/* evaluate the infinite sum first */

		ps = q - (int) q;
		if (ps == 0.0)
			ps = 1.0;
		xb = p * log(y) - lbeta(ps, p) - log(p);
		betai = 0.0;
		if (xb >= alnsml) {
			betai = exp(xb);
			term = betai * p;
			if (ps != 1.0) {
				n = fmax2(alneps / log(y), 4.0);
				for (i = 1; i <= n; i++) {
					term = term * (i - ps) * y / i;
					betai = betai + term / (p + i);
				}
			}
		}
		/* now evaluate the finite sum, maybe */

		if (q > 1.0) {
			xb = p * log(y) + q * log(1.0 - y) - lbeta(p, q) - log(q);
			ib = fmax2(xb / alnsml, 0.0);
			term = exp(xb - ib * alnsml);
			c = 1.0 / (1.0 - y);
			p1 = q * c / (p + q - 1.);
			finsum = 0.0;
			n = q;
			if (q == n)
				n = n - 1;
			for (i = 1; i <= n; i++) {
				if (p1 <= 1.0 && term / eps <= finsum)
					break;
				term = (q - (i - 1)) * c * term / (p + q - i);
				if (term > 1.0)
					ib = ib - 1;
				if (term > 1.0)
					term = term * sml;
				if (ib == 0)
					finsum = finsum + term;
			}
			betai = betai + finsum;
		}
		if (y != x || p != pin)
			betai = 1.0 - betai;
		betai = fmax2(fmin2(betai, 1.0), 0.0);
	}
	return betai;
}
