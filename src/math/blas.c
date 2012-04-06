#include "Mathlib.h"

extern double square(double);
extern double fsign(double, double);

double dasum(int n, double *dx, int incx)
{
	int i, ix;
	double sum;

	--dx;

	if (n <= 0) return 0.0;

	sum = 0.0;
	ix = 1;
	if (incx < 0)
		ix = (- n + 1) * incx + 1;

	for (i = 1; i <= n; ++i) {
		sum += fabs(dx[ix]);
		ix += incx;
	}
	return sum;
}

void daxpy(int n, double da, double *dx, int incx, double *dy, int incy)
{
	int i, ix, iy;

	--dy;
	--dx;

	if (n <= 0 || da == 0.0) return;

	ix = 1;
	iy = 1;
	if (incx < 0)
		ix = (- n + 1) * incx + 1;
	if (incy < 0)
		iy = (- n + 1) * incy + 1;
	for (i = 1; i <= n; ++i) {
		dy[iy] += da * dx[ix];
		ix += incx;
		iy += incy;
	}
}

void dcopy(int n, double *dx, int incx, double *dy, int incy)
{
	int i, ix, iy;

	--dy;
	--dx;

	if (n <= 0) return;

	ix = 1;
	iy = 1;
	if (incx < 0)
		ix = (- n + 1) * incx + 1;
	if (incy < 0)
		iy = (- n + 1) * incy + 1;

	for (i = 1; i <= n; ++i) {
		dy[iy] = dx[ix];
		ix += incx;
		iy += incy;
	}
}

double ddot(int n, double *dx, int incx, double *dy, int incy)
{
	int i, ix, iy;
	double sum;

	--dy;
	--dx;

	if (n <= 0) return 0.0;

	sum = 0.0;
	ix = 1;
	iy = 1;
	if (incx < 0)
		ix = (- n + 1) * incx + 1;
	if (incy < 0)
		iy = (- n + 1) * incy + 1;

	for (i = 1; i <= n; ++i) {
		sum += dx[ix] * dy[iy];
		ix += incx;
		iy += incy;
	}
	return sum;
}

double dmach(int job)
{
	double eps, tiny, huge, s;
	eps = 1.0;
	do {
		eps = eps / 2.0;
		s = 1.0 + eps;
	}
	while(s > 1.0);
	eps = 2.0 * eps;
	s = 1.0;
	do {
		tiny = s;
		s = s / 16.0;
	}
	while(s * 1.0 != 0.0);
	tiny = (tiny / eps) * 100.0;
	huge = 1.0 / tiny;

	if (job==1) return eps;
	if (job==2) return tiny;
	if (job==3) return huge;
	return 0.0;
}

double dnrm2(int n, double *x, int incx)
{
	double norm, scale, absxi;
	int ix, ixlim;
	double ssq;

	--x;

	/* Function Body */
	if (n < 1 || incx < 1) {
		norm = 0.0;
	} else if (n == 1) {
		norm = fabs(x[1]);
	} else {
		scale = 0.0;
		ssq = 1.0;
		ixlim = (n - 1) * incx + 1;
		ix = 1;
		for(;;) {
			if (x[ix] != 0.0) {
				absxi = fabs(x[ix]);
				if (scale < absxi) {
					ssq = ssq * square(scale / absxi) + 1.0;
					scale = absxi;
				} else {
					ssq += square(absxi / scale);
				}
			}
			if(ix == ixlim) break;
			ix += incx;
		}
		norm = scale * sqrt(ssq);
	}
	return norm;
}

void drot(int n, double *dx, int incx, double *dy, int incy, double c, double s)
{
	int i, ix, iy;
	double dtemp;

	--dy;
	--dx;

	if (n <= 0) return;

	ix = 1;
	iy = 1;
	if (incx < 0)
		ix = (- n + 1) * incx + 1;
	if (incy < 0)
		iy = (- n + 1) * incy + 1;

	for (i = 1; i <= n; ++i) {
		dtemp  = c * dx[ix] + s * dy[iy];
		dy[iy] = c * dy[iy] - s * dx[ix];
		dx[ix] = dtemp;
		ix += incx;
		iy += incy;
	}
}

drotg(double *da, double *db, double *c, double *s)
{
	double r, scale, z, roe;

	roe = *db;
	if (fabs(*da) > fabs(*db)) {
		roe = *da;
	}
	scale = fabs(*da) + fabs(*db);
	if (scale == 0.0) {
		*c = 1.;
		*s = 0.;
		r = 0.;
	}
	else {
		r = scale * sqrt(square(*da / scale) + square(*db / scale));
		r = fsign(1.0, roe) * r;
		*c = *da / r;
		*s = *db / r;
	}
	z = *s;
	if (fabs(*c) > 0.0 && fabs(*c) <= *s) {
		z = 1.0 / *c;
	}
	*da = r;
	*db = z;
}

void dscal(int n, double da, double *dx, int incx)
{
	int i, ix;

	--dx;

	if (n <= 0) return;

	ix = 1;
	if (incx < 0)
		ix = (- n + 1) * incx + 1;
	for (i = 1; i <= n; ++i) {
		dx[ix] = da * dx[ix];
		ix += incx;
	}
}

void dswap(int n, double *dx, int incx, double *dy, int incy)
{
	int i, ix, iy;
	double dtemp;

	--dy;
	--dx;

	if (n <= 0) return;

	ix = 1;
	iy = 1;
	if (incx < 0)
		ix = (- n + 1) * incx + 1;
	if (incy < 0)
		iy = (- n + 1) * incy + 1;
	for (i = 1; i <= n; ++i) {
		dtemp = dx[ix];
		dx[ix] = dy[iy];
		dy[iy] = dtemp;
		ix += incx;
		iy += incy;
	}
}

int idamax(int n, double *dx, int incx)
{
	int i, imax, ix;
	double xmax;

	--dx;

	if (n < 1) return 0;
	if (n == 1) return 1;

	imax = 1;
	ix = 1;
	if (incx < 0)
		ix = (- n + 1) * incx + 1;
	xmax = fabs(dx[ix]);
	ix += incx;
	for (i = 2; i <= n; ++i) {
		if (fabs(dx[ix]) > xmax) {
			imax = i;
			xmax = fabs(dx[ix]);
		}
		ix += incx;
	}
	return imax;
}
