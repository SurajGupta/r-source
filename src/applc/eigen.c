/*
 *  R : A Computer Langage for Statistical Data Analysis
 *  Copyright (C) 1996  Robert Gentleman and Ross Ihaka
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

#include "Blas.h"

/*
 *     this subroutine is a translation of the algol procedure tred2,
 *     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
 *     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
 *
 *     this subroutine reduces a real symmetric matrix to a
 *     symmetric tridiagonal matrix using and accumulating
 *     orthogonal similarity transformations.
 *
 *     on input
 *
 *        nm must be set to the row dimension of two-dimensional
 *          array parameters as declared in the calling program
 *          dimension statement.
 *
 *        n is the order of the matrix.
 *
 *        a contains the real symmetric input matrix.  only the
 *          lower triangle of the matrix need be supplied.
 *
 *     on output
 *
 *        d contains the diagonal elements of the tridiagonal matrix.
 *
 *        e contains the subdiagonal elements of the tridiagonal
 *          matrix in its last n-1 positions.  e(1) is set to zero.
 *
 *        z contains the orthogonal transformation matrix
 *          produced in the reduction.
 *
 *        a and z may coincide.  if distinct, a is unaltered.
 *
 *     questions and comments should be directed to burton s. garbow,
 *     mathematics and computer science div, argonne national laboratory
 *
 *     this version dated august 1983.
 */

static void tred2(int nm, int n, double *a, double *d, double *e, double *z)
{
	int i, j, k, l, ii, jp1;
	double f, g, h, hh, scale;

	a -= (nm + 1);
	d -= 1;
	e -= 1;
	z -= (nm + 1);

	for(i=1;i<=n;i++) {
		for(j=i;j<=n;j++)
			z[j+i*nm] = a[j+i*nm];
		d[i] = a[n+i*nm];
	}
	if (n!=1) {
		for(ii=2;ii<=n;ii++) {
			i = n+2-ii;
			l = i-1;
			h = 0.0;
			scale = 0.0;
			if (l>=2) {
				for(k=1;k<=l;k++)
					scale = scale+fabs(d[k]);
				if (scale!=0.0) {
					for(k=1;k<=l;k++) {
						d[k] = d[k]/scale;
						h = h+d[k]*d[k];
					}
					f = d[l];
					g = -fsign(sqrt(h), f);
					e[i] = scale*g;
					h = h-f*g;
					d[l] = f-g;
					for(j=1;j<=l;j++)
						e[j] = 0.0;
					for(j=1;j<=l;j++) {
						f = d[j];
						z[j+i*nm] = f;
						g = e[j]+z[j+j*nm]*f;
						jp1 = j+1;
						if (l>=jp1)
							for(k=jp1;k<=l;k++) {
								g = g+z[k+j*nm]*d[k];
								e[k] = e[k]+z[k+j*nm]*f;
							}
						e[j] = g;
					}
					f = 0.0;
					for(j=1;j<=l;j++) {
						e[j] = e[j]/h;
						f = f+e[j]*d[j];
					}
					hh = f/(h+h);
					for(j=1;j<=l;j++)
						e[j] = e[j]-hh*d[j];
					for(j=1;j<=l;j++) {
						f = d[j];
						g = e[j];
						for(k=j;k<=l;k++)
							z[k+j*nm] = z[k+j*nm]-f*e[k]-g*d[k];
						d[j] = z[l+j*nm];
						z[i+j*nm] = 0.0;
					}
					goto lab10;
				}
			}
			e[i] = d[l];
			for(j=1;j<=l;j++) {
				d[j] = z[l+j*nm];
				z[i+j*nm] = 0.0;
				z[j+i*nm] = 0.0;
			}
		lab10:	d[i] = h;
		}
		for(i=2;i<=n;i++) {
			l = i-1;
			z[n+l*nm] = z[l+l*nm];
			z[l+l*nm] = 1.0;
			h = d[i];
			if (h!=0.0) {
				for(k=1;k<=l;k++)
					d[k] = z[k+i*nm]/h;
				for(j=1;j<=l;j++) {
					g = 0.0;
					for(k=1;k<=l;k++)
						g = g+z[k+i*nm]*z[k+j*nm];
					for(k=1;k<=l;k++)
						z[k+j*nm] = z[k+j*nm]-g*d[k];
				}
			}
			for(k=1;k<=l;k++)
				z[k+i*nm] = 0.0;
		}
	}
	for(i=1;i<=n;i++) {
		d[i] = z[n+i*nm];
		z[n+i*nm] = 0.0;
	}
	z[n+n*nm] = 1.0;
	e[1] = 0.0;
}

/* finds dsqrt(a**2+b**2) without overflow or destructive underflow */

static double pythag(double a, double b)
{
	double p, r, s, t, tmp, u;

	p = fmax2(fabs(a), fabs(b));
	if (p != 0.0) {
		
		/* r = (fmin2(fabs(a), fabs(b))/p)**2 */

		tmp = fmin2(fabs(a), fabs(b))/p;
		r = tmp * tmp;
		for(;;) {
			t = 4.0 + r;
			if (t == 4.0)
				break;
			s = r / t;
			u = 1.0 + 2.0 * s;
			p = u * p;

			/* r = (s / u)**2 * r */

			tmp = (s / u);
			r = tmp * tmp * r;
		}
	}
	return p;
}

/*
 *     this subroutine is a translation of the algol procedure tql2,
 *     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
 *     wilkinson.
 *     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
 *
 *     this subroutine finds the eigenvalues and eigenvectors
 *     of a symmetric tridiagonal matrix by the ql method.
 *     the eigenvectors of a full symmetric matrix can also
 *     be found if  tred2  has been used to reduce this
 *     full matrix to tridiagonal form.
 *
 *     on input
 *
 *        nm must be set to the row dimension of two-dimensional
 *          array parameters as declared in the calling program
 *          dimension statement.
 *
 *        n is the order of the matrix.
 *
 *        d contains the diagonal elements of the input matrix.
 *
 *        e contains the subdiagonal elements of the input matrix
 *          in its last n-1 positions.  e(1) is arbitrary.
 *
 *        z contains the transformation matrix produced in the
 *          reduction by  tred2, if performed.  if the eigenvectors
 *          of the tridiagonal matrix are desired, z must contain
 *          the identity matrix.
 *
 *      on output
 *
 *        d contains the eigenvalues in ascending order.  if an
 *          error exit is made, the eigenvalues are correct but
 *          unordered for indices 1,2,...,ierr-1.
 *
 *        e has been destroyed.
 *
 *        z contains orthonormal eigenvectors of the symmetric
 *          tridiagonal (or full) matrix.  if an error exit is made,
 *          z contains the eigenvectors associated with the stored
 *          eigenvalues.
 *
 *        ierr is set to
 *          zero       for normal return,
 *          j          if the j-th eigenvalue has not been
 *                     determined after 30 iterations.
 *
 *     calls pythag for  dsqrt(a*a + b*b) .
 *
 *     questions and comments should be directed to burton s. garbow,
 *     mathematics and computer science div, argonne national laboratory
 *
 *     this version dated august 1983.
 */


void tql2(int nm, int n, double *d, double *e, double *z, int *ierr)
{
	int i, j, k, l, m, ii, l1, l2, mml;
	double c, c2, c3, dl1, el1, f, g, h, p, r, s, s2, tst1, tst2;

	d -= 1;
	e -= 1;
	z -= (nm + 1);

	*ierr = 0;
	if( n!=1 ) {
		for(i= 2;i<= n;i++)
			e[i-1] = e[i];
		f = 0.0;
		tst1 = 0.0;
		e[n] = 0.0;
		for(l= 1;l<= n;l++) {
			j = 0;
			h = fabs(d[l])+fabs(e[l]);
			if (tst1<h)
				tst1 = h;
			for(m= l;m<= n;m++) {
				tst2 = tst1+fabs(e[m]);
				if (tst2==tst1)
					break;
			}
			if (m!=l)
				do {
					if (j == 30) {
						*ierr = l;
						return;
					}
					j = j+1;
					l1 = l+1;
					l2 = l1+1;
					g = d[l];
					p = (d[l1]-g)/(2.0*e[l]);
					r = pythag(p, 1.0);
					d[l] = e[l]/(p+fsign(r, p));
					d[l1] = e[l]*(p+fsign(r, p));
					dl1 = d[l1];
					h = g-d[l];
					if (l2<=n)
						for(i= l2;i<= n;i++)
							d[i] = d[i]-h;
					f = f+h;
					p = d[m];
					c = 1.0;
					c2 = c;
					el1 = e[l1];
					s = 0.0;
					mml = m-l;
					for(ii= 1;ii<= mml;ii++) {
						c3 = c2;
						c2 = c;
						s2 = s;
						i = m-ii;
						g = c*e[i];
						h = c*p;
						r = hypot(p, e[i]);
						e[i+1] = s*r;
						s = e[i]/r;
						c = p/r;
						p = c*d[i]-s*g;
						d[i+1] = h+s*(c*g+s*d[i]);
						for(k= 1;k<= n;k++) {
							h = z[k+(i+1)*nm];
							z[k+(i+1)*nm] = s*z[k+i*nm]+c*h;
							z[k+i*nm] = c*z[k+i*nm]-s*h;
						}
					}
					p = -s*s2*c3*el1*e[l]/dl1;
					e[l] = s*p;
					d[l] = c*p;
					tst2 = tst1+fabs(e[l]);
				}
					while(tst2>tst1);
			d[l] = d[l]+f;
		}
		for(ii= 2;ii<= n;ii++) {
			i = ii-1;
			k = i;
			p = d[i];
			for(j= ii;j<= n;j++)
				if (d[j]<p) {
					k = j;
					p = d[j];
				}
			if (k!=i) {
				d[k] = d[i];
				d[i] = p;
				for(j= 1;j<= n;j++) {
					p = z[j+i*nm];
					z[j+i*nm] = z[j+k*nm];
					z[j+k*nm] = p;
				}
			}
		}
	}
}

int eigen(int *p, double *c, double *vecs, double *vals,
	double *wrk, int *ierr)
{
	tred2(*p, *p, c, vals, wrk, vecs);
	tql2(*p, *p, vals, wrk, vecs, ierr);
}
