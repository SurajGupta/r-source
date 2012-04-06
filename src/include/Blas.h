#ifndef BLAS_H_
#define BLAS_H_

#include <Mathlib.h>

extern double	F77_SYMBOL(dasum)(int*, double*, int*);
extern int	F77_SYMBOL(daxpy)(int*, double*, double*, int*, double*, int*);
extern int	F77_SYMBOL(dcopy)(int*, double*, int*, double*, int*);
extern double	F77_SYMBOL(ddot)(int*, double*, int*, double*, int*);
extern double	F77_SYMBOL(dmach)(int*);
extern double	F77_SYMBOL(dnrm2)(int*, double*, int*);
extern int	F77_SYMBOL(drot)(int*, double*, int*, double*, int*, double*, double*);
extern int	F77_SYMBOL(drotg)(double*, double*, double*, double*);
extern int	F77_SYMBOL(dscal)(int*, double*, double*, int*);
extern int	F77_SYMBOL(dswap)(int*, double*, int*, double*, int*);
extern int	F77_SYMBOL(idamax)(int*, double*, int*);

#endif
