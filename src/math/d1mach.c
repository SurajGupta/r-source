#include "Mathlib.h"

double D1MACH(int i)
{
	switch(i) {
	case 1:
		return DBL_MIN;
	case 2:
		return DBL_MAX;
	case 3:
		return pow((double)I1MACH(10), -(double)I1MACH(14));
	case 4:
		return pow((double)I1MACH(10), 1-(double)I1MACH(14));
	case 5:
		return log10(2.0);
	}
}

double F77_SYMBOL(d1mach)(int *i)
{
	return D1MACH(*i);
}
