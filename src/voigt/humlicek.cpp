#include <math.h>
#include "cmplx.h"
//
_complex approx1(_complex), approx2(_complex, _complex), approx3(_complex), approx4(_complex, _complex);
//
double humlicek(double x, double y)
/*--------------------------------------------------------------------------------------------------
TASK:
	To calculate K(x, y) = (y/pi) * integral{exp(-t2) / (y2 + (x - t)) dt, t = -inf ... + inf} using
	algorithm described in [1]; see Sec.5 for original code source code and [2: Appendix 2] for a
	refactored version.
IN:
	x, y   d   arguments of K(x, y); x > 0, y > 0
OUT:
	z.x    d   real part of the integral K(x, y); estimated relative error < 1.E-4 [1: Sec.5].
NOTE:
	This subroutine is a revised (coments, code tabulation) version of that from [3].
	Structure '_complex' is defined in <math.h>. If not, it is also defined in "cmplx.h".
REFS:
	1. Humlicek J, 1982: "Optimized computation of the Voigt and complex probability functions",
	   JQSRT 27(4), pp.437-444. https://doi.org/10.1016/0022-4073(82)90078-4
	2. Schreier F, 1992: "The Voigt and complex error function: A comparison of computational methods",
       JQSRT 48(5–6), pp.743-762. https://doi.org/10.1016/0022-4073(92)90139-U.v
	3. Lyapustin A.I., 2003: "Interpolation and Profile Correction (IPC) Method for Shortwave
	   Radiative Transfer in Spectral Intervals of Gaseous Absorption", J. Atmos. Sci., v.60, pp.865–871,
	   https://journals.ametsoc.org/view/journals/atsc/60/6/1520-0469_2003_060_0865_iapcim_2.0.co_2.xml
--------------------------------------------------------------------------------------------------*/
{
	double
		s;
	_complex // defines real x and imaginary y parts of z = x + i*y (see math.h and cmplx.h)
		t, u, z;
//--------------------------------------------------------------------------------------------------
//
	t.x =  y;
	t.y = -x;
//
    if (y >= 15.0)					// all points are in region I
		z = approx1(t);
    else if (y < 15.0 && y >= 5.5)	// points are in region I or II
	{
		s = fabs(x) + y;
		if(s >= 15.0)
			z = approx1(t);
		else
		{
			u = t*t;
			z = approx2(t, u);
		}
	}
	else if (y < 5.5 && y > 0.75)
	{
		s = fabs(x) + y;
		if(s >= 15.0)
			z = approx1(t);
		else if (s < 5.5)
			z = approx3(t);
		else
		{
			u = t*t;
			z = approx2(t, u);
		}
	}
	else
	{
		s = fabs(x) + y;
		if (s >= 15.0)
			z = approx1(t);
		else if (s < 15.0 && s >= 5.5)
		{
			u = t*t;
			z = approx2(t, u);
		}
		else if (s < 5.5 && y >= (0.195 * fabs(x) - 0.176) ) // region III
			z = approx3(t);
		else
		{
			u = t*t;
			s = exp(u.x);
			z.x = s*cos(u.y);
			z.y = s*sin(u.y);
			z = z - approx4(t, u);
		}
	 }
	return(z.x);
}
//
_complex approx1(_complex t)
{
	return( (t*0.5641896)/(0.5 + t*t) );
}
//
_complex approx2(_complex t, _complex u)
{
	return( (t*(1.410474 + u*0.5641896))/(0.75 + u*(3.0 + u)) );
}
//
_complex approx3(_complex t)
{
	return( (16.4955 + t * (20.20933 + t*(11.96482 + t*
			(3.778987 + t*0.5642236)))) / (16.4955 + t*(38.82363 + t*
			(39.27121 + t*(21.69274 + t*(6.699398 + t)))))  );
 }
//
_complex approx4(_complex t, _complex u)
{
 return( (t * (36183.31 - u*(3321.99 - u*(1540.787 - u*
		 (219.031 - u*(35.7668 - u*(1.320522 -u*0.56419))))))) /
         (32066.6 - u*(24322.8 - u*(9022.23 - u*(2186.18 - u*(364.219 -
         u*(61.5704 - u*(1.84144 - u)))))))  );
 }
/*--------------------------------------------------------------------------------------------------
20230924: revised, commented, minor changes in code structure.
--------------------------------------------------------------------------------------------------*/