/*
prototype:
    double simpson(double const *, int const, double const);
*/
double simpson(double const *f, int const nx, double const dx) {
/*--------------------------------------------------------------------------------------------------
TASK:
	To integrate f(x) numerically using Simpson's method [1].
IN:
	f   d[nx]   f(x)
	nx  i       total number of points, must be ODD & .ge. 3: nx = 3, 5, 7, 9, ...
	dx  d       equidistant integration step
OUT:
	z   d       result of integration
NOTE:
    f = f[0], f[1], ... f[n-1], f[n] - here, n must be EVEN [1], but nx=n+1, is ODD.
	
	Adopted from A.Lyapustin's spectroscopy package IPC [2] with modifications (see original
	subroutine below).
	
	The Simpson and trapezoidal error estimations are dx**4 and dx**2, respectively. Thus "Simpson's
	rule typically gives far superior results for numerical integration, compared to the trapezoidal
	rule, with basically no additional computational cost" [3].

	In case EVEN nx, first or last interval is added using e.g. trapz [4].
REFS:
	1.https://en.wikipedia.org/wiki/Simpson%27s_rule
	2.Lyapustin A, 2003: J. Atmospheric. Sci, 60, 865
	3.https://stackoverflow.com/questions/44915116/how-to-decide-between-scipy-integrate-simps-or-numpy-trapz
	4.https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.simps.html
--------------------------------------------------------------------------------------------------*/
	int
		ix;
	double
		z, c24, sign;
//--------------------------------------------------------------------------------------------------
	sign = 1.0;
	c24 = 4.0;
	z = f[0] + 4.0*f[1] + f[nx-1];
	for(ix = 2; ix < nx-1; ix++) {
		sign = -sign;
		c24 = c24 + sign*2.0;
		z += c24*f[ix];
	} // for ix
	return z*dx/3.0;
} // simpson
/*--------------------------------------------------------------------------------------------------
2020-01-15:
	First created and tested vs original subroutine (below) for nx = 11, dx = 1 & 0.1, 
	x = ix*dx, f=+/-x, x2, ...x5, sin(x). Also tested for x < 0, f < 0, |f|, nx = 101, etc.
	Absolute difference == 0.0e+00.
----------------------------------------------------------------------------------------------------
Original subroutine from gauss.cpp
/* Program calculates the integral using n-point Simpson quadrature.
	n - number of equidistant points. Must BE ODD (n>=3).
	rab[0 : n-1] - integrand function
	dx - distance between points (delta x)
*/
/*
double simpson(double *rab, int n, double dx)
{
  int j;
  double x, z;

	x = rab[0] + 4*rab[1] + rab[n-1];

	for(j=2; j<n-1; j++)
	{
		if(j%2 == 0)	z = 2;	// even number
		else			z = 4;	// odd number
		x += z*rab[j];
	}
	return x*dx/3;
}
*/