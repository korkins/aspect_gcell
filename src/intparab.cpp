#include <math.h> /* exp, log */
/*
prototype:
	double intparab(double const, double const,
	                double const, double const, double const,
	                double const, double const, double const);
*/
double intparab(double const xa, double const xb,
	            double const x1, double const x2, double const x3,
	            double const y1, double const y2, double const y3) {
/*--------------------------------------------------------------------------------------------------
TASK:
	To approximate the 3 points, [xi, yi], with a parbola and to integarte analyticaly from xa to xb
IN:
    xa, xb     d   interval of integration, xa < xb
	[xi, yi]   d   3 points
OUT:
	intparab   d   result of integration of parabola
NOTE:
	x must be within x1 & x2.
	Cramer's rule [1] is used to draw parabola.
REFS:
	1. https://en.wikipedia.org/wiki/Cramer%27s_rule
--------------------------------------------------------------------------------------------------*/
	double
		det, det1, det2, det3, a, b, c;
//--------------------------------------------------------------------------------------------------

	det = x1*x1*x2 + x1*x3*x3 + x2*x2*x3 - x2*x3*x3 - x1*x2*x2 - x1*x1*x3;
	det1 = y1*x2 + x1*y3 + y2*x3 - x2*y3 - x1*y2 - y1*x3;
	det2 = x1*x1*y2 + y1*x3*x3 + x2*x2*y3 - y2*x3*x3 - y1*x2*x2 - x1*x1*y3;
	det3 = x1*x1*x2*y3 + x1*y2*x3*x3 + y1*x2*x2*x3 - y1*x2*x3*x3 - x1*x2*x2*y3 - x1*x1*y2*x3;
	a = det1/det;
	b = det2/det;
	c = det3/det;
	//return  a*(x2*x2*x2 - x*x*x)/3.0 + b*(x2*x2 - x*x)/2.0 + c*(x2 - x);
	return  a*(xb*xb*xb - xa*xa*xa)/3.0 + b*(xb*xb - xa*xa)/2.0 + c*(xb - xa);
} // intparab
/*--------------------------------------------------------------------------------------------------
2020-02-09:
    Modification: 'integarte from x to x2' ->  'integarte from xa to xb'
2020-02-08:
	First created and tested
		a) test the parabola itself: x = [1, 3, 4]; y(x) = 3x2 + 1.5x + 2 - a,b,c are ok.
		b) xo = 2; int(y(x), x=2:3) = 24.75 - ok
		c) xo=x1=1: int(y(x), x=1:3)=36.0 - ok
		d) xo=x2=3: int(y(x), x=3:3)=0.0 - ok
--------------------------------------------------------------------------------------------------*/