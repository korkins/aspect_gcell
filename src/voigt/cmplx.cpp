/***********************************************************************
*		Define overload operators +, -, *, / for structure _complex     * 
************************************************************************/
#include <math.h>
#include "cmplx.h" // Korkin_2020-05-02: "Cmplx.h"->"cmplx.h"

_complex operator+ (_complex a, _complex b)
{	_complex c;
	c.x = a.x + b.x;
	c.y = a.y + b.y;
	return c;
}
_complex operator+ (double a, _complex b)
{	_complex c;
	c.x = a + b.x;
	c.y = b.y;
	return c;
}
_complex operator+ (_complex a, double b)
{	_complex c;
	c.x = a.x + b;
	c.y = a.y;
	return c;
}

_complex operator- (_complex a)
{	_complex c;
	c.x = -a.x;
	c.y = -a.y;
	return c;
}
_complex operator- (_complex a, _complex b)
{	_complex c;
	c.x = a.x - b.x;
	c.y = a.y - b.y;
	return c;
}
_complex operator- (double a, _complex b)
{	_complex c;
	c.x = a - b.x;
	c.y = -b.y;
	return c;
}
_complex operator- (_complex a, double b)
{	_complex c;
	c.x = a.x - b;
	c.y = a.y;
	return c;
}

_complex operator* (_complex a, _complex b)
{	_complex c;
	c.x = a.x*b.x - a.y*b.y;
	c.y = a.x*b.y + a.y*b.x;
	return c;
}
_complex operator* (double a, _complex b)
{	_complex c;
	c.x = b.x * a;
	c.y = b.y * a;
	return c;
}
_complex operator* (_complex a, double b)
{	_complex c;
	c.x = a.x * b;
	c.y = a.y * b;
	return c;
}

_complex operator/ (_complex a, _complex b)
{	_complex c;
	double v = b.x*b.x + b.y*b.y;
	c.x = (a.x*b.x + a.y*b.y)/v;
	c.y = (a.y*b.x - a.x*b.y)/v;
	return c;
}
_complex operator/ (double a, _complex b)
{	_complex c;
	double v = b.x*b.x + b.y*b.y;		
	c.x = a*b.x/v;
	c.y = -a*b.y/v;
	return c;
}
_complex operator/ (_complex a, double b)
{	_complex c;
	c.x = a.x/b;
	c.y = a.y/b;
	return c;
}


//	Calculate squared magnitude of a complex number (20230924-korkin: comment corrected)
double _SQ(_complex z)
{
	return(z.x*z.x + z.y*z.y);
}

//	Calculate complex conjugate of a complex number
_complex _CONJG(_complex z)
{
  _complex q;
	q.x = z.x;
	q.y = -z.y;
	return q;
}

double _RE(_complex z)
{
	return(z.x);
}

double _IM(_complex z)
{
	return(z.y);
}

double _CABS(_complex z)
{
	return sqrt( z.x*z.x + z.y*z.y );
}

