/* 20230924-korkin: minor changes in comments */
#ifndef _COMPLEX_DEFINED
struct _complex
{
    double x, y; // real and imaginary parts //
};
#define _COMPLEX_DEFINED
#endif
//
// +
_complex operator+ (_complex, _complex);
_complex operator+ (double, _complex);
_complex operator+ (_complex, double);
//
// -
_complex operator- (_complex);
_complex operator- (_complex, _complex);
_complex operator- (double, _complex);
_complex operator- (_complex, double);
//
// *
_complex operator* (_complex, _complex);
_complex operator* (double, _complex);
_complex operator* (_complex, double);
//
// /
_complex operator/ (_complex, _complex);
_complex operator/ (double, _complex);
_complex operator/ (_complex, double);
//
// Operators
_complex _CONJG(_complex);	// complex conjugate of complex number 
double _SQ(_complex);		// squared magnitude of a complex number
double _RE(_complex);		// real part of complex number
double _IM(_complex);		// imaginary part of complex number
double _CABS(_complex z);
