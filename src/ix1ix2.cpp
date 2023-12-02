/*
prototype:
    int ix1ix2(double const, double const, double const *, int const, int &, int &);
*/
int ix1ix2(double const x0, double const dx, double const *x, int const nx, // in
	       int &ix1, int &ix2) {                                            // out
/*--------------------------------------------------------------------------------------------------
TASK:
	To get boundary indices for elements of an array, x, that fall within +/- dx from a value x0.
	Elements of x are arranged in ascending order: x[i] < x[i+1].
	If successful, code=1 is returned, otherwise code = -1.
IN:
	x0   d       central point of interval [x0-dx, x0+dx]
	dx   d       half-width of the interval, dx > 0.0
	x    d[nx]   array of elements in ascending order
	nx   i       length of x
OUT:
	ix1, ix2   i   zero-offset left and right boundaries, respectively
NOTE:
	ix1 = ix2 = code = -1 are returned if no boundaries found.
	ix1 & ix2 are zero-offset (ready to use for C-arrays)
	ix1 = ix2 is allowed (only one value from 'x' is covered by x0+/-dx)
REFS:
	-
--------------------------------------------------------------------------------------------------*/
	int
		code;
	double
		x0_left, x0_right, x_min, x_max;
//--------------------------------------------------------------------------------------------------
	code = -1;
	ix1 = -1;
	ix2 = -1;
//
	x0_left  = x0 - dx;
	x0_right = x0 + dx;
	x_min = x[0];
	x_max = x[nx-1];
//  
	if (!(x0_left > x_max || x0_right < x_min)) // check for overlap; THINKME: simplify condition
	{
//		left index
		if (x0_left <= x_min)
			ix1 = 0;
		else
		{
			ix1 = 0;
			while (x[ix1] < x0_left && ix1 < nx)
				ix1 += 1;
		} // if x0_left <= x_min
//		right index
		if (x_max <= x0_right)
			ix2 = nx-1;
		else
		{
			ix2 = nx-1;
			while (x[ix2] > x0_right && ix2 >= ix1)
				ix2 -= 1;
		} // if xmax <= x0_right
	} // if (check_for_overlap)
//
	if (ix1 > -1 && ix2 > -1)
		code = 1;
	return(code);
} // inu1inu2
/*--------------------------------------------------------------------------------------------------
2020-05-15:
    minor changes in comments;
2020-01-08:
	First created and tested using nx=101, x=[0:1:100], x[ix]=1.0*ix
		Test 1 (no overlap):
			in  >> x0 = +/-200, dx=10.0
			out >> code=-1 - OK.
		Test 2 (complete overlap: x within x0+/-dx):
			in  >> x0 = 50, dx=100.0; x0=+/-200, dx=250
			out >> code=+1, ix1=0, ix2=100 - OK.
		Test 3 (x_min within, x_max not):
			in  >> x0 = 10, dx=20
			out >> code=+1, ix1=0, ix2=30 - OK.
		Test 4 (x_min is off, x_max within):
			in  >> x0 = 120, dx=30
			out >> code=+1, ix1=90, ix2=100 - OK.
		Test 5 (exact overlap: x == x0+/-dx):
			in  >> x0 = 50, dx=50
			out >> code=+1, ix1=0, ix2=100 - OK.
		Test 6 ("non-int" step):
			in  >> x0 = 50, dx=49.5
			out >> code=+1, ix1=1, ix2=99 - OK.
		Test 7 ( x0+/-dx within x):
			in  >> x0 = 45, dx=15.5
			out >> code=+1, ix1=30, ix2=60 - OK.
		Test 8 ( exactly one point):
			in  >> x0 = 13, dx=0
			out >> code=+1, ix1=13, ix2=13 - OK.
*/