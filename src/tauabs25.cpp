//#include <stdio.h>   /* printf, fprintf */
#include <math.h>      /* exp, log */
//
double simpson(double const *, int const, double const);
double intparab(double const, double const, double const, double const, double const,
	                                        double const, double const, double const);
/*
prototype:
    void tauabs25(double const *, double const *, double const *, int const, double *);
*/
void tauabs25(double const *ext_km, double const *zkm_mod, double const *ztau, int const ntau, // in
			  double *tau) {                                                                   // out
/*--------------------------------------------------------------------------------------------------
TASK:
	To integrate exteinction height profile from ztau < 25km to TOA and get tau_absorption
IN:
	ext_km   d[nz_mod]   extinction profile, 1/km, for MODTRAN atmosphere
	ztau     d[ntau]     array of heights, must not exceed 25 km
	ntau     i           len(ztau)
OUT:
	tau      d[ntau]     int{ext_km dz, z=ztau:z_toa}
NOTE:
	This subroutine is desinged for MODTRAN profiles [1].
REFS:
	1. http://modtran.spectral.com/modtran_faq
--------------------------------------------------------------------------------------------------*/
	int
		ix, jx, iz;
	double
		tauZ_25, tau25_120;
//--------------------------------------------------------------------------------------------------
//
	tau25_120 = simpson(&ext_km[25], 11, 2.5) + simpson(&ext_km[35], 15, 5.0); //intexp(&ext_km[35], &zkm_mod[35], 15);
	for (iz = 0; iz < ntau; iz++)
	{
		ix = (int) floor(ztau[iz]);   // left index of 1km-bin; WARNING: floor returns DOUBLE!
		if (ix > 24)
			tauZ_25 = 0.0;
		else if (ix == 24)      // bin=[24:25]km
			tauZ_25 = intparab(ztau[iz], zkm_mod[25], zkm_mod[24], zkm_mod[25], zkm_mod[26],
			                                          ext_km[24], ext_km[25], ext_km[26]);
		else if (ix == 23)      // bin=[23:24]km
			tauZ_25 = intparab(ztau[iz], zkm_mod[25], zkm_mod[23], zkm_mod[24], zkm_mod[25],
			                                          ext_km[23], ext_km[24], ext_km[25]);
		else
		{
			jx = ix % 2;        // jx = 0 for even ix
			tauZ_25 = intparab(ztau[iz], zkm_mod[ix+jx+1], zkm_mod[ix], zkm_mod[ix+1], zkm_mod[ix+2],
			                                               ext_km[ix], ext_km[ix+1], ext_km[ix+2]) +
                      simpson(&ext_km[ix+jx+1], 25-(ix+jx), 1.0);
		}
		tau[iz] = tauZ_25 + tau25_120;
	} // for iz=0:ntau
} // tauabs25
/*--------------------------------------------------------------------------------------------------
2020-05-15:
    minor changes in comments;
2020-02-25:
    replaced: intexp(&ext_km[35], &zkm_mod[35], 15) -> simpson(&ext_km[35], 15, 5.0) due to 'intexp'
	does not work for ext_km = const = 0.0 (no absorption)
2020-02-09:
	First created. Selection of the correct bin & integration interval tested for: ntau = 10,
	ztau[ntau] = {25.0, 24.9, 24.1, 23.5, 22.5, 21.5, 13.3, 2.9, 0.5, 0.0} - ok (manual test).
--------------------------------------------------------------------------------------------------*/