#include <stdio.h>        /* printf */
#include <time.h>         /* clock, CLOCKS_PER_SEC*/
#include <math.h>         /* sqrt, pow, exp, fabs */
#include "const_param.h"  /* int niso_max, double physical constants (e.g. k_noltzman),
//                           double math constants (e.g. sqrt_ln2) */
//
int hisotops(int const, double const *, int const, double *, double *, double *);
int ix1ix2(double const, double const, double const *, int const, int &, int &);
/*
	humlicek(x, y): .\src_c\alyapustin_IPC\ (+dependents) 
    So far, i dont know the details of humlicek(x, y)
*/
double humlicek(double, double);
/*
prototype:
    void kabs(int const,
	          double const *, int const,
	          int const, double const *, double const *, double const *,
	          int const, int *, double *, double *, double *,
	                            double *, double *, double *, double *,
	          double *);
*/
void kabs(int const molec_id,                                                                 // in
	      double const *nu, int const nnu,
	      int const nz, double const *Tkelv, double const *Patm, double const *Pgas,
	      int const nlines, int *isotop_id, double *nuij, double *Sij, double *gamma_air,
	                        double *gamma_self, double *Epp, double *n_air, double *delta_air,
	      double *knu) {                                                                      // out
/*--------------------------------------------------------------------------------------------------
TASK:
	To compute spectral absorption x-section per molecule using HITRAN data as a function of height.
IN:
	molec_id   i   as defiend in HITRAN
	nu         d[nu]   user wavenumber grid
	nnu        i       len(nu)
	nz         i       number of grid heights
	Tkelv      d[nz]   tempertuare profile, [K]
	Patm       d[nz]   total pressure profile, [atm]
	Pgas       d[nz]   partial pressure profile, [atm]
	nlines     i   number of lines to read (including iline0); iline0+nlines <= nlines_total
	* REQIRED HITRAN PARAMETERS: *
	isotop_id    i[nlines]   [X]                   isotop_id to get e.g. proper molar mass, gp, etc.
	nuij         d[nlines]   [cm-1]:               Transition wavenumber; '+1' to append '\0' = end of string
	Sij          d[nlines]   [cm-1/(molec.cm-2)]   Line intensity, scaled by isotopologue abundance, at T = 296K
	gamma_air    d[nlines]   [cm-1.atm-1]:         Air-broadened Lorentzian HWHM at p = 1 atm and T = 296K
	gamma_self   d[nlines]   [cm-1.atm-1]:         Self-broadened HWHM at 1 atm pressure and 296K
    Elower       d[nlines]   [cm-1]:               Lower-state energy = Epp = Ei
	n_air        d[nlines]   [X]:                  Temperature exponent for the air-broadened HWHM
	delta_air    d[nlines]   [cm-1.atm-1]:         Pressure shift induced by air, referred to p=1 atm
OUT:
    knu        d[nz*nnu]     [cm3/(cm.molec)]     Spectral absorption coeffiecient: x-section per molecule
NOTE:
	In knu[nz*nnu], 'nz'is the fast dimentsion.
	On input, knu[nz*nnu] is initilaized with 0.0.
REFS:
	1. -
--------------------------------------------------------------------------------------------------*/
	int
		iline, inu, inu1, inu2, iso_ix, iz, niso;
	double
		e11, e21, e12, e22, gam_lorentz, x;
	double
		time, time1, time2;
	double
		mmass_iso[niso_max], Ia_iso[niso_max];
	double
		*nuij_pshift, *Qratio, *alf_doppler, *SijT, *y;
//--------------------------------------------------------------------------------------------------
//
	time = 0.0;
//
	nuij_pshift = new double [nz];
	alf_doppler = new double [nz];
	SijT = new double [nz];
	Qratio = new double [niso_max*nz];
	y = new double [nz];
	niso = hisotops(molec_id, Tkelv, nz, Qratio, mmass_iso, Ia_iso); // Earth: Sij includes Ia_iso
	for (iline = 0; iline < nlines; iline++)
	{
		if (ix1ix2(nuij[iline], delta_nu, nu, nnu, inu1, inu2) > 0)
		{
			iso_ix = isotop_id[iline]-1;
			for (iz = 0; iz < nz; iz++) // precompute z-dependent parameters
			{
				nuij_pshift[iz] = nuij[iline] + delta_air[iline]*Patm[iz];
//              THINK ME:
				alf_doppler[iz] =
					nuij_pshift[iz]*sqrt(2.0*n_avogadro*k_boltzman*Tkelv[iz]*ln2/mmass_iso[iso_ix])/c_light;
				gam_lorentz =
					pow(T_ref/Tkelv[iz], n_air[iline])*
						(gamma_air[iline]*(Patm[iz] - Pgas[iz]) + gamma_self[iline]*Pgas[iz]);
				y[iz] = sqrt_ln2*gam_lorentz/alf_doppler[iz];
//				THINKME: combine e11 & e21 and save one exp call
				e11 = exp(-c2_rad*Epp[iline]/Tkelv[iz]);
				e21 = exp(-c2_rad*Epp[iline]/T_ref);
				e12 = 1.0 - exp(-c2_rad*nuij_pshift[iz]/Tkelv[iz]);
				e22 = 1.0 - exp(-c2_rad*nuij_pshift[iz]/T_ref);
//				THINK ME: Alexei if-checks for 1-exp(-big_number) = 1
				SijT[iz] = Sij[iline]*Qratio[iso_ix*nz+iz]*e11*e12/e21/e22; // in hitran, Sij accounts for terrestrial Ia
			} // for iz = 0:nz
			time1 = (double)clock() /(double)CLOCKS_PER_SEC;
			for (inu = inu1; inu < inu2+1; inu++) // inu2+1 to include inu2
			{
				for (iz = 0; iz < nz; iz++)
				{
//                  THINK ME: save operations, e.g. sqrt(ln2/sqrt_pi)
					x = sqrt_ln2*fabs(nu[inu] - nuij_pshift[iz])/alf_doppler[iz];
					knu[inu*nz+iz] += SijT[iz]*(sqrt_ln2/sqrt_pi)*humlicek(x, y[iz])/alf_doppler[iz]; // accumulate from all lines
				} // for iz = 0:nz
			} // for inu = inu1:inu2
			time2 = (double)clock() /(double)CLOCKS_PER_SEC;
			time += time2 - time1;
		} // if iline contributes to 'nu'
		if (!(iline%1000)) printf("\n-processed: iline=%d", iline);
	} // for iline
	printf("\n(in kabs.cpp): humlichek %6.2fs", time);
	delete[] nuij_pshift;
	delete[] alf_doppler;
	delete[] SijT;
	delete[] Qratio;
	delete[] y;
} // kabs
//--------------------------------------------------------------------------------------------------
/*
2020-05-15:
    minor changes in comments;
2020-05-12:
    Fixed: SijT[iz] = Ia_iso[iso_ix]*Sij[iline]*... -> SijT[iz] = Sij[iline]*...; Sij includes Ia 
2020-05-08:
    Minor corrections in comments;
2020-02-06:
    First created, tested as part of tau_abs_hprofile.
*/