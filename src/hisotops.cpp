#include <cstring>       /* strcpy, strcat */
#include <stdio.h>       /* FILE, printf, fopen, fscanf, fclose */
#include <math.h>        /* floor */
#include "paths.h"
#include "const_param.h"
/*
prototype:
	int hisotops(int const, double const *, int const,
	             double *, double *, double *);
*/
int hisotops(int const molec_id, double const *Tz, int const nz,                              //  in
	         double *Qratio, double *mmass_iso, double *Ia_iso) {                             // out
/*--------------------------------------------------------------------------------------------------
TASK:
	To compute the height-dependent total internal partition sum (TIPS) ratio, Q(296K)/Q(T), for a
	given temperature profile, T(z[km]), and return for each isotop the number of isotops, molar
	mass, and abundance.
IN:
	molec_id   i       as in HITRAN
	Tz         d[nz]   temperature, [K]
	nz         i       number of heights, nz = len(T_kelvin)
OUT:
	niso        i                number of isotops
	Qratio      d[nz*niso_max]   Q(296K)/Q(T): *** WARNING: Tref/Tuser - not vice versa! ***
	                             nz is the lead dimension;
	mmass_iso   d[niso_max]      molar mass of isotops [g/mol]
	Ia_iso      d[niso_max]      natural terrestrial isotopic abundances
NOTE:
    In terrestrial atmosphere, Ia_iso is already included in hitran's Sij.
    For output arrays,  only first 'niso' values are used.
	Q(T) is linear intepolated, LUT step over temperature = 1 Kelvin.
	Q(296K) comes from [1], not from [2] (less digits).
	See [3,4] for the definition of TIPS, and [5,6] for a Fortran code, TIPS.for.
REFS:
	1. https://hitran.org/docs/iso-meta/
	2. https://hitran.org/media/molparam.txt
	3. https://hitran.org/docs/definitions-and-units/
	4. Gamache RR, et al., 2017: JQSRT 203, 70-87 (doi: 10.1016/j.jqsrt.2017.03.045)
	5. https://hitran.org/supplementary/
	6. https://hitran.org/suppl/TIPS/
--------------------------------------------------------------------------------------------------*/
	FILE *fin;
	int const
		nlines = 240; // Tq = 150 + iline*dTkelv = [150:1:389] [K]
	double const
		Tq_min = 150.0, dTkelv = 1.0;
	char
		fpath[path_len_max], fname_iso[niso_max][fname_len_max];
	int
		ibin, iline, iso, iz, niso;
	double
		Qo, To;
	double
		Qref[niso_max], Qt[nlines], Tq[nlines];
//--------------------------------------------------------------------------------------------------
//
	niso = -999;
//
	switch (molec_id)
	{
		case 1: // h2o
			niso = niso_h2o;
			for (iso = 0; iso < niso; iso++)
			{
				Qref[iso] = Qref_h2o[iso];
				mmass_iso[iso] = molar_mass_h2o[iso];
				Ia_iso[iso] = Ia_iso_h2o[iso];
				strcpy(fname_iso[iso], fname_iso_h2o[iso]);
			} // for iso
			break;
		case 2: // co2
			niso = niso_co2;
			for (iso = 0; iso < niso; iso++)
			{
				Qref[iso] = Qref_co2[iso];
				mmass_iso[iso] = molar_mass_co2[iso];
				Ia_iso[iso] = Ia_iso_co2[iso];
				strcpy(fname_iso[iso], fname_iso_co2[iso]);
			} // for iso
			break;
		case 3: // o3
			niso = niso_o3;
			for (iso = 0; iso < niso; iso++)
			{
				Qref[iso] = Qref_o3[iso];
				mmass_iso[iso] = molar_mass_o3[iso];
				Ia_iso[iso] = Ia_iso_o3[iso];
				strcpy(fname_iso[iso], fname_iso_o3[iso]);
			} // for iso
			break;
		case 4: // n2o
			niso = niso_n2o;
			for (iso = 0; iso < niso; iso++)
			{
				Qref[iso] = Qref_n2o[iso];
				mmass_iso[iso] = molar_mass_n2o[iso];
				Ia_iso[iso] = Ia_iso_n2o[iso];
				strcpy(fname_iso[iso], fname_iso_n2o[iso]);
			} // for iso
			break;
		case 5: // co
			niso = niso_co;
			for (iso = 0; iso < niso; iso++)
			{
				Qref[iso] = Qref_co[iso];
				mmass_iso[iso] = molar_mass_co[iso];
				Ia_iso[iso] = Ia_iso_co[iso];
				strcpy(fname_iso[iso], fname_iso_co[iso]);
			} // for iso
			break;
		case 6: // ch4
			niso = niso_ch4;
			for (iso = 0; iso < niso; iso++)
			{
				Qref[iso] = Qref_ch4[iso];
				mmass_iso[iso] = molar_mass_ch4[iso];
				Ia_iso[iso] = Ia_iso_ch4[iso];
				strcpy(fname_iso[iso], fname_iso_ch4[iso]);
			} // for iso
			break;
		case 7: // o2
			niso = niso_o2;
			for (iso = 0; iso < niso; iso++)
			{
				Qref[iso] = Qref_o2[iso];
				mmass_iso[iso] = molar_mass_o2[iso];
				Ia_iso[iso] = Ia_iso_o2[iso];
				strcpy(fname_iso[iso], fname_iso_o2[iso]);
			} // for iso
			break;
		case 10: // no2
			niso = niso_no2;
			for (iso = 0; iso < niso; iso++)
			{
				Qref[iso] = Qref_no2[iso];
				mmass_iso[iso] = molar_mass_no2[iso];
				Ia_iso[iso] = Ia_iso_no2[iso];
				strcpy(fname_iso[iso], fname_iso_no2[iso]);
			} // for iso
			break;
		default:
			printf("\n(in hisotops.cpp) default is undefined");
	} // switch (molec_id)
//
	for (iso = 0; iso < niso; iso++)
	{
		strcpy(fpath, path_TIPS);
		strcat(fpath, fname_iso[iso]);
		fin = fopen(fpath, "r");
		fscanf(fin, "%lf %lf", &To, &Qo);
		while (To < Tq_min) fscanf(fin, "%lf %lf", &To, &Qo);
		Tq[0] = To;
		Qt[0] = Qo;
		for (iline = 1; iline < nlines; iline++) fscanf(fin, "%lf %lf", &Tq[iline], &Qt[iline]);
		for (iz = 0; iz < nz; iz++)
		{
			ibin = int(floor((Tz[iz] - To)/dTkelv)); // attn: floor returns double
			Qo = Qt[ibin] + (Tz[iz] - Tq[ibin])*(Qt[ibin+1] - Qt[ibin])/(Tq[ibin+1] - Tq[ibin]);
			Qratio[iso*nz+iz] = Qref[iso]/Qo;
		} // for iz = 0:nz
		fclose(fin);
		printf("\n(in isotops.cpp) isotop file processed: %s", fpath);
	} // for iso
//
	return(niso);
} // hisotops
/*--------------------------------------------------------------------------------------------------
2020-05-15:
    minor changes in comments;
2020-05-08:
    ibin = int(floor((Tz[iz] - To)/dTkelv)); note: int()
2020-02-12:
	H2O added, tested with main program with modtran profiles
2020-02-05:
	First created based on isotops.cpp, compiled and visually tested for one point:
    molec_id = 6, T=288.2K, isotop #1 and for several temperatures for the last isotop from the
	calling subroutine (it also makes sure the output array is formed properly):
	T = 288.20, Qratio = 1.042462e+00 vs python = 1.0424620540539196
	T = 300.00, Qratio = 9.792296e-01 vs exact LUT values:
	    Qref(9599.15882000)/Q300K(9802.76635000) = 0.97922958451417...
*/