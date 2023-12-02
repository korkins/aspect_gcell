#include <cstring>       /* strcpy, strcat */
#include <stdio.h>       /* FILE, printf, fopen, fscanf, fclose */
#include "paths.h"
#include "const_param.h"
//
int isotops(int const molec_id, double const T_kelvin,
	        double *Qratio, double *mmass_iso, double *Ia_iso) {
/*--------------------------------------------------------------------------------------------------
TASK:
	To compute the total internal partition sum (TIPS) ratio, Q(296K)/Q(T), for a given temperature
	T(K), and return for each isotope the number of isotopes, molar mass, and abundance.
IN:
	molec_id   i   as in HITRAN
	T_kelvin   d   temperature (K)
OUT:
	niso        i             number of isotopes
	Qratio      d[niso_max]   Q(296K)/Q(T): *** WARNING: reference/arbitrary - not vice versa! ***
	mmass_iso   d[niso_max]   molar mass of isotopes (g/mol)
	Ia_iso      d[niso_max]   natural terrestrial isotopic abundances
NOTE:
    For output arrays,  only first 'niso' values are used: niso <= niso_max.
	Q(T) is linear interpolated.
	Q(296K) comes from [1], not from [2] (the latter indicates less digits).
	See [3,4] for the definition of TIPS, and [5,6] for a Fortran code, TIPS.for.
REFS:
	1. https://hitran.org/docs/iso-meta/
	2. https://hitran.org/media/molparam.txt
	3. https://hitran.org/docs/definitions-and-units/
	4. Gamache RR, et al., 2017: JQSRT 203, 70-87 (doi: 10.1016/j.jqsrt.2017.03.045)
	5. https://hitran.org/supplementary/
	6. https://hitran.org/suppl/TIPS/
--------------------------------------------------------------------------------------------------*/
	FILE
		*fin;
	char
		fpath[path_len_max], fname_iso[niso_max][fname_len_max];
	int
		iso, niso;
	double
		Q, Q1, Q2, T1, T2;
	double
		Qref[niso_max];
//--------------------------------------------------------------------------------------------------
//
	niso = -999;
	T1 = 0.0;
	Q1 = 0.0;
//
	switch (molec_id)
	{
		case 1: // H2O
			niso = niso_h2o;
			for (iso = 0; iso < niso; iso++)
			{
				Qref[iso] = Qref_h2o[iso];
				mmass_iso[iso] = molar_mass_h2o[iso];
				Ia_iso[iso] = Ia_iso_h2o[iso];
				strcpy(fname_iso[iso], fname_iso_h2o[iso]);
			} // for iso
			break;
		case 6: // CH4
			niso = niso_ch4;
			for (iso = 0; iso < niso; iso++)
			{
				Qref[iso] = Qref_ch4[iso];
				mmass_iso[iso] = molar_mass_ch4[iso];
				Ia_iso[iso] = Ia_iso_ch4[iso];
				strcpy(fname_iso[iso], fname_iso_ch4[iso]);
			} // for iso
			break;
		case 7: // O2
			niso = niso_o2;
			for (iso = 0; iso < niso; iso++)
			{
				Qref[iso] = Qref_o2[iso];
				mmass_iso[iso] = molar_mass_o2[iso];
				Ia_iso[iso] = Ia_iso_o2[iso];
				strcpy(fname_iso[iso], fname_iso_o2[iso]);
			} // for iso
			break;
		case 10: // NO2
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
			printf("\n(in isotops.cpp) default is undefined");
	} // switch (molec_id)
//
	for (iso = 0; iso < niso; iso++)
	{
		strcpy(fpath, path_TIPS);
		strcat(fpath, fname_iso[iso]);
		fin = fopen(fpath, "r");
		fscanf(fin, "%lf %lf", &T2, &Q2);
		while (T2 < T_kelvin)
		{
//          THINKME: potentially infinite loop if T_kelvin > T_grid_max
//          THINKME: skip T from 1 to 190 as in Alexei's code? Cut lut from top
			T1 = T2;
			Q1 = Q2;
			fscanf(fin, "%lf %lf", &T2, &Q2);
		} // while T2 < T_kelvin
		Q = Q1 + (T_kelvin - T1)*(Q2 - Q1)/(T2 - T1);
		Qratio[iso] = Qref[iso]/Q;
		fclose(fin);
		printf("\n(in isotops.cpp) isotop file processed: %s", fpath);
	} // for iso
//
	return(niso);
} // isotops
/*--------------------------------------------------------------------------------------------------
2022-05-28:
    again, minor changes in comments;
2020-05-15:
    minor changes in comments;
2020-04-15:
	Replaced: case 6: CH4, niso = 4 -> niso = niso_ch4;
	Added: case 7: O2
2020-01-09:
	Renamed: molar_mass_iso -> mmass_iso
	Added:  Ia_iso
	Minor:  printf removed/modified
	Tested: compiled - ok.
2020-01-07:
    First created and tested for molec_id=6(ch4), T=Tref=296K -> Qratio[:]=1.0
	T=Tref+0.5K -> Q=0.5*(Q1+Q2) - ok
	niso(out) = 4 - ok
	mmass_iso = {16.0313, 17.034655, 17.037475, 18.04083} - ok
*/