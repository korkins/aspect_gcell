#include <cstring>       /* strcpy, strcat */
#include <stdio.h>       /* FILE */
#include <stdlib.h>      /* atof */
#include "const_param.h" /* char path_hitdb[], char fname_hitdb[], int path_len_max */
/*
protitype:
	void read_hitran100(int, int, int,
							int *, double *, double *, double *,
								double *, double *, double *, double *);
*/
void read_hitran100(int const molec_id, int const iline0, int const nlines,                   //  in
	                int *isotop_id, double *nuij, double *Sij, double *gamma_air,             // out
	                double *gamma_self, double *Elower, double *n_air, double *delta_air) {
/*--------------------------------------------------------------------------------------------------
TASK:
	To read 'nlines' starting from 'iline0' from a hitran 100-ascii file for a molecule 'molec_id'.
IN:
	molec_id   i   as defiend in HITRAN
	iline0     i   starting line, iline0 = 0 (top line), 1, 2, 3 ...
	nlines     i   number of lines to read (including iline0); iline0+nlines <= nlines_total
OUT:
	isotop_id    i[nlines]   [X]                  isotop_id to get e.g. proper molar mass, gp, etc.
	nuij         d[nlines]   [cm-1]:              Transition wavenumber; '+1' to append '\0' = end of string
	Sij          d[nlines]   [cm-1/(molec.cm-2)]: Line intensity, scaled by isotopologue abundance, at T = 296K
	gamma_air    d[nlines]   [cm-1.atm-1]:        Air-broadened Lorentzian HWHM at p = 1 atm and T = 296K
	gamma_self   d[nlines]   [cm-1.atm-1]:        Self-broadened HWHM at 1 atm pressure and 296K
    Elower       d[nlines]   [cm-1]:              Lower-state energy = Epp = Ei
	n_air        d[nlines]   [X]:                 Temperature exponent for the air-broadened HWHM
	delta_air    d[nlines]   [cm-1.atm-1]:        Pressure shift induced by air, referred to p=1 atm
NOTE:
    HWHM = 'half-width at half-maximum'
	i) The subroutine ASSUMES the HITRAN database file exisits and is not empty.
	   iline=0 corresponds to the top line in the HITRAN database file.
	   Ref.[1] gives a conveneint way to read a fixed-format ASCII file.
	ii) It is recommended [1] to check fscanf for output, e.g.
		    while (scanf(" %10c %20c", first, second) == 2)
	    However, this subroutines assumes the par-file is correct.
	iii) The spaces in " %10c %20c" correspond to space in data.
	iv) This subroutine does NOT check for iline0+nlines <= nlines_total and
	    does not check for end-of-file.
REFS:
	1. https://stackoverflow.com/questions/34759398/read-only-a-fixed-column-width-in-a-text-file-in-c (2019-12-26)
--------------------------------------------------------------------------------------------------*/
	FILE *fin;
	char
		fpath[path_len_max],
//      hitran 100 ascii par file: '+1' to append '\0' = end of string
		str_full[100],
	    str_molec_id[2],
		str_isotop_id[1+1],
		str_nuij[12+1],
		str_Sij[10+1],
		str_R[10],
		str_gamma_air[5+1],
		str_gamma_self[5+1],
        str_Elower[10+1],
		str_n_air[4+1],
		str_delta_air[8+1],
		str_tail[33],        // skip unused parameters
		chr_end_of_line[1];
	int
		iline;
//--------------------------------------------------------------------------------------------------
//
	strcpy(fpath, path_hitdb);
	strcat(fpath, fname_hitdb[molec_id-1]);
	fin = fopen(fpath, "r");
	printf("\n(in read_hitran160.cpp) opened: %s", fpath);
//
	for (iline = 0; iline < iline0; iline++)
		fscanf(fin, "%100c%c", str_full, chr_end_of_line);
//
	for (iline = 0; iline < nlines; iline++) {
		fscanf(fin, "%2c%1c%12c%10c%10c%5c%5c%10c%4c%8c%33c%c", // see (iii) in Comments
					str_molec_id, str_isotop_id,
						str_nuij, str_Sij, str_R,
							str_gamma_air, str_gamma_self, str_Elower,
								str_n_air, str_delta_air, str_tail,
									chr_end_of_line);
//
		str_isotop_id[1] = '\0';
		str_nuij[12] = '\0';
		str_Sij[10] = '\0';
		str_gamma_air[5] = '\0';
		str_gamma_self[5] = '\0';
        str_Elower[10] = '\0';
		str_n_air[4] = '\0';
		str_delta_air[8] = '\0';
//
		isotop_id[iline] = atoi(str_isotop_id);
		nuij[iline] = atof(str_nuij);
		Sij[iline] = atof(str_Sij);
		gamma_air[iline] = atof(str_gamma_air);
		gamma_self[iline] = atof(str_gamma_self);
		Elower[iline] = atof(str_Elower);
		n_air[iline] = atof(str_n_air);
		delta_air[iline] = atof(str_delta_air);
	} // for iline
//
	fclose(fin);
	printf("\n(in read_hitran100.cpp) closed: %s", fpath);
} // read_hitran100
//--------------------------------------------------------------------------------------------------
/*
2020-05-05: created from read_hitran160.cpp and tested as part of the whole code for tau_atm(AO2).
*/