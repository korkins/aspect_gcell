#include <cstring>       /* strcpy, strcat */
#include <stdio.h>       /* FILE */
#include <stdlib.h>      /* atof */
#include "paths.h"
#include "const_param.h"
//
void read_hitran160(int const molec_id, int const iline0, int const nlines,                   //  in
	                int *isotop_id, double *nuij, double *Sij, double *gamma_air,             // out
	                double *gamma_self, double *Elower, double *n_air, double *delta_air) {
/*--------------------------------------------------------------------------------------------------
TASK:
	To read 'nlines' starting from 'iline0' from a hitran 160-ascii file for a molecule 'molec_id'.
IN:
	molec_id   i   as defiend in HITRAN
	iline0     i   starting line, iline0 = 0 (top line), 1, 2, 3 ...
	nlines     i   number of lines to read (including iline0); iline0+nlines <= nlines_total
OUT:
	isotop_id    i[nlines]   [X]                  isotop_id to get e.g. proper molar mass, gp, etc.
	nuij         d[nlines]   [cm-1]:              transition wavenumber in vacuum
	Sij          d[nlines]   [cm-1/(molec.cm-2)]: line intensity, *scaled* by isotope abundance Ia, at T = 296K
	gamma_air    d[nlines]   [cm-1.atm-1]:        air-broadened Lorentzian HWHM at p = 1 atm and T = 296K
	gamma_self   d[nlines]   [cm-1.atm-1]:        self-broadened HWHM at 1 atm pressure and 296K
    Elower       d[nlines]   [cm-1]:              lower-state energy = Epp = Ei
	n_air        d[nlines]   [X]:                 temperature exponent for the air-broadened HWHM
	delta_air    d[nlines]   [cm-1.atm-1]:        pressure shift induced by air, referred to p=1 atm
NOTE:
    HITRAN line intensity, Sij, already accounts for terrestrial abundance.
    HWHM = 'half-width at half-maximum'
	i) The subroutine ASSUMES the HITRAN database file exisits and is not empty.
	   iline=0 corresponds to the top line in the HITRAN database file.
	   Ref.[1] gives a convenient way to read a fixed-format ASCII file.
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
//      hitran 160 ascii par file: '+1' to append '\0' = end of string
		str_full[160],
	    str_molec_id[2],
		str_isotop_id[1+1],
		str_nuij[12+1],
		str_Sij[10+1],
		str_Aij[10],
		str_gamma_air[5+1],
		str_gamma_self[5+1],
        str_Elower[10+1],
		str_n_air[4+1],
		str_delta_air[8+1],
		str_tail[93],        // skip unused parameters
		chr_end_of_line[1];
	int
		iline;
/*------------------------------------------------------------------------------------------------*/
//
	strcpy(fpath, path_hitdb);
	strcat(fpath, fname_hitdb[molec_id-1]);
	fin = fopen(fpath, "r");
	printf("\n(in read_hitran160.cpp) opened: %s", fpath);
//
	for (iline = 0; iline < iline0; iline++)
		fscanf(fin, "%160c%c", str_full, chr_end_of_line);
//
	for (iline = 0; iline < nlines; iline++)
	{
		fscanf(fin, "%2c%1c%12c%10c%10c%5c%5c%10c%4c%8c%93c%c", // see (iii) in Comments
					str_molec_id, str_isotop_id,
						str_nuij, str_Sij, str_Aij,
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
		if (molec_id == 2) // CO2 has 12 isotops. By definition: ID = 0, A, B -> isotop_id = 10, 11, 12
			if (str_isotop_id[0] == '0')
				isotop_id[iline] = 10;
			else if (str_isotop_id[0] == 'A')
				isotop_id[iline] = 11;
			else if (str_isotop_id[0] == 'B')
				isotop_id[iline] = 12;
			else
				isotop_id[iline] = atoi(str_isotop_id);
		else // all other molecules
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
	printf("\n(in read_hitran160.cpp) closed: %s", fpath);
} // read_hitran160
/*------------------------------------------------------------------------------------------------*/
/*
2022-05-27:
    Cosmetic changes in comments.
2021-04-20:
	CO2 has 12 isotops: ID=0, A, B -> isotop_id = 10, 11, 12. New piece of code:
	>>
		if (molec_id == 2) // CO2 has 12 isotops. By definition: ID = 0, A, B -> isotop_id = 10, 11, 12
			if (str_isotop_id[0] == '0')
				isotop_id[iline] = 10;
			else if (str_isotop_id[0] == 'A')
				isotop_id[iline] = 11;
			else if (str_isotop_id[0] == 'B')
				isotop_id[iline] = 12;
			else
				isotop_id[iline] = atoi(str_isotop_id);
		else // all other molecules
			isotop_id[iline] = atoi(str_isotop_id);
	<<
2020-05-15:
    minor changes in comments;
2020-01-13:
    Aij, gp, gpp are removed from output (not used); code modified: str_skip -> str_tail...
			Checked vs original version.
2019-12-26: created and tested using simplified 06_hitran.par, nlines_total = 37 (ilines = 0:36)

Case A: test 'fscanf' and 'atof' (first 3 lines), note nuij in the 2nd line is changed
 61 4000.000000 6.769E-27 1.237E-04.04930.065 1095.61940.62-.007545    1 0 0 1 1F2    0 0 0 0 1A1   13F1 62        14F2  2     233333332329 1 1 1    81.0   87.0
 92-.123456E-01 6.769E-27 1.237E-04.04930.065 1095.61940.62-.007545    1 0 0 1 1F2    0 0 0 0 1A1   13F1 62        14F2  2     233333332329 1 1 1    81.0   87.0
 61 4001.000000 3.759E-27 3.629E-03.04850.059 1976.57460.60-.007545    1 0 0 1 1F2    0 0 0 0 1A1   18F2 98        19F1  4     233333332329 1 1 1   111.0  117.0

Case B: result is visualy checked in a reduced par file (see below)
	 B1: iline = 0, nlines = 13 to test how the first line is read - OK.
	 B2: iline = 13, nlines = 21 to test a genral case: start/stop somwhere in the file - OK.
	 B3: iline = 1, nlines = 1 to test a single line reading - OK.
	 B4: iline = 34, nlines = 3 to test how the last line in the file is read - OK.
	 B5: iline = 9, nlines = 37 to test reading the whole file - OK.

Case C:  tested using real 06_hitdb.par:
     C1: iline0 = 194430, nlines = 6 (end of file including last line) - OK. Runtime ~2s. (debug)
	 C2: iline=0, nlines = 194436 (all) - OK (1st & last lines printed out). Runtime ~4s (debug). RAM ~ 18M per 194436 lines

Test.par:
 61 4000.000000 6.769E-27 1.237E-04.04930.065 1095.61940.62-.007545    1 0 0 1 1F2    0 0 0 0 1A1   13F1 62        14F2  2     233333332329 1 1 1    81.0   87.0
 92-.123456E-01 6.769E-27 1.237E-04.04930.065 1095.61940.62-.007545    1 0 0 1 1F2    0 0 0 0 1A1   13F1 62        14F2  2     233333332329 1 1 1    81.0   87.0
 61 4001.000000 3.759E-27 3.629E-03.04850.059 1976.57460.60-.007545    1 0 0 1 1F2    0 0 0 0 1A1   18F2 98        19F1  4     233333332329 1 1 1   111.0  117.0
 61 4002.000000 2.538E-24 2.284E-02.05250.067  949.84310.63-.007546    0 1 0 2 1F2    0 0 0 0 1A1   13F1 36        13F2  1     363333332329 1 1 1    81.0   81.0
 61 4003.000000 7.612E-26 1.281E-02.04780.061 1593.56200.61-.007546    0 1 0 2 1F1    0 0 0 0 1A1   16F1 79        17F2  2     233333332329 1 1 1    99.0  105.0
 61 4004.000000 1.509E-25 7.214E-03.04780.063 1251.30790.62-.007546    0 0 0 3 1F2    0 0 0 0 1A1   16E  11        15E   1     233333332329 1 1 1    66.0   62.0
 61 4005.000000 7.437E-26 6.118E-04.05230.068  814.99330.64-.007546    0 1 0 2 1F1    0 0 0 0 1A1   11E  32        12E   2     233333332329 1 1 1    46.0   50.0
 61 4006.000000 1.380E-25 4.509E-04.05370.070  689.70540.64-.007546    0 1 0 2 1F1    0 0 0 0 1A1   10F2 42        11F1  1     233333332329 1 1 1    63.0   69.0
 61 4007.000000 3.130E-27 4.973E-04.04780.061 1593.87850.61-.007546    0 0 0 3 1F1    0 0 0 0 1A1   17F1 46        17F2  3     233333332329 1 1 1   105.0  105.0
 61 4008.000000 7.502E-26 3.224E-03.04780.063 1417.86520.62-.007546    0 0 0 3 1F1    0 0 0 0 1A1   16A2 15        16A1  2     233333332329 1 1 1   165.0  165.0
 61 4009.000000 1.296E-26 1.932E-02.04700.058 2181.85130.60-.007546    1 0 0 1 1F2    0 0 0 0 1A1   19A1 37        20A2  1     133333332329 1 1 1   195.0  205.0
 61 4010.000000 5.068E-27 7.362E-03.04850.059 1977.19780.60-.007546    1 0 0 1 1F2    0 0 0 0 1A1   18E  65        19E   3     233333332329 1 1 1    74.0   78.0
 61 4011.000000 7.900E-27 5.622E-04.04780.063 1416.55430.62-.007546    0 0 0 3 1F1    0 0 0 0 1A1   16F2 43        16F1  1     233333332329 1 1 1    99.0   99.0
 61 4012.000000 4.573E-27 8.359E-05.04930.065 1095.63210.62-.007546    1 0 0 1 1F2    0 0 0 0 1A1   13F2 61        14F1  1     333333332329 1 1 1    81.0   87.0
 61 4013.000000 1.491E-24 1.342E-02.05250.067  949.84170.63-.007546    0 1 0 2 1F2    0 0 0 0 1A1   13F2 35        13F1  1     233333332329 1 1 1    81.0   81.0
 61 4014.000000 1.072E-27 2.663E-03.04700.058 2181.82130.60-.007546    0 0 1 1 1F1    0 0 0 0 1A1   19F1113        20F2  2     133333332329 1 1 1   117.0  123.0
 61 4015.000000 3.011E-25 3.396E-04.05500.073  470.86510.65-.007546    0 0 0 3 1F2    0 0 0 0 1A1   10F1 15         9F2  2     333333332329 1 1 1    63.0   57.0
 63 4016.000000 3.131E-25 1.191E+00.05270.050  621.46990.65-.008000       V4+V6 E          GROUND 10  8  A2      11  9  A1     454432554528 7 1 7   504.0  552.0
 61 4017.000000 9.889E-28 3.685E-03.04700.058 2181.80790.60-.007546    0 0 1 1 1F1    0 0 0 0 1A1   19E  74        20E   2     133333332329 1 1 1    78.0   82.0
 61 4018.000000 4.374E-27 4.225E-03.04850.059 1976.65080.60-.007546    0 0 1 1 1F2    0 0 0 0 1A1   18F1 95        19F2  4     233333332329 1 1 1   111.0  117.0
 61 4019.000000 5.641E-27 8.981E-05.04880.065 1095.61940.62-.007546    0 0 0 3 1F2    0 0 0 0 1A1   15F1 18        14F2  2     233333332329 1 1 1    93.0   87.0
 61 4020.000000 1.898E-25 1.116E-01.04900.060 1779.65140.61-.007546    1 0 0 1 1F2    0 0 0 0 1A1   17E  58        18E   2     133333332329 1 1 1    70.0   74.0
 61 4021.000000 3.875E-24 2.266E-02.05250.067  950.38530.63-.007546    0 1 0 2 1F2    0 0 0 0 1A1   12A1 19        13A2  1     343333332329 1 1 1   125.0  135.0
 61 4022.000000 2.461E-27 1.668E-02.04600.058 2398.56690.59-.007546    0 0 0 3 2F2    0 0 0 0 1A1   20F1123        21F2  3     233333332329 1 1 1   123.0  129.0
 61 4023.000000 1.740E-24 6.322E-04.06480.078  104.77470.72-.007546    0 0 0 3 1F1    0 0 0 0 1A1    5F2 10         4F1  1     363333332329 1 1 1    33.0   27.0
 61 4024.000000 2.501E-27 4.572E-05.04930.065 1095.61940.62-.007546    1 0 0 1 1F2    0 0 0 0 1A1   13F1 63        14F2  2     333333332329 1 1 1    81.0   87.0
 61 4025.000000 3.705E-26 3.572E-02.04850.059 1976.24230.60-.007546    1 0 0 1 1F2    0 0 0 0 1A1   18F2 98        19F1  3     233333332329 1 1 1   111.0  117.0
 61 4026.000000 2.005E-26 7.835E-03.04900.060 1779.03060.61-.007546    1 0 0 1 1F2    0 0 0 0 1A1   17F2 87        18F1  1     133333332329 1 1 1   105.0  111.0
 61 4027.000000 4.471E-27 3.004E-02.04600.058 2396.79420.59-.007546    1 0 0 1 1F2    0 0 0 0 1A1   20F2125        21F1  2     133333332329 1 1 1   123.0  129.0
 61 4028.000000 2.942E-25 1.344E-02.04780.063 1417.57930.62-.007546    0 1 0 2 1F1    0 0 0 0 1A1   15A1 24        16A2  1     233333332329 1 1 1   155.0  165.0
 61 4029.000000 3.624E-25 2.925E-04.06340.077  219.91350.72-.007546    0 0 0 3 2F2    0 0 0 0 1A1    6E  10         6E   1     333333332329 1 1 1    26.0   26.0
 61 4030.000000 4.179E-26 1.422E-03.04880.065 1251.80730.62-.007546    0 0 0 3 1F1    0 0 0 0 1A1   15F2 40        15F1  3     233333332329 1 1 1    93.0   93.0
 61 4031.000000 3.852E-26 3.164E-04.05230.068  814.64900.64-.007546    0 1 0 2 1F1    0 0 0 0 1A1   11E  32        12E   1     233333332329 1 1 1    46.0   50.0
 63 4032.000000 3.543E-26 6.142E-02.06230.069  218.48940.72-.008000      2V3+V5 E          GROUND  6  1  E        7  1  E      454432554528 7 1 7   156.0  180.0
 61 4033.000000 6.343E-27 6.125E-03.04850.059 1976.57460.60-.007546    0 0 1 1 1F2    0 0 0 0 1A1   18F2 99        19F1  4     233333332329 1 1 1   111.0  117.0
 61 4034.000000 3.905E-26 6.233E-04.04880.065 1096.13310.62-.007546    0 0 0 3 2F2    0 0 0 0 1A1   15F2 17        14F1  3     233333332329 1 1 1    93.0   87.0
 61 4035.000000 1.518E-27 2.408E-04.04780.061 1593.56200.61-.007546    0 0 0 3 1F1    0 0 0 0 1A1   17F1 46        17F2  2     233333332329 1 1 1   105.0  105.0
*/