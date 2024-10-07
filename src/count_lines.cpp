#include <cstring>       /* strcpy, strcat */
#include <stdio.h>       /* FILE, fopen, printf, fscanf, EOF, fclose */
#include <stdlib.h>      /* atof */
#include "paths.h"
#include "const_param.h"
/*
prototype:
    void count_lines(int, double, double, int &, int &);
*/
void count_lines(int const molec_id, double const nu_min, double const nu_max, // in
	             int &iline0, int &nlines) {                                   // out
/*--------------------------------------------------------------------------------------------------
TASK:
	To find the number of lines, nlines, in [nu_min:nu_max] and index of the first line in the
	interval, iline0.
IN:
	molec_id   i   as defiend in HITRAN
	nu_min     d   wavenumber interval lower boundary
	nu_max     d   wavenumber interval upper boundary
OUT:
	iline0     i   index of the first line in the interval: iline=-1 by default
	nlines     i   number of lines in the interval: nlines=0 by default
NOTE:
	The subroutine assumes the par-file exisits and is not empty.
	iline0 = 0 corresponds to the top line in the HITRAN database file.
REFS:
	-
--------------------------------------------------------------------------------------------------*/
	FILE *fin;
	char
		fpath[path_len_max],
	    str_molec_iso_id[3], // molec_is & isotop_id are not used
		str_nuij[12+1],      // '+1' to append '\0' = end of string
        str_tail[145],       // 145 = 160 - (12_nuij + 3_molec_iso_id)
		chr_end_of_line[1];
	bool
		end_of_file;
	int
		iline;
	double
		nuij;
//--------------------------------------------------------------------------------------------------
//
	iline0 = -1;
	nlines = 0;
	end_of_file = false;
//
	strcpy(fpath, path_hitdb);
	strcat(fpath, fname_hitdb[molec_id-1]);
	fin = fopen(fpath, "r");
	printf("\n(in count_lines.cpp) opened: %s", fpath);
//
//  top line exists but may be the only line
	iline = 0;
	fscanf(fin, "%3c%12c%145c%c", str_molec_iso_id, str_nuij, str_tail, chr_end_of_line);
	str_nuij[12] = '\0';
	nuij = atof(str_nuij);
	while (nuij < nu_min)
	{
		iline += 1;
		if (fscanf(fin, "%3c%12c%145c%c", str_molec_iso_id, str_nuij, str_tail, chr_end_of_line) != EOF)
		{
			str_nuij[12] = '\0';
			nuij = atof(str_nuij);
		} // if not EOF
		else
		{
			end_of_file = true;
			break;
		} // else: if not EOF
	}
//
	if (!end_of_file)
		while (nuij <= nu_max)
		{
			nlines += 1;
			iline += 1;
			if (fscanf(fin, "%3c%12c%145c%c", str_molec_iso_id, str_nuij, str_tail, chr_end_of_line) != EOF)
			{
				str_nuij[12] = '\0';
				nuij = atof(str_nuij);
			} // if not EOF
			else
				break;
		} // while nuij <= nu_max
//
	if (nlines > 0)
		iline0 = iline-nlines;
//
	fclose(fin);
	printf("\n(in count_lines.cpp) closed: %s", fpath);
} // count_lines
/*--------------------------------------------------------------------------------------------------
2020-05-15:
            minor changes in comments;
2020-04-11: 'delta_nu' removed from the subroutine (it was included twice: in the calling program &
            in this code). Now the subroutine simply counts lines within [nu_min:nu_max] as defined
			on input. Visually tested using O2 within 13012.0(-25.0):13022.0(+25.0) - ok.
2019-12-24: created and tested using HITRAN format 160-par file, delta_nu = 2.5 cm-1,
            nu_par = 4000.0:1.0:4035; nlines_par = 36; molec_id = 6 (not important)

            Case (a): {nu_min   nu_max}   [nu_LUT_min   nu_LUT_max]
			           2000.0   3000.0     4000.0       4035.0 
  			           iline0 = -1, nlines = 0                     OK.
			Case (b): {nu_min   [nu_LUT_min   nu_max}  nu_LUT_max]
			           2000.0    4000.0       4010.0   4035.0
			           iline0 = 0, nlines = 13                     OK.
			Case (c): [nu_LUT_min {nu_min   nu_max}  nu_LUT_max]
			           4000.0      4010.0   4020.0   4035.0
			           iline0 = 8, nlines = 15                      OK.
		    Case (d): [nu_LUT_min {nu_min   nu_LUT_max]   nu_max}
			           4000.0      4010.0   4035.0        4050.0
			           iline = 8 (as in (c)!), nlines = 28         OK.
			Case (e): [nu_LUT_min   nu_LUT_max] {nu_min   nu_max}   
			           4000.0       4035.0       5000.0   6000.0      
  			           iline0 = -1, nlines = 0                     OK.
			Case (f):  {nu_min  [nu_LUT_min   nu_LUT_max] nu_max}   
			            2000.0   4000.0       4035.0       5000.0     
  			           iline0 = 0, nlines = 36                     OK.

 61 4000.000000 6.769E-27 1.237E-04.04930.065 1095.61940.62-.007545    1 0 0 1 1F2    0 0 0 0 1A1   13F1 62        14F2  2     233333332329 1 1 1    81.0   87.0
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