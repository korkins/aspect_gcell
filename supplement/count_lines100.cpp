#include <cstring>       /* strcpy, strcat */
#include <stdio.h>       /* FILE */
#include <stdlib.h>      /* atof */
#include "const_param.h" /* char path_hitdb[], char fname_hitdb[], int path_len_max */
//
void count_lines100(int const molec_id, double const nu_min, double const nu_max, // in
	             int &iline0, int &nlines) {                                   // out
/*--------------------------------------------------------------------------------------------------
PURPOSE:
	To find the number of lines, nlines, in [nu_min:nu_max] and index of the first line in the
	interval, iline0.
INPUT:
	molec_id   i   as defiend in HITRAN
	nu_min     d   wavenumber interval lower boundary
	nu_max     d   wavenumber interval upper boundary
OUTPUT:
	iline0     i   index of the first line in the interval: iline=-1 by default
	nlines     i   number of lines in the interval: nlines=0 by default
COMMENT:
	The subroutine assumes the par-file exisits and is not empty.
	iline0 = 0 corresponds to the top line in the HITRAN database file.
REFERENCESS:
	-
PROTOTYPE:
	void count_lines100(int, double, double, int &, int &);
--------------------------------------------------------------------------------------------------*/
	FILE *fin;
	char
		fpath[path_len_max],
	    str_molec_iso_id[3], // molec_is & isotop_id are not used
		str_nuij[12+1],      // '+1' to append '\0' = end of string
        str_tail[85],        // 85 = 100 - (12_nuij + 3_mlec_iso_id)
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
	printf("\n(in count_lines100.cpp) opened: %s", fpath);
//
//  top line exists but may be the only line
	iline = 0;
	fscanf(fin, "%3c%12c%85c%c", str_molec_iso_id, str_nuij, str_tail, chr_end_of_line);
	str_nuij[12] = '\0';
	nuij = atof(str_nuij);
	while (nuij < nu_min)
	{
		iline += 1;
		if (fscanf(fin, "%3c%12c%85c%c", str_molec_iso_id, str_nuij, str_tail, chr_end_of_line) != EOF)
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
			if (fscanf(fin, "%3c%12c%85c%c", str_molec_iso_id, str_nuij, str_tail, chr_end_of_line) != EOF)
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
	printf("\n(in count_lines100.cpp) closed: %s", fpath);
} // count_lines100
/*--------------------------------------------------------------------------------------------------
2020-05-05: created from count_lines.cpp to use with legacy 100-char hitran files; tested for AO2.
*/