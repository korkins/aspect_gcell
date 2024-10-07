#include <stdio.h>  /* printf, fprintf */
#include <cstring>  /* strcpy, strcat */
#include <iostream> /* system("pause") */
#include <math.h>   /* pow */
#include <time.h>
#include "paths.h"
#include "const_param.h"
//
int isotops(int const, double const, double *, double *, double *);
void count_lines(int const, double const, double const, int &, int &);
void read_hitran160(int const, int const, int const, int *, double *, double *, double *,
	                double *, double *, double *, double *);
int ix1ix2(double const, double const, double const *, int const, int &, int &);
/*
	humlicek(x, y):
		- location: src_c\alyapustin_IPC\
		- requires: cmplx.cpp, cmplx.h
*/
double humlicek(double, double);
//
int main(int argc, char *argv[])
{
	int const
		nmolec = 10;
	double const
		fill_value = -999.0;
	char
		molec_name[nmolec][3+1] = {"H2O", "CO2", "O3", "N2O", "CO", "CH4", "O2",  "NA", "NA", "NO2"};
//
	FILE
		*pFile;
	bool
		flag_input_file, flag_molec_id;
	int
		code, molec_id, ichar, iline, iline0, nlines, niso, iso_ix, nnu, inu, inu1, inu2; // inu_lo, inu_hi;
	double
		time_start, time_end,
		nu_hit_min, nu_hit_max, T_kelvin, p_atm, l_cm, nu_usr_min, nu_usr_max, dnu,
		e11, e21, e12, e22, SijT, nuij_pshift, alf_doppler, gam_lorentz, x, y, n_column;
	int
		*isotop_id;
	double
		*nuij, *Sij, *gamma_air, *gamma_self, *Epp,
		*n_air, *delta_air, *nu, *k_abs;
	double
		Qratio[niso_max], mmass_iso[niso_max], Ia_iso[niso_max];
	char 
		fname_out[fname_len_max];
/*------------------------------------------------------------------------------------------------*/
//
	time_start = (double)clock()/(double)CLOCKS_PER_SEC;
//
//  Dafault scenario: O2A
	molec_id = 7;
	T_kelvin = 296.0;
	p_atm = 1.0;
	l_cm = 100.0;
	nu_usr_min = 13050.0;
	nu_usr_max = 13160.0;
	dnu = 0.01;
	strcpy(fname_out, molec_name[molec_id-1]);
	strcat(fname_out, ".txt");
	flag_input_file = false;
	flag_molec_id = false;
//
	printf(" program: %s", argv[0]);
	if (argc == 1)
	{
		printf("\n WARNING: no input file -- use O2A as default");
		flag_input_file = true;
	}
	else
	{
		printf("\n input file: %s", argv[1]);
		pFile = fopen(argv[1], "r");
		if (pFile == NULL)
		{
			printf("\n ERROR: missing file - execution terminated\n");
			flag_input_file = false;
		}
		else
		{ // read input file: data in the first line
			fscanf(pFile, "%d %lf %lf %lf %lf %lf %lf %s", &molec_id,
				&nu_usr_min, &nu_usr_max, &dnu,
					&l_cm, &T_kelvin, &p_atm, fname_out);
			fclose(pFile);
			printf("\n closed: %s", argv[1]);
			flag_input_file = true;
		} // if pFile is missing
	}
//
	if (molec_id == 1 || molec_id == 2 || molec_id == 3 || molec_id == 4 ||
		molec_id == 5 || molec_id == 6 || molec_id == 7 || molec_id == 10)
		flag_molec_id = true;
	else
		printf("\n ERROR: molec_id = %i is not included - execution terminated", molec_id);
//
	if (flag_input_file && flag_molec_id)
	{
		nu_hit_min = nu_usr_min - delta_nu;
		nu_hit_max = nu_usr_max + delta_nu;
//
//		User-defined spectral grid, nu, and spectral absorption coefficient, k_abs
		nnu = ceil((nu_usr_max - nu_usr_min)/dnu);
		nu = new double [nnu];
		k_abs = new double [nnu];
		for (inu = 0; inu < nnu; inu++)
		{
			nu[inu] = nu_usr_min + inu*dnu;
			k_abs[inu] = 0.0;
		} // inu = 0:nnu
//
//		Read HITRAN data, compute TIPS_ratio
		count_lines(molec_id, nu_hit_min, nu_hit_max, iline0, nlines);
		printf("\n iline0 = %i, nlines = %d", iline0, nlines);
//
		isotop_id = new int[nlines];
		nuij = new double[nlines];
		Sij = new double[nlines];        // already weighted with abundance!
		gamma_air = new double[nlines];  // not used for gas cell (one gas)
		gamma_self = new double[nlines];
		Epp = new double[nlines];
		n_air = new double[nlines];
		delta_air = new double[nlines];
//
		read_hitran160(molec_id, iline0, nlines, isotop_id, nuij, Sij,
						   gamma_air, gamma_self, Epp, n_air, delta_air);
//
		niso = isotops(molec_id, T_kelvin, Qratio, mmass_iso, Ia_iso);
//
		for (iline = 0; iline < nlines; iline++)
		{
			if (ix1ix2(nuij[iline], delta_nu, nu, nnu, inu1, inu2) > 0)
			{
				iso_ix = isotop_id[iline]-1;
				nuij_pshift = nuij[iline] + delta_air[iline]*p_atm; // Alexei checks for non-zero pressure shift delta_air > tiny
				alf_doppler = nuij_pshift*sqrt(2.0*n_avogadro*k_boltzman*T_kelvin*ln2/mmass_iso[iso_ix])/c_light; // save many operations
				gam_lorentz = pow(T_ref/T_kelvin, n_air[iline])*gamma_self[iline]*p_atm; // 1 gas: p_atm == p_self, gamma_air is not used
				y = sqrt_ln2*gam_lorentz/alf_doppler;
//				Temperature-corrected line intensity
//				THINKME: combine e11 & e21 and save one exp call
				e11 = exp(-c2_rad*Epp[iline]/T_kelvin);
				e21 = exp(-c2_rad*Epp[iline]/T_ref);
				e12 = 1.0 - exp(-c2_rad*nuij_pshift/T_kelvin); // nuij_pressure_shift ??
				e22 = 1.0 - exp(-c2_rad*nuij_pshift/T_ref);    // nuij_pressure_pshift ??
//				THINK ME: Alexei if-checks for 1-exp(-big_number) = 1
				SijT = Sij[iline]*Qratio[iso_ix]*e11*e12/e21/e22; // In HITRAN, Sij already accounts for isotop abundance, Ia
				for (inu = inu1; inu < inu2+1; inu++) // note: inu2 is included
				{
					x = sqrt_ln2*fabs(nu[inu] - nuij_pshift)/alf_doppler; // save operation sqrt(ln2)/gam_D
					k_abs[inu] += SijT*(sqrt_ln2/sqrt_pi)*humlicek(x, y)/alf_doppler; // in hum(x, y) make sure variables are in the right order!!
				} // inu = inu1 : inu2
			} // if iline contributes to 'nu'
			if (!(iline%1000)) printf("\n-processed: iline=%d", iline);
		} // for iline
//
//		Column (areal) *number* density: *number* of particles per unit area, [cm-2]
//		[k.T] = [erg]; 1erg = 1 g.cm2/s2: https://en.wikipedia.org/wiki/Erg 
		n_column = l_cm*atm_to_cm_g_s*p_atm/(k_boltzman*T_kelvin); // n_column = l_cm * n_volume
		printf("\n column number density, n (molec/cm2): %9.3e", n_column);
		pFile = fopen(fname_out, "w");
		fprintf(pFile, "# top line  : fill_value    column number density n_column (molec/cm2)\n");
		fprintf(pFile, "# next lines: wavenumber nu (cm-1)    optical thickness tau_abs(nu)\n");
		fprintf(pFile, "%9.3f  %12.6e\n", fill_value, n_column);
		for (inu = 0; inu < nnu; inu++)
			fprintf(pFile, "%9.3f  %12.6e\n", nu[inu], k_abs[inu]*n_column);
		fclose(pFile);
//
		delete[] k_abs;
		delete[] nu;
		delete[] isotop_id;
		delete[] nuij;
		delete[] Sij;
		delete[] gamma_air;
		delete[] gamma_self;
		delete[] Epp;
		delete[] n_air;
		delete[] delta_air;
//
		time_end = (double)clock() /(double)CLOCKS_PER_SEC;
		printf("\n output file: %s", fname_out);
		printf("\n cpu total time %6.2fs\n", time_end - time_start);
//		system("\n pause");
	} // if (flag_input_file && flag_molec_id)
}
/*------------------------------------------------------------------------------
2020-05-14:
    Fixed:  SijT = Ia_iso[iso_ix]*Sij[iline]* ...; -> SijT = Sij[iline]*...;
	Ia_iso[iso_ix] is included in HITRAN for the Earth atmopshere.
	Integarted values ratio: dagr_signal/my_signal=0.9994 - note three 9-s.
2020-03-15:
    First created and tested vs DAGR data for CH4.
	Integarted values ratio: dagr_signal/my_signal=0.9934
------------------------------------------------------------------------------*/