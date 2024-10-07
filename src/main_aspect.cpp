/*
NOTES:
    Command line arguments in MSVS: Project Properties\Configuration Properties\Debugging\Command Arguments

	Input:
	-molec_id:: 1:H2O, 2:CO2, 3:O3, 4:N2O, 5:CO, 6:CH4; 7:O2, 10:NO2
	-typical bands::
	--O2A: 13050-13160; BO2=14490-14600;
	--DAGR: 4081.901: 4505.699 (dnu = 0.002)
	--PACE: 335-1000 (nm) -> 10000-29900 (cm-1)
*/
#include <stdio.h>  /* printf, fprintf */
#include <math.h>   /* ceil */
#include <iostream> /* system("pause") for Win */
#include <cstring>  /* strcpy, strcat */
#include <time.h>
#include "const_param.h"
#include "hprofiles.h"
//
//int hisotops(int const, double const *, int const, double *, double *, double *);
double simpson(double const *, int const, double const);
double intparab(double const, double const, double const, double const, double const,
	                                        double const, double const, double const);
void count_lines(int const, double const, double const, int &, int &);
void read_hitran160(int const, int const, int const, int *, double *, double *, double *,
	                double *, double *, double *, double *);
void kabs(int const,
	      double const *, int const,
	      int const, double const *, double const *, double const *,
	      int const, int *, double *, double *, double *,
	                        double *, double *, double *, double *,
	      double *);
void tauabs25(double const *, double const *, double const *, int const, double *);
/*
http://modtran.spectral.com/modtran_faq
*/
int main(int argc, char *argv[])
{
	FILE
		*pFile;
	bool
		flag_input_file, flag_molec_id = false;
	int
		molec_id, iatm, iz, iz_lowest, iline0, nlines, nnu, inu, ik, nzkm, count;
	double
		time_start, time_end, hint0_24km, hint24_25km, hint25_50km, hint50_120km,
		nmolec, vol_m3, vol_cm3, scol_cm2, hcol_cm,
		nu_usr_min, nu_usr_max, dnu, nmolec_all, column_amount_usr, column_amount_mod, scalef,
		nu_hit_min, nu_hit_max, WV_mass_density, hcol_water_cm, c_stp,
		t1_sec, t2_sec, megabytes;
	int
		*isotop_id, *nuix;
	float
		*tau_float; // to minimize output binary file
	double
		*Patm, *Tkelv, *conc_cm3, *Pgas;
	double
		*nuij, *Sij, *gamma_air, *gamma_self, *Epp, *n_air, *delta_air, *nu, *knu, *gas_ratio,
		*ext_km, *tau_abs, *ztau, *conc_all, *kcont, *atmcm_km;
	char
		txt_or_bin[3+1],
		fname_out_base[fname_len_max],
		fname_out[fname_len_max];
/*------------------------------------------------------------------------------------------------*/
//
	time_start = (double)clock()/(double)CLOCKS_PER_SEC;
//
//  Dafault scenario: AO2
	molec_id = 7;
	iatm = 6;
	column_amount_usr = -1.0;
	nu_usr_min = 13050.0;
	nu_usr_max = 13160.0;
	dnu = 0.01;
	nzkm = 3;
	ztau = new double [nzkm];
	ztau[0] = 0.0; ztau[1] = 2.5; ztau[2] = 5.0;
	strcpy(txt_or_bin, "txt");
	strcpy(fname_out_base, molec_name[molec_id-1]);
//
	flag_input_file = false;
	flag_molec_id = false;
//
	printf(" program: %s", argv[0]);
	if (argc == 1)
	{
		printf("\n WARNING: no input file -- use AO2 as default");
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
			fscanf(pFile, "%d %d %lf %lf %lf %lf %d", &molec_id, &iatm, &column_amount_usr,
				&nu_usr_min, &nu_usr_max, &dnu, &nzkm);
			if (nzkm != 5) // default value
			{
				delete[] ztau;
				ztau = new double [nzkm];
			}
			for (iz = 0; iz < nzkm; iz++)
				fscanf(pFile, "%lf", &ztau[iz]);
			fscanf(pFile, "%s", fname_out_base);
			fscanf(pFile, "%s", txt_or_bin);
			fclose(pFile);
			printf("\n closed: %s\n", argv[1]);
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
		printf("\n atmosphere: %i - %s", iatm, atmos_name[iatm-1]);
		printf("\n molecule: %i - %s", molec_id, molec_name[molec_id-1]);
		printf("\n number of altitudes: %i", nzkm);
		printf("\n WARNING: z(km) must not exceed 25km");
		printf("\n altitudes, z(km):");
		for (iz = 0; iz < nzkm; iz++)
			printf("%8.3f", ztau[iz]);
		printf("\n spectral range: %9.2f - %9.2f (cm-1)", nu_usr_min, nu_usr_max);
		printf("\n output spectral resolution, dnu: %6.4f (cm-1)", dnu);
		printf("\n line wings range, delta_nu (in const_param.h): %5.2f (cm-1)", delta_nu);
//
//		Find index of the lowest altitude - usually 0, but not necessary
//      This is to check if tau_abs[iz_lowest] > tau_min to save/skip in bin file
		iz_lowest = 0;
		for (iz = 1; iz < nzkm; iz++)
			if (ztau[iz] < ztau[iz_lowest])
				iz_lowest = iz;
//
//		User-defined spectral grid, nu, and spectral absorption coefficient, k_abs
		nnu = int(ceil((nu_usr_max - nu_usr_min)/dnu))+1; // n_points = n_intervals+1
		nu = new double [nnu];
		for (inu = 0; inu < nnu; inu++)
			nu[inu] = nu_usr_min + inu*dnu;
		printf("\n nu[1] ... nu[nnu]: %12.4f - %12.4f (cm-1)", nu[0], nu[nnu-1]);
//
		gas_ratio = new double [nz_mod];
		switch(molec_id) {
			case 1:
				for (iz = 0; iz < nz_mod; iz++)
					gas_ratio[iz] = H2O_ppmv[iatm-1][iz]/1.0e6;
				break;
			case 2:
				for (iz = 0; iz < nz_mod; iz++)
					gas_ratio[iz] = CO2_ppmv[iz]/1.0e6;
				break;
			case 3:
				for (iz = 0; iz < nz_mod; iz++)
					gas_ratio[iz] = O3_ppmv[iatm-1][iz]/1.0e6;
				break;
			case 4:
				for (iz = 0; iz < nz_mod; iz++)
					gas_ratio[iz] = N2O_ppmv[iatm-1][iz]/1.0e6;
				break;
			case 5:
				for (iz = 0; iz < nz_mod; iz++)
					gas_ratio[iz] = CO_ppmv[iatm-1][iz]/1.0e6;
				break;
			case 6:
				for (iz = 0; iz < nz_mod; iz++)
					gas_ratio[iz] = CH4_ppmv[iatm-1][iz]/1.0e6;
				break;
			case 7:
				for (iz = 0; iz < nz_mod; iz++)
					gas_ratio[iz] = O2_ppmv[iz]/1.0e6;
				break;
			case 10:
				for (iz = 0; iz < nz_mod; iz++)
					gas_ratio[iz] = NO2_ppmv[iz]/1.0e6;
				break;
		} // switch(molec_id)
//
		Patm = new double[nz_mod];
		Tkelv = new double[nz_mod];
		conc_cm3 = new double[nz_mod];
		conc_all = new double[nz_mod];
		Pgas = new double [nz_mod];
		for (iz = 0; iz < nz_mod; iz++)
		{
			Tkelv[iz] = Tkelv_mod[iatm-1][iz];
			Patm[iz] = Pmbar_mod[iatm-1][iz]*mbar_to_atm;
			Pgas[iz] = Pmbar_mod[iatm-1][iz]*mbar_to_atm*gas_ratio[iz]; // below, it is scaled by ratio = content/standard_content
			conc_cm3[iz] = Dcm3_mod[iatm-1][iz]*gas_ratio[iz];
			conc_all[iz] = Dcm3_mod[iatm-1][iz];
		}
//
		hint0_24km = simpson(conc_cm3, 25, 1.0);
//	    hint24_25km = 0.5*(conc_cm3[24] + conc_cm3[25]);
		hint24_25km = intparab(zkm_mod[24], zkm_mod[25], zkm_mod[23], zkm_mod[24], zkm_mod[25],
														 conc_cm3[23], conc_cm3[24], conc_cm3[25]);
		hint25_50km = simpson(&conc_cm3[25], 11, 2.5);
		hint50_120km = simpson(&conc_cm3[35], 15, 5.0);
		nmolec = (hint0_24km + hint24_25km + hint25_50km + hint50_120km)*km_to_cm;
		vol_m3 = nmolec*k_boltzman*To_stp/Po_stp * 1.0e-7; //1.0e-7 for k_B in SI
		vol_cm3 = vol_m3*1.0e6;
		scol_cm2 = 1.0;
		hcol_cm = vol_cm3/scol_cm2;
		printf("\n column height: %10.4e (atm-cm)", hcol_cm);
//
//======================================================================================================
//  Exercise (not used in the code):
//      Convert density profile (MODTRAN) to atm-cm/km profile (Alexei):
//          atm-cm/km = atm-cm/(1e-5.cm) = 1e5.atm-cm/cm
//      where atm-cm (numerator) is the gas column height at the zero level with
//      z0, To_stp, Po_stp and cm(denominator) is the one at z, Tz, Pz.
//      Using the known equation for ideal gas (n is concentration in molec/m3)
//              p = nkT
//      one gets ratio of 2 gas column heights, h(z0)/h(z):
//              atm-cm/cm = [To_stp/Po_stp] * [Pz/Tz]
//      where
//              Pz/Tz = nk
//      Note, in SI [n] = molec/m3 = molec/(1e-6.cm3) = 1e6 cm-3
//      Combining together, one gets
//              atm-cm/km = 1e5.k.To_stp/Po_stp 1e6.conc_cm3
//      where the Boltzmann constant, k, must be in SI and both sides depend on z.
		c_stp = 1.0e5 * (k_boltzman * 1.0e-7) * To_stp / Po_stp * 1.0e6;
		atmcm_km = new double [nz_mod];
		for (iz = 0; iz < nz_mod; iz++) atmcm_km[iz] = c_stp*conc_cm3[iz]; // compare vs Alexei's profiles
//======================================================================================================
//
		if (molec_id == 1) // H2O in wv_cm
		{
			WV_mass_density = h2o_molar_mass*n_loschmidt/n_avogadro;
			printf("\n WV mass density: %10.4e (g/cm3)", WV_mass_density);
			hcol_water_cm = hcol_cm*WV_mass_density/water_mass_density;
			printf("\n MODTRAN precipitated water column amount: %6.2f (cm)", hcol_water_cm);
			if (column_amount_usr > 0.0)
			{
				printf("\n user precipitated water column amount: %6.2f (cm)", column_amount_usr);
				scalef = column_amount_usr/hcol_water_cm;
				for (iz = 0; iz < nz_mod; iz++)
				{ // gas concentration and its partial pressure grow with total column amount 
					conc_cm3[iz] *= scalef;
					Pgas[iz] *= scalef;
				} // for iz
//
				hint0_24km = simpson(conc_cm3, 25, 1.0);
//				hint24_25km = 0.5*(conc_cm3[24] + conc_cm3[25]);
				hint24_25km = intparab(zkm_mod[24], zkm_mod[25], zkm_mod[23], zkm_mod[24], zkm_mod[25],
														 conc_cm3[23], conc_cm3[24], conc_cm3[25]);
				hint25_50km = simpson(&conc_cm3[25], 11, 2.5);
				hint50_120km = simpson(&conc_cm3[35], 15, 5.0);
				nmolec = (hint0_24km + hint24_25km + hint25_50km + hint50_120km)*km_to_cm;
				vol_m3 = nmolec*k_boltzman*To_stp/Po_stp * 1.0e-7; //1.0e-7 for k_B in SI
				vol_cm3 = vol_m3*1.0e6;
				scol_cm2 = 1.0;
				hcol_cm = vol_cm3/scol_cm2;
				hcol_water_cm = hcol_cm*WV_mass_density/water_mass_density;
				printf("\n check column amount: %6.2f (cm)", hcol_water_cm);
			} // if (column_amount_usr > 0.0)
		} // if (molec_id == 1)
		else // species other than H2O
		{
			hint0_24km = simpson(conc_all, 25, 1.0);
//			hint24_25km = 0.5*(conc_all[24] + conc_all[25]);
			hint24_25km = intparab(zkm_mod[24], zkm_mod[25], zkm_mod[23], zkm_mod[24], zkm_mod[25],
														 conc_all[23], conc_all[24], conc_all[25]);
			hint25_50km = simpson(&conc_all[25], 11, 2.5);
			hint50_120km = simpson(&conc_all[35], 15, 5.0);
			nmolec_all = (hint0_24km + hint24_25km + hint25_50km + hint50_120km)*km_to_cm;
//
			column_amount_mod = nmolec*1.0e6/nmolec_all;
			printf("\n MODTRAN column amount: %8.2e (ppmv)", column_amount_mod);
			if (column_amount_usr > 0.0)
			{
				printf("\n user column amount: %8.2e (ppmv)", column_amount_usr);
				scalef = column_amount_usr/column_amount_mod;
				for (iz = 0; iz < nz_mod; iz++)
				{ // gas concentration and its partial pressure grow with total column amount
					conc_cm3[iz] *= scalef;
					Pgas[iz] *= scalef;
				} // for iz
//
				hint0_24km = simpson(conc_cm3, 25, 1.0);
//				hint24_25km = 0.5*(conc_cm3[24] + conc_cm3[25]);
				hint24_25km = intparab(zkm_mod[24], zkm_mod[25], zkm_mod[23], zkm_mod[24], zkm_mod[25],
														 conc_cm3[23], conc_cm3[24], conc_cm3[25]);
				hint25_50km = simpson(&conc_cm3[25], 11, 2.5);
				hint50_120km = simpson(&conc_cm3[35], 15, 5.0);
				nmolec = (hint0_24km + hint24_25km + hint25_50km + hint50_120km) * km_to_cm;
				column_amount_mod = nmolec*1.0e6/nmolec_all; // 1e6 for ppmv
				printf("\n check column amount: %8.2e (ppmv)", column_amount_mod);
			} // if (column_amount_usr > 0.0)
		} // else if molec_id == 1
//
//		read hitran
		nu_hit_min = nu_usr_min - delta_nu;
		nu_hit_max = nu_usr_max + delta_nu;
		printf("\n hitran spectral range: %9.2f - %9.2f (cm-1)", nu_hit_min, nu_hit_max);
//
//		run subroutine that takes hitran on input and spits k on output
		t1_sec = (double)clock()/(double)CLOCKS_PER_SEC;
		count_lines(molec_id, nu_hit_min, nu_hit_max, iline0, nlines);
		printf("\n iline0 = %i, nlines = %d", iline0, nlines);
		t2_sec = (double)clock() /(double)CLOCKS_PER_SEC;
		printf("\n runtime for 'count_lines': %6.2fs", t2_sec - t1_sec);
//
		isotop_id = new int[nlines];
		nuij = new double[nlines];
		Sij = new double[nlines];
		gamma_air = new double[nlines];
		gamma_self = new double[nlines];
		Epp = new double[nlines];
		n_air = new double[nlines];
		delta_air = new double[nlines];
//
		t1_sec = (double)clock() /(double)CLOCKS_PER_SEC;
		read_hitran160(molec_id, iline0, nlines, isotop_id, nuij, Sij,
						   gamma_air, gamma_self, Epp, n_air, delta_air);
		t2_sec = (double)clock() /(double)CLOCKS_PER_SEC;
		printf("\n runtime for 'read_hitran160': %6.2fs", t2_sec - t1_sec);
//
		megabytes = sizeof(double)/1024.0*nnu/1024.0*nz_mod;
		printf("\n knu[nnu*nz_mod] needs ~ %.1f (Mb) of RAM", megabytes);
		knu = new (std::nothrow) double [nnu*nz_mod];
		if (knu == NULL)
			printf("\n ERROR allocating knu -- execution terminated\n");
		else
		{
			t1_sec = (double)clock() /(double)CLOCKS_PER_SEC;
			for (ik = 0; ik < nnu*nz_mod; ik++)
				knu[ik] = 0.0;
			kabs(molec_id, nu, nnu, nz_mod, Tkelv, Patm, Pgas, nlines, isotop_id, nuij, Sij, gamma_air,
					  gamma_self, Epp, n_air, delta_air, knu);
			t2_sec = (double)clock() /(double)CLOCKS_PER_SEC;
			printf("\n runtime for 'kabs': %6.2fs", t2_sec - t1_sec);
//
			t1_sec = (double)clock() /(double)CLOCKS_PER_SEC;
//
			ext_km = new double [nz_mod];
			tau_abs = new double [nzkm];
//
			if (strcmp(txt_or_bin, "txt") == 0)
			{
				printf("\n Save result to ASCII file");
				strcpy(fname_out, fname_out_base);
				strcat(fname_out, ".txt");
				pFile = fopen(fname_out, "w");
				fprintf(pFile, "# columns: [1] index inu, [2] nu (cm-1), [3:] tau from TOA to zkm = ");
				for (iz = 0; iz < nzkm; iz++)
					fprintf(pFile, "%8.3f", ztau[iz]);
				for (inu = 0; inu < nnu; inu++)
				{
					fprintf(pFile, "\n%d", inu);
					fprintf(pFile, "%12.4f", nu[inu]);
					for (iz = 0; iz < nz_mod; iz++)
						ext_km[iz] = knu[inu*nz_mod+iz]*conc_cm3[iz]*km_to_cm;
					tauabs25(ext_km, zkm_mod, ztau, nzkm, tau_abs);
//				    for (iz = 0; iz < nzkm; iz++) fprintf(pFile, "%16.6e", knu[inu*nz_mod+iz]); // <-to test k_abs apart from z-integration
					for (iz = 0; iz < nzkm; iz++) fprintf(pFile, "%16.6e", tau_abs[iz]);
				} // for inu
				fclose(pFile);
			} // if ascii
			else
			{ // binary
				printf("\n Save result to BIN file: tau[inu][iz], iz is the lead dimension");
				strcpy(fname_out, fname_out_base);
				strcat(fname_out, ".bin");
				tau_float = new float [nzkm];
				pFile = fopen(fname_out, "wb");
//
				count = 0;
				nuix = new int [nnu];
				for (inu = 0; inu < nnu; inu++)
				{
					for (iz = 0; iz < nz_mod; iz++)
						ext_km[iz] = knu[inu*nz_mod+iz]*conc_cm3[iz]*km_to_cm;
					tauabs25(ext_km, zkm_mod, ztau, nzkm, tau_abs);
					if (tau_abs[iz_lowest] > tau_min)
					{
						for (iz = 0; iz < nzkm; iz++)
							tau_float[iz] = (float)tau_abs[iz];
						fwrite(tau_float, nzkm*sizeof(float), 1, pFile);
						nuix[count] = inu;
						count += 1;
					}
				} // for inu
				fclose(pFile);
				printf("\n saved: tau in *.bin, dtype=float32");
//
				strcpy(fname_out, fname_out_base);
				strcat(fname_out, "_inu.bin");
				pFile = fopen(fname_out, "wb");
				fwrite(nuix, count*sizeof(int), 1, pFile);
				fclose(pFile);
				printf("\n saved: inu.bin, dtype=int");
//
				strcpy(fname_out, fname_out_base);
				strcat(fname_out, "_dat.txt");
				pFile = fopen(fname_out, "wt");
				fprintf(pFile, " spectral interval: left bound {nu0}, resolution {dnu}, nu[inu] = nu0 + inu*dnu");
				fprintf(pFile, "\n %10.4f  %6.4f", nu_usr_min, dnu);
				fprintf(pFile, "\n data shape: number of records {nnu}, number of heights {nzkm}:");
				fprintf(pFile, "\n %i  %i", count, nzkm);
				fprintf(pFile, "\n py: data = np.fromfile('*.bin', dtype=np.float32)");
				fprintf(pFile, "\n py: tau[0:nnu, 0:nzkm] = np.reshape(data, (nnu, nzkm))");
				fprintf(pFile, "\n dtype(*_inu.bin) = int_32bit");
				fprintf(pFile, "\n dtype(*.bin) = float_32bit");
				fclose(pFile);
			} // if ascii or binary
//
			t2_sec = (double)clock()/(double)CLOCKS_PER_SEC;
			printf("\n runtime for 'tauabs25' + write to file: %6.2fs", t2_sec - t1_sec);
//
			time_end = (double)clock()/(double)CLOCKS_PER_SEC;
			printf("\n cpu total time %6.2fs", time_end - time_start);
			printf("\n file = %s\n", fname_out);
		} // if (knu = new double [nnu*nz_mod]) else
	} // if flag = true - run the code
//	system("pause");
}
/*------------------------------------------------------------------------------*/