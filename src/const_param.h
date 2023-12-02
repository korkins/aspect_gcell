int const
//	path_len_max = 256,
	fname_len_max = 24,
	nmolec_id = 10,
	fpar_len_str = 160+2,
	niso_max = 12;
double const
	tiny_double = 1.0e-12,
	pi = 3.1415926535897932384626433832795,
	sqrt_pi = 1.772453850905516027298167483,
	ln2 = 0.69314718055994530941723212145818,
	sqrt_ln2 = 0.8325546111576977563531646448952,
//
	tau_min = 1.0e-4, //1.0e-8,			 // smallest total columnar tau to be saved in bin-file
//
	n_avogadro = 6.02214129e+23, // [molec/mol=mol-1]; not to be confused with the Avogadro constant: same value, but unitless
	n_loschmidt = 2.6867811e+19, // [molec/cm3] number concentration of an ideal gas at std. conditions, To_stp & Po_stp
	k_boltzman = 1.3806488e-16,  // [erg/K]
	h_planck = 6.62606957e-27,   // [erg.s]
	c_light = 2.99792458e+10,    // [cm/s]
	c2_rad = 1.4387770,          // [cm.K], c2 = h*c/k
	T_ref = 296.0,               // reference temperature, Tref[K]
	To_stp = 273.15,             // standard temperature, [T]=K, & pressure, [p]=[Pa] ... 
	Po_stp = 101325.0,           // ...https://en.wikipedia.org/wiki/Standard_conditions_for_temperature_and_pressure (01-30-2020)
	water_mass_density = 1.0,    // [g/cm3] to convert WV column (atm-cm) ...
	h2o_molar_mass = 18.01528,   // [g/mol] ... to precipitated water equivalent (cm)
//
//  1[Pa] = 10[g/(cm.s)]; 1[atm] = 1.01325e+5[Pa] = 1.01325e+6[g/(cm.s)]
	atm_to_cm_g_s = 1.01325e+6,
//  1[atm] = 1013.25 mbar
    mbar_to_atm = 1.0/1013.25,
//
	km_to_cm = 1.0e5,
//
	delta_nu = 25.0; // [cm-1], half interval to account for transition line wings
char const fname_hitdb[nmolec_id][fname_len_max] = {"01_hitdb.par", "02_hitdb.par", "03_hitdb.par", "04_hitdb.par",
                                                    "05_hitdb.par", "06_hitdb.par", "07_hitdb.par",
                                                    "nill", "nill",
                                                    "10_hitdb.par"};
//
// (1) H2O:
int const niso_h2o = 7;
char const fname_iso_h2o[niso_h2o][fname_len_max] = {"q1.txt", "q2.txt", "q3.txt", "q4.txt", "q5.txt", "q6.txt", "q129.txt"};
double const
	Qref_h2o[niso_h2o] = {174.58130000, 176.05243000, 1052.14455000, 864.74240000, 875.57259000, 5226.79494000, 1027.787680},
    molar_mass_h2o[niso_h2o] = {18.010565, 20.014811, 19.01478, 19.01674, 21.020985, 20.020956, 20.022915},
	Ia_iso_h2o[niso_h2o] = {0.997317, 0.002000, 3.718840e-4, 3.106930e-4, 6.230030e-7, 1.158530e-7, 2.419700e-8};
//
// (2) CO2:
int const niso_co2 = 12;
char const fname_iso_co2[niso_co2][fname_len_max] = {"q7.txt", "q8.txt", "q9.txt", "q10.txt",
                                                     "q11.txt", "q12.txt", "q13.txt", "q14.txt",
                                                     "q121.txt", "q15.txt", "q120.txt", "q122.txt"};
double const
	Qref_co2[niso_co2] = {  286.09382000,  576.64381000,  607.80771000,  3542.61190000,
                           1225.46188000, 7141.31955000,  323.42391000,  3766.57655000,
                          10971.56690000,  652.24137000, 7595.03572000, 22120.465100},
    molar_mass_co2[niso_co2] = {43.989830, 44.993185, 45.994076, 44.994045,
	                            46.997431, 45.997400, 47.998322, 46.998291,
	                            45.998262, 49.001675, 48.001646, 47.0016182378},
	Ia_iso_co2[niso_co2] = {9.84204e-1, 1.10574e-2, 3.94707e-3, 7.33989e-4,
	                        4.43446e-5, 8.24623e-6, 3.95734e-6, 1.47180e-6,
	                        1.36847e-7, 4.44600e-8, 1.65354e-8, 1.537500e-9};
//
// (3) O3:
int const niso_o3 = 5;
char const fname_iso_o3[niso_o3][fname_len_max] = {"q16.txt", "q17.txt", "q18.txt", "q19.txt", "q20.txt"};;
double const
	Qref_o3[niso_o3] = {3483.71020000, 7465.67515000, 3647.08033000, 43330.85026000, 21404.96421000},
    molar_mass_o3[niso_o3] = {47.984745, 49.988991, 49.988991, 48.988960, 48.988960},
	Ia_iso_o3[niso_o3] = {9.92901e-1, 3.98194e-3, 1.99097e-3, 7.40475e-4, 3.70237e-4};
//
// (4) N2O:
int const niso_n2o = 5;
char const fname_iso_n2o[niso_n2o][fname_len_max] = {"q21.txt", "q22.txt", "q23.txt", "q24.txt", "q25.txt"};
double const
	Qref_n2o[niso_n2o] = {4984.89635000, 3362.01274000, 3458.58243000, 5314.73671000, 30971.79391000},
    molar_mass_n2o[niso_n2o] = {44.001062, 44.998096, 44.998096, 46.005308, 45.005278},
	Ia_iso_n2o[niso_n2o] = {9.90333e-1, 3.64093e-3, 3.64093e-3, 1.98582e-3, 3.69280e-4};
//
// (5) CO:
int const niso_co = 6;
char const fname_iso_co[niso_co][fname_len_max] = {"q26.txt", "q27.txt", "q28.txt", "q29.txt", "q30.txt", "q31.txt"};
double const
	Qref_co[niso_co] = {107.41982000, 224.69441000, 112.77498000, 661.17319000, 236.44257000, 1384.66172000},
    molar_mass_co[niso_co] = {27.994915, 28.998270, 29.999161, 28.999130, 31.002516, 30.002485},
	Ia_iso_co[niso_co] = {9.86544e-1, 1.10836e-2, 1.97822e-3, 3.67867e-4, 2.22250e-5, 4.13292e-6};
//
// (6) CH4:
int const niso_ch4 = 4;
char const fname_iso_ch4[niso_ch4][fname_len_max] = {"q32.txt", "q33.txt", "q34.txt", "q35.txt"};
double const
	Qref_ch4[niso_ch4] = {590.47834000, 1180.82268000, 4794.72925000, 9599.15882000},
    molar_mass_ch4[niso_ch4] = {16.0313, 17.034655, 17.037475, 18.04083},
	Ia_iso_ch4[niso_ch4] = {9.88274e-1, 1.11031e-2, 6.15751E-4, 6.91785E-6};
//
// (7) O2:
int const niso_o2 = 3;
char const fname_iso_o2[niso_o2][fname_len_max] = {"q36.txt", "q37.txt", "q38.txt"};
double const
	Qref_o2[niso_o2] = {215.73450400, 455.22995200, 2658.12071500},
    molar_mass_o2[niso_o2] = {31.989830, 33.994076, 32.994045},
	Ia_iso_o2[niso_o2] = {9.95262e-1, 3.99141e-3, 7.42235e-4};
//
// (10) NO2:
int const niso_no2 = 1;
char const fname_iso_no2[niso_no2][fname_len_max] = {"q44.txt"};
double const
	Qref_no2[niso_no2] = {13577.48029000},
    molar_mass_no2[niso_no2] = {45.992904},
	Ia_iso_no2[niso_no2] = {9.91616e-1};
//
/*------------------------------------------------------------------------------------------------*/