# Summary
This repository contains sources for two codes explained in S. Korkin, A. M. Sayer, A. Ibrahim, and A. Lyapustin, "A practical guide to coding line-by-line trace gas absorption in Earth’s atmosphere", *Journal of Quantitative Spectroscopy and Radiative Transfer*, vol. 37: 109345, 2025. doi: https://doi.org/10.1016/j.jqsrt.2025.109345 (also note the corrigendum: https://doi.org/10.1016/j.jqsrt.2025.109713).

The first code is 'gcell' simulating line-by-line (LBL, monochromatic) absorption in a gas cell. Another one is 'aspect' for atmopsheric absorption spectroscopy simualtion. The two share many files - it is threfore recoemmended to look into 'gcell' before styding 'aspect'. The line parameters for both codes come from HITRAN (v.2020), hard-coded profiles of atmopsheric temperture, pressure, and gas concentration are from AFGL/MODTRAN.

# Quick installation
1. Download to ./ (local folder with sources - any name is fine, paths to HITRAN database are hardcoded but relative)
2. Unzip *.par files in ./hitran/
3. PDFs comming with the source - paper and corrigendum - can be deleted
4. In ./src/paths.h, lines 6-7, check paths and update if necessary
5. Run 'make'. Two executables should be compiled shortly: ./gcell and ./aspect - see below for further instructions.

# Gas cell mode - code GCELL
## Input file format
Set of parametrs (space or tab dilimited)  
```
    molec_id    nu_usr_min    nu_usr_max    dnu    lcm    T_kelv    p_atm    fname  
```
- molec_id (integer) 1(H2O) 2(CO2) 6(CH4) 7(O2) 10(NO2)  
- nu_usr_min, nu_usr_max, dnu (floats) define grid of wavenumbers, ν (cm-1); the minimum (left) and maximum (right) points in the grid are 'nu_usr_min' and 'nu_usr_max', respectively; dnu is step (often 0.01 cm-1)  
- lcm (float) gas cell length l (cm)  
- T_kelv (float) temperature in degrees Kelvin (K)  
- p_atm (float) pressure in atmospheres (atm)  
- fname (string) defines output file name (optionally, with path - up to 256 symbols long)  

Example:  
```
    6    4081.901    4505.699    0.002    8.0    296.0    1.0    gcell-ch4.txt
```
## Output file format
top line, left number: fill_value=-999; right number: column number density, n_column (molec/cm2)  
next lines, left number: wavenumber v (cm-1); right numbers: optical thickness tau(v)  

Example:  
-999.000  3.621485e+14
500.000  1.613564e-17
500.010  1.301552e-17
...
If needed, calculate absorption cross-section per molecule k(v) = tau / n_column (cm2/molec)

## Code structure (Tree & LOC)
The lines of code (LOC) count excludes headers
```
main_gcell (102)                                 # reads input file, makes calculations, prints output
         |
         +-paths (header)                        # defines hardcoded paths to HITRAN and TIPS files
         |
         +-const_param (header)                  # defines physical constants, file names, accuracy parameters
         |
         +-count_lines (34)                      # counts HITRAN records within user-defined spectral band
         |           |
         |           +-paths (header)
         |           |
         |           +-const_param (header)
         |
         +-read_hitran160 (38)                   # reads HITRAN data from a 160 symbols per record *.par ASCII file
         |              |
         |              +-paths (header)
         |              |
         |              +-const_param (header)
         |
         +-isotops (59)                          # calculates TIPS ratio and returns parameters of isotopes
         |       |
         |       +-paths (header)
         |       |
         |       +-const_param (header)
         |
         +-ix1ix2 (28)                           # for x0 and dx, finds indices in an array x[] so that x0-dx <= x[ix1] < x[ix2] <= x0+dx
         |
         +-humlicek (42)                         # calculates Voigt profile from Lorentz and Doppler line shape parameters (slowest part of the code)
                  |
                  +-cmplx.h (header)
                  |
                  +-cmplx (61)

LOC (excluding humlicek) = 102 + 34 + 38 + 59 + 28 = 261
```

# Typos
1. Corrigendum: https://doi.org/10.1016/j.jqsrt.2025.109713
2. Page 5, Fig.1(a): ix1ix2() function block, comment (in green): "... x0-dx <= x[ix1] < x[ix**2**] <= x0+dx ..." (note the bold **2**)
