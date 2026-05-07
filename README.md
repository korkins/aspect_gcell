# Summary

This repository contains the source codes described in S. Korkin, A. M. Sayer, A. Ibrahim, and A. Lyapustin, "A practical guide to coding line-by-line trace gas absorption in Earth’s atmosphere", *Journal of Quantitative Spectroscopy and Radiative Transfer*, vol. 337: 109345, 2025. doi: https://doi.org/10.1016/j.jqsrt.2025.109345 (see also the corrigendum: https://doi.org/10.1016/j.jqsrt.2025.109713).

The first code, 'gcell', simulates line-by-line (LBL, monochromatic) absorption in a gas cell. The second one, 'aspect', is intended for atmospheric absorption spectroscopy simulations. The two codes share many files; therefore, it is recommended to study 'gcell' before exploring 'aspect'.

For both codes, spectral line parameters are taken from HITRAN (v.2020), while the hard-coded atmospheric profiles of temperature, pressure, and gas concentration are based on AFGL/MODTRAN models.

For another type of absorption, see S. Korkin, A. M. Sayer, A. Ibrahim, and A. Lyapustin, "A practical guide to simulating continuum absorption in Earth’s atmosphere", *Journal of Quantitative Spectroscopy and Radiative Transfer*, vol. TBD: TBD, 202X (currently under second revision).

# Quick installation

1. Download into `./` (a local source folder; any name is fine, since paths to the HITRAN database are hard-coded as relative paths).
2. Unzip the `*.par` files into `./hitran/`.
3. The PDFs included with the source code — the paper and corrigendum — may be deleted.
4. In `./src/paths.h` (lines 6-7), check the paths and update them if necessary.
5. Run `make`. Two executables should be compiled shortly: `./gcell` and `./aspect` — see below for further instructions.

# Troubleshooting

- Both codes were originally developed on a Windows machine. The HITRAN `*.par` files were first downloaded to Windows and occasionally opened during debugging or inspection. The codes were then uploaded to GitHub. Although we tested the GitHub version on Linux, differences in line endings between Windows and Linux in the HITRAN `*.par` files may still remain. This discrepancy can cause strange errors, such as zeros in the output (when some absorption should be present) or segmentation faults. Arguably the easiest fix is to open the ASCII file on a Linux machine, add and immediately remove a space anywhere in the file, save it, and run the code again. Alternatively, download fresh `*.par` files directly from HITRAN (https://hitran.org/lbl/), while keeping in mind the version difference: the paper uses HITRAN 2020, whereas newly downloaded files may correspond to HITRAN 2024 or later.

- Depending on the machine, the user may encounter an endianness discrepancy (https://en.wikipedia.org/wiki/Endianness) while comparing their `*.bin` output files (from `aspect` only) against ours. So far, we have never encountered this problem ourselves.

- If unsuccessful, seek help from `sergey.v.korkin@nasa.gov` or `korkins@gmail.com` (sharing the input files is helpful).

# Gas cell mode - code GCELL

## How to compile and run

Assuming `./` is the source directory (with `./src`, `./hitran`, etc.), run the following commands:

```
# Compile only the 'gcell' code
$ make -f makefile_g

# Run the code without an input file (quick test: the O2 A-band runs by default)
$ ./gcell

# Run the code with a user-defined input file (see example below)
$ ./gcell gcell-ch4.inp
```

The latter command should create the file `gcell-ch4.txt` (see next section). The file can be compared against our output file `./check/gcell-ch4_check.txt` (keeping possible HITRAN version differences in mind).

## Input file format

Set of parameters (space- or tab-delimited)

```
molec_id    nu_usr_min    nu_usr_max    dnu    lcm    T_kelv    p_atm    fname
```

- `molec_id` (integer): 1 (H2O), 2 (CO2), 6 (CH4), 7 (O2), 10 (NO2)
- `nu_usr_min`, `nu_usr_max`, `dnu` (floats): define the grid of wavenumbers, ν (cm-1); the minimum (left) and maximum (right) grid points are `nu_usr_min` and `nu_usr_max`, respectively; `dnu` is the grid step (often 0.01 cm-1)
- `lcm` (float): gas cell length, `l` (cm)
- `T_kelv` (float): temperature in Kelvin (K)
- `p_atm` (float): pressure in atmospheres (atm)
- `fname` (string): output file name (optionally with path; up to 256 characters long)

Example:

```
6    4081.901    4505.699    0.002    8.0    296.0    1.0    gcell-ch4.txt
```

## Output file format

- Top line: left value is `fill_value = -999`; right value is column number density, `n_column` (molec/cm2)
- Subsequent lines: left value is wavenumber ν (cm-1); right value is optical thickness τ(ν)

Example:

```
-999.000  1.983498e+20
 4081.901  1.629810e-02
 4081.903  1.615791e-02
......................
 4505.695  4.662289e-02
 4505.697  4.712665e-02
```
Absorption cross-section per molecule is computed as: k(ν) = τ(ν) / n_column (cm2/molec)

## Code structure (Tree & LOC)

The lines of code (LOC) count excludes headers.

```
main_gcell (102)                                 # reads input file, performs calculations, prints output
         |
         +-paths.h (header)                      # defines hardcoded paths to HITRAN and TIPS files
         |
         +-const_param.h (header)                # defines physical constants, file names, accuracy parameters
         |
         +-count_lines (34)                      # counts HITRAN records within user-defined spectral band
         |           |
         |           +-paths.h
         |           |
         |           +-const_param.h
         |
         +-read_hitran160 (38)                   # reads HITRAN data from a 160-character-per-record *.par ASCII file
         |              |
         |              +-paths.h
         |              |
         |              +-const_param.h
         |
         +-isotops (59)                          # calculates TIPS ratios and returns isotope parameters
         |       |
         |       +-paths.h
         |       |
         |       +-const_param.h
         |
         +-ix1ix2 (28)                           # finds indices in array x[] such that x0-dx <= x[ix1] < x[ix2] <= x0+dx
         |
         +-humlicek (42)                         # computes Voigt profile from Lorentz and Doppler parameters (slowest part of the code)
                  |
                  +-cmplx.h (header)
                  |
                  +-cmplx (61)                   # set of complex arithmetic functions for humlicek(...)
```

LOC (excluding humlicek) = 102 + 34 + 38 + 59 + 28 = 261  

Code `aspect` does not use `main_gcell` and `isotops` functions


# Typos
1. Corrigendum: https://doi.org/10.1016/j.jqsrt.2025.109713
2. Page 5, Fig.1(a): ix1ix2() function block, comment (in green): "... x0-dx <= x[ix1] < x[ix**2**] <= x0+dx ..." (note the bold **2**)
