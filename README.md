# Summary
This repository contains sources for two codes explained in S. Korkin, A. M. Sayer, A. Ibrahim, and A. Lyapustin, "A practical guide to coding line-by-line trace gas absorption in Earth’s atmosphere", Journal of Quantitative Spectroscopy and Radiative Transfer, vol. 37: 109345, 2025. doi: https://doi.org/10.1016/j.jqsrt.2025.109345 (also note the corrigendum: https://doi.org/10.1016/j.jqsrt.2025.109713).

The first code is 'gcell' simulating line-by-line (LBL, monochromatic) absorption in a gas cell. Another one is 'aspect' for atmopsheric absorption spectroscopy simualtion. The two share many files - it is threfore recoemmended to look into 'gcell' before styding 'aspect'. The line parameters for both codes come from HITRAN (v.2020/24), hard-coded profiles of atmopsheric temperture, pressure, and gas concentration are from AFGL/MODTRAN.

# Gas cell absorption - code GCELL:
## Code structure (Tree & LOC)
The lines of code (LOC) count excludes headers
```
main_gcell (102)
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
