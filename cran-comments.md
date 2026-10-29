# BFpack 1.6.3

* Date: 2026-10-29.

This is a resubmission/update shortly after 1.6.2 in response to the email
from Prof. Ripley (flang warnings on r-devel-linux-x86_64-fedora-clang,
deadline 2026-11-19).

Changes:
- Replaced the implicit 'external :: dpotrf' with an explicit interface
  consistent with the one used elsewhere in bct_mixedordinal.f90.
- Made the real constants pi (bct_mixedordinal.f90) and aln4
  (bct_prior.f90) double-precision literals.
- Removed compiled Fortran module files from src/.
   
## Test environments
* Local OS X Tahoe 26.5.1, R 4.5.3, clang-1400.0.29.202, GNU Fortran 12.2.0
* rhub check: atlas, c23, gcc13, intel, linux (R-devel)

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs

