# Version 0.3.1
* New extension for BF for rma.uni object (metafor package) for Bayesian meta-analysis
* BF for lmerMod object now also works on unbalanced data
* vignette was added
* Bugs fixes
* Minor changes to output format

## Test environments
* Local OS X 10.14.6, R 3.6.1
* rhub check: Ubuntu Linux 16.04 LTS, R-release
* rhub check: Fedora Linux, R-devel
* rhub check: Debian Linux, R-devel, GCC ASAN/UBSAN
  + PREPERROR: Dependency 'openssl' is not available. Bug reported on R-hub GitHub page https://github.com/r-hub/sysreqsdb/issues/77
* rhub check: Windows Server 2008 R2 SP1, R-devel, 32/64 bit
  + PREPERROR
* winbuilder: R 3.6.3 (2020-02-29)
* winbuilder: R Under development (unstable) (2020-05-08)

## R CMD check results
There were no ERRORs or WARNINGs, one NOTE
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit
checking for non-standard things in the check directory ... NOTE
Found the following files/directories:
  'BFpack-Ex_i386.Rout' 'BFpack-Ex_x64.Rout' 'examples_i386'
  'examples_x64' 'tests_i386' 'tests_x64'
