# Version 0.2.1
* Fixes bugs
* New function cor_test included for Bayesian correlation analysis
* BF for lmerMod object now also works on unbalanced data

## Test environments
* Local OS X 10.14.6, R 3.5.2
* Local Windows 10 (x64 and x32, build 15063), R 3.6.1
* Travis CI OS X 10.13.3, R 3.6.1
* Travis CI Ubuntu 16.04.6 LTS, R 3.6.1
* Travis CI Ubuntu 16.04.6 LTS, R 3.5.3
* Travis CI Ubuntu 16.04.6 LTS, R Under development (unstable) (2019-10-25 r77332)
* rhub check: Windows Server 2008 R2 SP1, R-devel, 32/64 bit
* rhub check: Ubuntu Linux 16.04 LTS, R-release
* rhub check: Fedora Linux, R-devel
* rhub check: Debian Linux, R-devel, GCC ASAN/UBSAN
  + PREPERROR: Dependency 'openssl' is not available. Bug reported on R-hub GitHub page https://github.com/r-hub/sysreqsdb/issues/77
* winbuilder: R 3.5.3 (2019-03-11)
* winbuilder: R 3.6.1 (2019-07-05)
 

## R CMD check results
There were no ERRORs or WARNINGs, one NOTE
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit
checking for non-standard things in the check directory ... NOTE
Found the following files/directories:
  'BFpack-Ex_i386.Rout' 'BFpack-Ex_x64.Rout' 'examples_i386'
  'examples_x64' 'tests_i386' 'tests_x64'
