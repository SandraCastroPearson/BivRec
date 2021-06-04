## Package update - version 1.2.1 as of 06-04-2021

## Test environments
* Local: Windows OS x86_64-w64-mingw32/x64 (64-bit), R 4.0.2
* On travis-ci: Linux Ubuntu Trusty 14.04.5 x86_64-pc-linux-gnu, R 3.5.0 and 4.0.2

## R CMD check results

*Local results:*   
Duration: 2m 13.6s
0 errors | 0 warnings | 0 notes

*Travis-ci results (https://travis-ci.org/SandraCastroPearson/BivRec):*   
Done. Your build exited with 0.

*check_win_devel() and check_win_release() results*
* DONE
Status: OK

*check_rhub() results:*

- BivRec 1.2.1: OK for all platforms except Windows Server 2008 R2 SP1

- PREPERROR for	Windows Server 2008 R2 SP1, R-devel, 32/64 bit
  Error: Bioconductor does not yet build and check packages for R version 4.2; see
  https://bioconductor.org/install
  Execution halted

  Found a reference to this issue at: https://github.com/r-hub/rhub/issues/471. No       information on how to fix this as it doesn't apear to be a package issue but a         testing set-up issue.
