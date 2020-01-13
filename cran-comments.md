## Package update - version 1.2.0
## New Release - previous 11/17/19 release was rejected. Issues with timing fixed.

## Test environments
* Local: Windows OS x86_64-w64-mingw32/x64 (64-bit), R 3.6.1
* On travis-ci: Linux Ubuntu Trusty 14.04.5 x86_64-pc-linux-gnu, R 3.5.1

## R CMD check results

*Local results:*   
0 errors | 0 warnings | 0 notes  

*Travis-ci results (https://travis-ci.org/SandraCastroPearson/BivRec):*   
Done. Your build exited with 0.

*check_rhub() results:*
Has 1 note when it runs under Windows Server 2008 R2 SP1, R-devel, 32/64 bit:  
* checking for non-standard things in the check directory ... NOTE
Found the following files/directories:
  'BivRec-Ex_i386.Rout' 'BivRec-Ex_x64.Rout' 'examples_i386'
  'examples_x64' 'tests_i386' 'tests_x64'


No errors, warnings or notes when using other servers. 
  
*check_win_devel() results*
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Sandra Castro-Pearson <cast0135@umn.edu>'

Cannot find any resources to fix this.

## Reverse dependencies
Previous functions deprecated but in working order.
