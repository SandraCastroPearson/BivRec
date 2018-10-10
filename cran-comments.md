## Test environments
* local windows OS  x86_64-w64-mingw32, R 3.5.0
* Linux Ubuntu Trusty 14.04.5 (on travis-ci), R 3.5.0

## R CMD check results
*Local results:*
0 errors | 0 warnings | 0 notes 

*Travis-ci results:*
NOTE:
plot.joint.cdf: no visible binding for global variable ‘heat.colors’
Undefined global functions or variables:
  heat.colors

This note does not affect program since heat.colors is part of grDevices for base R


