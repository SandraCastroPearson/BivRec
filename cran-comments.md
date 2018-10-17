## Test environments
* Local: Windows OS  x86_64-w64-mingw32, R 3.5.1
* On travis-ci: Linux Ubuntu Trusty 14.04.5 x86_64-pc-linux-gnu, R 3.5.1

## R CMD check results

*Local results:*  
0 errors | 0 warnings | 1 note   
checking R code for possible problems ... NOTE  
plot.joint.cdf: no visible binding for global variable 'heat.colors'  
Undefined global functions or variables:  
  heat.colors  

*Travis-ci results:*  
NOTE:  
plot.joint.cdf: no visible binding for global variable ‘heat.colors’  
Undefined global functions or variables:  
  heat.colors  

This note does not affect package since heat.colors is part of grDevices for base R.  

## Reverse dependencies

This is a new release, so there are no reverse dependencies.
