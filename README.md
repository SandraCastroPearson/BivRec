Bivariate Alternating Recurrent Event Data Analysis (BivRec)
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

Alternating recurrent event data arise frequently in biomedical and
social sciences where two types of events such as hospital admissions
and discharge occur alternatively over time. BivRec implements a
collection of non-parametric and semiparametric methods to analyze such
data.

The main functions are:  
\- bivrecReg: Use for the estimation of covariate effects on the two
alternating event gap times (Xij and Yij) using semiparametric methods.
The method options are “Lee.et.al” and “Chang”.  
\- bivrecNP: Use for the estimation of the joint cumulative distribution
funtion (cdf) for the two alternating events gap times (Xij and Yij) as
well as the marginal survival function for type I gap times (Xij) and
the conditional cdf of the type II gap times (Yij) given an interval of
type I gap times (Xij) in a non-parametric fashion.

The package also provides options to simulate and visualize the data and
results of analysis.

## Installation

BivRec depends on the following system requirements:  
\- Rtools. Download Rtools 35 from
<https://cran.r-project.org/bin/windows/Rtools/>

Once those requirements are met you can install BivRec from github as
follows:

``` r
#Installation requires devtools package.
#install.packages("devtools")
library(devtools)
install_github("SandraCastroPearson/BivRec")
```

## Example

This is an example using a simulated data set.

``` r
# Simulate bivariate alternating recurrent event data
library(BivRec)
#> Loading required package: survival
#> Warning: package 'survival' was built under R version 3.5.3
#> 
#> Attaching package: 'BivRec'
#> The following object is masked from 'package:stats':
#> 
#>     simulate
set.seed(1234)
simdata <- simulate(nsize=150, beta1=c(0.5,0.5), beta2=c(0,-0.5), tau_c=63, set=1.1)
head(simdata)
#>   id epi      xij      yij       ci d1 d2 a1        a2
#> 1  1   1 2.411938 1.608223 40.41911  1  1  0 0.4390421
#> 2  1   2 1.405158 1.358592 40.41911  1  1  0 0.4390421
#> 3  1   3 2.188781 1.633935 40.41911  1  1  0 0.4390421
#> 4  1   4 2.045351 1.071826 40.41911  1  1  0 0.4390421
#> 5  1   5 5.047795 2.175306 40.41911  1  1  0 0.4390421
#> 6  1   6 2.503392 1.324126 40.41911  1  1  0 0.4390421

# Create a bivrecSurv object
bivrec_data <- with(simdata, bivrecSurv(id, epi, xij, yij, d1, d2))
# Plot gap times
plot(bivrec_data)
```

![](man/figures/README-BivRecExample-1.png)<!-- -->

``` r

# Apply the non-parametric method of Huang and Wang (2005) and visualize marginal and conditional results

# To save plots in a pdf file un-comment the following line of code: 
# pdf("nonparamplots.pdf")
# 
# nonpar.result <- biv.rec.np(formula = id + epi + xij + yij + d1 + d2 ~ 1,
#            data=biv.rec.data, ai=1, u1 = seq(2, 25, 1), u2 = seq(1, 20, 1),
#            conditional = TRUE, given.interval=c(0, 10), jointplot=TRUE,
#            marginalplot = TRUE, condiplot = TRUE)

# To close the pdf file with the saved plots un-comment the following line of code
# dev.off()
# 
# head(nonpar.result$joint.cdf)
# head(nonpar.result$marginal.survival)
# head(nonpar.result$conditional.cdf)

# Apply Lee C, Huang CY, Xu G, Luo X (2017) method using multiple covariates
fitlee <- bivrecReg(bivrecSurv(id, epi, xij, yij, d1, d2) ~ a1 + a2,
                     data = simdata, method="Lee.et.al")
#> [1] "Fitting model with covariates: a1,a2"
#> [1] "Estimating standard errors"
summary(fitlee)
#> 
#> Call:
#> bivrecReg(formula = bivrecSurv(id, epi, xij, yij, d1, d2) ~ a1 + 
#>     a2, data = simdata, method = "Lee.et.al")
#> 
#> 
#> Number of Subjects:
#> 150
#> 
#> 
#> Coefficients:
#>                      Estimates        SE         z Pr(>|z|)    
#> a1 type 1 gap time  0.5744414 0.1306845  4.395635  0.00001 ***
#> a2 type 1 gap time  0.5128054 0.2707497  1.894020  0.02911   *
#> a1 type 2 gap time  0.2888306 0.1985364  1.454799  0.07286   .
#> a2 type 2 gap time -0.6207422 0.3842713 -1.615375  0.05311   .
#> 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#>                    exp(coef) exp(-coef) lower .95 upper .95
#> a1 type 1 gap time 1.7761381  0.5630193 1.3747948  2.294645
#> a2 type 1 gap time 1.6699695  0.5988133 0.9823042  2.839037
#> a1 type 2 gap time 1.3348656  0.7491391 0.9045718  1.969845
#> a2 type 2 gap time 0.5375453  1.8603083 0.2531178  1.141583

# To apply Chang (2004) method use method="Chang".
# biv.rec.chang<- biv.rec.fit(formula = id + epi + xij + yij + d1 + d2 ~ a1 + a2, 
# data=biv.rec.data, method="Chang", CI=0.99)
```
