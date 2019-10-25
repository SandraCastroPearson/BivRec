Bivariate Alternating Recurrent Event Data Analysis (BivRec)
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

Alternating recurrent event data arise frequently in biomedical and
social sciences where two types of events such as hospital admissions
and discharges occur alternatively over time. BivRec implements a
collection of nonparametric and semiparametric methods to analyze such
data.

The main functions are:  
\- bivrecReg: Use for the estimation of covariate effects on the two
alternating event gap times (Xij and Yij) using semiparametric methods.
The method options are “Lee.et.al” and “Chang”.  
\- bivrecNP: Use for the estimation of the joint cumulative distribution
funtion (cdf) for the two alternating events gap times (Xij and Yij) as
well as the marginal survival function for type I gap times (Xij) and
the conditional cdf of the type II gap times (Yij) given an interval of
type I gap times (Xij) in a nonparametric fashion.

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
#install_github("SandraCastroPearson/BivRec")
```

## Example

This is an example using a simulated data set.

``` r
# Simulate bivariate alternating recurrent event data
library(BivRec)
#> Registered S3 method overwritten by 'BivRec':
#>   method       from    
#>   plot.formula graphics
set.seed(1234)
sim_data <- simBivRec(nsize=150, beta1=c(0.5,0.5), beta2=c(0,-0.5), tau_c=63, set=1.1)
head(sim_data)
#>   id epi      xij      yij       ci d1 d2 a1        a2
#> 1  1   1 2.411938 1.608223 40.41911  1  1  0 0.4390421
#> 2  1   2 1.405158 1.358592 40.41911  1  1  0 0.4390421
#> 3  1   3 2.188781 1.633935 40.41911  1  1  0 0.4390421
#> 4  1   4 2.045351 1.071826 40.41911  1  1  0 0.4390421
#> 5  1   5 5.047795 2.175306 40.41911  1  1  0 0.4390421
#> 6  1   6 2.503392 1.324126 40.41911  1  1  0 0.4390421

# Create a bivrecSurv object
bivrec_object <- with(sim_data, bivrecSurv(id, epi, xij, yij, d1, d2))
# Plot gap times
plot(bivrec_object, main="Example", type = c("Type I", "Type II"))
```

![](man/figures/README-BivRecExample-1.png)<!-- -->

Nonparametric
Analysis

``` r
# Apply the nonparametric method of Huang and Wang (2005) and visualize joint, marginal and conditional results

library(BivRec)
npresult <- bivrecNP(response = bivrec_object, ai=1,
                u1 = seq(2, 25, 1), u2 = seq(1, 20, 1), conditional = TRUE,
                given.interval = c(0, 10), level = 0.99)
#> [1] "Estimating joint CDF and marginal survival"
#> [1] "Estimating conditional cdf with 99% confidence interval using 200 bootstrap samples"
head(npresult)
#> 
#> Joint CDF:
#>    x y Joint Probability         SE  Lower .99 Upper .99
#> 1 2 1        0.07765854 0.01916368 0.02829617 0.1270209
#> 2 2 2        0.12390735 0.02288661 0.06495534 0.1828594
#> 3 2 3        0.13381917 0.02345172 0.07341153 0.1942268
#> 4 2 4        0.13694252 0.02348318 0.07645387 0.1974312
#> 5 2 5        0.13928503 0.02373745 0.07814140 0.2004287
#> 6 2 6        0.13928503 0.02373745 0.07814140 0.2004287
#> 
#> Marginal Survival:
#>         Time Marginal Survival           SE Lower .99 Upper .99
#> 1 0.2697386         0.9998851 9.383939e-06 0.9998609 0.9999092
#> 2 0.2759313         0.9997701 1.149293e-04 0.9994741 1.0000000
#> 3 0.2912704         0.9996552 2.292833e-04 0.9990646 1.0000000
#> 4 0.3055601         0.9995402 3.437648e-04 0.9986548 1.0000000
#> 5 0.3203637         0.9994253 4.582784e-04 0.9982448 1.0000000
#> 6 0.3229318         0.9993103 5.728047e-04 0.9978349 1.0000000
#> 
#> Conditional CDF:
#>      Time Conditional Probability  Bootstrap SE Bootstrap Lower .99
#> 1 0.0000                  0.0000        0.0000                   0
#> 2 0.0744                  0.0000        0.0000                   0
#> 3 0.1487                  0.0055        0.0044                   0
#> 4 0.2231                  0.0133        0.0094                   0
#> 5 0.2975                  0.0194        0.0117                   0
#> 6 0.3718                  0.0247        0.0134                   0
#>   Bootstrap Upper .99
#> 1                0.00
#> 2                0.00
#> 3                0.02
#> 4                0.04
#> 5                0.05
#> 6                0.06
plot(npresult)
```

![](man/figures/README-BivRecExample2-1.png)<!-- -->![](man/figures/README-BivRecExample2-2.png)<!-- -->

``` r

# To save individual plots in a pdf file un-comment the following line of code: 
# pdf("nonparam_jointcdfplot.pdf")
# plotJoint(npresult)
# dev.off()
# pdf("nonparam_marginalplot.pdf")
# plotMarg(npresult)
# dev.off()
# pdf("nonparam_conditionaplot.pdf")
# plotCond(npresult)
# dev.off()
```

Semiparametric Regression
Analysis

``` r
#Explore how the response changes by levels of a categorical covariate using a plot.
plot(x = bivrecSurv(id, epi, xij, yij, d1, d2) ~ a1 + a2, data = sim_data,
    type = c("Type I", "Type II"))
#> [1] "a2 not used - either continuous or had more than 6 levels."
#> [1] "Original number of subjects: 150. Subjects for plots: 150"
```

![](man/figures/README-BivRecExample3-1.png)<!-- -->

``` r

# Apply Lee, Huang, Xu, Luo (2018) method using multiple covariates.
lee_fit <- bivrecReg(bivrecSurv(id, epi, xij, yij, d1, d2) ~ a1 + a2,
                    data= sim_data, "Lee.et.al")
#> Warning in if (length(d2check) == 1 & d2check == 0) {: the condition has
#> length > 1 and only the first element will be used
#> [1] "Fitting model with covariates: a1, a2"
#> [1] "Estimating standard errors"

summary(lee_fit)
#> 
#> Call:
#> bivrecReg(formula = bivrecSurv(id, epi, xij, yij, d1, d2) ~ a1 + 
#>     a2, data = sim_data, method = "Lee.et.al")
#> 
#> Number of Subjects:
#> 150
#> 
#> Coefficients:
#>         Estimates     SE        z  Pr(>|z|)    
#> xij a1   0.57444  0.13068  4.3956    1e-05 ***
#> xij a2   0.51281  0.27075  1.8940  0.02911 *  
#> yij a1   0.28883  0.19854  1.4548  0.07286 .  
#> yij a2  -0.62074  0.38427 -1.6154  0.05311 .  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> exp(coefficients):
#>         exp(coeff) Lower .95 Upper .95
#> xij a1    1.77614   1.37479    2.2946
#> xij a2    1.66997   0.98230    2.8390
#> yij a1    1.33487   0.90457    1.9698
#> yij a2    0.53755   0.25312    1.1416

# To apply Chang (2004) method use method="Chang".
# chang_fit <- bivrecReg(bivrecSurv(id, epi, xij, yij, d1, d2) ~ a1 + a2,
#                    data= sim_data, "Chang")
```
