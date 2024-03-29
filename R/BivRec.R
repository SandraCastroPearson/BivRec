#' Bivariate Alternating Recurrent Event Data Analysis
#'
#' @description Alternating recurrent event data arise frequently in biomedical and social sciences where two types of events such as hospital admissions and discharges occur alternatively over time.
#' As such we implement a collection of nonparametric and semiparametric methods to analyze this type of data.
#' The main functions are \verb{bivrecReg} and \verb{bivrecNP}. Use \verb{bivrecReg} for the estimation of covariate effects on the two alternating event gap times (xij and yij) using semiparametric methods. The method options are "Lee.et.al" and "Chang".
#' Use \verb{bivrecNP} for the estimation of the joint cumulative distribution function (cdf) for the two alternating events gap times (xij and yij) as well as the marginal survival function for the Type I gap times (xij) and the conditional cdf of the Type II gap times (yij) given an interval of the Type I gap times (xij) in a nonparametric fashion.
#' The package also provides options to simulate and visualize the data and the results of analysis.
#'
#' @aliases BivRec-package
#' @references
#' \enumerate{
#' \item Chang S-H. (2004). Estimating marginal effects in accelerated failure time models for serial sojourn times among repeated events. Lifetime Data Analysis, 10: 175-190.
#' \doi{10.1023/B:LIDA.0000030202.20842.c9}
#'
#' \item Huang CY, Wang MC. (2005). Nonparametric estimation of the bivariate recurrence time distribution. Biometrics, 61: 392-402.
#' \doi{10.1111/j.1541-0420.2005.00328.x}
#'
#' \item Lee CH, Huang CY, Xu G, Luo X. (2018). Semiparametric regression analysis for alternating recurrent event data. Statistics in Medicine, 37: 996-1008.
#' \doi{10.1002/sim.7563}
#'
#' \item Parzen MI, Wei LJ, Ying Z. (1994). A resampling method based on pivotal estimating functions. Biometrika, 81: 341-350.
#' \url{http://www.people.fas.harvard.edu/~mparzen/published/parzen1.pdf}
#'
#' }
#'
#' @docType package
#' @author Sandra Castro-Pearson, Aparajita Sur, Chi Hyun Lee, Chiung-Yu Huang, Xianghua Luo
#' @import survival
#' @import graphics
#' @importFrom utils head
#' @useDynLib BivRec
"_PACKAGE"
NULL
