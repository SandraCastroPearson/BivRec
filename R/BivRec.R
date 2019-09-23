#' BivRec
#'
#' @description Alternating recurrent event data arise frequently in biomedical and social sciences where 2 types of events such as hospital admissions and discharge occur alternatively over time.
#' As such we implement a collection of non-parametric and semiparametric methods to analyze such data.
#' main functions are \code{bivrecReg}() and \code{bivrecNP}(). Use \code{bivrecReg}() for estimation of covariate effects on the two alternating event gap times (xij and yij) using semiparametric methods. The method options are "Lee.et.al" and "Chang".
#' Use \code{bivrecNP}() for estimation of the joint cumulative distribution function (cdf) for the two alternating events gap times (xij and yij) as well as the marginal survival function for type I gap times (xij) and the conditional cdf of the type II gap times (yij) given an interval of type I gap times (xij) in a non-parametric fashion.
#' The package also provides options to simulate and visualize the data and results of analysis.
#'
#' @docType package
#' @author Sandra Castro-Pearson
#' @import Rcpp
#' @import survival
#' @import graphics
#' @useDynLib BivRec "bivrecur", .registration = TRUE
#' @useDynLib BivRec "mprovar", .registration = TRUE
#' @useDynLib BivRec "onesamp", .registration = TRUE
#' @useDynLib BivRec "xmproee", .registration = TRUE
#' @useDynLib BivRec "ymproee", .registration = TRUE
#' @name BivRec
