#' Deprecated: Use bivrecReg
#'
#' @description
#' Deprecated function from the previous version. Use \verb{bivrecReg}.
#'
#' @param formula A formula with six variables indicating the bivariate alternating gap time
#'
#' response on the left of the ~ operator and the covariates on the right.
#'
#' The six variables on the left must have the same length and be given as
#'
#' \verb{ID + episode + xij + yij + d1 + d2 ~ covariates}, where:
#'
#' \itemize{
#'  \item \verb{id}: Vector of subject's unique identifier (i).
#'  \item \verb{episode}: Vector indicating the observation or episode (j) for a subject (i). This will determine order of events for each subject.
#'  \item \verb{xij}: Vector with the lengths of time spent in event of Type I for individual i in episode j.
#'  \item \verb{yij}: Vector with the lengths of time spent in event of Type II for individual i in episode j.
#'  \item \verb{d1}: Vector of censoring indicator corresponding to Type I gap times (xij): = 1 for uncensored, and = 0 for censored gap times.
#'  \item \verb{d2}: Vector of censoring indicator corresponding to Type II gap times (yij): = 1 for uncensored, and = 0 for censored gap times.
#'  \item covariates: the names of the covariates in the form covariate_1 + ... + covariate_p.
#' }
#' @param data A data frame that includes all the vectors/covariates listed in the formula above.
#' @param method A string indicating which method to use to estimate effects of the covariates. See details.
#' @param CI The level to be used for confidence intervals. Must be between 0.50 and 0.99. The default is 0.95.
#'
#' @return See \verb{bivrecReg}
#'
#' @importFrom stats as.formula
#'
#' @details
#' Two different estimation methods are available:
#' \itemize{
#' \item  method = "Lee.et.al" (default) is a U-statistics-based smooth estimating function approach. See Lee, Huang, Xu, Luo (2018) for further details.
#' \item  method = "Chang" is a rank-based estimating function approach.  See Chang (2004) for further details.
#' Note that following the Chang method, the variances of the estimated regression coefficients are approximated using the resampling techniques developed by Parzen, Wei and Ying (1994).
#' This approximation requires extensive computing time for a relatively small sample size. In addition, using the Chang method does not guarantee convergence for the estimation of the coefficients.
#' }
#'
#' @references
#' \enumerate{
#' \item Chang S-H. (2004). Estimating marginal effects in accelerated failure time models for serial sojourn times among repeated events. Lifetime Data Analysis, 10: 175-190.
#' \doi{10.1023/B:LIDA.0000030202.20842.c9}
#'
#' \item Lee CH, Huang CY, Xu G, Luo X. (2018). Semiparametric regression analysis for alternating recurrent event data. Statistics in Medicine, 37: 996-1008.
#' \doi{10.1002/sim.7563}
#'
#' \item Parzen MI, Wei LJ, Ying Z. (1994). A resampling method based on pivotal estimating functions. Biometrika, 81: 341-350.
#' \url{http://www.people.fas.harvard.edu/~mparzen/published/parzen1.pdf}
#' }
#'
#' @export

biv.rec.fit <- function(formula, data, method, CI){

  .Deprecated("bivrecReg")

    resp <- paste(all.vars(formula[[2]]), collapse=" + ")
    part1 <- paste("bivrecSurv(", resp,")", sep="")
    part2 <- paste(all.vars(formula[[3]]), collapse=" + ")
    newformula <- as.formula(paste(part1, " ~ ", part2, sep=" "))

    ans <- bivrecReg(formula = newformula, data = data, method=method)

    return(ans)

}

