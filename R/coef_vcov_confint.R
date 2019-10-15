########################    COEF     ########################

#' Extract the Coefficients of a Semiparametric Regression Fit from a \code{bivrecReg} Object
#'
#' @description This function extracts the coefficients of a semiparametric regression fit object.
#' @param object A bivrecReg object.
#' @param ... Additional parameters if needed.
#'
#' @importFrom stats printCoefmat
#'
#' @export

coef.bivrecReg <- function(object, ...) {
  #add chang
  if (!inherits(object, "bivrecReg")) stop("Must be a bivrecReg")

  if (object$method=="Chang") {
    coeffs <- object$chang_fit$fit
  } else {coeffs <- object$leefit$fit}

  coeffs<- cbind(coeffs[,1:2], coeffs[,1] / coeffs[,2],
                 rep(0, nrow(coeffs)))
  for (i in 1:nrow(coeffs)) {
    coeffs[i,4] <- round(pnorm(abs(coeffs[i,3]), lower.tail = FALSE), digits=5)
    #coeffs_df[i,5] <- significance(coeffs_df[i,4])
  }

  colnames(coeffs) <- c("Estimates", "SE", "z", "Pr(>|z|)")
  printCoefmat(coeffs, digits = max(3, getOption("digits") - 2),
               signif.stars=TRUE, P.values=TRUE, has.Pvalue=TRUE)
}

########################    VCOV     ########################

#' Extract the Variance-Covariance Matrix from a \code{bivrecReg} object
#'
#' @description This function extracts the variance-covariance matrix from the fit of a semiparametric regression analysis.
#'
#' @param object A \code{bivrecReg} object.
#' @param ... Additional parameters if needed.
#'
#' @export

vcov.bivrecReg <- function(object, ...) {
  if (!inherits(object, "bivrecReg")) stop("Must be a bivrecReg")

  if (object$method=="Chang") {
    vcovmatrix <- object$chang_fit$vcovmat
    covnames <- rownames(object$chang_fit$fit)
  } else {
    vcovmatrix <- object$leefit$vcovmat
    covnames <- rownames(object$leefit$fit)}

  rownames(vcovmatrix) = colnames(vcovmatrix) =covnames
  vcovmatrix

}

########################    confint     ########################
#' Obtain the Confidence Interval for the Coefficients of a Semiparametric Regression Fit from a \code{bivrecReg} object
#'
#' @description This function obtains the confidence interval for the coefficients of a semiparametric regression fit object.
#' @importFrom stats pnorm
#' @importFrom stringr str_extract
#' @param object A bivrecReg object.
#' @param parm The parameters for which to run confidence interval. Default gets CI for all the covariates in the model.
#' @param level Significance level. Default is 0.95.
#' @param ... Additional parameters if needed.
#'
#' @export

confint.bivrecReg <- function(object, parm, level, ...) {

  if (!inherits(object, "bivrecReg")) stop("Must be a bivrecReg")

  if (object$method=="Chang") {
    coeffs <- object$chang_fit$fit
  } else {coeffs <- object$leefit$fit}

  if (missing(level)) {level = 0.95}

  conf_lev = 1 - ((1-level)/2)

  CIcalc <- t(apply(coeffs, 1, function(x) c(x[1]+qnorm(1-conf_lev)*x[2], x[1]+qnorm(conf_lev)*x[2])))
  lowstring <- paste((1 - conf_lev)*100, "%", sep="")
  upstring <- paste(conf_lev*100, "%", sep="")
  colnames(CIcalc) <- c(lowstring, upstring)

  if (missing(parm)) {
    parm = rownames(coeffs)
    rownames(CIcalc) <- parm
    ans <- CIcalc} else {
      parm_res <- str_extract(rownames(coeffs), parm)
      ans <- CIcalc[-which(is.na(parm_res)),]
      rownames(ans) <- rownames(coeffs)[-which(is.na(parm_res))]
    }

  ans
}
