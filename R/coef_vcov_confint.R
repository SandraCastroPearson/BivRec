########################    COEF     ########################

#' Obtain coefficients from semi-parametric regression fit using bivrecReg
#'
#' @title coef
#' @param object A bivrecReg object obtained by using bivrecReg() function
#' @export
coef <- function(object) {UseMethod("coef")}

coef.default <- function(object) stats::coef(object)

coef.bivrecReg <- function(object) {
  #add chang
  if (!inherits(object, "bivrecReg")) stop("Must be a bivrecReg")
  coeffs <- object$leefit$fit
  coeffs <- as.data.frame(cbind(coeffs[,1:2], coeffs[,1] / coeffs[,2],
                                rep(0, nrow(coeffs)), rep(0, nrow(coeffs))))
  for (i in 1:nrow(coeffs)) {
    coeffs[i,4] <- round(pnorm(abs(coeffs[i,3]), lower.tail = FALSE), digits=5)
    coeffs[i,5] <- significance(coeffs[i,4])
  }
  colnames(coeffs) <- c("Estimates", "SE", "z", "Pr(>|z|)", "")
  coeffs
}

########################    VCOV     ########################

#' Obtain variance-covariance matrix from semi-parametric regression fit using bivrecReg
#'
#' @title vcov
#' @param object A bivrecReg object obtained by using bivrecReg() function
#'
#' @export
vcov <- function(object) UseMethod("vcov")

vcov.default <- function(object) stats::vcov(object)

vcov.bivrecReg <- function(object) {
  if (!inherits(object, "bivrecReg")) stop("Must be a bivrecReg")
  #do for Chang
  vcovmatrix <- object$leefit$vcovmat
  covnames <- rownames(object$leefit$fit)
  rownames(vcovmatrix) = covnames
  colnames(vcovmatrix) = covnames
  vcovmatrix

}

########################    confint     ########################
#' Obtain confidence interval for exponentiated coefficients of semi-parametric regression fit using bivrecReg
#'
#' @title confint
#' @importFrom stats pnorm
#' @param object A bivrecReg object obtained by using bivrecReg() function
#' @param parm The parameters for which to run confidence interval. Default is giving CI for all the covariates in the model.
#' @param level Significance level. Example: 0.99 for a 99\% confidence interval. Default is 0.95.
#'
#' @export
confint <- function(object, parm, level) UseMethod("confint")

confint.default <- function(object, parm, level) stats::confint(object)

confint.bivrecReg <- function(object, parm, level) {
  if (!inherits(object, "bivrecReg")) stop("Must be a bivrecReg")
  #do for Chang
  coeffs <- object$leefit$fit
  if (missing(level)) {level = 0.95}
  if (missing(parm)) {parm = rownames(coeffs)}
  conf_lev = 1 - ((1-level)/2)
  CIcalc <- t(apply(coeffs, 1, function(x) c(x[1]+qnorm(1-conf_lev)*x[2], x[1]+qnorm(conf_lev)*x[2])))
  ans  <- cbind(coeffs, CIcalc)
  lowstring <- paste((1 - conf_lev), "%", sep="")
  upstring <- paste(conf_lev, "%", sep="")
  colnames(ans) <- c("Estimate", "SE", lowstring, upstring)
  rownames(ans) <- rownames(coeffs)
  ans
}
