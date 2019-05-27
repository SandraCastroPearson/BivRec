significance <- function(pval) {
  sigcode <- ifelse(pval < 0.001, "***",
                    ifelse(pval < 0.01, "**",
                           ifelse(pval < 0.05, "*",
                                  ifelse(pval < 0.1, ".", " ")
                           )
                    )
  )
  return(sigcode)
}

#' @importFrom stats qnorm
#' @importFrom stats pnorm
#' @export
#'

print.bivrecSurv <- function(x) {
  #Perry
}

#' @export
print.bivrecNP <- function(x) {
  #Perry
}

#' @export
summary.bivrecReg <- function(object) {
  coeffs <- object$leefit$fit
  coeffs <- as.data.frame(cbind(coeffs[,1:2], coeffs[,1] / coeffs[,2],
                                rep(0, nrow(coeffs)), rep(0, nrow(coeffs))))
  for (i in 1:nrow(coeffs)) {
    coeffs[i,4] <- round(pnorm(abs(coeffs[i,3]), lower.tail = FALSE), digits=5)
    coeffs[i,5] <- significance(coeffs[i,4])
  }
  colnames(coeffs) <- c("Estimates", "SE", "z", "Pr(>|z|)", "")
  conf_lev = 1 - ((1-0.95)/2)
  CIcalc <- t(apply(coeffs[,1:2], 1, function (x) c(x[1]+qnorm(1-conf_lev)*x[2], x[1]+qnorm(conf_lev)*x[2])))

  expcoeffs <- data.frame(exp(coeffs[,1]), exp(-coeffs[,1]), exp(CIcalc))
  colnames(expcoeffs) <- c("exp(coef)", "exp(-coef)", "lower .95", "upper .95")

  ans <- list(call = object$call, n=object$data$response$n,
              coefficients = coeffs,
              signifcodes = "Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1",
              expcoeffs = expcoeffs)

  class(ans) <- "summary.bivrecReg"
  ans
}

#' @export
print.summary.bivrecReg <- function(x) {
  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat("\nNumber of Subjects:\n",
      paste(x$n, sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat("\nCoefficients:\n", " ", sep = "")

  print(x$coefficients)

  cat("\n---\n",
      paste(x$signifcodes, sep = "\n", collapse = "\n"), "\n\n", sep = "")

  print(x$expcoeffs)

}

#' @export
coef.bivrecReg <- function(object, ...) {
  if (!is.bivrecReg(object)) stop("Must be a bivrecReg")
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

#' @export
print.bivrecReg <- function(x) {
  if (!is.bivrecReg(x)) stop("Must be a bivrecReg")
  coeffs1 <- coef.bivrecReg(x)
  print(coeffs1)
}

#' @export
vcov.bivrecReg <- function(object) {
  if (!is.bivrecReg(object)) stop("Must be a bivrecReg")
  vcovmatrix <- object$leefit$vcovmat
  covnames <- rownames(object$leefit$fit)
  rownames(vcovmatrix) = covnames
  colnames(vcovmatrix) = covnames
  vcovmatrix

}

#' @export
confint.bivrecReg <- function(object, parm, level) {
  if (!is.bivrecReg(object)) stop("Must be a bivrecReg")
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
