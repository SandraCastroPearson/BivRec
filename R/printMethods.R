#########################################################################
#' Print an Object of Class \code{bivrecReg}
#'
#' @description This function prints an object of class \code{bivrecReg}.
#' @param x An object of class \code{bivrecReg}.
#' @param ... Additional parameters if needed.
#'
#' @importFrom stats printCoefmat
#'
#' @export
#'

print.bivrecReg <- function(x, ...) {

  object <- x

  if (!inherits(object, "bivrecReg")) stop("Must be a bivrecReg object")

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

#########################################################################
#' Print an Object of Class \code{bivrecNP}
#'
#' @description This function prints an object of class \code{bivrecNP}.
#' @param x An object of class \code{bivrecNP}.
#' @param ... Additional parameters if needed.
#'
#' @importFrom stats printCoefmat
#'
#' @export
#'

print.bivrecNP <- function(x, ...) {

  object <- x

  if (!inherits(object, "bivrecNP")) stop("Must be a bivrecNP object")

  cat("\nJoint CDF:\n", " ", sep = "")

  printCoefmat(x$joint_cdf, digits = max(3, getOption("digits") - 2),
               signif.stars=FALSE, P.values=FALSE, has.Pvalue=FALSE)

  cat("\nMarginal Survival:\n", " ", sep = "")

  printCoefmat(x$marginal_survival, digits = max(3, getOption("digits") - 2),
               signif.stars=FALSE, P.values=FALSE, has.Pvalue=FALSE)

  if (x$conditional==TRUE) {
    cat("\nConditional CDF:\n", " ", sep = "")

    printCoefmat(x$conditional_cdf$conditional, digits = max(3, getOption("digits") - 2),
                 signif.stars=FALSE, P.values=FALSE, has.Pvalue=FALSE)
  }

}

#########################################################################
#' Print a \code{summary.bivrecReg} Object
#'
#' @description This function prints an object of class \code{summary.bivrecReg}.
#' @param x A \code{summary.bivrecReg} object.
#' @param ... Additional parameters if needed.
#'
#' @importFrom stats printCoefmat
#'
#' @export
#'
print.summary.bivrecReg <- function(x, ...) {

  if (!inherits(x, "summary.bivrecReg")) stop("Must be a bivrecReg summary object")

  cat("\nCall:\n")
  dput(x$call)

  cat("\nNumber of Subjects:\n")
  dput(as.double(x$n))

  cat("\nCoefficients:\n", " ", sep = "")
  printCoefmat(x$coefficients, digits = max(3, getOption("digits") - 2),
               signif.stars=TRUE, P.values=TRUE, has.Pvalue=TRUE)

  cat("\nexp(coefficients):\n", " ", sep = "")
  printCoefmat(x$expcoeffs, digits = max(3, getOption("digits") - 2),
               signif.stars=FALSE, P.values=FALSE, has.Pvalue=FALSE)
}
