#' Print the summary of a bivrecReg object.
#'
#' @param x a summary.bivrecReg object
#' @param ... additional parameters if needed
#' @importFrom stats printCoefmat
#'
#' @export
#'
print.summary.bivrecReg <- function(x, ...) {
  if (!inherits(x, "summary.bivrecReg")) stop("Must be a bivrecReg summary object")

  cat("\nCall:\n")
  dput(x$call)

  cat("\nNumber of Subjects:\n")
  dput(x$n)

  cat("\nCoefficients:\n", " ", sep = "")
  printCoefmat(x$coefficients, digits = max(3, getOption("digits") - 2),
               signif.stars=TRUE, P.values=TRUE, has.Pvalue=TRUE)

  cat("\nOdd Ratios:\n", " ", sep = "")
  printCoefmat(x$OddRatios, digits = max(3, getOption("digits") - 2),
               signif.stars=FALSE, P.values=FALSE, has.Pvalue=FALSE)
}


#' Print bivrecNP object
#' @title print
#' @param object of a bivrecNP object obtained by using the bivrecNP() function
#' @export

print.bivrecNP <- function(object){
  if (!inherits(object, "bivrecNP")) stop("Must be a bivrecNP object")
  head(object$joint_cdf)
  head(object$marginal_survival)
  if (object$conditional==TRUE) {
    head(object$conditional_cdf$conditional)
  }
  print(paste("Confidence Interval:",npresult2$CI))
}

#' Summary of a bivrecReg object
#'
#' @param object a bivrecReg object
#' @param ... additional parameters if needed
#'
#' @export

summary.bivrecReg <- function(object, ...){

  if (!inherits(object, "bivrecReg")) stop("Must be a bivrecReg object")

  #Make summary for Chang

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
  conf_lev = 1 - ((1-0.95)/2)
  CIcalc <- t(apply(coeffs[,1:2], 1, function (x) c(x[1]+qnorm(1-conf_lev)*x[2], x[1]+qnorm(conf_lev)*x[2])))

  expcoeffs <- data.frame(exp(coeffs[,1]), exp(-coeffs[,1]), exp(CIcalc))
  colnames(expcoeffs) <- c("Odds Ratio", "Inverse OR", "lower .95", "upper .95")

  ans <- list(call = object$call, n=object$data$response$n,
              coefficients = coeffs,
              OddRatios = expcoeffs)

  class(ans) <- "summary.bivrecReg"

  return(ans)
}

