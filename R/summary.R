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

#' Print the summary of a bivrecReg object.
#'
#' @param x a summary.bivrecReg object
#' @param ... additional parameters if needed
#'
#' @export
#'
print.summary.bivrecReg <- function(x, ...) {
  if (!inherits(x, "summary.bivrecReg")) stop("Must be a bivrecReg summary object")
  call <- x$call
  n <- x$n
  coefficients <- x$coefficients
  signifcodes <- x$signifcodes
  OddRatios <- x$OddRatios

  cat("\nCall:\n",
      paste(deparse(call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat("\nNumber of Subjects:\n",
      paste(deparse(n), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat("\nCoefficients:\n", " ", sep = "")
  print(coefficients)

  cat("\n---\n",
      paste(signifcodes, sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat("\nOdd Ratios:\n", " ", sep = "")
  print(OddRatios)

}

#' Summary of a bivrecReg object.
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

  coeffs_df <- as.data.frame(cbind(coeffs[,1:2], coeffs[,1] / coeffs[,2],
                                rep(0, nrow(coeffs)), rep(0, nrow(coeffs))))
  for (i in 1:nrow(coeffs_df)) {
    coeffs_df[i,4] <- round(pnorm(abs(coeffs_df[i,3]), lower.tail = FALSE), digits=5)
    coeffs_df[i,5] <- significance(coeffs_df[i,4])
  }
  colnames(coeffs_df) <- c("Estimates", "SE", "z", "Pr(>|z|)", "")
  conf_lev = 1 - ((1-0.95)/2)
  CIcalc <- t(apply(coeffs_df[,1:2], 1, function (x) c(x[1]+qnorm(1-conf_lev)*x[2], x[1]+qnorm(conf_lev)*x[2])))

  expcoeffs <- data.frame(exp(coeffs_df[,1]), exp(-coeffs_df[,1]), exp(CIcalc))
  colnames(expcoeffs) <- c("exp(coef)", "exp(-coef)", "lower .95", "upper .95")

  ans <- list(call = object$call, n=object$data$response$n,
              coefficients = coeffs_df,
              signifcodes = "Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1",
              OddRatios = expcoeffs)

  class(ans) <- "summary.bivrecReg"

  return(ans)
}

