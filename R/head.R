########################    Head     ########################

#' Print top elements of the Joint CDF, Marginal Survival, and Conditional CDF from a bivrecNP object
#'
#' @param x A bivrecNP object
#' @param ... additional parameters if needed
#'
#' @noRd
#' @export

head.bivrecNP <- function(x, ...) {

  if (!inherits(x, "bivrecNP")) stop("Must be a bivrecNP")

  joint_cdf <- x$joint_cdf
  marg_suvr <- x$marginal_survival

  cat("\nJoint CDF:\n", " ", sep = "")

  print(joint_cdf[1:6, ])

  cat("\nMarginal Survival:\n", " ", sep = "")

  print(marg_suvr[1:6, ])

  if (x$conditional==TRUE) {
    cond_cdf <- x$conditional_cdf

    cat("\nConditional CDF:\n", " ", sep = "")

    print(cond_cdf[1:6, ])
  }

}
