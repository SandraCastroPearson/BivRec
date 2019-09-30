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

  cat("\nJoint CDF:\n", " ", sep = "")

  print(x$joint_cdf[1:6, ])

  cat("\nMarginal Survival:\n", " ", sep = "")

  print(x$marginal_survival[1:6, ])

  if (x$conditional==TRUE) {

    cat("\nConditional CDF:\n", " ", sep = "")

    print(x$conditional_cdf[1:6, ])
  }

}
