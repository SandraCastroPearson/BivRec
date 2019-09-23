########################    Head     ########################

#' Print first five elements of the Joint CDF, Marginal Survival, and Conditional CDF from a bivrecNP object
#'
#' @param x A bivrecNP object
#' @param ... additional parameters if needed
#' @importFrom utils head
#'
#' @export

head.bivrecNP <- function(x, ...) {

  if (!inherits(x, "bivrecNP")) stop("Must be a bivrecNP")

  joint_cdf <- x$joint_cdf
  marg_suvr <- x$marginal_survival

  if (x$conditional==TRUE) {
    cond_cdf = x$conditional_cdf
  }

  cat("\nJoint CDF:\n", " ", sep = "")

  print(joint_cdf[1:5, ])

  cat("\nMarginal Survival:\n", " ", sep = "")

  print(marg_suvr[1:5, ])

  cat("\nConditional CDF:\n", " ", sep = "")

  print(cond_cdf[1:5, ])

}
