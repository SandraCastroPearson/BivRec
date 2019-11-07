#' Create a bivrecFormula object.
#' @param x A formula with a bivrecSurv object on the left of a '~' operator as response, and the covariate(s) on the right.
#' @return A bivrecFormula object.
#'
#' @export
bivrecFormula <- function(x) {
  if (inherits(x, "formula") & "bivrecSurv" %in% as.character(x[[2]])) {
    class(x) <- "bivrecFormula"
  }
}

#' Check if object is bivrecFormula
#' @param x bivrecFormula or formula objects.
#' @return Logical.
#' @noRd
#' @keywords internal
#'
is.bivrecFormula <- function(x) {
  if(inherits(x, "bivrecFormula")==TRUE) {
    return(TRUE)} else {
      inherits(x, "formula") & "bivrecSurv" %in% as.character(x[[2]])}
}

