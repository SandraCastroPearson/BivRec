#' @title print
#' @param object An object to print
#' @keywords internal
#'
print <- function(object, ...){
  UseMethod("print")
}

print.bivrecSurv <- function(object) {
      #Perry
}

print.bivrecNP <- function(object) {
      #Perry
}

print.bivrecReg <- function(object) {
      if (!is_bivrecReg(object)) stop("Must be a bivrecReg object")
      coeffs1 <- coef.bivrecReg(object)
      print(coeffs1)
    }

print.default <- function(object) {base::print(object)}
