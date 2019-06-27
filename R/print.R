#' @title print
#' @param object An object to print
#' @keywords internal
#'
print <- function(object, ...){
  UseMethod("print")
}

print.default <- function(object) {base::print(object)}

print.bivrecReg <- function(object) {
      if (!inherits(object, "bivrecReg")) stop("Must be a bivrecReg object")
      coeffs1 <- coef.bivrecReg(object)
      print(coeffs1)
    }


