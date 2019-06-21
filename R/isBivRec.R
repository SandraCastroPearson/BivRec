#' Function to check object class bivrecSurv
#'
#' @keywords internal
is.bivrecSurv <- function(x) {
  res <- inherits(x, "bivrecSurv")
  return(res)}

#' Function to check object class bivrecReg
#'
#' @keywords internal
#'
is.bivrecReg <- function(x) {
  res <- inherits(x, "bivrecReg")
  return(res)}

#' Function to check object class bivrecNP
#'
#' @keywords internal
is.bivrecNP <- function(x) {
  res <- inherits(x, "bivrecNP")
  return(res)}

#' Function to check object class wholenumber
#'
#' @keywords internal
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  {abs(x - round(x)) < tol}
