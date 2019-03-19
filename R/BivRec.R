#################### CREATE A BIVREC OBJECT ######################

#####
#' A function to create a biv.rec object
#'
#' @description
#' This function creates  BivRec response object.
#'
#' @importFrom stats na.omit
#' @importFrom survival survfit
#'
#' @param identifier Vector of subject's unique identifier (i).
#' @param xij Vector with the lengths of time spent in event of type X for individual i in episode j.
#' @param yij vector with the lengths of time spent in event of type Y for individual i in episode j.
#' @param Y_censoring Vector with values of 0 for the last episode for subject i or 1 otherwise. A subject with only one episode will have one 0.
#' @param X_censoring Vector with values of 0 if the last episode for subject i occurred for event of type X or 1 otherwise. A subject with only one episode could have either one 1 (if he was censored at event Y) or one 0 (if he was censored at event X). A subject with censoring in event Y will have a vector of 1's.
#' @param episode Vector indicating the observation or episode (j) for a subject (i). This will determine order of events.
#' @param method A string for method to be used. Passed from biv.rec.fit()
#'
#' @return a BivRec repsonse object ready to put in a formula.
#' @seealso \code{\link{BivRec.fit}}
#'
#' @keywords export
