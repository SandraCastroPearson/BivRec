#' Non-Parametric Analysis of Bivariate Recurrent Event Data
#'
#' @description
#' This function allows the user to find the joint cdf, marginal survival and conditional cdf for two alternating states/gap times.
#'
#' @importFrom stats model.frame
#' @importFrom stats na.omit
#' @importFrom stats quantile
#'
#' @param formula Formula in the form: \strong{ID + xij + yij + episode + c_indicatorY + c_indicatorX ~ 1}.
#'  where
#'  \strong{ID:} is a vector of subjects' unique identifier.
#'  \strong{xij:} is a vector with the lengths of time spent in event of type X for individual i in episode j.
#'  \strong{yij:} is a vector with the lengths of time spent in event of type Y for individual i in episode j.
#'  \strong{episode:} is a vector indicating the observation or episode j for a subject i.
#'  \strong{c_indicatorY:} is an vector of indicators with values of 0 for the last episode for subject i and 1 otherwise. A subject with only one episode will have one 0.
#'  \strong{c_indicatorX:} is an optional vector of indicators with values of 0 if the last episode for subject i occurred for event of type X or 1 otherwise. A subject with only one episode could have only a 1 (if he was censored at event Y) or a 0 (if he was censored at event X).
#' @param data A data frame that includes all the vectors/covariates listed in the formula
#' @param ai A real non-negative function of censoring time to be used as weights in the Non-Parametric Method. Current options are 1 such that ai = f(censoring times) = 1 or the default of 2 such that ai = f(censoring times) = censoring times.
#' @param u1 A vector or single number to be used for Estimation of joint CDF P(xij <= u1, yij <= u2) in the Non-Parametric Method. (see vignette for further details)
#' @param u2 A vector or single number to be used for Estimation of joint CDF P(xij <= u1, yij <= u2) in the Non-Parametric Method. (see vignette for further details)
#' @param marginalplot Logical. If TRUE (default) this function will plot the marginal survival for the first gap time with 95\% confidence interval.
#' @param conditional Logical. If TRUE will calculate the conditional cdf for the second gap time given an interval for the first gap time. Must provide given.interval.
#' @param given.interval is a 1x2 vector indicating an interval for the first gap time to estimate the cdf of second gap time given this interval. If given.interval = c(v1, v2) the function calculates P(yij <= u2 | v1 <= xij <= v2).
#' The given values v1 and v2 must be in the range of gap times in the estimated marginal survival (valid times are given in the Time column of the results for marginal survival).
#'
#' @return A list with dataframes of results and some of the parameters used to produce them.
#' @export
#'
#' @keywords biv.rec.np
#'
#' @examples
#' \dontrun{
#' library(BivRec)
#' set.seed(1234)
#' sim.data <- data.sim(nsize=300, beta1=c(0.5,0.5), beta2=c(0,-0.5), cr=63, sg2=0.5, set=1.1)
#' nonpar.result <- biv.rec.fit(id + xij + yij + epi + d2 + d1 ~ 1,
#'           data=sim.data, ai=1)
#' nonpar.result$cdf
#' }
#'

biv.rec.np <- function(formula, data, ai, u1, u2, marginalplot, conditional, given.interval){

  if (missing(ai)) {ai<-1}
  if (missing(u1)) {u1 <- round(seq(quantile(xij, probs = 0.4), max(xij), length.out=5))}
  if (missing(u2)) {u2 <- round(seq(quantile(yij, probs = 0.4), max(yij), length.out=4))}
  if (missing(marginalplot)) {marginalplot <- TRUE}
  if (missing(conditional)) {conditional <- FALSE}

  ### PULL INFORMATION FROM PARAMETERS TO SEND TO REFORMAT
  variables <- all.vars(formula)

  ####Ensure unique identifiers are numeric
  iden <- eval(parse(text = paste("data$", variables[1], sep="")))
  iden.u <- unique(iden)
  new.id <- NULL
  if (class(iden)!="num") {
    if (class(iden)!="int") {
      for (i in 1:length(iden.u)){
        for (j in 1:length(iden)) {
          if (iden[j] == iden.u[i]){
            new.id=c(new.id,i)
          }
        }
      }
      data$new.id <- new.id
    }
  }
  data <- data[,-which(colnames(data)==variables[1])]
  colnames(data)[ncol(data)] = variables[1]

  ####extract vectors/data needed to send to biv.rec.reformat
  names <- paste("data$", variables, sep="")
  identifier <- eval(parse(text = names[1]))
  xij <- eval(parse(text = names[2]))
  yij <- eval(parse(text = names[3]))
  episode <- eval(parse(text = names[4]))
  c_indicatorY <- eval(parse(text = names[5]))
  ifelse(length(names)==6, c_indicatorX <- eval(parse(text = names[6])), c_indicatorX <- rep(1, length(xij)))
  covariates <- rep(1, length(identifier))
  cov_names <- "No Covariates"
  temp <- rep(u1, each = length(u2))
  temp2 <- rep(u2, length(u1))
  u <- cbind(u1=temp, u2=temp2)
  method <- "Non-Parametric"
  condgx <- FALSE

  ###Send to biv.rec.reformat and complete analysis
  new_data <- biv.rec.reformat(identifier, xij, yij, c_indicatorY, c_indicatorX, episode, covariates, method, ai, condgx, data)
  res1 <- nonparam.cdf(fit_data=new_data$forcdf, u, ai)
  res2 <- nonparam.marginal(new_data$formarg)

  if (marginalplot==TRUE) {marg.surv.plot(list(marginal.survival = res2, formula=formula, data=data))}
  if (conditional == FALSE) {
    final.result <- list(joint.cdf = res1, marginal.survival = res2, formula=formula, ai=ai)
  } else {
    if (missing(given.interval)) {
      print("Error for conditional calculation given.interval argument missing.")
      final.result <- list(cdf = res1, marginal.survival = res2, formula=formula, data = data, ai=ai)
    } else {
      partial.result <- list(cdf = res1, marginal.survival = res2, formula=formula, data = data, ai=ai, new_data=new_data)
      res3 <- nonparam.conditional(partial.result, given.interval)
      final.result <- list(joint.cdf = res1, marginal.survival = res2, conditional.cdf = res3, formula=formula, ai=ai)
    }
  }

  return(final.result)
}

