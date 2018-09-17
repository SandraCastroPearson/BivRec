#' Non-Parametric Analysis of Bivariate Alternating Recurrent Event Data
#'
#' @description
#' This function allows the user to apply a non-parametric method to estimate the joint cumulative distribution funtion (cdf) for the two alternating events gap times (xij and yij)
#' as well as the marginal survival function for type I gap times (xij) and the conditional cdf of the type II gap times (yij) given an interval of type I gap times (xij).
#' See Huang and Wang (2005) for more details.
#'
#' @importFrom stats model.frame
#' @importFrom stats na.omit
#' @importFrom stats quantile
#'
#' @param formula A formula with six variables indicating the bivariate alternating gap time response on the left of the ~ operator and a 1 on the right.
#' The six variables on the left must have the same length and be given as follows \strong{ID + episode +  xij + yij + delta_x + delta_y ~ 1}, where
#' \itemize{
#'   \item ID: A vector of subjects' unique identifier which can be numeric or character.
#'   \item episode: A vector indicating the episode of the bivariate alternating gap time pairs, e.g.: 1, 2, ..., \eqn{m_i} where \eqn{m_i} indicates the last episode for subject i.
#'   \item xij: A vector with the lengths of the type I gap times.
#'   \item yij: A vector with the lengths of the type II gap times.
#'   \item delta_x: An optional vector of indicators with values:
#'   \itemize{
#'       \item 0 for the last episode for subject i (\eqn{m_i}) if subject was censored during period xij.
#'       \item 1 otherwise.
#'      }
#'   A subject with only one episode (\eqn{m_i=1}) could have a 0 if he was censored during period xi1 or 1 if he was censored during period yi1. If delta_x is not provided estimation proceeds with the assumption that no subject was censored during period xij.
#'   \item delta_y: A vector of indicators with values:
#'   \itemize{
#'       \item 0 for the last episode of subject i (\eqn{m_i}).
#'       \item 1 otherwise.
#'      }
#'   A subject with only one episode (\eqn{m_i=1}) will have one 0.
#' }
#' @param data A data frame that includes all the vectors listed in the formula.
#' @param ai A real non-negative function of censoring time. See details.
#' @param u1 A vector or single number to be used for estimation of joint cdf \eqn{P(X0 \le u1, Y0 \le u2)} in the non-parametric method.
#' @param u2 A vector or single number to be used for estimation of joint cdf \eqn{P(X0 \le u1, Y0 \le u2)} in the non-parametric method.
#' @param conditional A logical value. If TRUE this function will calculate the conditional cdf for the type II gap time given an interval of the type I gap time and bootstrap standard error and confidence interval at the specified confidence level. Default is FALSE.
#' @param given.interval A 1x2 vector c(v1, v2) that must be specified if conditional = TRUE. The vector indicates an interval for the type I gap time to use for estimation of the cdf of the type II gap time given this interval.
#' If given.interval = c(v1, v2) the function calculates \eqn{P(Y0 \le y | v1 \le X0 \le v2)}. The given values v1 and v2 must be in the range of gap times in the estimated marginal survival.
#' Valid values for these times are given in the "Time" column of the marginal survival data frame that results from biv.rec.np().
#' @param jointplot A logical value. If TRUE (default) this function will create a 3D plot of the joing cdf for the two gap times with pointwise large sample confidence interval at the specified confidence level.
#' @param marginalplot A logical value. If TRUE (default) this function will plot the marginal survival fuction for the type I gap times with pointwise large sample confidence interval at the specified confidence level.
#' @param condiplot A logical value. Can only be TRUE if conditional=TRUE. If TRUE this function will plot the conditional cdf with bootstrap confidence interval at the specified confidence level. Default is FALSE.
#' @param CI The level for confidence intervals for joing cdf plot, marginal plot and conditional cdf. Must be between 0.50 and 0.99 where 0.99 would give 99\% CI. Default is 0.95.

#'
#' @return Plots as specified from jointplot, marginalplot, conditional and a BivRec list object containing:
#' \itemize{
#'   \item \strong{joint.cdf:} Data frame with joint cdf and standard error for the two alternating gap times.
#'   \item \strong{marginal.survival:} Data frame with marginnal survival for the first gap time and standard error.
#'   \item \strong{conditional.cdf:} Data frame with conditional cdf, bootstrap standard error and bootstrap confidence interval.
#'   \item \strong{formula:} The formula used to specify components of bivariate recurrent response.
#'   \item \strong{ai:} The function of censoring time used as weights.
#' }
#'
#' @details
#' \strong{ai} indicates a real non-negative function of censoring times to be used as weights in the non-parametric method. This variable can take on values of 1 or 2 which indicate:
#' \itemize{
#' \item 1: the weights are simply 1 for all subjects \eqn{a(C_i) = 1} (default).
#' \item 2: the weight for each subject is his/her censoring time \eqn{a(C_i) = C_i}.
#' }
#' For further information, see Huang and Wang (2005).
#'
#' @references
#' Huang CY, Wang MC (2005). Nonparametric estimation of the bivariate recurrence time distribution. Biometrics, 61: 392-402.
#' \url{doi.org/10.1111/j.1541-0420.2005.00328.x}
#'
#' @export
#' @keywords biv.rec.np
#'
#' @examples
#'\dontrun{
#' library(BivRec)
#'#Simulate bivariate alternating recurrent event data
#' set.seed(1234)
#' biv.rec.data <- biv.rec.sim(nsize=150, beta1=c(0.5,0.5), beta2=c(0,-0.5), tau_c=63, set=1.1)
#' #Apply the non-parametric method of Huang and Wang (2005) and
#' #visualize marginal and conditional results.
#' par(mfrow=c(1,2))
#' nonpar.result <- biv.rec.np(formula = id + epi + xij + yij + d1 + d2 ~ 1,
#'           data=biv.rec.data, ai=1, u1 = c(2, 5, 10, 20), u2 = c(1, 5, 10, 15),
#'           conditional = TRUE, given.interval=c(0, 10), jointplot=TRUE,
#'           marginalplot = TRUE, condiplot = TRUE)
#' head(nonpar.result$joint.cdf)
#' head(nonpar.result$marginal.survival)
#' head(nonpar.result$conditional.cdf)
#' }

biv.rec.np <- function(formula, data, CI, ai, u1, u2, conditional, given.interval,jointplot,marginalplot,condiplot){

  if (missing(ai)) {ai<-1}
  if (missing(jointplot)) {jointplot <- TRUE}
  if (missing(marginalplot)) {marginalplot <- TRUE}
  if (missing(conditional)) {conditional <- FALSE}
  if (missing(CI)) {CI <- 0.95}
  if (CI > 0.99) {
    print("Error CI > 0.99")
    stop()} else {
      if (CI<0.5) {
        print("Error CI < 0.5")
        stop()
      }
    }

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
  episode <- eval(parse(text = names[2]))
  xij <- eval(parse(text = names[3]))
  yij <- eval(parse(text = names[4]))
  if (length(names)==6) {
    c_indicatorX <- eval(parse(text = names[5]))
    c_indicatorY <- eval(parse(text = names[6]))
  } else {
    c_indicatorX <- rep(1, length(xij))
    c_indicatorY <- eval(parse(text = names[5]))
  }
  covariates <- rep(1, length(identifier))
  method <- "Non-Parametric"
  condgx <- FALSE

  if (missing(u1)) {u1 <- round(seq(quantile(xij, probs = 0.4), max(xij), length.out=5))}
  if (missing(u2)) {u2 <- round(seq(quantile(yij, probs = 0.4), max(yij), length.out=4))}
  temp <- rep(u1, each = length(u2))
  temp2 <- rep(u2, length(u1))
  u <- cbind(u1=temp, u2=temp2)

  ###Send to biv.rec.reformat and complete analysis
  new_data <- biv.rec.reformat(identifier, xij, yij, c_indicatorY, c_indicatorX, episode, covariates, method, ai, condgx, data)
  print("Estimating joint cdf and marginal survival")
  res1 <- nonparam.cdf(fit_data=new_data$forcdf, u, ai, CI)
  res2 <- nonparam.marginal(new_data$formarg, CI)

  if (jointplot==TRUE) {plot.joint.cdf(list(cdf = res1, formula=formula, data=data), CI)}
  if (marginalplot==TRUE) {marg.surv.plot(list(marginal.survival = res2, formula=formula, data=data), CI)}
  if (conditional == FALSE) {
    final.result <- list(joint.cdf = res1, marginal.survival = res2, formula=formula, ai=ai)
  } else {
    if (missing(given.interval)) {
      print("Error for conditional calculation given.interval argument missing.")
      final.result <- list(cdf = res1, marginal.survival = res2, formula=formula, data = data, ai=ai)
    } else {
      partial.result <- list(cdf = res1, marginal.survival = res2, formula=formula, data = data, ai=ai, new_data=new_data)
      res3 <- nonparam.conditional(partial.result, given.interval, CI, condiplot)
      final.result <- list(joint.cdf = res1, marginal.survival = res2, conditional.cdf = res3, formula=formula, ai=ai)
    }
  }

  return(final.result)
}

