#' Non-Parametric Accelerated Failure Time Analysis of Bivariate Alternating Recurrent Event Gap Time Data
#'
#' @description
#' This function allows the user to obtain the joint, conditional and marginal cumulative distribution functions.
#' See details for the estimation methods provided.
#'
#' @importFrom stats model.frame
#' @importFrom stats na.omit
#' @importFrom stats quantile
#' @importFrom stats model.matrix
#'
#' @param x A response object of the \code{bivrecSurv} class.
#' @param data A data frame that includes all the covariates listed in the formula.
#' @param CI The level for confidence intervals for joint cdf plot, marginal plot and conditional cdf. Must be between 0.50 and 0.99, where 0.99 would give 99\% CI. Default is 0.95.
#' @param ai A real non-negative function of censoring time. See details.
#' @param u1 A vector or single number to be used for estimation of joint cdf \eqn{P(X0 \le u1, Y0 \le u2)} in the non-parametric method.
#' @param u2 A vector or single number to be used for estimation of joint cdf \eqn{P(X0 \le u1, Y0 \le u2)} in the non-parametric method.
#' @param conditional A logical value. If TRUE, this function will calculate the conditional cdf for the type II gap time given an interval of the type I gap time and the bootstrap standard error and confidence interval at the specified confidence level. Default is FALSE.
#' @param given.interval A vector c(v1, v2) that must be specified if conditional = TRUE. The vector indicates an interval for the type I gap time to use for estimation of the cdf of the type II gap time given this interval.
#' If given.interval = c(v1, v2), the function calculates \eqn{P(Y0 \le y | v1 \le X0 \le v2)}. The given values v1 and v2 must be in the range of gap times in the estimated marginal survival.
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
#'
#' @examples
#' library(BivRec)
#'# Simulate bivariate alternating recurrent event data
#' set.seed(1234)
#' bivrec_data <- simulate(nsize=150, beta1=c(0.5,0.5), beta2=c(0,-0.5), tau_c=63, set=1.1)
#' bdat <- is.bivrecSurv(bivrec_data)
#' npresult <- bivrecNP(bdat,ai=1, u1 = c(2, 5, 10, 20), u2 = c(1, 5, 10, 15),conditional = FALSE, given.interval=c(0, 10))
#'
#' \dontrun{
#' #This is an example with longer runtime (it runs the conditional graph)
# npresult2 <- bivrecNP(bdat,ai=1, u1 = c(2, 5, 10, 20), u2 = c(1, 5, 10, 15),conditional = TRUE, given.interval=c(0, 10))
#'}
#'
#' @keywords bivrecNP

bivrecNP <- function(x, CI, ai, u1, u2, conditional, given.interval){
  if (!is.bivrecSurv(x)) stop("Response must be a bivrecSurv class")
  if (missing(ai)) {ai<-1}
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

  ### Dataframe with id, epi, xij, yij, d1, d2, zij and ci. This part will be moved to bivrecSurv
  data <- x$df

  iden <- data$id #move this to BivrecSurv object
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
  data <- data[,-which(colnames(data)=="id")]
  colnames(data)[ncol(data)] = "id" #this just moved the id column to the end to match up with new.id values

  if (missing(u1)) {u1 <- round(seq(quantile(data$xij, probs = 0.4), max(data$xij), length.out=5))}
  if (missing(u2)) {u2 <- round(seq(quantile(data$yij, probs = 0.4), max(data$yij), length.out=4))}
  temp <- rep(u1, each = length(u2))
  temp2 <- rep(u2, length(u1))
  u <- cbind(u1=temp, u2=temp2)

  print("Estimating joint cdf and marginal survival")
  if (ai==1) {
  cdf_res1 <- nonparam.cdf(x$dat4np1$forcdf, u, ai, CI) #result for joint cdf if ai=1
  marg_res1 <- nonparam.marginal(x$dat4np1$formarg, CI) #result for marginal if ai=1
  if (conditional == FALSE) {
    final.result <- list(joint.cdf = cdf_res1, marginal.survival = marg_res1, ai=ai)
  } else {
    if (missing(given.interval)) {
      print("Error for conditional calculation given.interval argument missing.")
      final.result <- list(cdf = cdf_res1, marginal.survival = marg_res1, ai=ai)
    } else {
      partial.result <- list(cdf = cdf_res1, marginal.survival = marg_res1, data = data, ai=ai, new_data=x$dat4np1) #took out formula as a param
      ccdf_res1 <- nonparam.conditional(partial.result, given.interval, CI,x$df$yij) #took out condiplot as a param
      final.result <- list(joint.cdf = cdf_res1, marginal.survival = marg_res1, conditional.cdf = ccdf_res1,ai=ai)
    }
  }
  }
  if (ai==2) { #ai=2
  cdf_res2 <- nonparam.cdf(x$dat4np2$forcdf, u, ai, CI) #result for joint cdf if ai=2
  marg_res2 <- nonparam.marginal(x$dat4np2$formarg, CI) #result for marginal if ai=2
  if (conditional == FALSE) {
    final.result <- list(joint.cdf = cdf_res2, marginal.survival = marg_res2, ai=ai)
  } else {
    if (missing(given.interval)) {
      print("Error for conditional calculation given.interval argument missing.")
      final.result <- list(joint.cdf = cdf_res2, marginal.survival = marg_res2, ai=ai)
    } else {
      partial.result <- list(cdf = cdf_res1, marginal.survival = marg_res2, data = data, ai=ai, new_data=x$dat4np2) #took out formula as a param
      ccdf_res2 <- nonparam.conditional(partial.result, given.interval, CI,x$df$yij) #took out condiplot as a param
      final.result <- list(joint.cdf = cdf_res2, marginal.survival = marg_res2, conditional.cdf = ccdf_res2,ai=ai)
    }
  }
  }
  class(final.result)<-"bivrecNP"
  final.result$CI <- CI
  final.result$given.interval<-given.interval
  final.result$conditional <- conditional #boolean indicator
  #final.result$ygrid <-ccdf_res2$ygrid #original response data from bivrecSurv object
  return(final.result) #Essentially the bivrecNP object provides the data (the new_id stuff), CI, results for all 3 (if conditional=true),
  #the conditional indicator
}

is.bivrecNP <- function(x) inherits(x, "bivrecNP")
