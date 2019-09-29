#' Nonparametric Analysis of Bivariate Alternating Recurrent Event Gap Time Data
#'
#' @description
#' This function allows the user to obtain the joint and conditional cumulative distributions and the marginal survival functions of the two types of gap times ($X_{ij}$ and $Y_{ij}$, referred to as Type I and Type I gap times) of a bivariate alternating recurrent event process.
#' See details for the estimation methods provided.
#'
#' @importFrom stats model.frame
#' @importFrom stats na.omit
#' @importFrom stats quantile
#' @importFrom stats model.matrix
#'
#' @param response A response object of the \code{bivrecSurv} class.
#' @param level The level for confidence intervals for joint cdf plot, marginal plot and conditional cdf. Must be between 0.50 and 0.99. Default is 0.95.
#' @param ai See details.
#' @param u1 A vector or single number to be used for estimation of joint cdf \eqn{P(Type I gap times \le u1, Type II gap times \le u2)} in the non-parametric method.
#' @param u2 A vector or single number to be used for estimation of joint cdf \eqn{P(Type I gap times \le u1, Type II gap times \le u2)} in the non-parametric method.
#' @param conditional A logical value. If TRUE, this function will calculate the conditional cdf for the type II gap time given an interval of the type I gap time and the bootstrap standard error and confidence interval at the specified confidence level. Default is FALSE.
#' @param given.interval A vector c(v1, v2) that must be specified if conditional = TRUE. The vector indicates an interval for the type I gap time to use for estimation of the cdf of the type II gap time given this interval.
#' If given.interval = c(v1, v2), the function calculates \eqn{P(Type II gap times \le y | v1 \le Type I gap times \le v2)}. The given values v1 and v2 must be in the range of gap times in the estimated marginal survival.
#'
#' @details
#' \strong{ai} indicates a real non-negative function of censoring times to be used as weights in the nonparametric method. This variable can take on values of 1 or 2 which indicate:
#' \itemize{
#' \item 1: the weights are simply 1 for all subjects \eqn{a(Ci) = 1} (default).
#' \item 2: the weight for each subject is his/her censoring time \eqn{a(Ci) = Ci}.
#' }
#'
#' Related methods: \code{plot.bivrecNP}, \code{head.bivrecNP}, \code{print.bivrecNP}.
#'
#' @return An object of class bivrecNP containing:
#' \itemize{
#'   \item joint_cdf
#'   \item marginal_survival
#'   \item conditional_cdf (when conditional = TRUE)
#'   \item formula
#'   \item ai
#'   \item level
#'   \item given.interval (when conditional = TRUE)
#'   \item xij, yij
#'   \item new_data
#' }
#'
#' @export
#' @examples
#' library(BivRec)
#'
#' # Simulate bivariate alternating recurrent event data
#' set.seed(1234)
#' sim_data <- simBivRec(nsize=150, beta1=c(0.5,0.5), beta2=c(0,-0.5))
#' bivrecsurv_data <- with(sim_data, bivrecSurv(id, epi, xij, yij, d1, d2))
#' npresult <- bivrecNP(response = bivrecsurv_data, ai=1,
#'                      u1 = seq(2, 25, 1), u2 = seq(1, 20, 1), level=0.99)
#' head(npresult)
#' plot(npresult)
#'
#' \dontrun{
#' #This is an example with longer runtime (it runs the conditional graph)
#'  npresult2 <- bivrecNP(response = bivrecsurv_data, ai=1,
#'                u1 = seq(2, 25, 1), u2 = seq(1, 20, 1), conditional = TRUE,
#'                given.interval=c(0, 10), level=0.99)
#'  head(npresult2)
#'  plot(npresult2)
#' }

bivrecNP <- function(response, ai, u1, u2, level, conditional, given.interval){

  x <- response

  if (!inherits(x, "bivrecSurv")) stop("Response must be a bivrecSurv class")
  if (missing(ai)) {ai <- 1}
  if (missing(conditional)) {conditional <- FALSE}
  if (missing(level)) {CI <- 0.95} else {CI = level}

  if (CI > 0.99) {stop("Error: level is higher than 0.99")} else {
    if (CI<0.5) {stop("Error: level is less than 0.5")}
  }

  xij <- x$data4Creg$xij
  yij <- x$data4Creg$yij

  if (missing(u1)) {u1 <- round(seq(quantile(xij, probs = 0.4), max(xij), length.out=5))}
  if (missing(u2)) {u2 <- round(seq(quantile(yij, probs = 0.4), max(yij), length.out=4))}
  temp <- rep(u1, each = length(u2))
  temp2 <- rep(u2, length(u1))
  u <- cbind(u1=temp, u2=temp2)

  print("Estimating joint CDF and marginal survival")

  if (ai==1) {
    new_data = x$dat4np1
    forcdf <- new_data$forcdf
    formarg <- new_data$formarg
    }

  if (ai==2) {
    new_data = x$dat4np2
    forcdf <- new_data$forcdf
    formarg <- new_data$formarg
  }

  cdf_res <- nonparam_cdf(forcdf, u, ai, CI)
  marg_res <- nonparam_marginal(formarg, CI)

  if (conditional==FALSE) {

    final_result <- list(joint_cdf = cdf_res, marginal_survival = marg_res, ai=ai,
                         xij=xij, yij=yij, new_data=new_data)

  } else {

    if (missing(given.interval)) {

      final_result <- list(joint_cdf = cdf_res, marginal_survival = marg_res, ai=ai)

    } else {

      print("Estimating conditional distribution")

      partial_result <- list(cdf = cdf_res, marginal_survival = marg_res,
                             ai=ai, new_data=new_data)

      ccdf_res <- nonparam_conditional(res=partial_result, given.interval, CI, yij)

      final_result <- list(joint_cdf = cdf_res, marginal_survival = marg_res,
                           conditional_cdf = ccdf_res, ai=ai, xij=xij, yij=yij, new_data=new_data)

      final_result$given.interval <-given.interval
    }
  }

  final_result$level <- CI
  final_result$conditional <- conditional

  class(final_result) <- "bivrecNP"
  return(final_result)
}
