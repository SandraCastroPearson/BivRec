#' Semiparametric Regression Analysis of Bivariate Alternating Recurrent Event Gap Time Data
#'
#' @description
#' This function allows the user to evaluate covariate effects on two alternating recurrent events gap times (referred as Type I and Type II gap times)
#' under the assumption that the two gap times follow accelerated failure time (AFT) models. See details for the estimation methods provided.
#'
#' @importFrom stats model.frame
#' @importFrom stats na.omit
#' @importFrom stats model.matrix
#' @importFrom dplyr filter
#'
#' @param formula A formula with a \code{bivrecSurv} object as response.
#' @param data A data frame that includes the vectors needed for the \code{bivrecSurv} response and the covariates in the formula.
#' @param method A string indicating which method to use to estimate effects of the covariates. See details.
#'
#' @return A bivrecReg object that contains:
#' \itemize{
#'   \item \verb{call}
#'   \item \verb{lee_fit} or \verb{chang_fit}
#'   \item \verb{formula}
#'   \item \verb{data}
#' }
#'
#' @details
#' Two different estimation methods are available:
#' \itemize{
#' \item  method = "Lee.et.al" (default) is a U-statistics-based smooth estimating function approach. See Lee, Huang, Xu, Luo (2018) for further details.
#' \item  method = "Chang" is a rank-based estimating function approach.  See Chang (2004) for further details.
#' Note that following the Chang method, the variances of the estimated regression coefficients are approximated using the resampling techniques developed by Parzen, Wei, Ying (1994).
#' This approximation requires extensive computing time for a relatively small sample size. In addition, using the Chang method does not guarantee convergence for the estimation of the coefficients and user may get the message, "Error: Max Iterations reached. Did not converge.".
#' }
#'
#' Related methods: \code{coef.bivrecReg}, \code{confint.bivrecReg}, \code{plot.bivrecReg}, \code{print.bivrecReg}, \code{summary.bivrecReg}, \code{vcov.bivrecReg}.
#' @references
#' \enumerate{
#' \item Chang S-H. (2004). Estimating marginal effects in accelerated failure time models for serial sojourn times among repeated events. Lifetime Data Analysis, 10: 175-190.
#' \url{https://doi.org/10.1023/B:LIDA.0000030202.20842.c9}
#'
#' \item Lee CH, Huang CY, Xu G, Luo X. (2018). Semiparametric regression analysis for alternating recurrent event data. Statistics in Medicine, 37: 996-1008.
#' \url{https://doi.org/10.1002/sim.7563}
#'
#' \item Parzen MI, Wei LJ, Ying Z. (1994). A resampling method based on pivotal estimating functions. Biometrika, 81: 341-350.
#' \url{https://doi.org/10.1093/biomet/81.2.341}
#' }
#'
#' @export
#'
#' @examples
#' library(BivRec)
#'# Simulate bivariate alternating recurrent event data
#' set.seed(1234)
#' bivrec_data <- simBivRec(nsize=150, beta1=c(0.5,0.5), beta2=c(0,-0.5),
#'                tau_c=63, set=1.1)
#' # Apply Lee, Huang, Xu, Luo (2018) method using two covariates
#' lee_reg <- bivrecReg(bivrecSurv(id, epi, xij, yij, d1, d2) ~ a1 + a2,
#'                      bivrec_data, "Lee.et.al")
#' summary(lee_reg)
#' confint(lee_reg, level=0.99)
#' vcov(lee_reg)
#'
#' \dontrun{
#' # Chang (2004) method using method="Chang". This is an example with longer runtime.
#' chang_reg <- bivrecReg(bivrecSurv(id, epi, xij, yij, d1, d2) ~ a1 + a2,
#'                       bivrec_data, method = "Chang")
#' summary(chang_reg)
#' confint(chang_reg, level=0.99)
#' vcov(chang_reg)
#'}
#'
#' @keywords bivrecReg

bivrecReg <- function(formula, data, method) {

  call = match.call()

  #Manage missing information by method
  formula_ref = formula
  if (missing(method)) {method <- "Lee.et.al"}
  if (missing(data)) {stop("data argument missing")}

  variables = all.vars(formula)
  ref_data <- data
  data = na.omit(data[ , colnames(data) %in% variables])

  resp <- eval(formula[[2]], data)
  if (inherits(resp, "bivrecSurv")==FALSE) stop("Response must be a bivrecSurv object")
  formula[[2]] <- NULL

  #Check if there are any covariates
  if (ncol(model.matrix(formula, data)) == 1) {stop("This is a non-parametric analysis use non-parametric functions")}

  #Lee et all Method
  if (method == "Lee.et.al") {
    predictors <- data.frame(id = resp$data4Creg$id,
                             model.matrix(formula, data)[,-1])
    colnames(predictors) <-  c("id", colnames(model.matrix(formula, data))[-1])
    cov_names <- colnames(predictors)[-1]
    amat = dplyr::filter(predictors, !(duplicated(predictors$id)))
    amat = as.matrix(amat)[,-1]

    if (ncol(amat)==1) {
      results <- list(call = call,
                      leefit = leeall_univariate(response=resp$data4Lreg, amat, cov_names, SE=TRUE),
                      formula = formula_ref, method="Lee.et.al",
                      data = list(response=resp$data4Lreg, predictors = amat, original = ref_data))
    } else {
      results <- list(call = call,
                      leefit = leeall_multivariate(response=resp$data4Lreg, amat, cov_names, SE=TRUE),
                      formula = formula_ref, method="Lee.et.al",
                      data = list(response=resp$data4Lreg, predictors = amat, original = ref_data))}


  }

  ### Chang Method
  if (method == "Chang") {

    #predictor portion
    predictors <- data.frame(id = resp$data4Creg$id, epi = resp$data4Creg$epi,
                             model.matrix(formula, data)[,-1])
    colnames(predictors) <-  c("id", "epi", colnames(model.matrix(formula, data))[-1])
    cov_names <- colnames(predictors)[-c(1,2)]
    new_data <- merge(resp$data4Creg, predictors, c("id","epi"))
    new_data <- new_data[order(new_data$id, decreasing = FALSE),]

    if (length(cov_names)==1) {
      results <- list(call = call,
                      chang_fit = chang_univariate(new_data, cov_names, SE=TRUE),
                      formula = formula_ref, method = "Chang", data = list(new_data, original = ref_data))
    } else {
      results <- list(call = call,
                      chang_fit = chang_multivariate(new_data, cov_names, SE=TRUE),
                      formula=formula_ref, method="Chang",
                      data = list(new_data = new_data, original = ref_data))}
  }

  class(results) <- "bivrecReg"
  return(results)
}
