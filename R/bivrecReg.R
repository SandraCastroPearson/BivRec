#' Semi-Parametric Accelerated Failure Time Analysis of Bivariate Alternating Recurrent Event Gap Time Data
#'
#' @description
#' This function allows the user to evaluate covariate effects on two alternating recurrent events gap times (referred as type I and type II gap times)
#' under the assumption that the two gap times follow accelerated failure time (AFT) models. See details for the estimation methods provided.
#'
#' @importFrom stats model.frame
#' @importFrom stats na.omit
#' @importFrom stats quantile
#'
#' @param formula A formula with a bivrecSurv object as response.
#' @param data A data frame that includes all the covariates listed in the formula.
#' @param method A string indicating which method to use to estimate effects of the covariates. See details.
#'
#' @return A bivrecReg object containing:
#' \itemize{
#'   \item \strong{covariate.effects:} A data frame summarizing effects of the covariates including the point estimate, standard error and confidence interval.
#'   \item \strong{formula:} The formula used to specify the model.
#' }
#'
#' @details
#' Two different estimation methods are available:
#' \itemize{
#' \item  method = "Lee.et.al" (default) is a U-statistics-based smooth estimating function approach. See Lee CH, Huang C-Y, Xu G, Luo X (2017) for further details.
#' \item  method = "Chang" is a rank-based estimating function approach.  See Chang (2004) for further details.
#' Note that following the Chang method, the variances of the estimated regression coefficients are approximated using the resampling techniques developed by Parzen, Wei and Ying (1994).
#' This approximation requires extensive computing time for a relatively small sample size. In addition, using the Chang method does not guarantee convergence for the estimation of the coefficients.
#' }
#'
#' @references
#' \enumerate{
#' \item Chang S-H. (2004). Estimating marginal effects in accelerated failure time models for serial sojourn times among repeated events. Lifetime Data Analysis, 10: 175-190.
#' \url{https://doi.org/10.1023/B:LIDA.0000030202.20842.c9}
#'
#' \item Lee C, Huang CY, Xu G, Luo X (2017). Semiparametric regression analysis for alternating recurrent event data. Statistics in Medicine, 37: 996-1008.
#' \url{https://doi.org/10.1002/sim.7563}
#'
#' \item Parzen MI, Wei LJ, Ying Z (1994). A resampling method based on pivotal estimating functions. Biometrika, 81: 341-350.
#' \url{https://doi.org/10.1093/biomet/81.2.341}
#' }
#'
#' @examples
#' library(BivRec)
#'# Simulate bivariate alternating recurrent event data
#' set.seed(1234)
#' #bivrec_data <- biv.rec.sim(nsize=150, beta1=c(0.5,0.5), beta2=c(0,-0.5), tau_c=63, set=1.1)
#'
#' bivrec_data <- simulate(nsize=150, beta1=c(0.5,0.5), beta2=c(0,-0.5), tau_c=63, set=1.1)
#' # Apply Lee C, Huang CY, Xu G, Luo X (2017) method using one covariate
#' fit_lee <- bivrecReg(bivrecSurv(id, epi, xij, yij, d1, d2) ~ a1 + a2,
#'                     data = bivrec_data, method="Lee.et.al")
#' fit_lee$covariate.effects
#' \dontrun{
#'
#' #This is an example with longer runtime.
#'
#' library(BivRec)
#'# Simulate bivariate alternating recurrent event data
#' set.seed(1234)
#' bivrec_data <- simulate(nsize=150, beta1=c(0.5,0.5), beta2=c(0,-0.5), tau_c=63, set=1.1)
#'
#' # Apply Lee C, Huang CY, Xu G, Luo X (2017) method using multiple covariates
#' # and 99% confidence intervals.
#'fit_chang <- bivrecReg(bivrecSurv(id, epi, xij, yij, d1, d2) ~ a1 + a2,
#'                     data = bivrec_data, method = "Chang")
#'
#'# To apply Chang (2004) method use method="Chang"
#'}

#'
#' @keywords bivrecReg

bivrecReg =function(formula, data, method){

  #Manage missing information by method
  formula_ref = formula
  if (missing(method)) {method <- "Lee.et.al"}
  if (missing(data)) response <- eval(formula[[2]], parent.frame())
  if (!missing(data)) response <- eval(formula[[2]], data)
  if (!is.bivrecSurv(response)) stop("Response must be a bivrecSurv object")
  formula[[2]] <- NULL

  if (ncol(model.matrix(formula, data)) > 1) {
      predictors <- data.frame(id = response$id_ref, model.matrix(formula, data)[,-1])
      colnames(predictors) <-  c("id", colnames(model.matrix(formula, data))[-1])
      cov_names <- colnames(predictors)[-1]
      amat = NULL
      for (i in unique(predictors$id)) {
        tmp=predictors[predictors$id==i, ]
        amat=rbind(amat,tmp[1, 2:ncol(predictors)])
      }
      amat = as.matrix(amat)
  } else {
    stop("This is a non-parametric analysis use non-parametric functions")
  }

  #Lee et all Method
   if (method == "Lee.et.al") {
      ### note issue - must figure out what to do with subjects that are missing predictors
        if (ncol(amat)==1) {
          results <- list(
                  leeall_univariate(response$dat4Lreg, amat, cov_names, SE="Y"),
                  formula=formula_ref, data = list(response=response$dat4Lreg, predictors = amat))
        } else {
          results <- list(
                  leeall_multivariate(response$dat4Lreg, amat, cov_names, SE="Y"),
                  formula=formula_ref, data = list(response=response$dat4Lreg, predictors = amat))}
  #Chang Method
    } else {
      if (ncol(amat)==1) {
        results <- list(
          chang_univariate(response$dat4Creg, amat, cov_names, SE="Y"),
          formula=formula_ref, data = list(response=response$dat4Creg, predictors = amat))
      } else {
        results <- list(
          chang_multivariate(response$dat4Creg, amat, cov_names, SE="Y"),
          formula=formula_ref, data = list(response=response$dat4Creg, predictors = amat))}
    }
  return(results)
}

#' @rdname bivrecReg
#' @export

is.bivrecReg <- function(x) inherits(x, "bivrecReg")

#Essentially the person has to give either just the dataframe and the formula
#or the exact vectors (we have to make a choice). Also they will have to use 2
#functions. One to get the y response variable and the other to put the y response variable
#into the biv.rec.fit function. Want to use S3 classes so then they can do summary(model),
#coeff(model) etc.
