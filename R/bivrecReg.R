#' Semi-Parametric Accelerated Failure Time Analysis of Bivariate Alternating Recurrent Event Gap Time Data
#'
#' @description
#' This function allows the user to evaluate covariate effects on two alternating recurrent events gap times (referred as type I and type II gap times)
#' under the assumption that the two gap times follow accelerated failure time (AFT) models. See details for the estimation methods provided.
#'
#' @importFrom stats model.frame
#' @importFrom stats na.omit
#' @importFrom stats quantile
#' @importFrom stats model.matrix
#'
#' @param formula A formula with a bivrecSurv object as response.
#' @param data A data frame that includes all the covariates listed in the formula.
#' @param method A string indicating which method to use to estimate effects of the covariates. See details.
#'
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
#' @export
#'
#' @examples
#' library(BivRec)
#'# Simulate bivariate alternating recurrent event data
#' set.seed(1234)
#' bivrec_data <- simulate(nsize=150, beta1=c(0.5,0.5), beta2=c(0,-0.5), tau_c=63, set=1.1)
#' # Apply Lee C, Huang CY, Xu G, Luo X (2017) method using one covariate
#' lee_reg <- bivrecReg(formula = bivrecSurv(id, epi, xij, yij, d1, d2) ~ a1 + a2,
#'                     data = bivrec_data, method="Lee.et.al")
#' summary(lee_reg)
#'
#' \dontrun{
#' #This is an example with longer runtime.
#'
#' library(BivRec)
#'# Simulate bivariate alternating recurrent event data
#' set.seed(1234)
#' bivrec_data <- simulate(nsize=150, beta1=c(0.5,0.5), beta2=c(0,-0.5), tau_c=63, set=1.1)
#'
#' # Apply Lee C, Huang CY, Xu G, Luo X (2017) method using multiple covariates
#' # and 99% confidence intervals.
#'chang_reg <- bivrecReg(formula = bivrecSurv(id, epi, xij, yij, d1, d2) ~ a1 + a2,
#'                     data = bivrec_data, method = "Chang")
#'summary(chang_reg)
#'# To apply Chang (2004) method use method="Chang"
#'}
#'
#' @keywords bivrecReg

bivrecReg =function(formula, data, method){

  call = match.call()

  #Manage missing information by method
  formula_ref = formula_preds = formula

  if (missing(method)) {method <- "Lee.et.al"}

  if (missing(data)) {stop("data argument missing")}
  #else continue
    #manage missingness
      data_ref <- data
      formula_preds[[2]] = NULL
      mypreds <- as.matrix(model.frame(formula_preds, data, na.action=NULL))
      which_missing <- unique(unlist(apply(mypreds, 2, function(x) which(is.na(x)==TRUE))))
      if (length(which_missing)!=0) {
        data2 <- data[-which_missing, ]
      } else {data2=data}
    #Get response
    response <- eval(formula[[2]], data2)
    if (!inherits(response, "bivrecSurv")) stop("Response must be a bivrecSurv object")
    lee_response = response$data4Lreg
    chang_response = response$data4Creg

    #Check if there are any covariates
    formula[[2]] <- NULL
    if (ncol(model.matrix(formula, data)) == 1) {stop("This is a non-parametric analysis use non-parametric functions")}

#Lee et all Method
  if (method == "Lee.et.al") {

    #predictor portion
        predictors <- cbind(id = chang_response$id, epi = chang_response$epi, model.matrix(formula, data2)[,-1])
        cov_names <- colnames(predictors)[-c(1,2)]
        amat = NULL
        for (i in unique(predictors[,1])) {
          tmp=predictors[predictors[,1]==i, ]
              if (class(tmp)=="matrix") {
                amat=rbind(amat, tmp[1, 3:ncol(predictors)])
              } else {amat=rbind(amat, tmp[3:ncol(predictors)])}
        }

    #fit portion
      if (ncol(amat)==1) {
        results <- list(call = call,
                        leefit = leeall_univariate(lee_response, amat, cov_names, SE="TRUE"))
      } else {
        results <- list(call = call,
                        leefit = leeall_multivariate(lee_response, amat, cov_names, SE="TRUE"))
      }
      results$data = list(response=lee_response, predictors = amat, original = data_ref)

  }

### Chang Method
  if (method == "Change") {

    #predictor portion
    predictors <- data.frame(id = chang_response$id, epi = chang_response$epi, model.matrix(formula, data2)[,-1])
    colnames(predictors) <-  c("id", "epi", colnames(model.matrix(formula, data))[-1])
    cov_names <- colnames(predictors)[-c(1,2)]
    new_data <- merge(chang_response, predictors, c("id","epi"))
    new_data <- new_data[order(new_data$id, decreasing = FALSE),]

    #fit portion
    if (length(cov_names)==1) {
      results <- list(call = call,
                      chang_fit = chang_univariate(new_data, cov_names, SE="TRUE"))
    } else {
      results <- list(call = call,
                      chang_fit = chang_multivariate(new_data, cov_names, SE="TRUE"))
    }
    results$data = list(new_data=new_data, original = data_ref)
  }

  results$formula=formula_ref
  results$method=method

  class(results) <- "bivrecReg"

  return(results)
}

