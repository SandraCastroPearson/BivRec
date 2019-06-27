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
#' lee_reg <- bivrecReg(bivrecSurv(id, epi, xij, yij, d1, d2) ~ a1 + a2,
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
#'chang_reg <- bivrecReg(bivrecSurv(id, epi, xij, yij, d1, d2) ~ a1 + a2,
#'                     data = bivrec_data, method = "Chang")
#'summary(chang_reg)
#'# To apply Chang (2004) method use method="Chang"
#'}
#'
#' @keywords bivrecReg

bivrecReg =function(formula, data, method){

  call = match.call()

  #Manage missing information by method
  formula_ref = formula

  if (missing(method)) {method <- "Lee.et.al"}

  if (!missing(data)) {
    variables <- all.vars(formula_ref)
    names <- paste("data$", variables, sep="")
    numbercov <- length(variables)- 6
    dat2 <- data.frame(eval(parse(text = names[1])), eval(parse(text = names[2])))

    for (i in 1:numbercov) {
      dat2[, 2+i] <- eval(parse(text = names[6+i]))
    }
    which_missing <- unique(unlist(apply(dat2, 2, function(x) which(is.na(x)==TRUE))))
    data_ref <- data

    if (length(which_missing)!=0) {
      data <- data[-which_missing, ]
    }

    response <- eval(formula[[2]], data)
    formula[[2]] <- NULL

  } else {stop("data argument missing")}

  if (!inherits(response, "bivrecSurv")) stop("Response must be a bivrecSurv object")
  if (ncol(model.matrix(formula, data)) == 1) {stop("This is a non-parametric analysis use non-parametric functions")}

  #Lee et all Method
  if (method == "Lee.et.al") {
    predictors <- data.frame(id = response$id_ref, model.matrix(formula, data)[,-1])
    colnames(predictors) <-  c("id", colnames(model.matrix(formula, data))[-1])
    cov_names <- colnames(predictors)[-1]

    amat = NULL
    for (i in unique(predictors$id)) {
      tmp=predictors[predictors$id==i, ]
      amat=rbind(amat,tmp[1, 2:ncol(predictors)])
    }
    amat = as.matrix(amat)

    new_response = response$data4Lreg

    if (ncol(amat)==1) {
      results <- list(call = call,
                      leefit = leeall_univariate(response=new_response, amat, cov_names, SE="TRUE"),
                      formula=formula_ref, method="Lee.et.al",
                      data = list(response=new_response, predictors = amat, original = data_ref))
    } else {
      results <- list(call = call,
                      leefit = leeall_multivariate(response=new_response, amat, cov_names, SE="TRUE"),
                      formula = formula_ref, method="Lee.et.al",
                      data = list(response=new_response, predictors = amat, original = data_ref))}

    ### Chang Method
  } else {
    predictors <- data.frame(id = response$id_ref, epi = na.omit(dat2)[,2], model.matrix(formula, data)[,-1])
    colnames(predictors) <-  c("id", "epi", colnames(model.matrix(formula, data))[-1])
    cov_names <- colnames(predictors)[-c(1,2)]
    new_data <- merge(response$data4Creg, predictors, c("id","epi"))
    new_data <- new_data[order(new_data$id, decreasing = FALSE),]

    if (length(cov_names)==1) {
      results <- list(call = call,
                      chang_fit = chang_univariate(new_data, cov_names, SE="TRUE"),
                      formula = formula_ref, method = "Chang",
                      data = list(new_data=new_data, original = data_ref))
    } else {
      results <- list(call = call,
                      chang_fit = chang_multivariate(new_data, cov_names, SE="TRUE"),
                      formula=formula_ref, method="Chang",
                      data = list(new_data=new_data, original = data_ref))}
  }

  class(results) <- "bivrecReg"
  return(results)
}

