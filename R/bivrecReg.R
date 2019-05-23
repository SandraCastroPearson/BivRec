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
#' @rdname bivrecReg
#' @export
#'
#' @examples
#' library(BivRec)
#'# Simulate bivariate alternating recurrent event data
#' set.seed(1234)
#' #bivrec_data <- biv.rec.sim(nsize=150, beta1=c(0.5,0.5), beta2=c(0,-0.5), tau_c=63, set=1.1)
#'
#' bivrec_data <- simulate(nsize=150, beta1=c(0.5,0.5), beta2=c(0,-0.5), tau_c=63, set=1.1)
#' # Apply Lee C, Huang CY, Xu G, Luo X (2017) method using one covariate
#' lee_reg <- bivrecReg(bivrecSurv(id, epi, xij, yij, d1, d2) ~ a1 + a2,
#'                     data = bivrec_data, method="Lee.et.al")
#' summary(lee_reg)
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

  call = match.call()

  #Manage missing information by method
  formula_ref = formula
  if (missing(method)) {method <- "Lee.et.al"}
  if (!missing(data)) {response <- eval(formula[[2]], data)
      } else {stop("data argument missing")}
  if (!is.bivrecSurv(response)) stop("Response must be a bivrecSurv object")
  formula[[2]] <- NULL

  if (ncol(model.matrix(formula, data)) == 1) {stop("This is a non-parametric analysis use non-parametric functions")}

  #Lee et all Method
   if (method == "Lee.et.al") {
     predictors <- data.frame(id = response$id_ref, model.matrix(formula, data)[,-1])
     colnames(predictors) <-  c("id", colnames(model.matrix(formula, data))[-1])
     cov_names <- colnames(predictors)[-1]
     missing <- unlist(apply(predictors[,-1], 2, function(x) which(is.na(x))))

     if (length(missing)!=0) {    #need to test missing fix
       predictors <- predictors[-missing,]
       temp1 = list(xmat = response$dat4Lreg$xmat, ymat = response$data4Lreg$ymat, zmat = response$dat4Lreg$zmat,
                    delta1 = response$dat4Lreg$delta1, delta2 = response$data4Lreg$delta2,
                    g1mat = response$dat4Lre$g1mat, g2mat = response$data4Lre$g2mat,
                    l1mat = response$dat4Lreg$l1mat, l2mat = response$data4Lreg$l2mat)
       new_response = lapply(temp1, function(x) x = x[-missing,])
       new_response$mstar = response$data4Lreg$mstar[-missing]
       new_response$ctime = response$data4Lreg$ctime[-missing]
       new_response$n = response$dat4Lreg$n - length(unique(missing))
       new_response$mc = response$data4Lreg$mc
       new_response$l1 = response$data4Lreg$l1
       new_response$l2 = response$data4Lreg$l2

     } else {new_response = response$data4Lreg}

     amat = NULL
     for (i in unique(predictors$id)) {
       tmp=predictors[predictors$id==i, ]
       amat=rbind(amat,tmp[1, 2:ncol(predictors)])
     }
     amat = as.matrix(amat)

      if (ncol(amat)==1) {
          results <- list(call = call,
                  leefit = leeall_univariate(response=new_response, amat, cov_names, SE="TRUE"),
                  formula=formula_ref, method="Lee.et.al",
                  data = list(response=new_response, predictors = amat, original = data))
        } else {
          results <- list(call = call,
                  leefit = leeall_multivariate(response=new_response, amat, cov_names, SE="TRUE"),
                  formula = formula_ref, method="Lee.et.al",
                  data = list(response=new_response, predictors = amat, original = data))}

    ### Chang Method
    } else {
      predictors <- data.frame(id = response$id_ref, model.matrix(formula, data)[,-1])
      colnames(predictors) <-  c("id", colnames(model.matrix(formula, data))[-1])
      cov_names = colnames(predictors)[-1]
      new_data <- rbind(as.data.frame(response$dat4Creg), predictors, by="id")
      orig_num <- length(unique(new_data$id))
      new_data <- na.omit(new_data)
      new_num <-length(unique(new_data$id))
      message <- paste("Original number of subjects: ", orig_num, ". Subjects for Chang Analysis: ", new_num, sep="")
      print(message)

      if (length(cov_names)==1) {
        results <- list(call = call,
          chang_fit = chang_univariate(new_data, cov_names, SE="TRUE"),
          formula = formula_ref, method = "Chang", data = list(new_data, original = data))
      } else {
        results <- list(call = call,
          chang_fit = chang_multivariate(new_data, cov_names, SE="TRUE"),
          formula=formula_ref, method="Chang", data = list(new_data, original = data))}
    }

  class(results) <- "bivrecReg"
  return(results)
}

is.bivrecReg <- function(x) inherits(x, "bivrecReg")

#Essentially the person has to give either just the dataframe and the formula
#or the exact vectors (we have to make a choice). Also they will have to use 2
#functions. One to get the y response variable and the other to put the y response variable
#into the biv.rec.fit function. Want to use S3 classes so then they can do summary(model),
#coeff(model) etc.
