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
#' @param formula A formula with six variables indicating the bivariate alternating gap time response on the left of the ~ operator and the covariates on the right.
#' The six variables on the left must have the same length and be given as \strong{ID + episode +  xij + yij + delta_x + delta_y ~ covariates}, where
#' \itemize{
#'   \item ID: A vector of subjects' unique identifier which can be numeric or character.
#'   \item episode: A vector indicating the episode of the bivariate alternating gap time pairs, e.g.: 1, 2, ..., m_i where m_i indicates the last episode for subject i.
#'   \item xij: A vector with the lengths of the type I gap times.
#'   \item yij: A vector with the lengths of the type II gap times.
#'   \item delta_x: A vector of indicators with values
#'   \itemize{
#'       \item 0 for the last episode for subject i (m_i) if subject was censored during period xij.
#'       \item 1 otherwise.
#'      }
#'   A subject with only one episode (m_i=1) could have a 0 if he was censored during period xi1 or 1 if he was censored during period yi1. If delta_x is not provided estimation proceeds with the assumption that no subject was censored during period xij.
#'   \item delta_y: A vector of indicators with values
#'   \itemize{
#'       \item 0 for the last episode of subject i (m_i).
#'       \item 1 otherwise.
#'      }
#'   A subject with only one episode (m_i=1) will have one 0.
#'   \item covariates: the names of the covariates in the form covariate_1 + ... + covariate_p.
#' }
#' @param data A data frame that includes all the vectors/covariates listed in the formula above.
#' @param method A string indicating which method to use to estimate effects of the covariates. See details.
#' @param CI The level to be used for confidence intervals. Must be between 0.50 and 0.99, where 0.99 would give 99\% CI. The default is 0.95. CI=NULL gives point estimates without confidence intervals.
#'
#' @return A BivRec list object containing:
#' \itemize{
#'   \item \strong{covariate.effects:} A data frame summarizing effects of the covariates including the point estimate, standard error and confidence interval.
#'   \item \strong{formula:} The formula used to specify components of bivariate recurrent response and covariates.
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
#' biv.rec.data <- biv.rec.sim(nsize=100, beta1=c(0.5,0.5), beta2=c(0,-0.5), tau_c=63, set=1.1)
#' # Apply Lee C, Huang CY, Xu G, Luo X (2017) method using one covariate
#' fit.lee <- biv.rec.fit(formula = id + epi + xij + yij + d1 + d2 ~ a1,
#'                 data=biv.rec.data, method="Lee.et.al", CI=NULL)
#' fit.lee$covariate.effects
#' \dontrun{
#'
#' #This is an example with longer runtime.
#'
#' library(BivRec)
#'# Simulate bivariate alternating recurrent event data
#' set.seed(1234)
#' biv.rec.data <- biv.rec.sim(nsize=150, beta1=c(0.5,0.5), beta2=c(0,-0.5), tau_c=63, set=1.1)
#'
#' # Apply Lee C, Huang CY, Xu G, Luo X (2017) method using multiple covariates
#' # and 99% confidence intervals.
#' fit.lee <- biv.rec.fit(formula = id + epi + xij + yij + d1 + d2 ~ a1 + a2,
#'                 data=biv.rec.data, method="Lee.et.al", CI=0.99)
#' fit.lee$covariate.effects
#'
#' }

#'# To apply Chang (2004) method use method="Chang"
#'
#'
#' @export
#'
#' @keywords biv.rec.fit
biv.rec.fit <- function(formula, data, method, CI){

  #Manage missing information by method

  if(missing(method)) {method <- "Lee.et.all"}
  if(missing(CI)) {CI <- 0.95}

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
      covariates <- model.frame(formula, data, na.action = NULL)[-1]
      if (ncol(covariates)==0) {
        print(paste("Error: Covariates needed for", method, "method", sep = " "))
        return()
      }
      cov_names <- colnames(covariates)
      removable_from_all <- seq(length(variables), length(variables) - length(cov_names) +1, -1)
      names <- paste("data$", variables[-removable_from_all], sep="")
      identifier <- eval(parse(text = names[1]))
      episode <- eval(parse(text = names[2]))
      xij <- eval(parse(text = names[3]))
      yij <- eval(parse(text = names[4]))
      c_indicatorX <- eval(parse(text = names[5]))
      c_indicatorY <- eval(parse(text = names[6]))
      covariates <- as.matrix(covariates)

  #### Manage data vectors, fit correct model, return results ###
  new_data <- biv.rec.reformat(identifier, xij, yij, c_indicatorY, c_indicatorX, episode, covariates, method, ai=2, condgx=FALSE, data)
    #Proposed Method
    if (method == "Lee.et.al") {
        if (length(cov_names)==1) {results <- semi.param.univariate(new_data, cov_names, CI)
        } else {results <- semi.param.multivariate(new_data, cov_names, CI)}
    #Chang Method
    } else {
        if (length(cov_names)==1) {results <- chang.univariate(new_data, cov_names, CI)
        } else {results <- chang.multivariate(new_data, cov_names, CI)}
    }
  return(list(covariate.effects = results, formula=formula))
}

