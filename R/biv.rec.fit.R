#' Semi-Parametric Accelerated Failure Time Analysis of Bivariate Alternating Recurrent Event Data
#'
#' @description
#' This function allows the user to evaluate covariate effects on two alternating recurrent events gap times (referred as type I and type II gap times)
#' under the assumption that the two gap times follow accelerated failure time (AFT) models. Two estimation methods are provided.
#'
#' @importFrom stats model.frame
#' @importFrom stats na.omit
#' @importFrom stats quantile
#'
#' @param formula Formula with six variables indicating the bivariate alternating gap time response on the left and the ~ operator and the covariates on the right.
#' The six variables on the left must have the same length and be given as follows \strong{ID + episode +  xij + yij + delta_x + delta_y ~ covariates} where
#' \itemize{
#'   \item ID: a vector of subjects' unique identifier which can be numeric or character.
#'   \item episode: is a vector indicating the episode of the bivariate alternating gap time pairs.
#'   \item xij: is a vector with the lengths of the type I gap times.
#'   \item yij: is a vector with the lengths of the type II gap times.
#'   \item delta_x: a vector of indicators with values
#'   \itemize{
#'       \item 0 for the last episode for subject i (j=m_i) if subject was censored during period xij.
#'       \item 1 otherwise.
#'      }
#'   A subject with only one episode could have a 0 if he was censored during period xi1 or 1 if he was censored during period yi1.
#'   \item delta_y: a vector of indicators with values
#'   \itemize{
#'       \item 0 for the last episode of subject i (j=m_i).
#'       \item 1 otherwise.
#'      }
#'   A subject with only one episode will have one 0.
#'   \item covariates: the names of the covariates in the form covariate_1 + ... + covariate_p.
#' }
#' @param data A data frame that includes all the vectors/covariates listed in the formula above.
#' @param method A string indicating which method to use to estimate effects of the covariates. See details.
#' @param CI Level for confidence interval must be between 0.50 and 0.99, 0.99 would give 99\% CI. The default is 0.95.
#'
#' @return A BivRec list object containing:
#' \itemize{
#'   \item \strong{covariate.effects:} Data Frame summarizing effects of the covariates including point estimate, standard error and confidence interval.
#'   \item \strong{formula:} The formula used to specify components of bivariate recurrent response and covariates.
#' }
#'
#' @details
#' Two different estimation methods are available:
#' \itemize{
#' \item  method = "Lee.et.all" (default) is a U-statistics-based smooth estimating function approach. See Lee CH, Huang C-Y, Xu G, Luo X (2017) reference for further details.
#' \item  method = "Chang" is a rank-based estimating function approach.  See Chang (2004) reference for further details.
#' Note that following the Chang method, the variances of the estimated regression coefficients are approximated using the resampling techniques developed by Parzen, Wei and Ying (1994).
#' This approximation requires extensive computing time for a relatively small sample size. In addition, using the Chang method does not guarantee convergence for the estimation of the coefficients.
#' }
#'
#' @export
#'
#' @keywords biv.rec.fit
#'
#' @references
#' \enumerate{
#' \item Chang, Shu-Hui. Estimating Marginal Effects in Accelerated Failure Time Models for Serial Sojourn Times Among Repeated Events. Lifetime data analysis. 2004;
#' \url{https://doi.org/10.1023/B:LIDA.0000030202.20842.c9}
#' \item Lee CH, Huang C-Y, Xu G., Luo X. Semiparametric regression analysis for alternating recurrent event data. Statistics in Medicine. 2017.
#' \url{https://doi.org/10.1002/sim.7563}
#' \item Parzen MI, Wei LJ, Ying Z. A resampling method based on pivotal estimating functions. Biometrika. Volume 81. Issue 2. June 1994. Pages 341-350.
#' \url{https://doi.org/10.1093/biomet/81.2.341}
#' }
#'
#' @examples
#' library(BivRec)
#'#Simulate Bivariate Alternating Recurrent Event Data
#' set.seed(1234)
#' sim.data <- biv.rec.sim(nsize=150, beta1=c(0.5,0.5), beta2=c(0,-0.5), cr=63, sg2=0.5, set=1.1)
#' #Apply Lee et All (2017) method using multiple covariates.
#' multiple.lee <- biv.rec.fit(formula = id + epi + xij + yij + d1 + d2 ~ a1 + a2,
#'                 data=sim.data, method="Lee.et.all", CI=0.99)
#' multiple.lee
#'

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
      if (length(names)==6) {
        c_indicatorX <- eval(parse(text = names[5]))
        c_indicatorY <- eval(parse(text = names[6]))
      } else {
        c_indicatorX <- rep(1, length(xij))
        c_indicatorY <- eval(parse(text = names[5]))
      }
      covariates <- as.matrix(covariates)

  #### Manage data vectors, fit correct model, return results ###
  new_data <- biv.rec.reformat(identifier, xij, yij, c_indicatorY, c_indicatorX, episode, covariates, method, ai=2, condgx=FALSE, data)
    #Proposed Method
    if (method == "Lee.et.all") {
        if (length(cov_names)==1) {results <- semi.param.univariate(new_data, cov_names, CI)
        } else {results <- semi.param.multivariate(new_data, cov_names, CI)}
    #Chang Method
    } else {
        if (length(cov_names)==1) {results <- chang.univariate(new_data, cov_names, CI)
        } else {results <- chang.multivariate(new_data, cov_names, CI)}
    }
  return(list(covariate.effects = results, formula=formula))
}

