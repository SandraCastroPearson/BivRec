#' Semi-Parametric Accelerated Failure Time Analysis of Bivariate Recurrent Data
#'
#' @description
#' This function allows the user to evaluate covariate effects on two alternating recurrent events/gap times using the methods outlined in: Lee CH, Huang C-Y, Xu G, Luo X. Semiparametric regression analysis for alternating recurrent event data. Statistics in Medicine. 2018;37:996-1008. \url{https://doi.org/10.1002/sim.7563}
#'
#' @importFrom stats model.frame
#' @importFrom stats na.omit
#' @importFrom stats quantile
#'
#' @param formula Formula in the form \strong{ID + xij + yij + episode + c_indicatorY + c_indicatorX ~ covariates} indicating model to fit.
#'  where
#' \itemize{
#'   \item ID: is a vector of subjects' unique identifier.
#'   \item xij: is a vector with the lengths of time spent in event of type X for individual i in episode j.
#'   \item yij: is a vector with the lengths of time spent in event of type Y for individual i in episode j.
#'   \item episode: is a vector indicating the observation or episode j for a subject i.
#'   \item c_indicatorY: is an vector of indicators with values of 0 for the last episode for subject i and 1 otherwise. A subject with only one episode will have one 0.
#'   \item c_indicatorX: is an optional vector of indicators with values of 0 if the last episode for subject i occurred for event of type X or 1 otherwise. A subject with only one episode could have only a 1 (if he was censored at event Y) or a 0 (if he was censored at event X).
#'   \item covariates: the names of the covariates in the form covariate_1 + ... + covariate_N. If using Non-parametric method replaces the names with a 1.
#' }
#' @param data A data frame that includes all the vectors/covariates listed in the formula
#' @param method A string indicating which method to use to estimate effect of covariates (and corresponding 95\% CI when specified). Options are: "Lee.all" (default) or "Chang". For details on methods see vignette.
#' @param CI Level for confidence interval between 0.50 and 0.99, 0.99 would give 99\% CI. Default is 0.95.
#'
#' @return A BivRec list object containing:
#' \itemize{
#'   \item \strong{covariate.effects:} Data Frame summarizing effects of the covariates including point estimate, standard error and confidence interval.
#'   \item \strong{formula:} The formula used to specify components of bivariate recurrent response and covariates.
#' }
#' @export
#'
#' @keywords biv.rec.fit
#'
#' @examples
#' \dontrun{
#' library(BivRec)
#' set.seed(1234)
#' sim.data <- data.sim(nsize=150, beta1=c(0.5,0.5), beta2=c(0,-0.5), cr=63, sg2=0.5, set=1.1)
#' multiple.lee <- biv.rec.fit(formula = id + xij + yij + epi + d2 + d1 ~ a1 + a2,
#'                 data=sim.data, method="Lee.all", CI=TRUE)
#' multiple.lee
#' }
#'

biv.rec.fit <- function(formula, data, method, CI){

  #Manage missing information by method

  if(missing(method)) {method <- "Lee.all"}
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
      xij <- eval(parse(text = names[2]))
      yij <- eval(parse(text = names[3]))
      episode <- eval(parse(text = names[4]))
      c_indicatorY <- eval(parse(text = names[5]))
      ifelse(length(names)==6, c_indicatorX <- eval(parse(text = names[6])), c_indicatorX <- rep(1, length(xij)))
      cov_names <- colnames(covariates)
      covariates <- as.matrix(covariates)

  #### Manage data vectors, fit correct model, return results ###
  new_data <- biv.rec.reformat(identifier, xij, yij, c_indicatorY, c_indicatorX, episode, covariates, method, ai=2, condgx=FALSE, data)
    #Proposed Method
    if (method == "Lee.all") {
        if (length(cov_names)==1) {results <- semi.param.univariate(new_data, cov_names, CI)
        } else {results <- semi.param.multivariate(new_data, cov_names, CI)}
    #Chang Method
    } else {
        if (length(cov_names)==1) {results <- chang.univariate(new_data, cov_names, CI)
        } else {results <- chang.multivariate(new_data, cov_names, CI)}
    }
  return(list(covariate.effects = results, formula=formula))
}

