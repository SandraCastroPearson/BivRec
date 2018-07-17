#' Analysis of Bivariate Recurrent Data
#'
#' @description
#' This function allows the user to:
#' 1) Evaluate covariate effects on two alternating states/gap times.
#' 2) Find the joint cdf, marginal survival and conditional cdf for two alternating states/gap times.
#'
#' @importFrom stats model.frame
#' @importFrom stats na.omit
#' @importFrom stats quantile
#'
#' @param formula Formula indicating model to fit in the form: \strong{ID + xij + yij + episode + c_indicatorY + c_indicatorX ~ covariates}.
#'  where
#'  \strong{ID:} is a vector of subjects' unique identifier.
#'  \strong{xij:} is a vector with the lengths of time spent in event of type X for individual i in episode j.
#'  \strong{yij:} is a vector with the lengths of time spent in event of type Y for individual i in episode j.
#'  \strong{episode:} is a vector indicating the observation or episode j for a subject i.
#'  \strong{c_indicatorY:} is an vector of indicators with values of 0 for the last episode for subject i and 1 otherwise. A subject with only one episode will have one 0.
#'  \strong{c_indicatorX:} is an optional vector of indicators with values of 0 if the last episode for subject i occurred for event of type X or 1 otherwise. A subject with only one episode could have only a 1 (if he was censored at event Y) or a 0 (if he was censored at event X).
#'  \strong{covariates:} the names of the covariates in the form covariate_1 + ... + covariate_N. If using Non-parametric method replaces the names with a 1.
#' @param data A data frame that includes all the vectors/covariates listed in the formula
#' @param method A string indicating which method to use to estimate effect of covariates (and corresponding 95\% CI when specified) or cdf. Options are: "Lee.all" (default), "Chang", or "Non-Parametric". For details on methods see vignette.
#' @param CI A logical variable. If TRUE a 95\% Confidence Interval for the effect estimates is calculated. Default is FALSE.
#' @param ai A real non-negative function of censoring time to be used as weights in the Non-Parametric Method. Current options are 1 such that ai = f(censoring times) = 1 or the default of 2 such that ai = f(censoring times) = censoring times.
#' @param u1 A vector or single number to be used for Estimation of joint CDF P(xij <= u1, yij <= u2) in the Non-Parametric Method. (see vignette for further details)
#' @param u2 A vector or single number to be used for Estimation of joint CDF P(xij <= u1, yij <= u2) in the Non-Parametric Method. (see vignette for further details)
#' @param condgx Automatically passed - not for user
#'
#' @return A list with dataframes of results and some of the parameters used to produce them.
#' @export
#'
#' @keywords biv.rec.fit
#'
#' @examples
#' \dontrun{
#' library(BivRec)
#' set.seed(1234)
#' sim.data <- data.sim(nsize=300, beta1=c(0.5,0.5), beta2=c(0,-0.5), cr=63, sg2=0.5, set=1.1)
#' multiple.lee <- biv.rec.fit(id + xij + yij + epi + d2 + d1 ~ a1 + a2,
#'                 data=sim.data, method="Lee.all", CI=FALSE)
#' multiple.lee
#' nonpar <- biv.rec.fit(id + xij + yij + epi + d2 + d1 ~ 1,
#'           data=sim.data, method="Non-Parametric", ai=1)
#' nonpar
#' }
#'

biv.rec.fit <- function(formula, data, method, CI, ai, u1, u2, condgx){

  #Manage missing information by method
  if(missing(condgx)) {condgx <- FALSE}
  if(missing(method)) {method <- "Lee.all"}
  if(missing(CI)) {CI <- FALSE}

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
  if (method != "Non-Parametric") {
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
  } else {
      names <- paste("data$", variables, sep="")
      identifier <- eval(parse(text = names[1]))
      xij <- eval(parse(text = names[2]))
      yij <- eval(parse(text = names[3]))
      episode <- eval(parse(text = names[4]))
      c_indicatorY <- eval(parse(text = names[5]))
      ifelse(length(names)==6, c_indicatorX <- eval(parse(text = names[6])), c_indicatorX <- rep(1, length(xij)))
      covariates <- rep(1, length(identifier))
      cov_names <- "No Covariates"
  }

  #### Manage data vectors, fit correct model, return results ###
  #Proposed Method
  if (method == "Lee.all") {
    new_data <- biv.rec.reformat(identifier, xij, yij, c_indicatorY, c_indicatorX, episode, covariates, method, ai=2, condgx=FALSE, data)
      if (length(cov_names)==1) {
        results <- semi.param.univariate(new_data, cov_names, CI)
      } else {
        results <- semi.param.multivariate(new_data, cov_names, CI)}
    return(list(covariate.effects = results, formula=formula, data = data))
  #Chang Method
  } else {
    if (method == "Chang") {
    new_data <- biv.rec.reformat(identifier, xij, yij, c_indicatorY, c_indicatorX, episode, covariates, method, ai=2, condgx=FALSE, data)
        if (length(cov_names)==1) {
          results <- chang.univariate(new_data, cov_names, CI)
        } else {
          results <- chang.multivariate(new_data, cov_names, CI)
        }
    return(list(covariate.effects = results, formula=formula, data = data))
  #Non-Parametric
    } else {
        if (missing(ai)) {ai<-1}
        if (missing(u1)) {u1 <- round(seq(quantile(xij, probs = 0.4), max(xij), length.out=5))}
        if (missing(u2)) {u2 <- round(seq(quantile(yij, probs = 0.4), max(yij), length.out=4))}
        temp <- rep(u1, each = length(u2))
        temp2 <- rep(u2, length(u1))
        u <- cbind(u1=temp, u2=temp2)
        new_data <- biv.rec.reformat(identifier, xij, yij, c_indicatorY, c_indicatorX, episode, covariates, method, ai, condgx, data)
        res1 <- nonparam.cdf(fit_data=new_data$forcdf, u, ai)
        if (condgx==FALSE) {
          res2 <- nonparam.marginal(new_data$formarg)
          return(list(cdf = res1, marginal.survival = res2, formula=formula, data = data, ai=ai, new_data=new_data))
        } else {return(list(cdf = res1, formula=formula, data = data, ai=ai, new_data=new_data))}
    }
  }
}

