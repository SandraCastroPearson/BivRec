########################    plot.bivrecFormula     ########################
#' Plot Bivariate Alternating Recurrent Series by Categorical Covariates
#'
#' @description
#' This function plots bivariate recurrent event gap times by the categorical covariates (with up to 6 categories) indicated in a formula with a \verb{bivrecSurv} object as the response variable.
#'
#' @import graphics
#' @importFrom utils tail
#' @importFrom stats model.frame
#' @importFrom stats na.omit
#' @importFrom stats model.matrix
#'
#' @param x A bivrecFormula object.
#' @param y Either empty or NULL.
#' @param data Required argument when x is a formula. Should indicate the data frame that contains the vectors to create the response and the categorical covariates indicated in the given formula.
#' @param xlab Optional string with label for horizontal axis. Default is "Time".
#' @param ylab Optional string with label for vertical axis. Default is "Individual".
#' @param type Optional vector of strings to label Type I and Type II gap times. Default is c("Type I", "Type II").
#' @param ... Additional arguments to be passed to graphical methods if needed.
#'
#' @export
#'
#' @examples
#'# Simulate bivariate alternating recurrent event data
#' library(BivRec)
#' set.seed(1234)
#' bivrec_data <- simBivRec(nsize=150, beta1=c(0.5,0.5), beta2=c(0,-0.5),
#'                tau_c=63, set=1.1)
#' bivrecfoo <- bivrecFormula(bivrecSurv(id, epi, xij, yij, d1, d2) ~ a1 + a2)
#' plot(x = bivrecfoo, data = bivrec_data,
#'      type = c("In Hospital", "Out of Hospital"))
#'

plot.bivrecFormula <- function(x, y=NULL, data, type = NULL, xlab = NULL, ylab = NULL, ...) {

  if (is.bivrecFormula(x)==FALSE) {
    stop("Object must be a formula with a bivrecSurv object as response.")}

  formula <- as.formula(x)

  #check arguments for labels
  if (missing(type)) {type=c("Type I","Type II")}
  if (missing(xlab)) {xlab="Time"}
  if (missing(ylab)) {ylab="Individual"}
  main=""

  args = c(main, xlab, ylab, type)

  #pull objects out of formula
  formula_ref = formula
  if (!missing(data)) {response <- eval(formula[[2]], data)
  } else {stop("Data argument missing.")}
  if (!inherits(response, "bivrecSurv")) stop("Response must be a bivrecSurv object.")
  formula[[2]] <- NULL
  predictors <- data.frame(id = response$id_ref, model.matrix(formula, data)[,-1])
  colnames(predictors) <-  c("id", colnames(model.matrix(formula, data))[-1])
  df = as.data.frame(response$data4Creg)

  #number of levels for each predictor
  num_levs <- apply(predictors[,-1], 2, function(x) length(unique(x)))
  preds_to_delete <- which(num_levs > 6) + 1

  message1 <- paste(colnames(predictors)[preds_to_delete], " not used - either continuous or had more than 6 levels.", sep="")
  print(message1)

  predictors <- predictors[-preds_to_delete]

  if (ncol(predictors)==1) {stop("Cannot break by covariate. All covariates are continuous or have more than 6 levels.")}

  cov_names <- colnames(predictors)[-1]
  new_data <- cbind(df, predictors)
  orig_num <- length(unique(new_data$id))
  new_data <- na.omit(new_data)
  new_num <-length(unique(new_data$id))

  message <- paste("Original number of subjects: ", orig_num,
                   ". Subjects for plots: ", new_num, ".", sep="")
  print(message)

  #plot each predictor
  ploteach <- function(pred_levels, plotdat, predictor, rdim, cov_name, formula_ref) {
    par(mar=c(1,1,1,1))
    par(mfrow=c(rdim, 2))
    for (p in 1:length(pred_levels)) {
      index <- which(predictor[,2] == pred_levels[p])
      tempdat <- plotdat[index,]
      new_main = paste("For ", cov_name, " = ", pred_levels[p], sep="")
      bdat = eval(formula_ref[[2]], tempdat)

      ##EXTRACT VECTORS FOR PLOTTING FUNCTION
      parameters <- bdat$data4Creg[-(5:7)]
      colnames(parameters) <- c("id", "episode", "xij", "yij", "ci")
      ctimes <- bdat$data4Lreg$ctime
      nsubject <- bdat$data4Lreg$n

      #Plot for one level of covariate
      args2 = args
      args2[1] = new_main
      par(mar=c(5,4,4,2)+0.1)
      basicplot(parameters, ctimes, nsubject, temp=NULL, args = args2, a = 2/5, b=1/4)
    }
    if (p == length(pred_levels)) { par(mfrow=c(1, 1)) }
  }

  if (length(cov_names==1)) {
    pred_levels = unique(predictors[,2])
    plotdat = new_data[ , 1:8]
    predictor = new_data[, c(1,10)]
    rdim = ceiling(length(pred_levels) / 2)
    ploteach(pred_levels, plotdat, predictor, rdim, cov_name = cov_names[1], formula_ref)
  } else {
    for (k in 1:length(cov_names)) {
      pred_levels = unique(predictors[,k+1])
      plotdat = new_data[ , 1:8]
      predictor = new_data[, c(1,10)]
      rdim = ceiling(length(pred_levels) / 2)
      ploteach(pred_levels, plotdat, predictor, rdim, cov_name = cov_names[1], formula_ref)
    }

  }

}


#' Plot Bivariate Alternating Recurrent Series by up to six covariates after analysis using \verb{bivrecReg}.
#'
#' @description
#' This function plots bivariate recurrent event gap times from a bivrecReg object. For examples see \verb{bivrecReg}.
#'
#' @import graphics
#' @importFrom utils tail
#' @importFrom stats model.frame
#' @importFrom stats na.omit
#' @importFrom stats model.matrix
#'
#' @param x An object of class \verb{bivrecReg}.
#' @param y Either empty or NULL
#' @param main Optional string with plot title. Default is no title.
#' @param xlab Optional string with label for horizontal axis. Default is "Time".
#' @param ylab Optional string with label for vertical axis. Default is "Individual".
#' @param type Optional vector of strings to label Type I and Type II gap times. Default is c("Type 1", "Type 2").
#' @param ... arguments to be passed to graphical methods as needed.
#'
#' @noRd
#' @keywords internal
#'

plot.bivrecReg <- function (x, y = NULL, type = NULL, main = NULL, xlab = NULL, ylab = NULL, ...){

  if (!inherits(x, "bivrecReg")) stop("Object must be a bivrecReg class")
  object <- x
  formula = object$formula
  data = object$data$original
  if (missing(type)) {type=c("Type 1","Type 2")}
  if (missing(xlab)) {xlab="Time"}
  if (missing(ylab)) {ylab="Individual"}
  if (missing(main)) {main=""}

  args = c(main, xlab, ylab, type)

  #pull objects out of formula
  formula_ref = formula
  if (!missing(data)) {response <- eval(formula[[2]], data)
  } else {stop("data argument missing")}
  if (!inherits(response, "bivrecSurv")) stop("Response must be a bivrecSurv object")
  formula[[2]] <- NULL
  predictors <- data.frame(id = response$id_ref, model.matrix(formula, data)[,-1])
  colnames(predictors) <-  c("id", colnames(model.matrix(formula, data))[-1])
  df = as.data.frame(response$data4Creg)

  #number of levels for each predictor
  num_levs <- apply(predictors[,-1], 2, function(x) length(unique(x)))
  preds_to_delete <- which(num_levs > 6) + 1

  message1 <- paste(colnames(predictors)[preds_to_delete], " not used - either continuous or had more than 6 levels.", sep="")
  print(message1)

  predictors <- predictors[-preds_to_delete]

  if (ncol(predictors)==1) {stop("Cannot break by covariate. All covariates are continuous or have more than 6 levels.")}

  cov_names <- colnames(predictors)[-1]
  new_data <- cbind(df, predictors)
  orig_num <- length(unique(new_data$id))
  new_data <- na.omit(new_data)
  new_num <-length(unique(new_data$id))

  message <- paste("Original number of subjects: ", orig_num,
                   ". Subjects for plots: ", new_num, sep="")
  print(message)

  #plot each predictor
  ploteach <- function(pred_levels, plotdat, predictor, rdim, cov_name, formula_ref) {
    par(mar=c(1,1,1,1))
    par(mfrow=c(rdim, 2))
    for (p in 1:length(pred_levels)) {
      index <- which(predictor[,2] == pred_levels[p])
      tempdat <- plotdat[index,]
      new_main = paste("For ", cov_name, " = ", pred_levels[p], sep="")
      bdat = eval(formula_ref[[2]], tempdat)

      ##EXTRACT VECTORS FOR PLOTTING FUNCTION
      parameters <- bdat$data4Creg[-(5:7)]
      colnames(parameters) <- c("id", "episode", "xij", "yij", "ci")
      ctimes <- bdat$data4Lreg$ctime
      nsubject <- bdat$data4Lreg$n

      #Plot for one level of covariate
      args2 = args
      args2[1] = new_main
      basicplot(parameters, ctimes, nsubject, temp=NULL, args = args2, a = 2/6, b=1/5, c=0.5)
    }
    if (p == length(pred_levels)) { par(mfrow=c(1, 1)) }
  }

  if (length(cov_names==1)) {
    pred_levels = unique(predictors[,2])
    plotdat = new_data[ , 1:8]
    predictor = new_data[, c(1,10)]
    rdim = ceiling(length(pred_levels) / 2)
    ploteach(pred_levels, plotdat, predictor, rdim, cov_name = cov_names[1], formula_ref)
  } else {
    for (k in 1:length(cov_names)) {
      pred_levels = unique(predictors[,k+1])
      plotdat = new_data[ , 1:8]
      predictor = new_data[, c(1,10)]
      rdim = ceiling(length(pred_levels) / 2)
      ploteach(pred_levels, plotdat, predictor, rdim, cov_name = cov_names[1], formula_ref)
    }

  }

}
