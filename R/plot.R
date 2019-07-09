########################    PLOT     ########################
#' Plot Bivariate Alternating Recurrent Series
#'
#' @description
#' This function plots bivariate recurrent event gap times in different format
#'
#' @import graphics
#' @importFrom utils tail
#' @importFrom stats model.frame
#' @importFrom stats na.omit
#' @importFrom stats model.matrix
#'
#' @param x either coordinates for a plot or an object of class \code{bivrecSurv}, \code{bivrecReg}, \code{bivrecNP} or a formula.
#' @param y the y coordinates of points for a plot if coordinates were used in x, not required otherwise
#' @param data argument only when x is a formula. Should inidicate the data frame that contains categorical covariates and/or vectors to create the response indicate in the given formula.
#' @param ... arguments to be passed to graphical methods as needed.
#' \itemize{
#'  \item main: String with plot title (default is no title)
#'  \item xlab: String with label for horizontal axis (default is "Gap Times")
#'  \item ylab: String with label for vertical axis (default is "Individual")
#'  \item type1: String to label type 1 gap times (default is "Type 1")
#'  \item type2: String to label type 2 gap times (default is "Type 2")
#'}
#' @export
#'
#' @examples
#'# Simulate bivariate alternating recurrent event data
#' library(BivRec)
#' set.seed(1234)
#' bivrec_data <- simulate(nsize=150, beta1=c(0.5,0.5), beta2=c(0,-0.5), tau_c=63, set=1.1)
#'
#' plot(with(bivrec_data, bivrecSurv(id, epi, xij, yij, d1, d2)), main="Example")
#' plot(bivrecSurv(id, epi, xij, yij, d1, d2) ~ a1 + a2, data = bivrec_data, main="Example")
#'

plot <- function(x, y, data, ...) {UseMethod("plot")}

plot.default <- function(x, y, data,...) {
  data=NULL
  graphics::plot(x, y,...)}

plot.formula <- function(x, y=NULL, data,...) {
  if (!inherits(x, "formula")) stop("Object must be a formula")
  formula <- x

  #check arguments for labels

  if (missing(xlab)) {xlab="Gap Times"}
  if (missing(ylab)) {ylab="Individual"}
  if (missing(main)) {main=""}
  if (missing(type1)) {type1="Type 1"}
  if (missing(type2)) {type2="Type 2"}

  args = c(main, xlab, ylab, type1, type2)

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

plot.bivrecReg <- function(x, y=NULL, data=NULL,...) {

  if (!inherits(x, "bivrecReg")) stop("Object must be a bivrecReg class")
  object <- x
  formula = object$formula
  data = object$data$original
  if (missing(xlab)) {xlab="Gap Times"}
  if (missing(ylab)) {ylab="Individual"}
  if (missing(main)) {main=""}
  if (missing(type1)) {type1="Type 1"}
  if (missing(type2)) {type2="Type 2"}

  args = c(main, xlab, ylab, type1, type2)

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

plot.bivrecSurv <- function(x, y=NULL,...){
  if (!inherits(x, "bivrecSurv")) stop("Object must be a bivrecSurv class")
  object <- x

  #check arguments for labels
  if (missing(main)) {main=""}
  if (missing(xlab)) {xlab="Gap Times"}
  if (missing(ylab)) {ylab="Individual"}
  if (missing(type1)) {type1="Type 1"}
  if (missing(type2)) {type2="Type 2"}

  args = c(main, xlab, ylab, type1, type2)

  ##EXTRACT VECTORS FOR PLOTTING FUNCTION
  parameters <- object$data4Creg[-(5:7)]
  colnames(parameters) <- c("id", "episode", "xij", "yij", "ci")
  ctimes <- object$data4Lreg$ctime
  nsubject <- object$data4Lreg$n
  temp <- NULL

  basicplot(parameters, ctimes, nsubject, temp, args, a = 1/2, b = 1/3)

}

plot.bivrecNP <-function(x, y=NULL,...) {
  if (!inherits(x, "bivrecNP")) stop("Object must be a bivrecNP class")
  cond=x$conditional #boolean saying if conditional is in bivrecNP object
  if (cond==FALSE){
    plotJoint(x)
    plotMarg(x)
  }
  else {
    plotJoint(x)
    #par(mfrow=c(1,2))
    plotMarg(x)
    plotCond(x)
  }
}
