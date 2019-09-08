########################    plot.formula     ########################
#' Plot Bivariate Alternating Recurrent Series from Formula
#'
#' @description
#' This function plots bivariate recurrent event gap times from a formula
#'
#' @import graphics
#' @importFrom utils tail
#' @importFrom stats model.frame
#' @importFrom stats na.omit
#' @importFrom stats model.matrix
#'
#' @param x either coordinates for a plot or an object of class formula.
#' @param y either empty or NULL
#' @param data required argument when x is a formula. Should indicate the data frame that contains categorical covariates and/or vectors to create the response indicate in the given formula.
#' @param main Optional string with plot title. Default is no title.
#' @param xlab Optional string with label for horizontal axis. Default is "Gap Times".
#' @param ylab Optional string with label for vertical axis. Default is "Individual".
#' @param type Optional vector of strings to label type 1and type 2 gap times. Default is c("Type 1", "Type 2").
#' @param ... arguments to be passed to graphical methods as needed.
#'
#' @export
#'
#' @examples
#'# Simulate bivariate alternating recurrent event data
#' library(BivRec)
#' set.seed(1234)
#' bivrec_data <- simulate(nsize=150, beta1=c(0.5,0.5), beta2=c(0,-0.5), tau_c=63, set=1.1)
#' plot(x = bivrecSurv(id, epi, xij, yij, d1, d2) ~ a1 + a2, data = bivrec_data,
#'      main="Example", xlab="type 1", ylab="type 2")
#'

plot.formula <- function(x, y=NULL, data, type = NULL,
                         main = NULL, xlab = NULL, ylab = NULL, ...) {
  if (!inherits(x, "formula")) stop("Object must be a formula")
  formula <- x

  #check arguments for labels
  if (missing(type)) {type=c("Type 1","Type 2")}
  if (missing(xlab)) {xlab="Gap Times"}
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
