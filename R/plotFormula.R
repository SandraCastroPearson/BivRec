#' #' Produce Bivariate Alternating Recurrent Series Plots by Categories
#'
#' @description
#' This function plots bivariate recurrent event gap times by levels of the categorical covariates as indicated in a formula.
#'
#' @import graphics
#' @importFrom utils tail
#' @importFrom stats model.frame
#' @importFrom stats na.omit
#' @importFrom stats model.matrix

#'
#' @param formula must be a formula that has a bivrecSurv class object as the response.
#' @param data a data frame that contains categorical covariates and/or vectors to create the response indicate in the given formula.
#' @param main optional string with plot title (default is no title)
#' @param xlab optional string with label for horizontal axis (default is "Gap Times")
#' @param ylab optional string with label for vertical axis (default is "Individual")
#' @param type1 optional string to label type 1 gap times (default is "Type 1")
#' @param type2 optional string to label type 2 gap times (default is "Type 2")
#'
#' @examples
#' library(BivRec)
#'
#' #Simulate data
#' set.seed(1234)
#' bivrec_data <- simulate(nsize=150, beta1=c(0.5,0.5), beta2=c(0,-0.5), tau_c=63, set=1.1)
#'
#' #Plot by the categorical covariates in a formula
#' plot(bivrecSurv(id, epi, xij, yij, d1, d2) ~ a1 + a2, data = bivrec_data)
#'
#' @export

plot.formula <- function(formula, data, main, xlab, ylab, type1, type2) {
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
  if (!is.bivrecSurv(response)) stop("Response must be a bivrecSurv object")
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
