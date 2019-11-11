#'
#' #' Plot Bivariate Alternating Recurrent Series by up to six covariates after analysis using \verb{bivrecReg}.
#' #'
#' #' @description
#' #' This function plots bivariate recurrent event gap times from a bivrecReg object. For examples see \verb{bivrecReg}.
#' #'
#' #' @import graphics
#' #' @importFrom utils tail
#' #' @importFrom stats model.frame
#' #' @importFrom stats na.omit
#' #' @importFrom stats model.matrix
#' #'
#' #' @param x An object of class \verb{bivrecReg}.
#' #' @param y Either empty or NULL
#' #' @param main Optional string with plot title. Default is no title.
#' #' @param xlab Optional string with label for horizontal axis. Default is "Time".
#' #' @param ylab Optional string with label for vertical axis. Default is "Individual".
#' #' @param type Optional vector of strings to label Type I and Type II gap times. Default is c("Type 1", "Type 2").
#' #' @param ... arguments to be passed to graphical methods as needed.
#' #'
#' #' @noRd
#' #' @keywords internal
#' #'
#' #
# plot.bivrecReg <- function (x, y = NULL, type = NULL, main = NULL, xlab = NULL, ylab = NULL, ...){
#
#   if (!inherits(x, "bivrecReg")) stop("Object must be a bivrecReg class")
#   object <- x
#   formula = object$formula
#   data = object$data$original
#   if (missing(type)) {type=c("Type 1","Type 2")}
#   if (missing(xlab)) {xlab="Time"}
#   if (missing(ylab)) {ylab="Individual"}
#   if (missing(main)) {main=""}
#
#   args = c(main, xlab, ylab, type)
#
#   #pull objects out of formula
#   formula_ref = formula
#   if (!missing(data)) {response <- eval(formula[[2]], data)
#   } else {stop("data argument missing")}
#   if (!inherits(response, "bivrecSurv")) stop("Response must be a bivrecSurv object")
#   formula[[2]] <- NULL
#   predictors <- data.frame(id = response$id_ref, model.matrix(formula, data)[,-1])
#   colnames(predictors) <-  c("id", colnames(model.matrix(formula, data))[-1])
#   df = as.data.frame(response$data4Creg)
#
#   #number of levels for each predictor
#   num_levs <- apply(predictors[,-1], 2, function(x) length(unique(x)))
#   preds_to_delete <- which(num_levs > 6) + 1
#
#   message1 <- paste(colnames(predictors)[preds_to_delete], " not used - either continuous or had more than 6 levels.", sep="")
#   print(message1)
#
#   predictors <- predictors[-preds_to_delete]
#
#   if (ncol(predictors)==1) {stop("Cannot break by covariate. All covariates are continuous or have more than 6 levels.")}
#
#   cov_names <- colnames(predictors)[-1]
#   new_data <- cbind(df, predictors)
#   orig_num <- length(unique(new_data$id))
#   new_data <- na.omit(new_data)
#   new_num <-length(unique(new_data$id))
#
#   message <- paste("Original number of subjects: ", orig_num,
#                    ". Subjects for plots: ", new_num, sep="")
#   print(message)
#
#   #plot each predictor
#   ploteach <- function(pred_levels, plotdat, predictor, rdim, cov_name, formula_ref) {
#     par(mar=c(1,1,1,1))
#     par(mfrow=c(rdim, 2))
#     for (p in 1:length(pred_levels)) {
#       index <- which(predictor[,2] == pred_levels[p])
#       tempdat <- plotdat[index,]
#       new_main = paste("For ", cov_name, " = ", pred_levels[p], sep="")
#       bdat = eval(formula_ref[[2]], tempdat)
#
#       ##EXTRACT VECTORS FOR PLOTTING FUNCTION
#       parameters <- bdat$data4Creg[-(5:7)]
#       colnames(parameters) <- c("id", "episode", "xij", "yij", "ci")
#       ctimes <- bdat$data4Lreg$ctime
#       nsubject <- bdat$data4Lreg$n
#
#       #Plot for one level of covariate
#       args2 = args
#       args2[1] = new_main
#       basicplot(parameters, ctimes, nsubject, temp=NULL, args = args2, a = 2/6, b=1/5, c=0.5)
#     }
#     if (p == length(pred_levels)) { par(mfrow=c(1, 1)) }
#   }
#
#   if (length(cov_names==1)) {
#     pred_levels = unique(predictors[,2])
#     plotdat = new_data[ , 1:8]
#     predictor = new_data[, c(1,10)]
#     rdim = ceiling(length(pred_levels) / 2)
#     ploteach(pred_levels, plotdat, predictor, rdim, cov_name = cov_names[1], formula_ref)
#   } else {
#     for (k in 1:length(cov_names)) {
#       pred_levels = unique(predictors[,k+1])
#       plotdat = new_data[ , 1:8]
#       predictor = new_data[, c(1,10)]
#       rdim = ceiling(length(pred_levels) / 2)
#       ploteach(pred_levels, plotdat, predictor, rdim, cov_name = cov_names[1], formula_ref)
#     }
#
#   }
#
# }

#' Plot by function
#'
#' @param df passed from plot.bivrecSurv
#' @param predictors passed from plot.bivrecSurv
#' @param args passed from plot.bivrecSurv
#'
#' @keywords internal
#' @noRd

plotBy <- function(df, args) {

  #number of levels for each predictor
  num_levs <- apply(df[, 6:ncol(df)], 2, function(x) length(unique(x)))
  to_delete <- which(num_levs > 6) + 5

  ploteach <- function(pred_levels, plotdat, rdim, cov_name, args) {
    plot.new()
    par(mar=c(1,1,1,1))
    par(mfrow=c(rdim, 2))
    for (p in 1:length(pred_levels)) {
      ##EXTRACT VECTORS FOR PLOTTING FUNCTION
      index <- which(plotdat[ ,6] == pred_levels[p])
      parameters <- plotdat[index, 1:5]
      new_main = paste("For ", cov_name, " = ", pred_levels[p], sep="")
      ctimes <- unique(parameters$ci)
      unik_ids <- unique(parameters$id)
      nsubject2 <- length(unik_ids)
      parameters$id2 <- NA
      for (i in 1:nsubject2){
          parameters$id2[parameters$id == unik_ids[i]]=i
      }
      parameters2 <- cbind(parameters[,6], parameters[,2:5])
      parameters <- data.frame(parameters2)
      colnames(parameters) <- c("id", "episode", "xij", "yij", "ci")

      #Plot for one level of covariate
      args2 = args
      args2[1] = new_main
      par(mar=c(5,4,4,2)+0.1)
      basicplot(parameters, ctimes, nsubject=nsubject2,
                temp=NULL, args = args2, c=0.6, cm=0.9, byp=TRUE)
    }
    if (p == length(pred_levels)) {
      par(mfrow=c(1, 1))}
  }

  message1 <- paste(colnames(df)[to_delete], " not used - either continuous or had more than 6 levels.", sep="")
  print(message1)
  df <- df[, -to_delete]

  if (ncol(df)==5) {stop("Cannot break by covariate. All covariates are continuous or have more than 6 levels.")}

  cov_names <- colnames(df)[6:ncol(df)]
  nsubject <-length(unique(df$id))

  message <- paste("Subjects for plots: ", nsubject, ".", sep="")
  print(message)

  if (length(cov_names)==1) {
    pred_levels = unique(df[,6])
    plotdat = df[ , 1:6]
    plotdat = na.omit(plotdat)
    rdim = ceiling(length(pred_levels) / 2)
    ploteach(pred_levels, plotdat, rdim, cov_name = cov_names, args)
  } else {
    for (k in 1:length(cov_names)) {
      pred_levels = unique(df[ ,5+k])
      plotdat = df[, c(1:5, 5+k)]
      plotdat = na.omit(plotdat)
      rdim = ceiling(length(pred_levels) / 2)
      ploteach(pred_levels, plotdat, rdim, cov_name = cov_names[1])
    }
  }

}


