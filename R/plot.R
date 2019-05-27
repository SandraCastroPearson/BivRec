basicplot <- function(parameters, ctimes, nsubject, temp, args, a, b) {

  #### Reformat data to start-stop times ########
  for (iter in 1:nsubject) {
    subject <- parameters[parameters$id==iter,] #pull dat for each subject
    #print(subject)
    if (nrow(subject)==1) {
      start_time <- c(0, subject$xij)
      stop_time <- c(subject$xij, subject$xij+subject$yij)
      temp <- rbind(temp, cbind(sub_id=rep(iter, 2), sub_epi=rep(1, 2),
                                start_time, stop_time, sub_ci=rep(subject$ci, 2)))
    } else {
      subject <- subject[order(subject$episode),]
      times <- c(rbind(subject$xij, subject$yij))
      start_time <- 0
      stop_time <- times[1]
      for (j in 2:length(times)) {
        start_time[j] <- stop_time[j-1]
        stop_time[j] <- start_time[j]+times[j]
      }
      sub_id <- rep(iter, length(times))
      sub_ci <- rep(subject$ci[1], length(times))
      sub_epi <- rep(subject$episode, each=2)
      temp <- rbind(temp, cbind(sub_id, sub_epi, start_time, stop_time, sub_ci))
    }
  }

  data4plot <- data.frame(temp)
  data4plot <- data4plot[-which(data4plot$start_time==data4plot$stop_time), ]
  colnames(data4plot) <- c("id", "epi", "start_time", "stop_time", "ci")
  data4plot  <- data4plot[order(data4plot$ci),]

  ###### PLOT ########
  # get the range for the x and y axis
  xrange <- c(0, max(ctimes)+ 10) #Time
  yrange <- c(0, nsubject) #subjects

  # set up the plot
  plot.new( )
  plot.window(xlim=xrange, ylim=yrange)
  axis(side=1)
  title(main=args[1], xlab=args[2], ylab=args[3])
  legendtext = c(args[4], args[5])

  legend(xrange[2]*a, yrange[2]*b, legend=legendtext, box.lty=0,
         col = c("blue", "red"), lty = 1, cex=0.6)

  # add line segments
  newid = 0
  for (c_iter in sort(ctimes)) {
    newid = newid + 1
    subject <- subset(data4plot, data4plot$ci == c_iter)
    subject$newid = newid
    if (nrow(subject)==1) {
      segments(subject$start_time[1], subject$newid,
               subject$stop_time[1], subject$newid, col="blue")
    } else {
      for (j in 1:length(subject$id)) {
        if (j %% 2 == 1) colors <- "blue" else colors <- "red"
        segments(subject$start_time[j], subject$newid,
                 subject$stop_time[j], subject$newid, col=colors)
      }
    }
  }

}

#' #' Produce Bivariate Alternating Recurrent Series Plot
#'
#' @description
#' This function plots bivariate recurrent event gap times for either a bivrec object or a function. If a function is used plots will be produced for levels of categorical covariates.
#'
#' @import graphics
#' @importFrom utils tail
#' @importFrom stats model.frame
#' @importFrom stats na.omit
#' @importFrom stats model.matrix

#'
#' @param object must be an object of \code{bivrecSurv}, \code{bivrecReg} or \code{bivrecNP} class or a formula that has a bivrecSurv class object as the response.
#' @param data a data frame that contains categorical covariates from a formula (if a formula is used for the object to be plotted).
#'
#' @examples
#' library(BivRec)
#'
#' #Simulate data
#' set.seed(1234)
#' bivrec_data <- simulate(nsize=150, beta1=c(0.5,0.5), beta2=c(0,-0.5), tau_c=63, set=1.1)
#'
#' #Plot a bivrecObject
#'
#' bivrecSurvObject <- with(bivrec_data, bivrecSurv(id, epi, xij, yij, d1, d2))
#' plot(bivrecSurvObject)
#'
#' #Plot a bivrecObject by categorical covariate from a formula
#' plot(bivrecSurv(id, epi, xij, yij, d1, d2) ~ a1 + a2, data = bivrec_data)
#'
#' @export
#'
#'

plot.bivrecSurv <- function(object, main, xlab, ylab, type1, type2){
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

#' @export
plot.formula <- function(object, data, main, xlab, ylab, type1, type2) {
  #check arguments for labels
  if (missing(xlab)) {xlab="Gap Times"}
  if (missing(ylab)) {ylab="Individual"}
  if (missing(main)) {main=""}
  if (missing(type1)) {type1="Type 1"}
  if (missing(type2)) {type2="Type 2"}

  args = c(main, xlab, ylab, type1, type2)

  #pull objects out of formula
  formula_ref = formula = object
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

#' @export
plot.bivrecReg <- function(object, main, xlab, ylab, type1, type2) {
    newobject = object$formula
    newdata = object$data$original
    plot.formula(newobject, newdata, main, xlab, ylab, type1, type2)
    }

