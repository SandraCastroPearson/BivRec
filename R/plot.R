#' #' Produce Bivariate Alternating Recurrent Series Plot
#'
#' @description
#' This function plots bivariate recurrent event gap times for either a bivrec object or a function. If a function is used plots will be produced for levels of categorical covariates.
#'
#' @importFrom graphics plot
#' @importFrom graphics segments
#' @importFrom graphics legend
#' @importFrom utils tail
#' @importFrom stats model.frame
#' @importFrom stats na.omit
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
#' plot(bivrecSurv(id, epi, xij, yij, d1, d2) ~ a1,
#'                     data = bivrec_data)
#'
#' @export
#'
#'

plot <- function(object, data, ...) UseMethod("plot")

#check arguments for labels
if (missing(xlab)) {xlab="Gap Times"}
if (missing(ylab)) {ylab="Individual"}
if (missing(main)) {main=""}
if (missing(type1)) {type1="Type 1"}
if (missing(type2)) {type2="Type 2"}


  if (is.bivrecSurv(object)) {plot.bivrecSurv(object, main, xlab, ylab, type1, type2)} else {
  if (inherits(object,"formula")) {plot.formula(object, data, main, xlab, ylab, type1, type2)} else {
  if (is.bivrecReg(object)) {plot.bivrecReg(object, data, main, xlab, ylab, type1, type2)} else {
  if (is.bivrecNP(object)) {plot.bivrecNP(object, data, main, xlab, ylab, type1, type2)}
    else {stop("Response must be of bivrecSurv, bivrecReg or bivrecNP class or a formula")
  }}}}

plot.bivrecSurv <- function(object, main, xlab, ylab, type1, type2){
  ##EXTRACT VECTORS FOR PLOTTING FUNCTION
  parameters <- object$data4Creg[-(5:7)]
  colnames(parameters) <- c("id", "episode", "xij", "yij", "ci")
  ctimes <- object$data4Lreg$ctime
  nsubject <- object$data4Lreg$n
  temp <- NULL

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
  xrange <- c(0, max(ctimes)+ 5) #Time
  yrange <- c(0, nsubject) #subjects

  # set up the plot

  plot(xrange, yrange, type="n", main=main, xlab=xlab, ylab = ylab, yaxt='n')
  legendtext = c(type1, type2)
  legend("bottomright", legend=legendtext, bty='n', inset = c(0,0),
         col = c("blue", "red"), lty = 1, cex=0.9)

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

plot.formula <- function(object, data, main, xlab, ylab, type1, type2) {

  formula_ref = formula
  if (!missing(data)) response <- eval(formula[[2]], data)
  if (!is.bivrecSurv(response)) stop("Response must be a bivrecSurv object")
  formula[[2]] <- NULL
  predictors <- data.frame(id = response$id_ref, model.matrix(formula, data)[,-1])
  colnames(predictors) <-  c("id", colnames(model.matrix(formula, data))[-1])
  cov_names <- colnames(predictors)[-1]

  new_data <- cbind(as.data.frame(response$dat4Creg), predictors)
  orig_num <- length(unique(new_data$id))
  new_data <- na.omit(new_data)
  new_num <-length(unique(new_data$id))

  #plotdat <- list()
  for (k in 2:(length(cov_names)+1)){
      #check variable is categorical and has less than 6 categories
      nlevels = length(levels(as.factor(predictors[,k])))
      if (nlevels <= 6) {
        pred_levels = unique(predictors[,k])
        rdim = ifelse(nlevels>4, 3, ifelse(nlevels<3, 1, 2))
        for (p in 1:pred_levels) {
            par(mfrow=c(rdim, 2))
            index <- which(predictors[,k] == pred_levels[p])
            tempdat <- new_data[index,]
            tempbdat <- eval(formula_ref[[2]], tempdat)
            main = paste("For ", cov_names[k-1], " = ", pred_levels[p], sep="")
            plot.bivrecSurv(tempbdat, main, xlab, ylab, type1, type2)
        }
      } else {
        message = paste("Can't plot by ", cov_names[k-1], "Too many levels.")
        print(message)
      }

  }


  message <- paste("Original number of subjects: ", orig_num,
                   ". Subjects for plots: ", new_num)
  print(message)




}

plot.bivrecReg <- function(object, data, main, xlab, ylab, type1, type2) {
    newobject = object$formula
    newdata = object$data$original
    plot.formula(newobject, newdata, main, xlab, ylab, type1, type2)}

# plot.bivrecNP <- function(object, data, main, xlab, ylab, type1, type2) {
#
# }

