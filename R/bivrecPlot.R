#' #' Produce Bivariate Alternating Recurrent Series Plot
#'
#' @description
#' This function plots bivariate recurrent event gap times.
#'
#' @importFrom graphics plot
#' @importFrom graphics segments
#' @importFrom graphics legend
#' @importFrom utils tail
#' @importFrom stats model.frame
#' @importFrom stats na.omit
#'
#' @param x an object of class \code{bivrecSurv}, or \code{bivrecReg}, or \code{bivrecNP}, usually returned by the \code{bivrecSurv} function.
#'
#' @examples
#' library(BivRec)
#' set.seed(1234)
#' dat <- simulate(nsize=150, beta1=c(0.5,0.5), beta2=c(0,-0.5), tau_c=63, set=1.1)
#' bdat <- with(dat, bivrecSurv(id, epi, xij, yij, d1, d2))
#' bivrecPlot(x)
#'
#' @export
#'
#'
#bivrecPlot <- function(x) UseMethod("bivrecPlot")

plot.bivrecSurv=function(x,main,xlab,ylab,type1,type2){
  #if (!is.bivrecSurv(x)) stop("Response must be a bivrecSurv class")
  #EXTRACT VECTORS FOR PLOTTING FUNCTION
  parameters <- x$df[-(5:7)]
  colnames(parameters) <- c("id", "episode", "xij", "yij", "ci")
  ctimes <- x$ctime
  nsubject <- x$n
  temp <- NULL

  #### Reformat data to start-stop times ########
  for (iter in 1:nsubject) {
      subject <- parameters[parameters$id==iter,] #pull dat for each subject
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
  #check arguments for labels 
  if (missing(xlab)) {xlab="Gap Times"}
  if (missing(ylab)) {ylab="Individual"}
  if (missing(main)) {main=""}
  if (missing(type1)) {type1="Type 1"}
  if (missing(type2)) {type2="Type 2"}
  # if (!is.numeric(xlab)&&!is.character(xlab)) {xlab="Gap Times"}
  # if (!is.numeric(ylab)&&!is.character(ylab)) {xlab="Individual"}
  # if (!is.numeric(main)&&!is.character(main)) {xlab="Bivariate Recurrent Gap Times"}
  # if (!is.numeric(type1)&&!is.character(type1)) {xlab="Type1"}
  # if (!is.numeric(type2)&&!is.character(type2)) {xlab="Type2"}
  
  plot(xrange, yrange, type="n", main=main, xlab=xlab, ylab = ylab, yaxt='n')
  legendtext = c(type1, type2)
  legend("bottomright", legend=legendtext, bty='n', inset = c(0,0),
         col = c("gray", "salmon"), lty = 1, cex=0.9)

  # add line segments
  newid = 0
  for (c_iter in sort(ctimes)) {
    newid = newid + 1
    subject <- subset(data4plot, data4plot$ci == c_iter)
    subject$newid = newid
    if (nrow(subject)==1) {
      segments(subject$start_time[1], subject$newid,
               subject$stop_time[1], subject$newid, col="gray")
    } else {
      for (j in 1:length(subject$id)) {
        if (j %% 2 == 1) colors <- "gray" else colors <- "salmon"
        segments(subject$start_time[j], subject$newid,
                 subject$stop_time[j], subject$newid, col=colors)
      }
    }
  }

}
