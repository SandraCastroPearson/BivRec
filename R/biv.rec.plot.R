#' Bivariate Recurrent Response Plotting
#'
#' @description
#' This function plots recurrent event data.
#' @param formula Formula of the form \strong{ID ~ xij + yij}.
#' \itemize{
#'   \item ID: is a vector of unique subject identifiers.
#'   \item xij: is a vector with the length of time subject i spends in event of type 1 during observation j.
#'   \item yij: is a vector with the length of time subject i spends in event of type 2 during observation j.
#' }
#' @param data A data frame that contains all the vectors listed in the formula
#'
#' @importFrom graphics plot
#' @importFrom graphics segments
#' @importFrom utils tail
#' @importFrom stats model.frame
#' @importFrom stats na.omit
#'
#' @examples
#' \dontrun{
#' library(BivRec)
#' set.seed(1234)
#' sim.biv.rec <- data.sim(300, beta1=c(0.5,0.5), beta2=c(0,-0.5), cr=63, sg2=0.5, set=1.1)
#' biv.rec.plot(formula = id ~ xij + yij, data = sim.biv.rec)
#' }
#'
#' @export
#'

biv.rec.plot <- function(formula, data) {

####### SET-UP DATA #######
  # PULL DATA FROM FORMULA
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

  #EXTRACT VECTORS FOR PLOTTING FUNCTION
  parameters <- model.frame(formula, data, na.action = NULL)
  colnames(parameters) <- c("id", "time1", "time2")
  nsubject <- max(parameters$id)
  temp <- NULL

    for (iter in 1:nsubject) {
      subject <- subset(parameters, parameters$id==iter)
      times <- c(rbind(subject$time1, subject$time2))
      start_time <- 0
      stop_time <- times[1]
      for (j in 2:length(times)) {
        start_time[j] <- stop_time[j-1]
        stop_time[j] <- start_time[j]+times[j]
      }
      sub_id <- rep(iter, length(times))
      temp <- rbind(temp, cbind(sub_id, start_time, stop_time))
    }
  data4plot <- data.frame(temp)
  data4plot <- data4plot[-which(data4plot$start_time==data4plot$stop_time), ]
  colnames(data4plot) <- c("id", "start_time", "stop_time")

###### PLOT ########
  # get the range for the x and y axis
  xrange <- range(c(data4plot$start_time, data4plot$stop_time))
  yrange <- range(data4plot$id)

  # set up the plot
  plot(xrange, yrange, type="n", xlab="Time", ylab="Subject ID" )

  #colors <- rainbow(nsubject)
  plotchar <- seq(18, 18+nsubject, 1)

  withgaps <- NULL

  # add lines
  for (iter in 1:nsubject) {
    subject <- subset(data4plot, data4plot$id == iter)
    if (nrow(subject)==1) {
      segments(subject$start_time[j], subject$id,
               subject$stop_time[j], subject$id, col="red")
    } else {
      check <- as.factor(subject$start_time[2:nrow(subject)]==subject$stop_time[1:(nrow(subject)-1)])
      if (length(levels(check))==1 && levels(check)=="TRUE") {
        for (j in 1:length(subject$id)) {
          if (j %% 2 == 1) colors <- "red" else colors <- "blue"
          segments(subject$start_time[j], subject$id,
                   subject$stop_time[j], subject$id, col=colors)
        }
      } else {withgaps<-c(withgaps, iter)}
    }
  }
  if (is.null(withgaps)==FALSE) {
    print(matrix(c("Warning: Data with time gaps not plotted for following subject IDs", withgaps, "consider removing these for analysis to avoid errors"), ncol=1))
  }
}
