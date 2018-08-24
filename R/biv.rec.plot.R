#' Bivariate Recurrent Response Plotting
#'
#' @description
#' This function plots bivariate recurrent event gap times.
#' @param formula Formula of the form \strong{ID + episode ~ xij + yij}.
#' \itemize{
#'   \item ID: a vector of subjects' unique identifier which can be numeric or character.
#'   \item episode: is a vector indicating the episode of the bivariate alternating gap time pairs.
#'   \item xij: is a vector with the lengths of the type I gap times.
#'   \item yij: is a vector with the lengths of the type II gap times.
#' }
#' @param data A data frame that contains all the vectors listed in the formula
#'
#' @importFrom graphics plot
#' @importFrom graphics segments
#' @importFrom graphics legend
#' @importFrom utils tail
#' @importFrom stats model.frame
#' @importFrom stats na.omit
#'
#' @examples
#' library(BivRec)
#' set.seed(1234)
#' sim.biv.rec <- data.sim(150, beta1=c(0.5,0.5), beta2=c(0,-0.5), cr=63, sg2=0.5, set=1.1)
#' biv.rec.plot(formula = id + epi ~ xij + yij, data = sim.biv.rec)
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
  names <- paste("data$", variables, sep="")
  id <- eval(parse(text = names[1]))
  episode <- eval(parse(text = names[2]))
  time1 <- eval(parse(text = names[3]))
  time2 <- eval(parse(text = names[4]))
  nsubject <- max(unique(id))
  parameters <- data.frame(id, episode, time1, time2)
  temp <- NULL

  for (iter in 1:nsubject) {
    subject <- parameters[parameters$id==iter,]
    if (length(which(duplicated(subject$episode)==TRUE))!=0) {
      print(paste("Error subject", iter, "has non-unique episodes", sep=" "))
      stop()
    }
    if (length(subject$episode)==1) {
      start_time <- c(0, subject$time1)
      stop_time <- c(subject$time1, subject$time1+subject$time2)
      temp <- rbind(temp, cbind(sub_id=rep(iter, 2), sub_epi=rep(1, 2), start_time, stop_time))
    } else {
        subject <- subject[sort(subject$episode),]
        times <- c(rbind(subject$time1, subject$time2))
        start_time <- 0
        stop_time <- times[1]
        for (j in 2:length(times)) {
          start_time[j] <- stop_time[j-1]
          stop_time[j] <- start_time[j]+times[j]
        }
        sub_id <- rep(iter, length(times))
        sub_epi <- rep(subject$episode, each=2)
        temp <- rbind(temp, cbind(sub_id, sub_epi, start_time, stop_time))
    }
  }
  data4plot <- data.frame(temp)
  data4plot <- data4plot[-which(data4plot$start_time==data4plot$stop_time), ]
  colnames(data4plot) <- c("id", "epi", "start_time", "stop_time")

  ###### PLOT ########
  # get the range for the x and y axis
  xrange <- c(0, max(data4plot$stop_time) + max(data4plot$stop_time)/4)
  yrange <- range(data4plot$id)

  # set up the plot
  plot(xrange, yrange, type="n", xlab="Time", ylab="Subject", yaxt='n')
  legendtext = c(variables[3], variables[4])
  legend("topright", legendtext, bty='n', inset = c(0,0), col = c("red", "blue"), lty = 1)
  # legend(list(x=max(data4plot$stop_time), y=max(data4plot$id)-max(data4plot$id)/2.5), legendtext, bty='n',
  #        inset = c(0,0),cex = 0.7, col = c("red", "blue"), lty = 1, lwd=0.5,
  #        y.intersp=0.25, x.intersp=0.3, xjust=0, yjust=0)

  #colors <- rainbow(nsubject)
  plotchar <- seq(18, 18+nsubject, 1)

  withgaps <- NULL

  # add lines
  for (iter in 1:nsubject) {
    subject <- subset(data4plot, data4plot$id == iter)
    if (nrow(subject)==1) {
      segments(subject$start_time[1], subject$id,
               subject$stop_time[1], subject$id, col="red")
    } else {
        for (j in 1:length(subject$id)) {
          if (j %% 2 == 1) colors <- "red" else colors <- "blue"
          segments(subject$start_time[j], subject$id,
                   subject$stop_time[j], subject$id, col=colors)
        }
      }
    }
}
