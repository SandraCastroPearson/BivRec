#FIRST FUNCTION CALLS ON COMPILED BIVRECUR FORTRAN CODE

#                                 FORTRAN CODE                                 #
#______________________________________________________________________________#
# Original by Shu-Hui Chang                                                    #
# Modified by Chiung-Yu for bivariate recurrent event - Fortran (Feb, 2001)    #
# Modified and compiled for package by Sandra Castro-Pearson (July, 2018)      #
# Received from Xianghua Luo (May, 2018)                                       #
#______________________________________________________________________________#

r_bivrecur <- function(n, gtime, ctime, mc, m,
                       cen, ucen, nd, udt, tot, gap, event,
                       r, d, sest, var, markvar1, markvar2,
                       mark1, mark2, u1, u2, Fest, tmpindex, prob, std){

  out1 <- .Fortran("bivrecur",
                   n=as.integer(n),
                   gtime=as.double(gtime),
                   ctime=as.double(ctime),
                   count=as.double(m),
                   mc=as.integer(mc),
                   m=as.integer(m),
                   cen=as.double(cen),
                   ucen=as.double(ucen),
                   nd=as.integer(nd),
                   udt=as.double(udt),
                   tot=as.integer(tot),
                   gap=as.double(gap),
                   event=as.double(event),
                   r=as.double(r),
                   d=as.double(d),
                   sest=as.double(sest),
                   var=as.double(var),
                   markvar1=as.double(markvar1),
                   markvar2=as.double(markvar2),
                   mark1=as.double(mark1),
                   mark2=as.double(mark2),
                   u1=as.double(u1),
                   u2=as.double(u2),
                   Fest=as.double(Fest),
                   tmpindex=as.integer(tmpindex),
                   prob=as.double(prob),
                   std= as.double(std))

  out2 <- c(prob = out1$prob, std = out1$std)

  return(out2)
}

###################################################################
#################### FUNCTION NOT FOR USER ########################
###################################################################
#' A Function for non-parametric analysis on a biv.rec object for joint cdf
#'
#' @description
#' This function calculates the joint CDF for bivariate recurrent events. Called from bivrecNP(). No user interface.
#' @param fit_data Passed from bivrecNP().
#' @param u Passed from bivrecNP().
#' @param ai Passed from bivrecNP().
#' @param CI Passed from bivrecNP().
#'
#' @return A dataframe with the joint CDF
#'
#' @useDynLib BivRec bivrecur
#'
#' @keywords internal
#'

nonparam_cdf <- function(fit_data, u, ai, CI) {

  n <- fit_data$n
  m <- fit_data$m
  mc <- fit_data$mc
  nd <- fit_data$nd
  tot <- fit_data$tot
  gap <- fit_data$gap
  event <- fit_data$event
  markvar1 <- fit_data$markvar1
  markvar2 <- fit_data$markvar2
  udt <- fit_data$udt
  ctime <- fit_data$ctime
  ucen <- fit_data$ucen
  r <- fit_data$r
  d <- fit_data$d
  sest <- fit_data$sest
  Fest <- fit_data$Fest
  var <- fit_data$var
  prob <- fit_data$prob
  std <- fit_data$std
  gtime <- fit_data$gtime
  cen <- fit_data$cen
  mark1 <- fit_data$mark1
  mark2 <- fit_data$mark2

  estcdf <- list()

  for (u_count in 1:nrow(u)) {
    u1 <- u[u_count, 1]
    u2 <- u[u_count, 2]


    tmpindex <-sum(as.integer(udt<=(u1+u2)))  ### index ORINALLY PART OF BIVGAP FUNCTION
    if (tmpindex==0) {
      temp <- data.frame(u1, u2, prob=0, std=0)
      rownames(temp) <- "1"
      estcdf[[u_count]] <- temp
    } else {
      estimates <- r_bivrecur(n, gtime, ctime, mc, m,
                              cen, ucen, nd, udt, tot, gap, event,
                              r, d, sest, var, markvar1, markvar2,
                              mark1, mark2, u1, u2, Fest, tmpindex, prob, std)
      estcdf[[u_count]] <- data.frame(u1, u2, prob=estimates[1], std=estimates[2])
    }
  }

  out1 <- data.frame(matrix(unlist(estcdf), nrow=nrow(u), byrow=T))

  conf_lev = 1 - ((1-CI)/2)
  out1$lower <- out1[,3] - qnorm(conf_lev)*out1[,4]
  out1$upper <- out1[,3] + qnorm(conf_lev)*out1[,4]
  out1$lower[which(out1$lower<0)] <- 0
  out1$upper[which(out1$upper>1)] <- 1

  lowstring <- paste((1 - conf_lev), "%", sep="")
  upstring <- paste(conf_lev, "%", sep="")
  colnames(out1) <- c("x", "y", "Joint_Probability", "SE", lowstring, upstring)

  return(cdf=out1)

}
