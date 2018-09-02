#FIRST FUNCTION CALLS ON COMPILED ONESAMP FORTRAN CODE

#                                 FORTRAN CODE                                 #
#______________________________________________________________________________#
# Original by Shu-Hui Chang                                                    #
# Modified by Chiung-Yu for bivariate recurrent event - Fortran (Feb, 2001)    #
# Modified and compiled for package by Sandra Castro-Pearson (July, 2018)      #
# Received from Xianghua Luo (May, 2018)                                       #
#______________________________________________________________________________#

r.onesamp <- function(n,gtime,ctime,mc,m,
                       cen,ucen,nd,udt,tot,gap,event,
                       r,d,sest,std){

  out1 <- .Fortran("onesamp",
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
                   std= as.double(std))

  out2 <- data.frame(time = out1$udt, surv = out1$sest, std = out1$std)

  return(out2)
}

###################################################################
#################### FUNCTION NOT FOR USER ########################
###################################################################
#' A Function for non-parametric analysis on a biv.rec object
#'
#' @description
#' This function calculates the marginal survival for bivariate recurrent events. Called from biv.rec.np(). No user interface.
#' @param fit_data An object that has been reformatted using the biv.rec.reformat() function. Passed from biv.rec.np().
#' @param CI Passed from biv.rec.np().
#'
#' @return A data frame with marginal survival
#'
#' @useDynLib BivRec onesamp
#'
#' @keywords internal
#'

nonparam.marginal <- function(fit_data, CI) {

  n <- fit_data$n
  m <- fit_data$m
  mc <- fit_data$mc
  nd <- fit_data$nm1
  tot <- fit_data$tot
  gap <- fit_data$markvar1
  event <- fit_data$event
  udt <- fit_data$umark1
  ctime <- fit_data$ctime
  ucen <- fit_data$ucen
  gtime <- fit_data$mark1
  cen <- fit_data$cen
  r = d = sest = std = rep(0, nd)


  surv <- r.onesamp(n,gtime,ctime,mc,m,
                    cen,ucen,nd,udt,tot,gap,event,
                    r,d,sest,std)

  conf.lev = 1 - ((1-CI)/2)
  surv$lower <- surv[,2] - qnorm(conf.lev)*surv[,3]
  surv$upper <- surv[,2] + qnorm(conf.lev)*surv[,3]
  surv$lower[which(surv$lower<0)] <- 0
  surv$upper[which(surv$upper>1)] <- 1

  low.string <- paste((1 - conf.lev), "%", sep="")
  up.string <- paste(conf.lev, "%", sep="")
  colnames(surv) <- c("Time", "Marginal.Survival", "SE", low.string, up.string)

  return(marg.survival = surv)

}
