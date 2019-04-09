
##-----reformat any dataset HIDDEN FROM USER

#                 m.dat, np.dat FUNCTIONS                             #
#_______________________________________________________________________________
# Original by Chihyun Lee (August, 2017)                                       #
# Last Modified by Sandra Castro-Pearson (March, 2019)                         #
# Received from Chihyun Lee (January, 2018)                                    #
#_______________________________________________________________________________



#####For non-parametric analysis
# dat : a data.frame including
#     1)id numbers, 2)orders of episodes, 3)first gap time, 4)second gap time
#     5)censoring times, 6) censoring indicators in each column
# ai: a non-negative function of censoring time

np.dat <- function(dat, ai) {

  id <- dat$id
  uid <- unique(id)   # vector of unique id's
  n.uid <- length(uid)   # scalar : number of unique IDs
  event <- dat$d2 #event indicator : must always be 0 for the last obs per ID and 1 otherwise
  markvar1 <- dat$vij
  markvar2 <- dat$wij
  gap <- markvar1 + markvar2

  m.uid <- as.integer(table(id))   # vector: number of observed pairs per id/subject (m)
  max.m <- max(m.uid, na.rm=T) # scalar : maximum number of repeated observations

  ifelse (ai == 1, weight <- rep(1, n.uid), weight <- dat$ci[which(dat$epi == 1)]) #Set weights

  tot <- length(gap) # total number of observations
  ugap <- sort(unique(gap[event == 1]))   # sorted unique uncensored X_0 gap times (support points for sum)
  n.ugap <- length(ugap)   # number of unique X_0 gap times (or support points for sum)

  umark1 <- sort(unique(markvar1[event == 1]))   # sorted unique uncensored V_0 times (support points for marginal)
  n.umark1 <- length(umark1) # number of unique V_0 gap times (or support points for marginal)

  # Space holders
  r <- sest <- Fest <- rep(0, n.ugap)
  d <- matrix(0, nrow = n.ugap, ncol = 2)
  prob <- var <- std <- 0
  gtime <- cen <- mark1 <- mark2 <- matrix(0, nrow = n.uid, ncol = max.m)

  out <- list(n = n.uid, m = m.uid, mc = max.m, nd = n.ugap, tot=tot,
              gap =gap, event = event, markvar1 = markvar1, markvar2 =markvar2,
              udt = ugap,  ctime = weight, ucen = m.uid-1,
              r = r, d=d, sest = sest, Fest = Fest, var = var,
              prob = prob, std = std, gtime = gtime, cen = cen,
              mark1 = mark1, mark2 = mark2, umark1=umark1, nm1 = n.umark1)

  return(out)
}

formarginal <- function(dat){

  mdata <- tmp<- NULL
  freq <-cumsum(c(0,table(dat[,1])))

  for (i in 1:(length(freq)-1)){
    tmp<- dat[(freq[i]+1):freq[i+1], -c(5,7)]
    if (nrow(tmp)==1){
      if(tmp$wij>0){
        mdata <- rbind(mdata,
                       c(id=tmp$id, vij=tmp$vij, wij=tmp$wij, d2=1, epi=tmp$epi, ci=tmp$ci),
                        c(id=tmp$id, vij=0, wij=0, d2=0, epi=2, ci=tmp$ci))
      } else{
        mdata<-rbind(mdata, c(id=tmp$id, vij=tmp$vij, wij=tmp$wij, d2=0, epi=tmp$epi, ci=tmp$ci))
      }
    } else{mdata<-rbind(mdata, tmp)}
  }
  return(mdata)
}

#####

#' A function to create a biv.rec object
#'
#' @description
#' This function takes a BivRec object and covariate information and reformats it for proper fit.
#'
#' @importFrom stats na.omit
#' @importFrom survival survfit
#'
#' @param identifier Vector of subject's unique identifier (i). Passed from biv.rec.fit()
#' @param xij Vector with the lengths of time spent in event of type X for individual i in episode j. Passed from biv.rec.fit()
#' @param yij vector with the lengths of time spent in event of type Y for individual i in episode j. Passed from biv.rec.fit()
#' @param c_indicatorY Vector with values of 0 for the last episode for subject i or 1 otherwise. A subject with only one episode will have one 0. Passed from biv.rec.fit()
#' @param c_indicatorX Optional vector with values of 0 if the last episode for subject i occurred for event of type X or 1 otherwise. A subject with only one episode could have either one 1 (if he was censored at event Y) or one 0 (if he was censored at event X). Passed from biv.rec.fit()
#' @param episode Vector indicating the observation or episode (j) for a subject (i).
#' @param method A string for method to be used. Passed from biv.rec.fit()
#' @param ai Passed from biv.rec.fit()
#' @param condgx Passed from biv.rec.fit()
#' @param data Passed from biv.rec.fit()
#'
#' @return a bivrec object ready for fitting.
#' @seealso \code{\link{biv.rec.fit}}
#'
#' @keywords internal

####MAIN FUNCTION

#### Check error - if yij =0 then both deltas are 0 at the end.

biv.rec.reformat <- function(identifier, xij, yij, c_indicatorY, c_indicatorX, episode, covariates, method, ai, condgx, data){

  if (condgx==FALSE){
  ### PUT POTENTIAL RESPONSE TOGETHER
  temp1 <- data.frame(identifier, xij, yij, c_indicatorY, c_indicatorX, episode, covariates)

<<<<<<< HEAD:R/biv.rec.reformat.R
  ####CHECK MISSINGNESS AND KEEP RECORD OF WHAT IS BEING OMITTED
  temp <- na.omit(temp1)
  n_missing <- nrow(temp1) - nrow(temp)

  ####CHECK xij, yij VARIABLES HAVE CORRECT VALUES

  invalid_xij <- which(temp$xij<=0)
  if (length(invalid_xij)>0) {
    print("Invalid values for length of time in event X for rows. All must be >0.")
    temp[invalid_xij,]
    stop()
  } else {
    invalid_yij <- which(temp$yij<0)
    if (length(invalid_yij)>0) {
      print("Invalid values for length of time in event Y for rows. All must be >=0.")
      temp[invalid_yij,]
      stop()
    }
      ####CHECK INDICATORS AND EPISODE SEQUENCES, yij VARIABLES HAVE CORRECT VALUES
    } else {

      #First check for indicators - all 0 and 1
      cx_check1 <- length(which(temp$c_indicatorX!=1&&temp$c_indicatorX!=0))
      cy_check1 <- length(which(temp$c_indicatorY!=1&&temp$c_indicatorY!=0))
      if (cx_check1!=0) {
        print("Invalid values for c_indicatorX")
        temp[which(c_indicatorX!=1&&c_indicatorX!=0),]
        stop()
      } else {
        if (cx_check1!=0) {
          print("Invalid values for c_indicatorY")
          temp[which(c_indicatorY!=1&&c_indicatorY!=0),]
          return("Error")
        } else {
=======
>>>>>>> 8c89f8021f00c92b91f53b5d49c67811755bed55:R/reformat.R
          ## Second check for indicators
          ## Indicators match episode (last is 0 or 1 for X an 0 for Y) and episode doesn't have gaps
          wrong_xind <- NULL
          wrong_yind <- NULL
          wrong_epi <- NULL
          unique_id <- unique(temp$identifier)
          for (i in 1:length(unique_id)) {
            sub_id <- unique_id[i]
            temp_by_subject <- subset(temp, temp$identifier==sub_id)
            temp_by_subject <- temp_by_subject[order(temp_by_subject$episode),]
            sub_n <- nrow(temp_by_subject)
            last_cx <- temp_by_subject$c_indicatorX[sub_n]
            rest_of_cx <- temp_by_subject$c_indicatorX[-sub_n]
            last_cy <- temp_by_subject$c_indicatorY[sub_n]
            rest_of_cy <- temp_by_subject$c_indicatorY[-sub_n]

            #Check for indicators (cx is all 1's or one zero at end, last cy is 0, all others are 1)
            if (sum(rest_of_cx)!=(sub_n-1)) {wrong_xind <- c(wrong_xind, sub_id)}
            if (sum(rest_of_cy)!=(sub_n-1)) {wrong_yind <- c(wrong_yind, sub_id)}
            if (last_cy!=0) {wrong_yind <- c(wrong_yind, sub_id)}
            if (last_cx==0) {if (last_cy==1) {wrong_yind <- c(wrong_yind, sub_id)}}

            #Check episodes don't have gaps
            for (j in 1:sub_n){
              if (temp_by_subject$episode[j]!=j) {
                wrong_epi <- c(wrong_epi, sub_id)}
            }
          }

          ##Print Errors and exit function with recommendations


            print(paste("Original number of observations:", nrow(temp1), "for", length(unique(temp1$identifier)), "individuals", sep=" "))
            print(paste("Observations to be used in analysis:", nrow(temp), "for", length(unique_id), "individuals",sep=" "))

            #Get ready to send to m.dat if needed
            covariate_indexes <- seq(7, ncol(temp), 1)
            response <- temp[, -covariate_indexes]
            #calculate censoring time
            ci=id=NULL
            j=1
            response$zij <- response$xij + response$yij
            for (i in unique(response$identifier)){
              tempi=response[response$identifier == i, ]
              if (nrow(tempi) == 1){
                ci=c(ci,tempi$zij)
                id=c(id,j)
              } else {
                ci=c(ci,rep(sum(tempi$zij),nrow(tempi)))
                id=c(id,rep(j,nrow(tempi)))
              }
              j=j+1
            }
            response <- cbind(id, response[-1], ci)

            if (method=="Non-Parametric") {
              colnames(response) <- c("id", "vij", "wij", "d2", "d1", "epi", "x0ij", "ci")
              my_data <- response
            } else {
              colnames(response) <- c("id", "xij", "yij", "d2", "d1", "epi", "zij", "ci")
              my_data <- cbind(response, temp[, covariate_indexes])
              cov_names <- colnames(temp)[7:ncol(temp)]
            }

            #### Finish reformat based on method ####

            ## Use m.dat for Lee.et.al
            if (method == "Lee.et.al") {
              if (length(cov_names)==1) {
                colnames(my_data)[ncol(my_data)] <- "a1"
                fit_data <- m.dat.1(my_data)
              } else {
                colnames(my_data)[(ncol(response)+1):ncol(my_data)] <- paste("a", seq(1,length(cov_names)), sep="")
                fit_data <- m.dat(my_data)
              }
            #Return data frame for Chang
            } else {if (method == "Chang") {
              if (length(cov_names)==1) {
                colnames(my_data)[ncol(my_data)] <- "a1"
                fit_data <- my_data
              } else {
                colnames(my_data)[(ncol(response)+1):ncol(my_data)] <- paste("a", seq(1,length(cov_names)), sep="")
                fit_data <- my_data
              }
            #Use np.dat for Non-Parametric
            } else {
              forcdf <- np.dat(dat=my_data, ai=ai)
              marg_dat <- formarginal(dat = my_data)
              formarg <- np.dat(dat=marg_dat, ai=ai)
              fit_data <- list(forcdf=forcdf, formarg=formarg, refdata = my_data) #my_data is just response object df
              }
          }
        }
      }
    }
   }
} else {
  my_data = na.omit(data)
  forcdf <- np.dat(dat=my_data, ai=ai)
  fit_data <- list(forcdf=identifierforcdf, refdata = my_data)
}

  return(fit_data)
}
