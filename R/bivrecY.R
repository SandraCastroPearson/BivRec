#bivrecY
#' @name bivrecY
#' @rdname bivrecY
#' @title Create a \code{bivrecFit} Object
#'
#' @description
#' This function creates a bivariate recurrent event response.
#' @importFrom stats na.omit
#' 
#' @param formula A formula with six variables indicating the bivariate alternating gap time response on the left of the ~ operator and the covariates on the right.
#' The six variables on the left must have the same length and be given as \strong{ID + episode +  xij + yij + delta_x + delta_y ~ covariates}, where
#' \itemize{
#'   \item ID: A vector of subjects' unique identifier which can be numeric or character.
#'   \item episode: A vector indicating the episode of the bivariate alternating gap time pairs, e.g.: 1, 2, ..., m_i where m_i indicates the last episode for subject i.
#'   \item xij: A vector with the lengths of the type I gap times.
#'   \item yij: A vector with the lengths of the type II gap times.
#'   \item delta_x: A vector of indicators with values 
#'   \itemize{
#'       \item 0 for the last episode for subject i (m_i) if subject was censored during period xij.
#'       \item 1 otherwise.
#'      }
#'   A subject with only one episode (m_i=1) could have a 0 if he was censored during period xi1 or 1 if he was censored during period yi1. 
#'   If delta_x is not provided estimation proceeds with the assumption that no subject was censored during period xij.
#'   \item delta_y: A vector of indicators with values
#'   \itemize{
#'       \item 0 for the last episode of subject i (m_i) if subject was censored during period yij.
#'       \item 1 otherwise.
#'      }
#'   A subject with only one episode (m_i=1) will have one 0.
#'   }
#' @param data A data frame that includes all the vectors/covariates listed in the formula above.
#' @return a bivrec object ready for fitting.
#'
#' @keywords bivrecY
bivrecY <- function(formula=ID+episode+xij+yij+delta_x+delta_y, data){
  
  ### PUT POTENTIAL RESPONSE TOGETHER
  temp1 <- data.frame(ID+episode+xij+yij+delta_x+delta_y)
  
  ####CHECK MISSINGNESS AND KEEP RECORD OF WHAT IS BEING OMITTED
  temp <- na.omit(temp1)
  n_missing <- nrow(temp1) - nrow(temp)
  
  ####CHECK xij, yij VARIABLES HAVE CORRECT VALUES
  invalid_xij <- which(temp$xij<=0)
  if (length(invalid_xij)>0) {
    print("Invalid values for length of time in event X for rows. All must be >0.")
    temp[invalid_xij,]
    stop()
  } 
  invalid_yij <- which(temp$yij<0)
    if (length(invalid_yij)>0) {
      print("Invalid values for length of time in event Y for rows. All must be >=0.")
      temp[invalid_yij,]
      stop()
    }
      ####CHECK INDICATORS AND EPISODE SEQUENCES
      #First check for indicators - all 0 and 1
      dx_check1 <- length(which(temp$delta_X!=1&&temp$delta_X!=0))
      dy_check1 <- length(which(temp$delta_Y!=1&&temp$delta_Y!=0))
      if (dx_check1!=0) {
        print("Invalid values for delta_X")
        temp[which(delta_X!=1&&delta_X!=0),]
        stop()
      } 
      if (dy_check1!=0) {
          print("Invalid values for delta_Y")
          temp[which(delta_Y!=1&&delta_Y!=0),]
          return("Error")
      }
      ## Second check for indicators
      ## Indicators match episode (last is 0 or 1 for X and 0 for Y) and episode doesn't have gaps
        wrong_xind <- NULL
        wrong_yind <- NULL
        wrong_epi <- NULL
        unique_id <- unique(temp$identifier)
        for (i in 1:length(unique_id)) {
          sub_id <- unique_id[i]
          temp_by_subject <- subset(temp, temp$identifier==sub_id)
          temp_by_subject <- temp_by_subject[order(temp_by_subject$episode),]
          sub_n <- nrow(temp_by_subject)
          last_dx <- temp_by_subject$delta_X[sub_n]
          rest_of_dx <- temp_by_subject$delta_X[-sub_n]
          last_dy <- temp_by_subject$delta_Y[sub_n]
          rest_of_dy <- temp_by_subject$delta_Y[-sub_n]
            
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
          
          if (length(wrong_epi)!=0 || length(wrong_xind)!=0 || length(wrong_yind)!=0) {
            if (length(wrong_epi)!=0) {
              wrong_epi <- unique(wrong_epi)
              print(paste("The subjects with following ID's have invalid episode sequences"))
              for (w in 1:length(wrong_epi)) {
                print(subset(temp, temp$identifier==wrong_epi[w]))
              }
              stop()
            } 
              if (length(wrong_xind)!=0 || length(wrong_yind)!=0) {
                wrong_ind <- c(unique(wrong_xind), unique(wrong_yind))
                print(paste("The subjects with following ID's have invalid censoring times"))
                for (w2 in 1:length(wrong_ind)) {
                  print(subset(temp, temp$identifier==wrong_xind[w2]))
                }
                stop()
              }
          } else {
            
            print(paste("Original number of observations:", nrow(temp1), "for", length(unique(temp1$identifier)), "individuals", sep=" "))
            print(paste("Observations to be used in analysis:", nrow(temp), "for", length(unique_id), "individuals",sep=" "))

          }
        #Creating response object and adding the censoring time "ci"
        response <- temp
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
      #class(response)<-"bivrec"
      return(response)
}