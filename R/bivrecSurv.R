#                 m.dat, np.dat FUNCTIONS                                      #
#_______________________________________________________________________________
# Original by Chihyun Lee (August, 2017)                                       #
# Last Modified by Sandra Castro-Pearson (March, 2019)                         #
# Received from Chihyun Lee (January, 2018)                                    #
#_______________________________________________________________________________
##-----reformat dataset
mdat=function(dat) {
  n=length(unique(dat$id))
  mc=max(dat$epi)-1
  
  g1dat=cbind(dat[dat$epi==1,]$xij,1-dat[dat$epi==1,]$d1) #a vector of event 1 gap time and 1-d1 for all the first episodes for each subject
  g2dat=cbind(dat[dat$epi==1,]$zij,1-dat[dat$epi==1,]$d2) #a vector of total gap time and 1-d2 for all the first episodes for each subject
  l1=max(g1dat[g1dat[,2]==0,1])-(1e-07) #finding maximum event 1 gap time for first episode
  l2=max(g2dat[g2dat[,2]==0,1])-(1e-07) #finding the maximum total time for first episode
  g1surv=survfit(Surv(g1dat[,1],g1dat[,2])~1) #outputs the number of censored (?) events in hospital and the median survival time
  g2surv=survfit(Surv(g2dat[,1],g2dat[,2])~1) #outputs the number of total censored (?) events and the median survival time 
  
  xmat=ymat=zmat=delta1=delta2=g1mat=g2mat=matrix(0,n,mc,byrow=TRUE)
  mstar=ctime=NULL #mstar is the number of episodes a subject has 
  for (i in 1:n) { #for each subject
    tmp=dat[dat$id==i,] #dataframe of all episodes for a subject
    tmp.mstar=ifelse(nrow(tmp)==1,1,nrow(tmp)-1) #if 1 episode mstar=1, otherwise 6 episodes before censored
    mstar=c(mstar,tmp.mstar) #compiles a vector of the number of episodes before censoring for each person
    ctime=c(ctime,tmp$ci[1]) #censoring time for each subject 
    
    xmat[i,1:tmp.mstar]=tmp$xij[1:tmp.mstar] #a matrix of all the in hospital times for each episode before censoring 
    ymat[i,1:tmp.mstar]=tmp$yij[1:tmp.mstar] #a matrix of all the out of hospital times for each episode before censoring
    zmat[i,1:tmp.mstar]=tmp$zij[1:tmp.mstar] #a matrix of all the combined times for each episode before censoring 
    delta1[i,1:tmp.mstar]=tmp$d1[1:tmp.mstar] #this should all be 1's always?
    delta2[i,1:tmp.mstar]=tmp$d2[1:tmp.mstar] #this should also always be 1's?
    
    g1mat[i,1:tmp.mstar]=sapply(xmat[i,1:tmp.mstar],function(x)summary(g1surv,times=min(x,l1))$surv) #something for survival curves 
    g2mat[i,1:tmp.mstar]=sapply(zmat[i,1:tmp.mstar],function(x)summary(g2surv,times=min(x,l2))$surv)
  }
  
  cumh1=cumsum(g1surv$n.event/g1surv$n.risk) #cumulative survival probabilities?
  cumh2=cumsum(g2surv$n.event/g2surv$n.risk)
  l1mat=cbind(g1surv$time,diff(c(cumh1,tail(cumh1,1))),g1surv$surv) #idk what this is tbh
  l2mat=cbind(g2surv$time,diff(c(cumh2,tail(cumh2,1))),g2surv$surv)
  
  out=list(n=n,mc=mc,xmat=xmat,ymat=ymat,zmat=zmat,delta1=delta1,delta2=delta2,g1mat=g1mat,g2mat=g2mat,l1=l1,l2=l2,l1mat=l1mat,l2mat=l2mat, mstar=mstar,ctime=ctime)
  return(out)
}

#################### CREATE A BIVREC RESPONSE OBJECT ######################

#####
#' A function to create a bivrecSurv object
#'
#' @description
#' This function creates a bivariate recurrent response survival object (bivrecSurv object).
#'
#' @importFrom stats na.omit
#' @importFrom survival survfit
#'
#' @param id Vector of subject's unique identifier (i).
#' @param episode Vector indicating the observation or episode (j) for a subject (i). This will determine order of events for each subject.
#' @param xij Vector with the lengths of time spent in event of type X for individual i in episode j.
#' @param yij vector with the lengths of time spent in event of type Y for individual i in episode j.
#' @param Ycind Vector of indicators, with values of 0 for the last episode for subject i or 1 otherwise. A subject with only one episode will have one 0.
#' @param Xcind Vector of indicators, with values of 0 if the last episode for subject i occurred for event of type X or 1 otherwise. A subject with only one episode could have either one 1 (if he was censored at event Y) or one 0 (if he was censored at event X). A subject with censoring in event Y will have a vector of 1's.
#'
#' @return a BivRec repsonse object ready to put in a formula.
#' @seealso \code{\link{BivRec.fit}}
#'
#' @rdname BivRec
#' @export
#' @examples
#' set.seed(1)
#' dat <- biv.rec.sim(nsize=150, beta1=c(0.5,0.5), beta2=c(0,-0.5))
#' bdat<-with(dat, bivrecSurv(id, epi, xij, yij, d1, d2))
#'
bivrecSurv <- function(id, episode, xij, yij, Xcind, Ycind) {
  
  #Check if anything is missing
  if (missing(xij)) stop("Missing - gap times for type 1 event (xij).")
  if (missing(yij)) stop("Missing - gap times for type 2 event (yij).")
  if (missing(id)) stop("Missing - subject identifiers (id).")
  if (missing(episode)) stop("Missing - episodes for each subject (episode).")
  if (missing(Ycind)) stop("Missing - censoring indicator for type 2 event (Ycind).")
  if (missing(Xcind)) stop("Missing - censoring indicator for type 1 event (Xcind).")
  
  #Check all vectors have same length
  all_lengths <- c(length(id),length(episode),length(xij),length(yij),length(Ycind),length(Xcind))
  if (length(unique(all_lengths)) != 1) top("One or more input vectors (id, episode, xij, yij, Ycind, Xcind) differs in length from the rest.")
  
  #Check xij > 0 and yij >=0 both numeric vectors
  if (!is.numeric(xij)) stop("Time arguments (xij and yij) must be numeric.")
  if (!is.numeric(yij)) stop("Time arguments (xij and yij) must be numeric.")
  if (any(xij <= 0)) stop("Time arguments for event type 1 (xij) must be positive.")
  if (any(yij < 0)) stop("Time arguments for event type 2 (yij) must be non-negative")
  
  #Check censoring indicators are made of only 0 or 1 values
  if (any(Xcind!=0 && Xcind!=1)) stop("Indicator vector for type 1 gap times (Xcind) must be made of 0 or 1 values only.")
  if (any(Ycind!=0 && Ycind!=1)) stop("Indicator vector for type 2 gap times (Ycind) must be made of 0 or 1 values only.")
  
  inputdf <- data.frame(id=id, epi=episode, xij=xij, yij=yij, d1=Xcind, d2=Ycind)
  #Checks for each subject
  err_xind = err_yind = err_epi = NULL
  unique_id <- unique(inputdf$id)
  for (i in 1:length(unique_id)) {
    sub_id <- unique_id[i]
    temp_by_subject <- subset(inputdf, inputdf$id==sub_id)
    temp_by_subject <- temp_by_subject[order(temp_by_subject$epi),]
    sub_n <- nrow(temp_by_subject)
    last_cx <- temp_by_subject$d1[sub_n]
    last_cy <- temp_by_subject$d2[sub_n]
    other_cx <- temp_by_subject$d1[-sub_n]
    other_cy <- temp_by_subject$d2[-sub_n]
    
    #Check indicators (cx is all 1's or one zero at end /  cy last is 0 and all others are 1)
    if (sum(other_cx)!=(sub_n-1)) {err_xind <- c(err_xind, sub_id)}
    if (sum(other_cy)!=(sub_n-1)) {err_yind <- c(err_yind, sub_id)}
    if (last_cy!=0) {err_yind <- c(err_yind, sub_id)}
    if (last_cx==0) {if (last_cy==1) {err_yind <- c(err_yind, sub_id)}}
    
    #Check episodes don't have gaps
    for (j in 1:sub_n){
      if (temp_by_subject$epi[j]!=j) {
        err_epi <- c(err_epi, sub_id)}
    }
  }
  error_subjects <- unique(c(err_xind, err_yind, err_epi))
  if (length(error_subjects>0)){
    errmsg <- paste(error_subjects, collapse = ", ")
    msg <- paste("Subjects with id", errmsg,
                 "removed becuase of gaps in episodes or incorrect values for Xcind, Ycind.",
                 sep=" ")
    print(msg)
    df4mdat <- inputdf[-which(inputdf$id %in% error_subjects), ]
  } else {df4mdat <- inputdf}
  
  #calculate censoring time
  ci=id2=NULL
  j=1
  df4mdat$zij <- df4mdat$xij + df4mdat$yij
  for (i in unique(df4mdat$id)){
    tempi=df4mdat[df4mdat$id == i,]
    if (nrow(tempi) == 1){
      ci=c(ci,tempi$zij)
      id2=c(id2,j)
    } else {
      ci=c(ci,rep(sum(tempi$zij),nrow(tempi)))
      id2=c(id2,rep(j,nrow(tempi)))
    }
    j=j+1
  }
  df4mdat <- cbind(id=id2, df4mdat[-1], ci)
  
  result <- mdat(dat=df4mdat)
  class(result) <- "bivrecSurv"
  return(result)
}

#' @rdname bivrecSurv
#' @export

is.bivrecSurv <- function(x) inherits(x, "bivrecSurv")
is.BivRec <- function(x) inherits(x, "BivRec")
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol