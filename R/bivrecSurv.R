#                 m.dat, np.dat FUNCTIONS                                      #
#_______________________________________________________________________________
# Original by Chihyun Lee (August, 2017)                                       #
# Last Modified by Sandra Castro-Pearson and Aparajita Sur (March, 2019)       #
# Received from Chihyun Lee (January, 2018)                                    #
#______________________________________________________________________________#


##-----reformat dataset for Lee Regression
mdat=function(dat) {
  n=length(unique(dat$id))
  mc=max(dat$epi)-1
  
  g1dat=cbind(dat[dat$epi==1,]$xij,1-dat[dat$epi==1,]$d1)
  g2dat=cbind(dat[dat$epi==1,]$zij,1-dat[dat$epi==1,]$d2)
  l1=max(g1dat[g1dat[,2]==0,1])-(1e-07)
  l2=max(g2dat[g2dat[,2]==0,1])-(1e-07)
  g1surv=survfit(Surv(g1dat[,1],g1dat[,2])~1)
  g2surv=survfit(Surv(g2dat[,1],g2dat[,2])~1)
  
  xmat=ymat=zmat=delta1=delta2=g1mat=g2mat=matrix(0,n,mc,byrow=TRUE)
  mstar=ctime=NULL
  for (i in 1:n) {
    tmp=dat[dat$id==i,]
    tmp.mstar=ifelse(nrow(tmp)==1,1,nrow(tmp)-1)
    mstar=c(mstar,tmp.mstar)
    ctime=c(ctime,tmp$ci[1])
    
    xmat[i,1:tmp.mstar]=tmp$xij[1:tmp.mstar]
    ymat[i,1:tmp.mstar]=tmp$yij[1:tmp.mstar]
    zmat[i,1:tmp.mstar]=tmp$zij[1:tmp.mstar]
    delta1[i,1:tmp.mstar]=tmp$d1[1:tmp.mstar]
    delta2[i,1:tmp.mstar]=tmp$d2[1:tmp.mstar]
    
    g1mat[i,1:tmp.mstar]=sapply(xmat[i,1:tmp.mstar],function(x)summary(g1surv,times=min(x,l1))$surv)
    g2mat[i,1:tmp.mstar]=sapply(zmat[i,1:tmp.mstar],function(x)summary(g2surv,times=min(x,l2))$surv)
  }
  
  cumh1=cumsum(g1surv$n.event/g1surv$n.risk)
  cumh2=cumsum(g2surv$n.event/g2surv$n.risk)
  l1mat=cbind(g1surv$time,diff(c(cumh1,tail(cumh1,1))),g1surv$surv)
  l2mat=cbind(g2surv$time,diff(c(cumh2,tail(cumh2,1))),g2surv$surv)
  
  out=list(n=n,mc=mc,xmat=xmat,ymat=ymat,zmat=zmat,delta1=delta1,delta2=delta2,g1mat=g1mat,g2mat=g2mat,l1=l1,l2=l2,l1mat=l1mat,l2mat=l2mat, mstar=mstar,ctime=ctime)
  return(out)
}

#####Reformat data set for non-parametric analysis

# dat : a data.frame including
#     1)id numbers, 2)orders of episodes, 3)first gap time, 4)second gap time
#     5)censoring times, 6) censoring indicators in each column
# ai: a non-negative function of censoring time

np.dat <- function(dat, ai) {
  
  id <- dat$id
  uid <- unique(id)   # vector of unique id's
  n.uid <- length(uid)   # scalar : number of unique IDs
  event <- dat$d2 #event indicator : must always be 0 for the last obs per ID and 1 otherwise
  markvar1 <- dat$vij #gap times of type 1
  markvar2 <- dat$wij #gap times of type 2
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
  
  mdata <- tmp <- NULL
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

#################### CREATE A BIVREC OBJECT ######################

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
#' set.seed(1234)
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
  
  result <- list()
  result$dat4Lreg <- mdat(dat=df4mdat) #data for Lee regression
  result$df <- df4mdat #data for Chang regression (this is also the df that is used in bivrecPlot)
  
  #####ADD data for cdf and marginal of NP model
  df4np <- df4mdat
  colnames(df4np)=c("id", "epi", "vij", "wij", "d1", "d2", "x0ij", "ci")
  df4np=df4np[,c("id","vij","wij","d2","d1","epi","x0ij","ci")] #change order of columns 
  forcdf1 <- np.dat(df4np, ai=1)
  forcdf2 <- np.dat(df4np, ai=2)
  marg1 <- formarginal(dat = df4np) #this is from the reformat code 
  marg2 <- formarginal(dat = df4np)
  formarg1 <- np.dat(marg1, ai=1)
  formarg2 <- np.dat(marg2, ai=2)
  #two np objects that have data for cdf and marg depending on ai
  result$dat4np1 <- list(forcdf=forcdf1, formarg=formarg1,refdata = df4np) #for ai=1
  result$dat4np2 <- list(forcdf=forcdf2, formarg=formarg2,refdata = df4np) #for ai=2
  
  class(result) <- "bivrecSurv"
  return(result)
}

#' @rdname bivrecSurv
#' @export

is.bivrecSurv <- function(x) inherits(x, "bivrecSurv")
is.bivrecReg <- function(x) inherits(x, "bivrecReg")
is.bivrecNP <- function(x) inherits(x, "bivrecNP")
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol