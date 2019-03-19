  ###########################################################################
  ############## FUNCTIONS FOR REFERENCE BY MAIN - NOT FOR USER #############
  ###########################################################################

  ###                 changmdat FORTRAN subroutine and R call                    #
  #______________________________________________________________________________#
  # By Sandra Castro-Pearson (Last modified July, 2018)                          #
  # Based on Chihyun Lee's R code                                                #
  #______________________________________________________________________________#

  r2f.changmdat <- function(tugap1r, tugap2r, tdatr, idcount,
                            maxidcount, n, tmpr){

    tmpstart = tmpend = tugapstart = tugapend = 0

    out <- .Fortran("changmdat",
                    tugap1=as.double(tugap1r),
                    tugap2=as.double(tugap2r),
                    tdat=as.double(tdatr),
                    tmpin=as.double(tmpr),

                    ugapcols=as.integer(ncol(tugap1r)),
                    tmpstart=as.integer(tmpstart),
                    tmpend=as.integer(tmpend),
                    tugapstart=as.integer(tmpend),
                    tugapend=as.integer(tugapend),

                    idcount=as.integer(idcount),
                    nrowtdat=as.integer(nrow(tdatr)),
                    nrowugap=as.integer(nrow(tugap1r)),
                    maxidcount=as.integer(maxidcount),

                    n=as.integer(n),
                    tmpcols=as.integer(ncol(tmpr)),
                    tdatcols=as.integer(ncol(tdatr))
    )

    mytugap1 <- data.frame(matrix(out$tugap1, ncol=ncol(tugap1r)))
    mytugap2 <- data.frame(matrix(out$tugap2, ncol=ncol(tugap2r)))
    colnames(mytugap1) <- colnames(tugap1r)
    colnames(mytugap2) <- colnames(tugap2r)

    return(list(tugap1 = mytugap1, tugap2 = mytugap2))
  }

#           m.dat.chang1, all RE, v.est1 and sd.estpar1 FUNCTIONS              #
#_______________________________________________________________________________
# Original by Chihyun Lee (August, 2017)                                       #
# Last Modified by Sandra Castro-Pearson (June, 2018)                          #
# Received from Chihyun Lee (January, 2018)                                    #
#_______________________________________________________________________________

######################
##-----reformat dataset
m.dat.chang1=function(dat,beta) {
  n = length(unique(dat$id))
  mc = max(dat$epi)-1
  beta1 = beta[1]
  beta2 = beta[2]
  maxb = apply(cbind(beta1,beta2), 1, max)
  amat_indexes <- which(substr(colnames(dat), 1,1)=="a")
  amat = dat$a1
  dat$txij = dat$xij*exp(-amat*beta1)
  dat$tzij = dat$txij + dat$yij*exp(-amat*beta2)
  dat$tci = dat$ci*exp(-amat*maxb)

  tdatr <- as.matrix(dat[,c(1:2, 11:13, amat_indexes)])
  tugap1r = tugap2r = matrix(rep(0, n*max(dat$epi)*(3+length(amat_indexes))), ncol=5)
  colnames(tugap1r) = c("tgtime", "delta", "mstar", "a1")
  colnames(tugap2r) = c("tgtime", "delta", "mstar", "a1")
  idcount = as.vector(table(dat$id))
  maxidcount = max(idcount)
  tmpr <- matrix(rep(0, maxidcount*(9 + length(amat_indexes))), ncol=9 + length(amat_indexes))
  colnames(tmpr) = c(colnames(tdatr)[1:5], "uxij" , "uzij", "udx", "udz", "a1")

  f.ugaps <- r2f.changmdat(tugap1r, tugap2r, tdatr, idcount,
                           maxidcount, n, tmpr)
  tugap1 <- f.ugaps$tugap1
  tugap2 <- f.ugaps$tugap2

  ugap1 = data.frame(tugap1[-which(tugap1$tgtime==0),])
  ugap2 = data.frame(tugap2[-which(tugap2$tgtime==0),])

  #order
  ugap1=ugap1[order(ugap1$tgtime,decreasing=TRUE),]
  ugap2=ugap2[order(ugap2$tgtime,decreasing=TRUE),]
  out=list(n=n,ugap1=ugap1,ugap2=ugap2)

  return(out)
}

##-----point estimation
##rev-biv
RE.biv1=function(beta,dat) {
  mdat=m.dat.chang1(dat,beta)
  n=mdat$n
  ugap1=mdat$ugap1
  ugap2=mdat$ugap2

  ss10=cumsum(1/ugap1$mstar/n)
  ss11=cumsum(ugap1$a1/ugap1$mstar/n)
  sub1=ugap1$delta*(ugap1$a1-ss11/ss10)/ugap1$mstar/sqrt(n)

  ss20=cumsum(1/ugap2$mstar/n)
  ss21=cumsum(ugap2$a1/ugap2$mstar/n)
  sub2=ugap2$delta*(ugap2$a1-ss21/ss20)/ugap2$mstar/sqrt(n)

  out=c(sum(sub1),sum(sub2))
  return(out)
}

RE.uf1=function(beta,dat) {
  tmp.out=RE.biv1(beta,dat)
  out=tmp.out%*%tmp.out
  return(out)
}

RE.uest1=function(init,dat) {
  res=optim(init, RE.uf1, dat=dat, control=list(maxit=20000))
  return(list(par=res$par,value=res$value,conv=res$convergence))
}

##-----variance estimation
############################################
#Zeng
############################################
v.est1=function(beta, dat, R)
  #----------------------------------------------------------------------------------------------------------------------------
# first step of variance estimate: estimate V by bootstrap, only need to evaluate the estimating function, no need to solve it
#----------------------------------------------------------------------------------------------------------------------------
{
  id=dat$id
  ids=unique(dat$id)
  n=length(ids)
  freq=table(dat$id)
  index=cumsum( c(0, freq[-n]) )
  p=length(beta)

  A=matrix(rep(NA,R*p),ncol=p)
  for (i in 1:R)
  {
    w=table(sample(ids,n,replace=TRUE))
    s=as.numeric(names(w)) # because of this line, id must be 1:n
    w=as.numeric(w)
    location=NULL
    newid=NULL
    for(ss in 1:length(s))
    {
      location=c(location, rep(index[s[ss]]+(1:freq[s[ss]]),times=w[ss]))
      # since the same subject may be drawn multiple times, new id need to be created to distinguish different duplicates
      # e.g., the first duplicate's id will be original id+1000, the next will be id+2000, etc.
      newid=c(newid,rep(id[index[s[ss]]+(1:freq[s[ss]])],times=w[ss])+rep(((1:w[ss])-1)*1000,each=freq[s[ss]]))
    }
    dat.boot=dat[location,]
    dat.boot$id=newid
    A[i,]=RE.biv1(dat=dat.boot, beta=beta)
  }
  v=cov(A)
  return(v)
}

############################################
#parzen
############################################
#should do v.est first
############################################

RE.bivR1=function(beta,dat,R) {
  mdat=m.dat.chang1(dat,beta)
  n=mdat$n
  p=length(beta)
  ugap1=mdat$ugap1
  ugap2=mdat$ugap2

  ss10=cumsum(1/ugap1$mstar/n)
  ss11=cumsum(ugap1$a1/ugap1$mstar/n)
  sub1=ugap1$delta*(ugap1$a1-ss11/ss10)/ugap1$mstar/sqrt(n)

  ss20=cumsum(1/ugap2$mstar/n)
  ss21=cumsum(ugap2$a1/ugap2$mstar/n)
  sub2=ugap2$delta*(ugap2$a1-ss21/ss20)/ugap2$mstar/sqrt(n)

  out1=sum(sub1)-R[1]
  out2=sum(sub2)-R[2]
  out=c(out1,out2)
  return(out)
}

RE.ufR1=function(beta,dat,R) {
  tmp.out = RE.bivR1(beta,dat,R)
  out = tmp.out%*%tmp.out
  return(out)
}

RE.uestR1=function(init,dat,R) {
  res=optim(init, RE.ufR1, dat=dat, R=R,control=list(maxit=20000))
  return(list(par=res$par,value=res$value,conv=res$convergence))
}

sd.estpar1=function(init,dat,v, B) {
  p=length(init)
  A=matrix(rep(NA,B*p),ncol=p)
  i=0
  while (i < B)
  {
    R=mvrnorm(1,rep(0,p),v)
    est.R=RE.uestR1(init,dat,R)
    if (est.R$conv!=0) next
    i=i+1
    A[i,]=est.R$par
  }
  var.est=cov(A,A) #cov compute the cov between columns
  out=sqrt(diag(var.est))
  return(out)
}

###################################################################
#################### FUNCTION NOT FOR USER ########################
###################################################################
#' A Function for univariate fits using semiparametric regression method on a biv.rec object
#'
#' @description
#' This function fits the model using Chang's Method given one covariate. Called from biv.rec.fit(). No user interface.
#' @param new_data An object that has been reformatted for fit using the biv.rec.reformat() function. Passed from biv.rec.fit().
#' @param cov_names A string with the name of the covariate. Passed from biv.rec.fit().
#' @param CI Passed from biv.rec.fit().
#'
#' @return A dataframe summarizing covariate effect estimate, SE and CI.
#' @seealso \code{\link{biv.rec.fit}}
#'
#' @importFrom stats na.omit
#' @importFrom stats optim
#' @importFrom stats optimize
#' @importFrom stats qnorm
#' @importFrom stats cov
#' @importFrom MASS mvrnorm
#'
#' @keywords internal


#multivariable regression analysis-Chang's method
chang.univariate <- function(new_data, cov_names, CI) {

  print(paste("Fitting model with covariate", cov_names))
  beta <- rep(0, 2)

  #solve first equation to get all  estimates
  chang1 <- RE.uest1(beta, new_data)

  if (chang1$conv!=0) {
    print("Error - No Convergence")
    stop()
  }

  if (is.null(CI)==TRUE) {
    chang.fit <- data.frame(chang1$par)
    colnames(chang.fit) <- c("Estimate")
    rownames(chang.fit) <- c(paste("xij", cov_names), paste("yij", cov_names))

  } else {

    print("Point Estimates complete. Estimating Standard Errors/Confidence Intervals.")

    #estimate covariance matrix / std. errors
    chang1.v <- v.est1(chang1$par,new_data, R=100)
    chang1.sd <- sd.estpar1(beta, new_data, chang1.v ,B=50)

    #calculate CIs, join all info, put in nice table
    conf.lev = 1 - ((1-CI)/2)
    CIlow <- chang1$par - qnorm(conf.lev)*chang1.sd
    CIup <- chang1$par + qnorm(conf.lev)*chang1.sd
    chang.fit <- data.frame(chang1$par, chang1.sd, CIlow, CIup)
    low.string <- paste((1 - conf.lev), "%", sep="")
    up.string <- paste(conf.lev, "%", sep="")
    colnames(chang.fit) <- c("Estimate", "SE", low.string, up.string)
    rownames(chang.fit) <- c(paste("xij", cov_names), paste("yij", cov_names))

  }

  return(chang.fit)
}
