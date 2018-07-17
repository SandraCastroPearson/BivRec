  ###########################################################################
  ############## FUNCTIONS FOR REFERENCE BY MAIN - NOT FOR USER #############
  ###########################################################################

#                 o.fun, all MPRO and MVAR FUNCTIONS                           #
#_______________________________________________________________________________
# Original by Chihyun Lee (August, 2017)                                       #
# Last Modified by Sandra Castro-Pearson (April, 2018)                         #
# Received from Chihyun Lee (January, 2018)                                    #
#_______________________________________________________________________________

##------symmetric O function
o.fun=function(t,s,L) {log(min(max(t,s),L))-log(L)}

##-----estimation functions
##proposed method
MPro.ee1=function(beta1,mdat) {
  n=mdat$n
  xmat=mdat$xmat
  #ymat=mdat$ymat
  #zmat=mdat$zmat
  delta1=mdat$delta1
  #delta2=mdat$delta2
  g1mat=mdat$g1mat
  #g2mat=mdat$g2mat
  l1=mdat$l1
  #l2=mdat$l2
  mstar=mdat$mstar
  amat=mdat$amat

  tmp.out=NULL
  for (i in 1:n) {
    A=t(t(amat)-amat[i,])
    expA=apply(A,1,function(x)exp(x%*%beta1))
    subsum=sapply(expA,function(x)mean(delta1[i,1:mstar[i]]*sapply(xmat[i,1:mstar[i]],function(t)o.fun(t,x*t,l1))/g1mat[i,1:mstar[i]]))
    tmp.out=rbind(tmp.out,apply(A*subsum,2,sum))
  }
  out=apply(tmp.out,2,sum)/(n^2)
  return(out)
}

MPro.uf1=function(beta1,mdat) {
  tmp.out=MPro.ee1(beta1,mdat)
  out=tmp.out%*%tmp.out
  return(out)
}

MPro.uest1=function(init,mdat) {
  res=optim(init, MPro.uf1, mdat=mdat,control=list(maxit=20000))
  return(list(par=res$par,value=res$value,conv=res$convergence))
}

MPro.ee2=function(beta2,beta1,mdat) {
  n=mdat$n
  xmat=mdat$xmat
  ymat=mdat$ymat
  #zmat=mdat$zmat
  #delta1=mdat$delta1
  delta2=mdat$delta2
  #g1mat=mdat$g1mat
  g2mat=mdat$g2mat
  #l1=mdat$l1
  l2=mdat$l2
  mstar=mdat$mstar
  amat=mdat$amat

  tmp.out=NULL
  for (i in 1:n) {
    A=t(t(amat)-amat[i,])
    expA1=apply(A,1,function(x)exp(x%*%beta1))
    expA2=apply(A,1,function(x)exp(x%*%beta2))
    expA=cbind(expA1,expA2)
    subsum=apply(expA,1,function(x)mean(delta2[i,1:mstar[i]]*apply(cbind(xmat[i,1:mstar[i]],ymat[i,1:mstar[i]]),1,function(t)o.fun(sum(t),x[1]*t[1]+x[2]*t[2],l2))/g2mat[i,1:mstar[i]]))
    tmp.out=rbind(tmp.out,apply(A*subsum,2,sum))
  }
  out=apply(tmp.out,2,sum)/(n^2)
  return(out)
}

MPro.uf2 <- function(beta2,beta1,mdat) {
  tmp.out <- MPro.ee2(beta2,beta1,mdat)
  out <- tmp.out%*%tmp.out
  return(out)
}

MPro.uest2 <- function(init,beta1,mdat) {
  res <- optim(init, MPro.uf2, beta1=beta1, mdat=mdat, control=list(maxit=20000))
  return(list(par=res$par,value=res$value,conv=res$convergence))
}

##variance estimation
Mvar.est=function(beta1,beta2,mdat) {
  n=mdat$n
  xmat=mdat$xmat
  ymat=mdat$ymat
  #zmat=mdat$zmat
  delta1=mdat$delta1
  delta2=mdat$delta2
  g1mat=mdat$g1mat
  g2mat=mdat$g2mat
  l1=mdat$l1
  l2=mdat$l2
  mstar=mdat$mstar
  amat=mdat$amat

  xi=matrix(0,length(c(beta1,beta2)),length(c(beta1,beta2)))
  gam1=gam21=gam22=rep(0,length(beta1))

  for (i in 1:n) {
    A=t(t(amat)-amat[i,])
    expA1=apply(A,1,function(x)exp(x%*%beta1))
    expA2=apply(A,1,function(x)exp(x%*%beta2))
    expA=cbind(expA1,expA2)

    sub1.xi1=sapply(expA1,function(x) mean(delta1[i,1:mstar[i]]*sapply(xmat[i,1:mstar[i]],function(t)o.fun(t,x*t,l1))/g1mat[i,1:mstar[i]]))
    sub1.xi2=apply(expA,1,function(x) mean(delta2[i,1:mstar[i]]*apply(cbind(xmat[i,1:mstar[i]],ymat[i,1:mstar[i]]),1,function(t)o.fun(sum(t),x[1]*t[1]+x[2]*t[2],l2))/g2mat[i,1:mstar[i]]))

    sub2.xi1=rep(0,n)
    sub2.xi2=rep(0,n)

    for (j in 1:n) {
      sub2.xi1[j]=mean(delta1[j,1:mstar[j]]*sapply(xmat[j,1:mstar[j]],function(t)o.fun(t,t/expA1[j],l1))/g1mat[j,1:mstar[j]])
      sub2.xi2[j]=mean(delta2[j,1:mstar[j]]*apply(cbind(xmat[j,1:mstar[j]],ymat[j,1:mstar[j]]),1,function(t)o.fun(sum(t),t[1]/expA1[j]+t[2]/expA2[j],l2))/g2mat[j,1:mstar[j]])
    }

    tmp.xi1=apply(A*(sub1.xi1-sub2.xi1),2,sum)/(n^(3/2))
    tmp.xi2=apply(A*(sub1.xi2-sub2.xi2),2,sum)/(n^(3/2))

    xi=xi+c(tmp.xi1,tmp.xi2)%o%c(tmp.xi1,tmp.xi2)

    Amat=apply(A,1,function(x) x%o%x)
    tmp.sub.gam1=apply(cbind(xmat[i,1],expA1*xmat[i,1],l1),1,function(x)(x[1]<=x[2])*(max(x[1],x[2])<=x[3]))
    tmp.sub.gam2=apply(cbind(xmat[i,1]+ymat[i,1],expA1*xmat[i,1]+expA2*ymat[i,1],l2),1,function(x)(x[1]<=x[2])*(max(x[1],x[2])<=x[3]))
    sub.gam1=t(Amat)*tmp.sub.gam1*mean(delta1[i,1:mstar[i]]/g1mat[i,1:mstar[i]])
    sub.gam21=t(Amat)*tmp.sub.gam2*apply(expA,1,function(x) mean(delta2[i,1:mstar[i]]*(x[1]*xmat[i,1:mstar[i]])/((x[1]*xmat[i,1:mstar[i]]+x[2]*ymat[i,1:mstar[i]])*g2mat[i,1:mstar[i]])))
    sub.gam22=t(Amat)*tmp.sub.gam2*apply(expA,1,function(x) mean(delta2[i,1:mstar[i]]*(x[2]*ymat[i,1:mstar[i]])/((x[1]*xmat[i,1:mstar[i]]+x[2]*ymat[i,1:mstar[i]])*g2mat[i,1:mstar[i]])))

    gam1=gam1+apply(sub.gam1,2,sum)/(n^2)
    gam21=gam21+apply(sub.gam21,2,sum)/(n^2)
    gam22=gam22+apply(sub.gam22,2,sum)/(n^2)
  }

  gam1=matrix(gam1,length(beta1),length(beta1))
  gam21=matrix(gam21,length(beta2),length(beta2))
  gam22=matrix(gam22,length(beta2),length(beta2))
  gamm=rbind(cbind(gam1,matrix(0,length(beta1),length(beta2))),cbind(gam21,gam22))

  mat=solve(gamm)%*%xi%*%t(solve(gamm))
  se = sqrt(diag(mat)/n)
  return(list(se, mat))
}

###################################################################
##################### FUNCTION NOT FOR USER #######################
###################################################################
#' A Function for multivariate fits using semiparametric regression method on a biv.rec object
#'
#' @description
#' This function fits the semiparametric model given multiple  covariates. Called from biv.rec.fit(). No user interface.
#' @param new_data An object that has been reformatted for fit using the biv.rec.reformat() function. Passed from biv.rec.fit().
#' @param cov_names A vector with the names of the covariates. Passed from biv.rec.fit().
#' @param CI Logical. Passed from biv.rec.fit().
#' @return A dataframe summarizing the estimates for effects of the covariates, their standard errors and 95% confidence intervals.
#' @seealso \code{\link{biv.rec.fit}}
#'
#' @importFrom stats na.omit
#' @importFrom stats optim
#' @importFrom stats optimize
#' @importFrom stats qnorm
#' @importFrom stringr str_c
#'
#' @export


#multivariable regression analysis
semi.param.multivariate <- function(new_data, cov_names, CI) {

  print(paste("fitting model with covariates:", str_c(cov_names, collapse = ","), sep=" "))
  n_params <- length(cov_names)

  #solve first equation to get beta1
  mpro1 <- MPro.uest1(init=rep(0, n_params),mdat=new_data)
  print(paste("First estimate complete.", str_c(mpro1$par, collapse = ","), "Second estimation in progress.", sep=" "))

  #solve second equation to get beta2
  mpro2 <- MPro.uest2(init=rep(0, n_params), beta1=mpro1$par, mdat=new_data)


  if (CI == TRUE) {
    print(paste("Second estimate complete.", str_c(mpro2$par, collapse = ","), "Estimating standard error", sep=" "))
    #estimate covariance matrix and get diagonal then std. errors
    se_est <- Mvar.est(beta1=mpro1$par, beta2=mpro2$par, mdat=new_data)
    print("Estimation of std. errors complete.")
    print(se_est[[1]])
    print("Calculating confidence intervals")

    #join all info and calculate CIs, put in nice table
    multi.fit <- data.frame(c(mpro1$par, mpro2$par), se_est[[1]])
    CI <- t(apply(multi.fit, 1, function (x) c(x[1]+qnorm(0.025)*x[2], x[1]+qnorm(0.975)*x[2])))
    multi.fit  <- cbind(multi.fit, CI)
    colnames(multi.fit) <- c("Estimate", "SE", "0.25%", "0.95%")

  } else {
    multi.fit <- data.frame(c(mpro1$par, mpro2$par))
    colnames(multi.fit) <- "Estimate"
  }

  rownames(multi.fit) <- c(paste("xij", cov_names), paste("yij", cov_names))
  return(multi.fit)
}
