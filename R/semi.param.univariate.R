  ###########################################################################
  ############## FUNCTIONS FOR REFERENCE BY MAIN - NOT FOR USER #############
  ###########################################################################

#                  o.fun, all PRO and var.est FUNCTIONS                        #
#_______________________________________________________________________________
# Original by Chihyun Lee (August, 2017)                                       #
# Last Modified by Sandra Castro-Pearson (April, 2018)                         #
# Received from Chihyun Lee (January, 2018)                                    #
#_______________________________________________________________________________

##------symmetric O function
o.fun=function(t,s,L) {log(min(max(t,s),L))-log(L)}

##-----estimation functions

##proposed method
Pro.ee1=function(beta1,mdat) {
  n=mdat$n
  xmat=mdat$xmat
  delta1=mdat$delta1
  g1mat=mdat$g1mat
  l1=mdat$l1
  mstar=mdat$mstar
  amat=mdat$amat

  tmp.out=NULL
  for (i in 1:n) {
    A=(amat)-amat[i]
    expA=exp(A*beta1)
    subsum=sapply(expA,function(x)mean(delta1[i,1:mstar[i]]*sapply(xmat[i,1:mstar[i]],function(t)o.fun(t,x*t,l1))/g1mat[i,1:mstar[i]]))
    tmp.out=c(tmp.out,sum(A*subsum))
  }
  out=sum(tmp.out)/(n^2)
  return(out)
}

Pro.uf1=function(beta1,mdat) {
  tmp.out=Pro.ee1(beta1,mdat)
  out=tmp.out%*%tmp.out
  return(out)
}

Pro.uest1=function(int,mdat) {
  res=optimize(Pro.uf1,interval=int,mdat=mdat)
  return(list(par=res$minimum,value=res$objective))
}

Pro.ee2=function(beta2,beta1,mdat) {
  n=mdat$n
  xmat=mdat$xmat
  ymat=mdat$ymat
  delta2=mdat$delta2
  g2mat=mdat$g2mat
  l2=mdat$l2
  mstar=mdat$mstar
  amat=mdat$amat

  tmp.out=NULL
  for (i in 1:n) {
    A=(amat)-amat[i]
    expA1=exp(A*beta1)
    expA2=exp(A*beta2)
    expA=cbind(expA1,expA2)
    subsum=apply(expA,1,function(x)mean(delta2[i,1:mstar[i]]*apply(cbind(xmat[i,1:mstar[i]],ymat[i,1:mstar[i]]),1,function(t)o.fun(sum(t),x[1]*t[1]+x[2]*t[2],l2))/g2mat[i,1:mstar[i]]))
    tmp.out=c(tmp.out,sum(A*subsum))
  }
  out=sum(tmp.out)/(n^2)
  return(out)
}

Pro.uf2=function(beta2,beta1,mdat) {
  tmp.out=Pro.ee2(beta2,beta1,mdat)
  out=tmp.out%*%tmp.out
  return(out)
}

Pro.uest2=function(int,beta1,mdat) {
  res=optimize(Pro.uf2,interval=int,beta1=beta1,mdat=mdat)
  return(list(par=res$minimum,value=res$objective))
}


##variance estimation
var.est=function(beta1,beta2,mdat) {
  n=mdat$n
  xmat=mdat$xmat
  ymat=mdat$ymat
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
    A=(amat)-amat[i]
    expA1=exp(A*beta1)
    expA2=exp(A*beta2)
    expA=cbind(expA1,expA2)

    sub1.xi1=sapply(expA1,function(x) mean(delta1[i,1:mstar[i]]*sapply(xmat[i,1:mstar[i]],function(t)o.fun(t,x*t,l1))/g1mat[i,1:mstar[i]]))
    sub2.xi1=rep(0,n)
    for (j in 1:n) {
      sub2.xi1[j]=mean(delta1[j,1:mstar[j]]*sapply(xmat[j,1:mstar[j]],function(t)o.fun(t,t/expA1[j],l1))/g1mat[j,1:mstar[j]])
    }
    tmp.xi1=sum(A*(sub1.xi1-sub2.xi1))/(n^(3/2))

    sub1.xi2=apply(expA,1,function(x) mean(delta2[i,1:mstar[i]]*apply(cbind(xmat[i,1:mstar[i]],ymat[i,1:mstar[i]]),1,function(t)o.fun(sum(t),x[1]*t[1]+x[2]*t[2],l2))/g2mat[i,1:mstar[i]]))
    sub2.xi2=rep(0,n)
    for (j in 1:n) {
      sub2.xi2[j]=mean(delta2[j,1:mstar[j]]*apply(cbind(xmat[j,1:mstar[j]],ymat[j,1:mstar[j]]),1,function(t)o.fun(sum(t),t[1]/expA1[j]+t[2]/expA2[j],l2))/g2mat[j,1:mstar[j]])
    }
    tmp.xi2=sum(A*(sub1.xi2-sub2.xi2))/(n^(3/2))

    xi=xi+c(tmp.xi1,tmp.xi2)%o%c(tmp.xi1,tmp.xi2)

    Amat=sapply(A,function(x) x%o%x)
    tmp.sub.gam1=apply(cbind(xmat[i,1],expA1*xmat[i,1],l1),1,function(x)(x[1]<=x[2])*(max(x[1],x[2])<=x[3]))
    tmp.sub.gam2=apply(cbind(xmat[i,1]+ymat[i,1],expA1*xmat[i,1]+expA2*ymat[i,1],l2),1,function(x)(x[1]<=x[2])*(max(x[1],x[2])<=x[3]))
    sub.gam1=(Amat)*tmp.sub.gam1*mean(delta1[i,1:mstar[i]]/g1mat[i,1:mstar[i]])
    sub.gam21=(Amat)*tmp.sub.gam2*apply(expA,1,function(x) mean(delta2[i,1:mstar[i]]*(x[1]*xmat[i,1:mstar[i]])/((x[1]*xmat[i,1:mstar[i]]+x[2]*ymat[i,1:mstar[i]])*g2mat[i,1:mstar[i]])))
    sub.gam22=(Amat)*tmp.sub.gam2*apply(expA,1,function(x) mean(delta2[i,1:mstar[i]]*(x[2]*ymat[i,1:mstar[i]])/((x[1]*xmat[i,1:mstar[i]]+x[2]*ymat[i,1:mstar[i]])*g2mat[i,1:mstar[i]])))

    gam1=gam1+sum(sub.gam1)/(n^2)
    gam21=gam21+sum(sub.gam21)/(n^2)
    gam22=gam22+sum(sub.gam22)/(n^2)
  }
  gam1=matrix(gam1,length(beta1),length(beta1))
  gam21=matrix(gam21,length(beta2),length(beta2))
  gam22=matrix(gam22,length(beta2),length(beta2))
  gamm=rbind(cbind(gam1,matrix(0,length(beta1),length(beta2))),cbind(gam21,gam22))

  mat=solve(gamm)%*%xi%*%t(solve(gamm))

  var1=sqrt(diag(mat)/n)[1]
  var2=sqrt(diag(mat)/n)[2]
  return(list(var1=var1,var2=var2))
}

###################################################################
######################## FUNCTION FOR USE ########################
###################################################################
#' A Function for univariate fits using semiparametric regression method on a biv.rec object
#'
#' @description
#' This function fits the semiparametric model given only one covariate. Called from biv.rec.fit(). No user interface.
#' @param new_data An object that has been reformatted for fit using the biv.rec.reformat() function. Passed from biv.rec.fit().
#' @param cov_names A vector with the names of the single covariate. Passed from biv.rec.fit().
#' @param CI Logical. Passed from biv.rec.fit().
#' @return A dataframe summarizing the estimates for effects of the covariate, their standard errors and 95% confidence intervals.
#' @seealso \code{\link{biv.rec.fit}}
#'
#' @importFrom stats na.omit
#' @importFrom stats optim
#' @importFrom stats optimize
#' @importFrom stats qnorm
#' @importFrom stats rbinom
#' @importFrom stats rgamma
#' @importFrom stats rnorm
#' @importFrom stats runif
#'
#' @keywords internal

#MAIN PROGRAM FOR univariate regression analysis
semi.param.univariate <- function(new_data, cov_names, CI) {

  print(paste("fitting model with covariate", cov_names))
  pro1 <- Pro.uest1(c(-2,2),new_data)[[1]]
  print("First estimate complete. Second estimation in progress.")
  pro2 <- Pro.uest2(c(-2,2), pro1, new_data)[[1]]

  if (CI == TRUE) {
    print("Second estimate complete. Estimating Variance")
    variance_est=var.est(pro1, pro2, new_data)
    print("Estimations complete. Calculating confidence intervals")
    univ_fits <- data.frame(c(pro1, pro2), c(variance_est$var1,variance_est$var2))
    CI <- t(apply(univ_fits, 1, function (x) c(x[1]+qnorm(0.025)*x[2], x[1]+qnorm(0.975)*x[2])))
    univ_fits <- cbind(univ_fits, CI)
    colnames(univ_fits) <- c("Estimate", "SE", "0.25%", "0.95%")

  } else {
    univ_fits <- data.frame(c(pro1, pro2))
    colnames(univ_fits) <- c("Estimate")
  }

  rownames(univ_fits) <- c("xij", "yij")
  return(univ_fits)
}


