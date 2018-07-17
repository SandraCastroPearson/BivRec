#' A Data Simulation Function
#'
#' @description This function simulates a series of alternating recurrent events based on a linear combination of two covariates.
#'
#' @param nsize sample size
#' @param beta1 true coefficients for first gap time
#' @param beta2 true coefficients for second gap time
#' @param cr maximum support of censoring time
#' @param sg2 variance of frailty
#' @param set optional simulation setting. Choose 1.1 (default) for $rho=1$ in covariance matrix, 1.2 for $rho=0.5$ in covariance matrix or 2.1 for $rho=0$ in covariance matrix.
#'
#' @return Data frame with alternating recurrent event data and two covariates
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats qnorm
#' @importFrom stats rbinom
#' @importFrom stats rgamma
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom utils tail
#'
#' @keywords data.sim, simulation
#'
#' @examples
#' data.sim(nsize=300, beta1=c(0.5,0.5), beta2=c(0,-0.5), cr=63, sg2=0.5, set=1.1)
#'
#' @export

##-----data generation
data.sim <- function(nsize,beta1,beta2,cr,sg2,set) {

  if (missing(set)) {set <- 1.1}
  id=1:nsize

  ##generate covariates (A1,A2)
  A1=rbinom(nsize,1,0.5)
  A2=runif(nsize)
  A=cbind(A1,A2)

  ##generate gamma
  if (set==1.1) {
    Sig=matrix(c(sg2,sqrt(sg2)*sqrt(sg2)*1,sqrt(sg2)*sqrt(sg2)*1,sg2),2,2)
    gamma=mvrnorm(nsize,c(1,1),Sig)
    gamma1=gamma[,1]
    gamma2=gamma[,2]
  }
  if (set==1.2) {
    Sig=matrix(c(sg2,sqrt(sg2)*sqrt(sg2)*0.5,sqrt(sg2)*sqrt(sg2)*0.5,sg2),2,2)
    gamma=mvrnorm(nsize,c(1,1),Sig)
    gamma1=gamma[,1]
    gamma2=gamma[,2]
  }
  if (set==2.1) {
    gamma1=rnorm(nsize,1,sqrt(sg2))
    gamma2=rnorm(nsize,1,sqrt(sg2))
  }
  if (set==2.2) {
    gamma1=rnorm(nsize,1,sqrt(sg2))
    gamma2=rgamma(nsize,shape=1/sg2,rate=1/sg2)
  }

  dat=NULL
  deltas=NULL
  for (i in id) {
    x.tmp=exp(gamma1[i]+A[i,]%*%beta1)
    y.tmp=exp(gamma2[i]+A[i,]%*%beta2)

    #ci=rexp(1,1/cr)
    ci=runif(1,0,cr)

    sum.z.tmp=0
    dat.tmp=NULL
    while(sum.z.tmp<ci) {
      #generate alternating gap times until C_i
      xij=x.tmp*exp(rnorm(1,0,sqrt(0.1)))
      yij=y.tmp*exp(rnorm(1,0,sqrt(0.1)))

      dat.tmp=rbind(dat.tmp,c(xij,yij,sum.z.tmp))
      sum.z.tmp=sum.z.tmp+(xij+yij)
    }
    epi=nrow(dat.tmp)
    #set censoring indicator and last censored gap time pair
    if (epi==1) {cen=cbind(1,0)} else {
      cen=cbind(rep(1,epi),c(rep(1,epi-1),0))}
    if (dat.tmp[epi,1]>(ci-dat.tmp[epi,3])) {
      dat.tmp[epi,1]=ci-dat.tmp[epi,3]
      dat.tmp[epi,2]=0
      cen[epi,1]=0
    } else {
      dat.tmp[epi,2]=ci-dat.tmp[epi,3]-dat.tmp[epi,1]
    }
    subdat=cbind(i,1:epi,dat.tmp[,1],dat.tmp[,2],dat.tmp[,1]+dat.tmp[,2],ci, cen,matrix(A[i,],epi,2,byrow=TRUE))
    dat=rbind(dat,subdat)
    deltas=rbind(deltas, cen)
  }
  dat=as.data.frame(dat)
  colnames(dat)=c("id","epi","xij","yij","zij","ci","d1", "d2","a1","a2")
  # deltas=as.data.frame(deltas)
  # colnames(deltas)=c("d1","d2")

  ##return data, censoring rate, average number of gap time pairs
  # cenrate=sum(table(dat$id)==1)/nsize
  # mbar=mean(table(dat$id))
  # maxm=max(table(dat$id))
  return(data=dat)
}
