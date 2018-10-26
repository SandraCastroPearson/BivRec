#' Bivariate Recurrent Response and Covariate Data Simulation
#'
#' @description This function simulates a series of alternating recurrent events based on simulations in Lee CH, Huang C-Y, Xu G, Luo X (2017).
#'
#' @param nsize sample size which refers to the number of subjects in the data set where each subject could have multiple episodes of events.
#' @param beta1 true coefficients for first gap time in the accelerated failure time model (AFT).
#' @param beta2 true coefficients for second gap time in the accelerated failure time model (AFT).
#' @param tau_c maximum support of censoring time. Can take values as follows:
#' \itemize{
#' \item tau_c=63: corresponds to cr=15\% and corresponding \ifelse{html}{\out{m_bar}}{\eqn{\bar{m}}} for each scenario in tables 1 and 2 of Lee CH, Huang C-Y, Xu G, Luo X (2017).
#' \item tau_c=30: corresponds to cr=30\% and corresponding \ifelse{html}{\out{m_bar}}{\eqn{\bar{m}}} for each scenario in tables 1 and 2 of Lee CH, Huang C-Y, Xu G, Luo X (2017).
#' }
#'
#' @param set Simulation setting based on scenarios outlined in tables 1 and 2 in Lee CH, Huang C-Y, Xu G, Luo X (2017). Choose 1.1 (default) for scenario 1 with \eqn{\rho=1} in the covariance matrix of the frailty vector, 1.2 for scenario 1 with \eqn{\rho=0.5}, 1.3 for scenario 1 with \eqn{\rho=0} and 2.0 for scenario 2.
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
#' @keywords biv.rec.sim, simulation
#'
#' @examples
#' biv.rec.sim(nsize=150, beta1=c(0.5,0.5), beta2=c(0,-0.5), tau_c=63, set=1.1)
#' @references
#' \enumerate{
#' \item Lee C, Huang CY, Xu G, Luo X (2017). Semiparametric regression analysis for alternating recurrent event data. Statistics in Medicine, 37: 996-1008.
#' \url{https://doi.org/10.1002/sim.7563}
#' }
#' @export

##-----data generation
biv.rec.sim <- function(nsize,beta1,beta2,tau_c,set) {

  if (missing(tau_c)) {tau_c <- 63}
  if (missing(set)) {set <- 1.1}
  sg2 <- 0.5

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
  if (set==1.3) {
    gamma1=rnorm(nsize,1,sqrt(sg2))
    gamma2=rnorm(nsize,1,sqrt(sg2))
  }
  if (set==2.0) {
    gamma1=rnorm(nsize,1,sqrt(sg2))
    gamma2=rgamma(nsize,shape=1/sg2,rate=1/sg2)
  }

  dat=NULL
  deltas=NULL
  for (i in id) {
    x.tmp=exp(gamma1[i]+A[i,]%*%beta1)
    y.tmp=exp(gamma2[i]+A[i,]%*%beta2)

    ci=runif(1,0,tau_c)

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
    subdat=cbind(i,1:epi,dat.tmp[,1],dat.tmp[,2],ci, cen,matrix(A[i,],epi,2,byrow=TRUE))
    dat=rbind(dat,subdat)
    deltas=rbind(deltas, cen)
  }
  dat=as.data.frame(dat)
  colnames(dat)=c("id","epi","xij","yij","ci","d1", "d2","a1","a2")

  ##return data, censoring rate, average number of gap time pairs
  # cenrate=sum(table(dat$id)==1)/nsize
  # mbar=mean(table(dat$id))
  # maxm=max(table(dat$id))
  return(data=dat)
}
