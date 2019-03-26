#bivrecNP function to obtain results for marginal survival, conditional cdf and joint cdf

#x is the bivrecSurv object 
#id + epi + xij + yij + d1 + d2
#####For non-parametric analysis
# dat : a data.frame including
#     1)id numbers, 2)orders of episodes, 3)first gap time, 4)second gap time
#     5)censoring times, 6) censoring indicators in each column
# ai: a non-negative function of censoring time
biv.rec.np <- function(x, CI, ai, u1, u2, conditional, given.interval){
  if (!is.bivrecSurv(x)) stop("Response must be a bivrecSurv class")
  if (missing(ai)) {ai<-1}
  if (missing(conditional)) {conditional <- FALSE}
  if (missing(CI)) {CI <- 0.95}
  if (CI > 0.99) {
    print("Error CI > 0.99")
    stop()} else {
      if (CI<0.5) {
        print("Error CI < 0.5")
        stop()
      }
    }
  
  ### Dataframe with id, epi, xij, yij, d1, d2, zij and ci 
  data <- x$df
  covariates <- rep(1, length(data$id))

  if (missing(u1)) {u1 <- round(seq(quantile(data$xij, probs = 0.4), max(data$xij), length.out=5))}
  if (missing(u2)) {u2 <- round(seq(quantile(data$yij, probs = 0.4), max(data$yij), length.out=4))}
  temp <- rep(u1, each = length(u2))
  temp2 <- rep(u2, length(u1))
  u <- cbind(u1=temp, u2=temp2)
  
  print("Estimating joint cdf and marginal survival")
  res1 <- nonparam.cdf(x$np$forcdf, u, ai, CI)
  res2 <- nonparam.marginal(x$np$formarg, CI)
  
  if (conditional == FALSE) {
    final.result <- list(joint.cdf = res1, marginal.survival = res2, formula=formula, ai=ai)
  } else {
    if (missing(given.interval)) {
      print("Error for conditional calculation given.interval argument missing.")
      final.result <- list(cdf = res1, marginal.survival = res2, formula=formula, data = data, ai=ai)
    } else {
      partial.result <- list(cdf = res1, marginal.survival = res2, formula=formula, data = data, ai=ai, new_data=new_data)
      res3 <- nonparam.conditional(partial.result, given.interval, CI) #took out condiplot as a param
      final.result <- list(joint.cdf = res1, marginal.survival = res2, conditional.cdf = res3, formula=formula, ai=ai)
    }
  }
  
  return(final.result)
}

