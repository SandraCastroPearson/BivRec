
#####bivrecNP function to obtain results for marginal survival, conditional cdf and joint cdf

#x is the bivrecSurv object

##Example
#npresult <- bivrecNP(bdat,ai=1, u1 = c(2, 5, 10, 20), u2 = c(1, 5, 10, 15),conditional = FALSE, given.interval=c(0, 10))
#npresult2 <- bivrecNP(bdat,ai=1, u1 = c(2, 5, 10, 20), u2 = c(1, 5, 10, 15),conditional = TRUE, given.interval=c(0, 10))

#' nonpar.result <- biv.rec.np(formula = id + epi + xij + yij + d1 + d2 ~ 1,
#'           data=dat, ai=1, u1 = c(2, 5, 10, 20), u2 = c(1, 5, 10, 15),
#'           conditional = FALSE, given.interval=c(0, 10), jointplot=FALSE,
#'           marginalplot = FALSE, condiplot = FALSE)
#'
#'           nonpar.result2 <- biv.rec.np(formula = id + epi + xij + yij + d1 + d2 ~ 1,
#'           data=dat, ai=1, u1 = c(2, 5, 10, 20), u2 = c(1, 5, 10, 15),
#'           conditional = TRUE, given.interval=c(0, 10), jointplot=FALSE,
#'           marginalplot = FALSE, condiplot = FALSE)
#'
bivrecNP <- function(x, CI, ai, u1, u2, conditional, given.interval){
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

  ### Dataframe with id, epi, xij, yij, d1, d2, zij and ci. This part will be moved to bivrecSurv
  data <- x$df
  covariates <- rep(1, length(data$id))

  iden <- data$id #move this to BivrecSurv object
  iden.u <- unique(iden)
  new.id <- NULL
  if (class(iden)!="num") {
    if (class(iden)!="int") {
      for (i in 1:length(iden.u)){
        for (j in 1:length(iden)) {
          if (iden[j] == iden.u[i]){
            new.id=c(new.id,i)
          }
        }
      }
      data$new.id <- new.id
    }
  }
  data <- data[,-which(colnames(data)=="id")]
  colnames(data)[ncol(data)] = "id" #this just moved the id column to the end to match up with new.id values

  if (missing(u1)) {u1 <- round(seq(quantile(data$xij, probs = 0.4), max(data$xij), length.out=5))}
  if (missing(u2)) {u2 <- round(seq(quantile(data$yij, probs = 0.4), max(data$yij), length.out=4))}
  temp <- rep(u1, each = length(u2))
  temp2 <- rep(u2, length(u1))
  u <- cbind(u1=temp, u2=temp2)

  print("Estimating joint cdf and marginal survival")
  if (ai==1) {
  cdf_res1 <- nonparam.cdf(x$dat4np1$forcdf, u, ai, CI) #result for joint cdf if ai=1
  marg_res1 <- nonparam.marginal(x$dat4np1$formarg, CI) #result for marginal if ai=1
  if (conditional == FALSE) {
    final.result <- list(joint.cdf = cdf_res1, marginal.survival = marg_res1, ai=ai)
  } else {
    if (missing(given.interval)) {
      print("Error for conditional calculation given.interval argument missing.")
      final.result <- list(cdf = cdf_res1, marginal.survival = marg_res1, ai=ai)
    } else {
      partial.result <- list(cdf = cdf_res1, marginal.survival = marg_res1, data = data, ai=ai, new_data=x$dat4np1) #took out formula as a param
      ccdf_res1 <- nonparam.conditional(partial.result, given.interval, CI,x$df$yij) #took out condiplot as a param
      final.result <- list(joint.cdf = cdf_res1, marginal.survival = marg_res1, conditional.cdf = ccdf_res1,ai=ai)
    }
  }
  }
  if (ai==2) { #ai=2
  cdf_res2 <- nonparam.cdf(x$dat4np2$forcdf, u, ai, CI) #result for joint cdf if ai=2
  marg_res2 <- nonparam.marginal(x$dat4np2$formarg, CI) #result for marginal if ai=2
  if (conditional == FALSE) {
    final.result <- list(joint.cdf = cdf_res2, marginal.survival = marg_res2, ai=ai)
  } else {
    if (missing(given.interval)) {
      print("Error for conditional calculation given.interval argument missing.")
      final.result <- list(joint.cdf = cdf_res2, marginal.survival = marg_res2, ai=ai)
    } else {
      partial.result <- list(cdf = cdf_res1, marginal.survival = marg_res2, data = data, ai=ai, new_data=x$dat4np2) #took out formula as a param
      ccdf_res2 <- nonparam.conditional(partial.result, given.interval, CI,x$df$yij) #took out condiplot as a param
      final.result <- list(joint.cdf = cdf_res2, marginal.survival = marg_res2, conditional.cdf = ccdf_res2,ai=ai)
    }
  }
  }
  class(final.result)<-"bivrecNP"
  final.result$CI <- CI
  final.result$given.interval<-given.interval
  final.result$conditional <- conditional #boolean indicator
  #final.result$ygrid <-ccdf_res2$ygrid #original response data from bivrecSurv object
  return(final.result) #Essentially the bivrecNP object provides the data (the new_id stuff), CI, results for all 3 (if conditional=true),
  #the conditional indicator
}

is.bivrecNP <- function(x) inherits(x, "bivrecNP")
