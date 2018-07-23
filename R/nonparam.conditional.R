np.fit4conditional <- function(formula, data, ai, u1, u2){

  ### PULL INFORMATION FROM PARAMETERS TO SEND TO REFORMAT
  identifier=xij=yij=c_indicatorY=c_indicatorX=episode=covariates=NULL
  method <- "Non-Parametric"
  condgx <- TRUE

  ###Send to biv.rec.reformat and complete analysis
  new_data <- biv.rec.reformat(identifier, xij, yij, c_indicatorY, c_indicatorX, episode, covariates, method, ai, condgx, data)
  temp <- rep(u1, each = length(u2))
  temp2 <- rep(u2, length(u1))
  u <- cbind(u1=temp, u2=temp2)
  res1 <- nonparam.cdf(fit_data=new_data$forcdf, u, ai)

  return(res1)
}

bstp <- function(seedi, ps1, ps2, x.grid, y.grid, n, refdata, ai, mintime) {
  set.seed(seedi)
  samp.id <- sample(1:n, n, replace = TRUE)
  boot.dat <- NULL
  for (j in 1:n){
    temp <- refdata[which(refdata$id==samp.id[j]), ]
    boot.dat <- rbind(boot.dat, cbind(id = j, temp[, -1]))
  }

  joint2 <- np.fit4conditional(formula = id + vij + wij + epi + d2 ~ 1, data=boot.dat,
                             ai=ai, u1=x.grid$Time[2], u2=y.grid)

  if (x.grid$Time[1] == mintime) {
    conditional <- joint2$`Prob(x_ij < u1,  y_ij < u2)` / (1-ps2)
  } else {
    print("esp")
    joint1 <- np.fit4conditional(formula = id + vij + wij + epi + d2 ~ 1, data=boot.dat,
                          ai=ai, u1=x.grid$Time[1], u2=y.grid)
    conditional <- (joint2$`Prob(x_ij < u1,  y_ij < u2)` - joint1$`Prob(x_ij < u1,  y_ij < u2)`)/
      (ps1 - ps2)
  }

  return(conditional)
}

###################################################################
#################### FUNCTION NOT FOR USER ########################
###################################################################
#' A Function for additional non-parametric analysis of bivariate recurrent event data.
#'
#' @description
#' This function calculates the conditional cdf after estimation of joint cdf and marginal survival.
#' @param bivrec.nonparam.result is a list obtained from running biv.rec.fit with method="Non-Parametric",
#' @param given.interval is a 1x2 vector indicating an interval for the first gap time to estimate the cdf of second gap time given this interval. If interval given is c(v1, v2) the  function calculates P(yij <= w | v1 <= xij <= v2).
#' The given values v1 and v2 must be in the range of gap times in the estimated marginal survival (valid times are given in the Time column of the results for marginal survival).
#'
#' @return A data frame with the conditional CDF for the given an interval of the first gap time and corresponding plot.
#'
#' @keywords internal
#'

nonparam.conditional <- function(bivrec.nonparam.result, given.interval) {

  ####Extract items from results
  marginal <- bivrec.nonparam.result$marginal.survival
  marginal$rounded <- round(marginal$`Marginal Survival`, digits=2)
  formula <- bivrec.nonparam.result$formula
  data <- bivrec.nonparam.result$data
  ai <- bivrec.nonparam.result$ai
  new_data <- bivrec.nonparam.result$new_data
  margdata <- new_data$formarg
  refdata <- new_data$refdata
  n <- margdata$n
  variables <- all.vars(formula)
  yij <- eval(parse(text =paste("data$", variables[3], sep="")))

  x.grid <- marginal[which(marginal$Time<=given.interval[2]), ]
  x.grid <- x.grid[which(x.grid$Time>=given.interval[1]), ]
  x.grid$diffs <- x.grid$rounded - x.grid$rounded[nrow(x.grid)]
  x.grid <- x.grid[c(1, max(which(x.grid$diffs >= 0.10)), nrow(x.grid)),]
  y.grid <- seq(min(yij), max(yij), length.out = 200)
  ps1 <- x.grid$`Marginal Survival`[1]
  ps2 <- x.grid$`Marginal Survival`[3]

  B = 100
  cond.prob <- matrix(rep(NA, length(y.grid)*B), ncol=B)
  colnames(cond.prob) = seq(1,B,1)
  print("Estimating Conditional Probability with 95% CI using 100 Sample Bootstrap")

  for (i in 1:B) {
    print(paste("Sample", i, sep = " "))
    cond.prob[,i] <- bstp(seedi=i, ps1, ps2, x.grid, y.grid, n, refdata, ai, mintime = min(marginal$Time))
  }

  bootstrapCIs <- apply(cond.prob, 1, function(x) c(mean(x), sort(x)[0.025*B], sort(x)[0.975*B]))
  cond <- round(data.frame(y.grid, bootstrapCIs[1,], bootstrapCIs[2,], bootstrapCIs[3,]), digits = 3)
  colnames(cond) <- c("Time", "Conditional.Probability" ,"lowerCI", "upperCI")
  flat.ind <- which(cond$upperCI>1.004)
  if (length(flat.ind)!=0) {cond[flat.ind, 2:4] <- cond[min(flat.ind)-1, 2:4]}

  mainlab <- paste("Conditional Probability given first gap time between", round(given.interval[1], digits=2),
                  "and", round(given.interval[2], digits=2), sep=" ")
  plot(cond$Time, cond$upperCI, type="l", lty = 2, xlab = "Time", ylab = "Conditional Probability",
       xlim=c(0, round(max(y.grid), digits=1)), ylim=c(0, round(max(cond$upperCI), digits=1)), main=mainlab)
  lines(cond$Time, cond$lowerCI, lty = 2)
  lines(cond$Time, cond$Conditional.Probability,lty = 1)

  return(conditional=cond)

}
