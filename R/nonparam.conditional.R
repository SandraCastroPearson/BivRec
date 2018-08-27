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
    conditional <- joint2[,3] / (1-ps2)
  } else {
    joint1 <- np.fit4conditional(formula = id + vij + wij + epi + d2 ~ 1, data=boot.dat,
                          ai=ai, u1=x.grid$Time[1], u2=y.grid)
    conditional <- (joint2[,3] - joint1[,3])/(ps1 - ps2)
  }

  return(conditional)
}

###################################################################
#################### FUNCTION NOT FOR USER ########################
###################################################################
#' A Function for additional non-parametric analysis of bivariate recurrent event data.
#'
#' @description
#' This function calculates the conditional cdf after estimation of joint cdf and marginal survival.  Called from biv.rec.np(). No user interface.
#' @param bivrec.nonparam.result List with joing.cdf and marginal.survival. Passed from biv.rec.np()
#' @param given.interval is a 1x2 vector indicating an interval for the first gap time to estimate the cdf of second gap time. Passed from biv.rec.np()
#' @param CI confidence level. Passed from biv.rec.np()
#' @param condiplot a logical value. Passed from biv.rec.np()
#'
#' @return A data frame with the conditional CDF for the given an interval of the first gap time and corresponding plot.
#' @importFrom stats sd
#' @keywords internal
#'

nonparam.conditional <- function(bivrec.nonparam.result, given.interval, CI, condiplot) {

  ####Extract items from results
  marginal <- bivrec.nonparam.result$marginal.survival
  marginal$rounded <- round(marginal[,2], digits=2)
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
  if (length(which(x.grid$diffs >= 0.1))==0) {
    print("Cannot estimate conditional cdf, given.interval is too narrow")
    stop()
  }
  x.grid <- x.grid[c(1, max(which(x.grid$diffs >= 0.05)), nrow(x.grid)),]
  y.grid <- seq(min(yij), max(yij), length.out = 200)
  ps1 <- x.grid[1,2]
  ps2 <- x.grid[3,2]


  B = ifelse(CI==0.99, 200, 100)
  cond.prob <- matrix(rep(NA, length(y.grid)*B), ncol=B)
  colnames(cond.prob) = seq(1,B,1)
  print(paste("Estimating Conditional CDF with ", CI*100, "% CI using ", B, " Sample Bootstrap", sep=""))

  for (i in 1:B) {
    #print(paste("Sample", i, sep = " "))
    cond.prob[,i] <- bstp(seedi=i, ps1, ps2, x.grid, y.grid, n, refdata, ai, mintime = min(marginal$Time))
  }

  conf.lev = 1 - ((1-CI)/2)
  bootstrapCIs <- apply(cond.prob, 1, function(x) c(mean(x), sd(x), sort(x)[(1-conf.lev)*B], sort(x)[conf.lev*B]))
  cond <- round(data.frame(y.grid, bootstrapCIs[1,], bootstrapCIs[2,], bootstrapCIs[3,], bootstrapCIs[4,]), digits = 4)
  low.string <- paste("Bootstrap ", (1 - conf.lev), "%", sep="")
  up.string <- paste("Bootstrap ", conf.lev, "%", sep="")
  colnames(cond) <- c("Time", "Conditional.Probability" , " Bootstrap SE", low.string, up.string)

  flat.ind <- which(cond[,5]>=1.001)
  if (length(flat.ind)!=0) {cond[flat.ind, 2:5] <- cond[(min(flat.ind)-1), 2:5]}
  if (condiplot == TRUE) {
    plot(cond$Time, cond[,5], type="l", lty = 2, xlab = "Type II Gap Times (y)", ylab = "Conditional Probability",
         xlim=c(0, round(max(y.grid), digits=1)), ylim=c(0, round(max(cond[,5]), digits=1)),
         main=substitute(
           paste("P(", Y^0 <= y, "|", X^0 %in% "[", gi1, ",", gi2, "])"),
           list(gi1 = given.interval[1], gi2 = given.interval[2]))
         )
    lines(cond$Time, cond[,4], lty = 2)
    lines(cond$Time, cond$Conditional.Probability,lty = 1)

  }

  cond[, 4:5] <- round(cond[,4:5], digits = 2)

  return(conditional=cond)

}
