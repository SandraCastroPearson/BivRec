#' Plot Marginal Survival after non-parametric analysis
#'
#' @description
#' This function plots marginal survival for recurrent event data.
#'
#' @param bivrec.nonparam.result is a list obtained from running biv.rec.fit with method="Non-Parametric"
#'
#' @importFrom graphics lines
#' @return A plot of marginal survival vs. first time with 95\% confidence interval.
#'
#' @examples
#' \dontrun{
#' library(BivRec)
#' set.seed(1234)
#' sim.data <- data.sim(nsize=300, beta1=c(0.5,0.5), beta2=c(0,-0.5), cr=63, sg2=0.5, set=1.1)
#' nonpar.example <- biv.rec.fit(id + xij + yij + epi + d2 + d1 ~ 1,
#'           data=sim.data, method="Non-Parametric", ai=1)
#' nonpar.example
#' plot.marg.surv(nonpar.example)
#' }
#'
#' @keywords plot marginal survival
#' @export

marg.surv.plot <- function(bivrec.nonparam.result) {

  forplot <- bivrec.nonparam.result$marginal.survival
  formula <- bivrec.nonparam.result$formula
  data <- bivrec.nonparam.result$data
  variables <- all.vars(formula)
  xij <- eval(parse(text =paste("data$", variables[2], sep="")))
  mx <- round(max(xij), digits = 0)
  str_mx <- substring(as.character(mx), 1, nchar(as.character(mx))-1)
  str_mx <- paste(as.numeric(str_mx)+1, 0, sep="")
  mx <- as.numeric(str_mx)
  forplot <- rbind(c(0, 1, 0), forplot, c(mx, 0, forplot[nrow(forplot),3]))

  #####95% Wald CI and plot
  forplot$lower <- forplot[,2] - 1.96*forplot[,3]
  forplot$upper <- forplot[,2] + 1.96*forplot[,3]
  index <- which(forplot$lower<0)
  forplot[index, -1] <- forplot[index[1]-1, -1]
  plot(forplot$Time, forplot$`Marginal Survival`, type = "l", xlab = "Time", ylab = "Marginal Survival",
       yaxp  = c(0, 1, 10), xaxp  = c(0, mx, 15))
  lines(forplot$Time, forplot$lower, lty = 2)
  lines(forplot$Time, forplot$upper, lty = 2)

}
