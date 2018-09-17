#' Plot Marginal Survival after non-parametric analysis
#'
#' @description
#' This function plots marginal survival for recurrent event data.  Called from biv.rec.np(). No user interface.
#'
#' @param bivrec.nonparam.result List with marginal.survival, formula and data. Passed from biv.rec.np()
#' @param CI Confidence level for CI. Passed from biv.rec.np()
#'
#' @return A plot of marginal survival vs. first gap time with confidence interval.
#'
#' @keywords internal
#'
marg.surv.plot <- function(bivrec.nonparam.result, CI) {

  forplot <- bivrec.nonparam.result$marginal.survival[1:3]
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
  conf.lev = 1 - ((1-CI)/2)
  forplot$lower <- forplot[,2] - qnorm(conf.lev)*forplot[,3]
  forplot$upper <- forplot[,2] + qnorm(conf.lev)*forplot[,3]
  index <- which(forplot$lower<0)
  forplot[index, -1] <- forplot[index[1]-1, -1]
  plot(forplot$Time, forplot$Marginal.Survival, type = "l", xlab = "Time", ylab = "Marginal Survival",
       yaxp  = c(0, 1, 10), xaxp  = c(0, mx, 15), main = expression(P(X^0 <= x)))
  graphics::lines(forplot$Time, forplot$lower, lty = 2)
  graphics::lines(forplot$Time, forplot$upper, lty = 2)

}
