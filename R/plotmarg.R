#' Plot Marginal Survival after non-parametric analysis
#'
#' @description
#' This function plots marginal survival for recurrent event data. For example see \code{bivrecNP}.
#'
#' @param x must be an object of \code{bivrecNP} class.
#' @param CI passed from object
#' @noRd
#' @return A plot of marginal survival vs. first gap time with confidence interval.


plotMarg <- function(x, CI) {
  if (!inherits(x, "bivrecNP")) stop("Object must be a bivrecNP class")
  xij <- x$
  forplot <- x$marginal_survival[1:3]
  #formula <- bivrec.nonparam.result$formula

  #variables <- all.vars(formula)
  xij <- x$data4Creg$xij
  mx <- round(max(xij), digits = 0)
  str_mx <- substring(as.character(mx), 1, nchar(as.character(mx))-1)
  str_mx <- paste(as.numeric(str_mx)+1, 0, sep="")
  mx <- round(as.numeric(str_mx), digits=1)
  forplot <- rbind(c(0, 1, 0), forplot, c(mx, 0, forplot[nrow(forplot),3]))

  ##### Wald CI and plot
  conf.lev = 1 - ((1-x$CI)/2)
  forplot$lower <- forplot[,2] - qnorm(conf.lev)*forplot[,3]
  forplot$upper <- forplot[,2] + qnorm(conf.lev)*forplot[,3]
  index <- which(forplot$lower<0)
  forplot[index, -1] <- forplot[index[1]-1, -1]

  plot(forplot$Time, forplot$Marginal_Survival, type = "l", xlab = "Type I Gap Times (x)",
       ylab = "Marginal Survival", yaxp  = c(0, 1, 10),
       xaxp  = round(c(0, mx, 10), digits=1), main = expression(1 - P(X^0 <= x))
       )
  graphics::lines(forplot$Time, forplot$lower, lty = 2)
  graphics::lines(forplot$Time, forplot$upper, lty = 2)

}
