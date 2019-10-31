context("Sample size 1")

library(BivRec)
n = round(runif(1, 80, 100), digits = 0)
bivrec_data <- simBivRec(nsize=n, beta1=c(0.5,0.5), beta2=c(0,-0.5),
                tau_c=63, set=1.1)

test_that("lee check 1", {
  lee_reg <- bivrecReg(bivrecSurv(id, epi, xij, yij, d1, d2) ~ a1,
                       bivrec_data2, "Lee.et.al")
  expect_is(lee_reg, "bivrecReg")
  lee_coeffs <- coef.bivrecReg(lee_reg)
  expect_is(lee_coeffs, "matrix")
  expect_is(lee_coeffs[,1], "numeric")
})


