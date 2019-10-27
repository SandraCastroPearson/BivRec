context("Sample size 1")

library(BivRec)
bivrec_data <- simBivRec(nsize=1, beta1=c(0.5,0.5), beta2=c(0,-0.5),
                tau_c=63, set=1.1)

test_that("data check 1", {
  expect_is(bivrec_data, "data.frame")
  expect_equal(unique(bivrec_data$id), 1)
})

check_reg_np <- function() {
  if (nrow(bivrec_data)==1) {
    skip("np check")
  }
}

### Lee and Chang methods may work or may lead to singular systems / no convergence

test_that("np check", {
  npresult <- bivrecNP(response = with(bivrec_data, bivrecSurv(id, epi, xij, yij, d1, d2)),
                       ai=1, u1 = seq(2, 15, 1), u2 = seq(1, 10, 1), conditional = TRUE,
                       given.interval = c(0, 10), level = 0.90)
  expect_is(npresult, "bivrecNP")
  expect_is(npresult$joint_cdf, "data.frame")
  expect_is(npresult$marginal_survival, "data.frame")
  expect_is(npresult$conditional_cdf$conditional, "data.frame")
})
