library(testthat)

if (!testthat:::on_cran()) {
  library(QuantRegGLasso)
  test_check("QuantRegGLasso")
}
