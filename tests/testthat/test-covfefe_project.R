context("covfefe_project")

#------------------------------------------------
# test assert_numeric
test_that("covfefe_project basic properties work", {
  p0 <- covfefe_project()
  expect_true(assert_covfefe_project(p0))
  
  capture.output(p1 <- print(p0))
  expect_equal(p1, p0)
  
  capture.output(p1 <- summary(p0))
  expect_equal(p1, p0)
  
  expect_true(is.covfefe_project(p0))
})
