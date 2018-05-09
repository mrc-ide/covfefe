context("assertions")

#------------------------------------------------
# test assert_covfefe_project

test_that("assert_covfefe_project works", {
  p <- covfefe_project()
  expect_true(assert_covfefe_project(p))
})
