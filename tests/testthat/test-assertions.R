context("assertions")

#------------------------------------------------
# test assert_pos
test_that("assert_pos works", {
  expect_true(assert_pos(seq(1, 5, 0.5)))
  expect_true(assert_pos(seq(0, 5, 0.5), zero_allowed = TRUE))
  
  expect_error(assert_pos(seq(-1, -5, -0.5)))
  expect_error(assert_pos(seq(0, 5, 0.5), zero_allowed = FALSE))
  expect_error(assert_pos(seq(-5, 5, 0.5), zero_allowed = TRUE))
  expect_error(assert_pos("foo"))
})

#------------------------------------------------
# test assert_int
test_that("assert_int works", {
  expect_true(assert_int(-5:5))
  
  expect_error(assert_int(0.5))
  expect_error(assert_int("foo"))
})

#------------------------------------------------
# test assert_pos_int
test_that("assert_pos_int works", {
  expect_true(assert_pos_int(1:5))
  expect_true(assert_pos_int(0:5, zero_allowed = TRUE))
  
  expect_error(assert_pos_int(-1:-5))
  expect_error(assert_pos_int(0:5, zero_allowed = FALSE))
  expect_error(assert_pos_int(-5:5, zero_allowed = TRUE))
  expect_error(assert_pos_int("foo"))
})

#------------------------------------------------
# test assert_bounded
test_that("assert_bounded works", {
  expect_true(assert_bounded(seq(0, 1, 0.1)))
  expect_true(assert_bounded(seq(-5, 5, 0.1), left = -5, right = 5))
  
  expect_error(assert_bounded(-5, left = -5, inclusive_left = FALSE))
  expect_error(assert_bounded(5, right = 5, inclusive_right = FALSE))
  expect_error(assert_bounded("foo"))
})

#------------------------------------------------
# test assert_covfefe_project
test_that("assert_covfefe_project works", {
  p <- covfefe_project()
  expect_true(assert_covfefe_project(p))
})