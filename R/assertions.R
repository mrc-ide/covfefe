
#------------------------------------------------
# is positive (with or without zero allowed)
assert_pos <- function(x, zero_allowed = TRUE) {
  assert_that(is.numeric(x))
  if (zero_allowed) {
    assert_that(all(x>=0))
  } else {
    assert_that(all(x>0))
  }
}

#------------------------------------------------
# is integer
assert_int <- function(x) {
  assert_that(is.numeric(x))
  assert_that(all.equal(x, as.integer(x)))
}

#------------------------------------------------
# is positive integer (with or without zero allowed)
assert_pos_int <- function(x, zero_allowed = TRUE) {
  assert_int(x)
  assert_pos(x, zero_allowed)
}

#------------------------------------------------
# is between bounds (inclusive or exclusive)
assert_bounded <- function(x, left = 0, right = 1, inclusive_left = TRUE, inclusive_right = TRUE) {
  assert_that(is.numeric(x))
  if (inclusive_left) {
    assert_that(all(x>=left))
  } else {
    assert_that(all(x>left))
  }
  if (inclusive_right) {
    assert_that(all(x<=right))
  } else {
    assert_that(all(x<right))
  }
}

#------------------------------------------------
# is class covfefe_project
assert_covfefe_project <- function(x) {
  is.covfefe_project(x)
}
