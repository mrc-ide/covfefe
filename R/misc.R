
#------------------------------------------------
#' @title Import file
#'
#' @description Import file from the inst/extdata folder of this package
#'
#' @details TODO
#' 
#' @param name name of file
#'
#' @export
#' @examples
#' # TODO

covfefe_file <- function(name) {

  # load file from inst/extdata folder
  name_full <- system.file("extdata/", name, package='covfefe', mustWork = TRUE)
  ret <- readRDS(name_full)

  # return
  return(ret)
}

# -----------------------------------
# convert matrix to list format for use within Rcpp code
# (not exported)
#' @noRd
mat_to_rcpp <- function(x) {
  return(split(x, f = 1:nrow(x)))
}

# -----------------------------------
# convert Rcpp list format to matrix
# (not exported)
#' @noRd
rcpp_to_mat <- function(x) {
  return(matrix(unlist(x), length(x), byrow = TRUE))
}

#------------------------------------------------
# convert 3-dimensional array to list format for use within Rcpp code
# (not exported)
#' @noRd
array_to_rcpp <- function(a) {
  s <- split(a, arrayInd(seq_along(a), dim(a))[,1])
  lapply(s, function(x){split(x, rep(1:dim(a)[2], each=dim(a)[3]))})
}

#------------------------------------------------
# replace NULL value with default
# (not exported)
#' @noRd
define_default <- function(x, default_value) {
  if (is.null(x)) {
    x <- default_value
  }
  return(x)
}

#------------------------------------------------
# force scalar to vector by repeating value 
# (not exported)
#' @noRd
force_vector <- function(x, l) {
  if (length(x)==1) {
    x <- rep(x, l)
  }
  return(x)
}

# -----------------------------------
# ask user a yes/no question. Return TRUE/FALSE
# (not exported)
#' @noRd
user_yes_no <- function(x = "continue? (Y/N): ") {
  userChoice <- NA
  while (!userChoice %in% c("Y", "y" ,"N", "n")) {
    userChoice <- readline(x)
  }
  return(userChoice %in% c("Y", "y"))
}

# -----------------------------------
# make a series of draws from a defined distribution
# (not exported)
#' @noRd
R_draws <- function(prob, n) {
  sample(0:(length(prob)-1), n, replace = TRUE, prob = prob)
}
