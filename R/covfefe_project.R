
#------------------------------------------------
#' @title Create covfefe project
#'
#' @description Create covfefe project
#'
#' @details TODO
#'
#' @export
#' @examples
#' # TODO

covfefe_project <- function() {

  # TODO - do we need parameters?
  parameters <- list(foo = -9)

  # create new project
  ret <- list(parameters = parameters,
              infection_history = NULL,
              pop_counts = NULL,
              samp_mat = NULL,
              samp_hosts = NULL,
              pruned_tree = NULL,
              distributions = NULL,
              genotypes = NULL)

  # create class
  class(ret) <- "covfefe_project"

  # return invisibly
  invisible(ret)
}

#------------------------------------------------
# overload print() function for covfefe_project
# (not exported)

print.covfefe_project <- function(x, ...) {

  # print without class name
  print(unclass(x))
  
  # return invisibly
  invisible(x)
}

#------------------------------------------------
# overload summary() function for covfefe_project
# (not exported)

summary.covfefe_project <- function(x, ...) {

  print("TODO")

  # return invisibly
  invisible(x)
}

#------------------------------------------------
# determine if object is of class covfefe_project
# (not exported)

is.covfefe_project <- function(x) {
  inherits(x, "covfefe_project")
}
