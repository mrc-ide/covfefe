
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

covfefe_project2 <- function() {

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

