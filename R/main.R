
#------------------------------------------------
# The following commands ensure that package dependencies are listed in the
# NAMESPACE file.

#' @useDynLib covfefe
#' @import assertthat
#' @importFrom Rcpp sourceCpp
#' @import stats
#' @import graphics
NULL

#------------------------------------------------
#' Simplify line list
#'
#' Simplify line list
#'
#' @param proj the current project
#'
#' @export
#' @examples
#' # TODO

simplify_line_list <- function(line_list, samp_mat, demes) {
  
  # split samp_mat into its elements
  samp_times <- unique(samp_mat[,1])
  samp_demes <- split(samp_mat[,2], samp_mat[,1])
  samp_num <- split(samp_mat[,3], samp_mat[,1])
  
  args <- list(samp_times = samp_times,
               samp_demes = samp_demes,
               samp_num = samp_num,
               demes = demes)
  
  t0 <- Sys.time()
  
  output_raw <- simplify_line_list_cpp(line_list, args)
  
  message(sprintf("completed in %s seconds", round(Sys.time() - t0, 2)))
  
  return(output_raw)
}

#------------------------------------------------
#' Simulate genotypes
#'
#' Simulate genotypes
#'
#' @param proj the current project
#'
#' @export
#' @examples
#' # TODO

sim_genotypes <- function(proj) {

  # check arguments
  assert_covfefe_project(proj)

  # get useful parameters
  demes <- length(proj$durations)
  max_time <- length(proj$durations[[1]])
  
  # define argument list
  args <- list(parameters = proj$parameters,
               durations = proj$durations,
               migrations = array_to_rcpp(proj$migrations),
               demes = demes,
               max_time = max_time
               )
  
  # run efficient C++ function
  output_raw <- sim_genotypes_cpp(args)
  
  return(output_raw)
}
