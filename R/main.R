
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
#' Bind mosquito lag durations to project
#'
#' Bind mosquito lag durations to project. Values of -1 indicate that this infected host should have genotypes generated de novo, rather than from another host in this simulation.
#'
#' @param proj the current project
#' @param lags lag durations
#'
#' @export
#' @examples
#' # TODO

bind_lags <- function(proj, lags) {
  
  # check inputs
  assert_covfefe_project(proj)
  
  # bind mosquito lag durations to project
  proj$lags <- lags
  
  # return project
  return(proj)
}

#------------------------------------------------
#' Bind infected durations to project
#'
#' Bind infected durations to project
#'
#' @param proj the current project
#' @param durations list of infection lengths
#'
#' @export
#' @examples
#' # TODO

bind_durations <- function(proj, durations) {

  # check inputs
  assert_covfefe_project(proj)

  # bind durations to project
  proj$durations <- durations

  # create infecteds matrix large enough to store longest infection in any deme
  #demes <- length(infected_durations)
  #max_time <- max(mapply(function(y){
  #  0:(length(y)-1) + mapply(function(x){max(0,x)}, y)
  #  }, infected_durations))
  #infecteds <- matrix(0, demes, max_time)

  # calculate total number of infecteds at each time point in each deme
  #for (k in 1:demes) {
  #  n <- length(infected_durations[[k]])
  #  start_times <- rep(0:(n-1), times = mapply(length, infected_durations[[k]]))
  #  end_times <- start_times + unlist(infected_durations[[k]])
  #  infecteds[k,] <- mapply(function(x){sum(start_times<=x & end_times>x)}, 0:(max_time-1))
  #}

  # binf total infecteds to project
  #proj$infecteds <- infecteds

  # return project
  return(proj)
}

#------------------------------------------------
#' Bind migration array to project
#'
#' Bind migration array to project
#'
#' @param proj the current project
#' @param migrations migration array
#'
#' @export
#' @examples
#' # TODO

bind_migrations <- function(proj, migrations) {

  # check inputs
  assert_covfefe_project(proj)

  # bind infected durations to project
  proj$migrations <- migrations

  # return project
  return(proj)
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
