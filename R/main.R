
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
#' Bind infected durations to project
#'
#' Bind infected durations to project
#'
#' @param proj the current project
#' @param infected_durations list of infection lengths
#'
#' @export
#' @examples
#' # TODO

bind_infected_durations <- function(proj, infected_durations) {

  # check inputs
  assert_covfefe_project(proj)

  # bind infected durations to project
  proj$infected_durations <- infected_durations

  # create infecteds matrix large enough to store longest infection in any deme
  demes <- length(infected_durations)
  max_time <- max(mapply(function(y){
    0:(length(y)-1) + mapply(function(x){max(0,x)}, y)
    }, infected_durations))
  infecteds <- matrix(0, demes, max_time)

  # calculate total number of infecteds at each time point in each deme
  for (k in 1:demes) {
    n <- length(infected_durations[[k]])
    start_times <- rep(0:(n-1), times = mapply(length, infected_durations[[k]]))
    end_times <- start_times + unlist(infected_durations[[k]])
    infecteds[k,] <- mapply(function(x){sum(start_times<=x & end_times>x)}, 0:(max_time-1))
  }

  # binf total infecteds to project
  proj$infecteds <- infecteds

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

  # check inputs
  assert_covfefe_project(proj)

  message("foobar")

}
