
#------------------------------------------------
#' @title Define empty covfefe project
#'
#' @description Define empty covfefe project
#'
#' @export

covfefe_project <- function() {
  
  # create null project
  project <- list(sim_parameters = NULL,
                  sim_output = NULL)
  class(project) <- "covfefe_project"
  
  # define default epi parameters
  project <- define_epi_parameters(project)
  
  # define default deme parameters
  project <- define_deme_parameters(project)
  
  # define default demography
  project <- define_demography(project)
  
  # define default migration parameters
  project <- define_migration(project)
  
  # return
  invisible(project)
}

#------------------------------------------------
# overload print() function for covfefe_project
#' @noRd
print.covfefe_project <- function(x, ...) {
  
  # print summary
  summary(x)
  
  # return invisibly
  invisible(x)
}

#------------------------------------------------
# overload summary() function for covfefe_project
#' @noRd
summary.covfefe_project <- function(x, ...) {
  
  message("TODO")
  
}
