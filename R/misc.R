
#------------------------------------------------
#' Import infected durations
#'
#' Import infected durations from the inst/extdata folder of this package
#'
#' @param name name of file
#'
#' @export
#' @examples
#' # TODO

covfefe_infected_durations <- function(name) {

  # load file from inst/extdata folder
  name_full <- system.file("extdata/", name, package='covfefe', mustWork = TRUE)
  ret <- readRDS(name_full)

  # return
  return(ret)
}
