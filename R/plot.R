
#------------------------------------------------
#' Plot total infecteds
#'
#' Plot total infecteds
#'
#' @param proj project
#'
#' @export
#' @examples
#' # TODO

plot_infecteds <- function(proj) {

  # check inputs
  assert_covfefe_project(proj)

  # produce plot
  barplot(proj$infecteds, space=0, border=NA)
}
