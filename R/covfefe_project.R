
#------------------------------------------------
#' Create covfefe project
#'
#' Create covfefe project
#'
#' @param recom_rate TODO
#' @param bloodstage_skew TODO
#' @param biting TODO
#' @param oocysts TODO
#' @param hepatocytes TODO
#'
#' @export
#' @examples
#' # TODO

covfefe_project <- function(recom_rate = 1e-2, bloodstage_skew = 1, biting = dpois(0:9, lambda = 1), oocysts = dpois(0:9, lambda = 1), hepatocytes = dpois(0:9, lambda = 1)) {

  # create parameters list
  parameters <- list(recom_rate = recom_rate,
                     bloodstage_skew = bloodstage_skew)

  # create distributions list
  distributions <- list(biting = biting,
                        oocysts = oocysts,
                        hepatocytes = hepatocytes)

  # create new project
  ret <- list(parameters = parameters,
              distributions = distributions,
              infected_durations = NULL,
              infecteds = NULL,
              genotypes = NULL)

  # create class
  class(ret) <- "covfefe_project"

  # return
  return(ret)
}

#------------------------------------------------
# overload print() function for covfefe_project
# (not exported)

print.covfefe_project <- function(x, ...) {

  # print selected elements
  print(unclass(x)[c("parameters", "distributions", "infected_durations", "infecteds", "genotypes")])

  invisible(x)
}

#------------------------------------------------
# overload summary() function for covfefe_project
# (not exported)

summary.covfefe_project <- function(x, ...) {

  print("TODO")

}

#------------------------------------------------
# determine if object is of class covfefe_project
# (not exported)

is.covfefe_project <- function(x) {
  inherits(x, "covfefe_project")
}
