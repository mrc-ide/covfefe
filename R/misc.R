
#------------------------------------------------
#' Import file
#'
#' Import file from the inst/extdata folder of this package
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

#------------------------------------------------
# convert 3-dimensional array to list of lists
# (not exported)

array_to_rcpp <- function(a) {
  s <- split(a, arrayInd(seq_along(a), dim(a))[,1])
  lapply(s, function(x){split(x, rep(1:dim(a)[2], each=dim(a)[3]))})
}