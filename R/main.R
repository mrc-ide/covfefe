
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
#' Simulate genotypes
#'
#' Prune infection tree and simulate genotypes
#'
#' @param line_list TODO
#' @param samp_mat TODO
#' @param demes TODO
#' @param loci TODO
#' @param recom_rate TODO
#'
#' @export
#' @examples
#' # TODO

sim_genotypes <- function(line_list, samp_mat, demes, loci, recom_rate) {
  
  # split samp_mat into its elements
  samp_times <- unique(samp_mat[,1])
  samp_demes <- split(samp_mat[,2], samp_mat[,1])
  samp_num <- split(samp_mat[,3], samp_mat[,1])
  
  args <- list(samp_times = samp_times,
               samp_demes = samp_demes,
               samp_num = samp_num,
               demes = demes,
               loci = loci,
               recom_rate = recom_rate)
  
  t0 <- Sys.time()
  
  output_raw <- sim_genotypes_cpp(line_list, args)
  
  message(sprintf("completed in %s seconds", round(Sys.time() - t0, 2)))
  
  return(output_raw)
}
