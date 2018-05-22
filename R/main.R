
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
#' @title Assign infection history list to covfefe project
#'
#' @description Assign infection history list to covfefe project
#'
#' @details TODO
#' 
#' @param proj current covfefe project
#' @param infection_history TODO
#'
#' @export
#' @examples
#' # TODO

assign_infection_history <- function(proj, infection_history) {
  
  # TODO - check inputs
  
  # TODO - more thorough checks on format of infection_history?
  
  # add to project
  proj$infection_history <- infection_history
  
  # return invisibly
  invisible(proj)
}

#------------------------------------------------
#' @title Sample hosts from infection history
#'
#' @description Sample hosts from infection history
#'
#' @details TODO
#' 
#' @param proj current covfefe project
#' @param samp_mat TODO
#' @param demes TODO
#' @param recover_pop_counts TODO
#' @param max_infections TODO
#'
#' @export
#' @examples
#' # TODO

sample_hosts <- function(proj, samp_mat, demes, recover_pop_counts = FALSE, max_infections = 5) {
  
  # TODO - check inputs
  assert_covfefe_project(proj)
  assert_that(!is.null(proj$infection_history))
  # check times increasing
  
  # get useful values
  max_time <- length(proj$infection_history)
  
  # add samp_mat to project
  proj$samp_mat <- samp_mat
  
  # split samp_mat into its elements
  samp_times <- unique(samp_mat[,1])
  samp_demes <- split(samp_mat[,2], samp_mat[,1])
  samp_num <- split(samp_mat[,3], samp_mat[,1])
  
  # define argument list
  args <- list(samp_times = samp_times,
               samp_demes = samp_demes,
               samp_num = samp_num,
               demes = demes,
               recover_pop_counts = recover_pop_counts,
               max_infections = max_infections)
  
  # start timer
  t0 <- Sys.time()
  
  # run efficient C++ function
  output_raw <- draw_hosts_cpp(proj$infection_history, args)
  
  # process population counts
  if (recover_pop_counts) {
    pop_counts <- list()
    for (k in 1:demes) {
      pop_counts[[k]] <- cbind(1:max_time, rcpp_to_mat(output_raw$pop_counts[[k]]))
      colnames(pop_counts[[k]]) <- c("time", paste0("Ih", 1:max_infections))
    }
    names(pop_counts) <- paste0("deme", 1:demes)
    
    # add to project
    proj$pop_counts <- pop_counts
  }
  
  # process samp_hosts
  samp_hosts <- list()
  for (i in 1:length(output_raw$samp_hosts)) {
    samp_hosts[[i]] <- list()
    for (j in 1:length(output_raw$samp_hosts[[i]])) {
      v <- output_raw$samp_hosts[[i]][[j]]
      samp_hosts[[i]][[j]] <- cbind(host_ID = as.integer(names(v)),
                                    infections = v)
    }
  }
  names(samp_hosts) <- paste0("deme", 1:demes)
  
  # add to project
  proj$samp_hosts <- samp_hosts
  
  # end timer
  message(sprintf("completed in %s seconds", round(Sys.time() - t0, 2)))
  
  # return invisibly
  invisible(proj)
}

#------------------------------------------------
#' @title Prune infection history
#'
#' @description Prune infection history
#'
#' @details TODO
#' 
#' @param proj current covfefe project
#'
#' @export
#' @examples
#' # TODO

prune <- function(proj) {
  
  # TODO - check inputs
  # check infection history present
  # check times increasing
  
  # extract sample times
  samp_times <- unique(proj$samp_mat[,1])
  
  # define argument list
  args <- list(samp_times = samp_times)
  
  # start timer
  t0 <- Sys.time()
  
  # run efficient C++ function
  output_raw <- prune_cpp(proj$infection_history,
                          proj$samp_hosts,
                          args)
  
  # add to project
  proj$pruned_tree <- output_raw$pruned
  
  # end timer
  message(sprintf("completed in %s seconds", round(Sys.time() - t0, 2)))
  
  # return invisibly
  invisible(proj)
}

#------------------------------------------------
#' @title Delete infection history list from covfefe project
#'
#' @description Delete infection history list from covfefe project
#'
#' @details TODO
#' 
#' @param proj current covfefe project
#'
#' @export
#' @examples
#' # TODO

delete_infection_history <- function(proj) {
  
  # TODO - check inputs
  
  # add to project
  proj["infection_history"] <- list(NULL)
  
  # return invisibly
  invisible(proj)
}

#------------------------------------------------
#' @title Define oocyst count distribution
#'
#' @description Define oocyst count distribution
#'
#' @details TODO
#' 
#' @param proj current covfefe project
#' @param oocysts TODO
#'
#' @export
#' @examples
#' # TODO

define_oocyst_distribution <- function(proj, oocysts) {
  
  # TODO - check inputs
  
  # add to project
  proj$distributions$oocysts <- oocysts/sum(oocysts)
  
  # TODO - option to plot?
  
  # return invisibly
  invisible(proj)
}

#------------------------------------------------
#' @title Define hepatocyte count distribution
#'
#' @description Define hepatocyte count distribution
#'
#' @details TODO
#' 
#' @param proj current covfefe project
#' @param hepatocytes TODO
#'
#' @export
#' @examples
#' # TODO

define_hepatocyte_distribution <- function(proj, hepatocytes) {
  
  # TODO - check inputs
  
  # add to project
  proj$distributions$hepatocytes <- hepatocytes/sum(hepatocytes)
  
  # TODO - option to plot?
  
  # return invisibly
  invisible(proj)
}

#------------------------------------------------
#' @title Simulate genotypes from pruned infection tree
#'
#' @description Simulate genotypes from pruned infection tree
#'
#' @details TODO
#' 
#' @param proj current covfefe project
#' @param loci TODO
#' @param recom_rate TODO
#'
#' @export
#' @examples
#' # TODO

sim_genotypes <- function(proj, loci, recom_rate) {
  
  # TODO - check inputs
  
  # extract sample times
  samp_times <- unique(proj$samp_mat[,1])
  
  # find null chromosomes
  loci_null <- rep(FALSE, length(loci))
  for (i in 1:length(loci)) {
    if (is.null(loci[[i]])) {
      loci_null[i] <- TRUE
      loci[[i]] <- -9
    }
  }
  
  # define argument list
  args <- list(samp_times = samp_times,
               pruned = proj$pruned_tree,
               dist_oocysts = proj$distributions$oocysts,
               dist_hepatocytes = proj$distributions$hepatocytes,
               loci = loci,
               loci_null = loci_null,
               recom_rate = recom_rate)
  
  # start timer
  t0 <- Sys.time()
  
  # run efficient C++ function
  output_raw <- sim_genotypes_cpp(proj$samp_hosts, args)
  
  # add to project
  proj$genotypes <- output_raw$genotypes
  
  # end timer
  message(sprintf("completed in %s seconds", round(Sys.time() - t0, 2)))
  
  # return invisibly
  invisible(proj)
}
