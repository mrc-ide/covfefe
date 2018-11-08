
#------------------------------------------------
# The following commands ensure that package dependencies are listed in the
# NAMESPACE file.

#' @useDynLib covfefe
#' @import ggplot2
#' @importFrom Rcpp sourceCpp
#' @import stats
#' @import graphics
NULL

#------------------------------------------------
#' @title Define simulation epidemiological parameters
#'
#' @description Define simulation epidemiological parameters
#'
#' @param project covfefe project.
#' @param a human blood feeding rate. The proportion of mosquitoes that feed on
#'   humans each day.
#' @param p mosquito probability of surviving one day.
#' @param mu mosquito instantaneous death rate. mu = -log(p) unless otherwise
#'   specified.
#' @param u intrinsic incubation period. The number of days from infection to
#'   blood-stage infection in a human host.
#' @param v extrinsic incubation period. The number of days from infection to
#'   becoming infectious in a mosquito.
#' @param g lag time between human blood-stage infection and production of
#'   gametocytes.
#' @param r daily recovery rate.
#' @param b probability a human becomes infected after being bitten by an
#'   infected mosquito.
#' @param c probability a mosquito becomes infected after biting an infected
#'   human.
#' @param max_innoculations maximum number of innoculations that an individual
#'   can hold simultaneously.
#'
#' @export

define_epi_parameters <- function(project,
                                  a = 0.3,
                                  p = 0.9,
                                  mu = -log(p),
                                  u = 10,
                                  v = 10,
                                  g = 10,
                                  r = 1/200,
                                  b = 1,
                                  c = 1,
                                  max_innoculations = 5) {
  
  # check inputs
  assert_custom_class(project, "covfefe_project")
  assert_single_pos(a)
  assert_bounded(a)
  assert_single_pos(p)
  assert_bounded(p)
  assert_single_pos(mu)
  assert_single_pos_int(u, zero_allowed = FALSE)
  assert_single_pos_int(v, zero_allowed = FALSE)
  assert_single_pos_int(g, zero_allowed = FALSE)
  assert_single_pos(r)
  assert_pos(b)
  assert_bounded(b)
  assert_single_pos(c)
  assert_bounded(c)
  assert_single_pos_int(max_innoculations, zero_allowed = FALSE)
  
  # modify project
  project$sim_parameters$epi_parameters <- list(a = a,
                                                p = p,
                                                mu = mu,
                                                u = u,
                                                v = v,
                                                g = g,
                                                r = r,
                                                b = b,
                                                c = c,
                                                max_innoculations = max_innoculations)
  
  # return
  invisible(project)
}

#------------------------------------------------
#' @title Define simulation deme parameters
#'
#' @description Define simulation deme parameters
#'
#' @param project covfefe project.
#' @param H vector specifying human population size in each deme.
#' @param Eh_init vector specifying the initial number of infected humans in
#'   each deme.
#' @param M vector specifying mosquito population size (strictly the number of
#'   adult female mosquitoes) in each deme.
#'
#' @export

define_deme_parameters <- function(project,
                                   H = 1000,
                                   seed_infections = 100,
                                   M = 1000) {
  
  # check inputs
  assert_custom_class(project, "covfefe_project")
  assert_pos_int(H)
  assert_pos_int(seed_infections)
  assert_pos_int(M)
  assert_same_length_multiple(H, seed_infections, M)
  assert_leq(seed_infections, H)
  
  # modify project
  project$sim_parameters$deme_parameters <- list(H = H,
                                                 seed_infections = seed_infections,
                                                 M = M)
  
  # return
  invisible(project)
}

#------------------------------------------------
#' @title Define simulation demographic parameters
#'
#' @description Define simulation demographic parameters
#'
#' @param project covfefe project.
#' @param demography vector specifying proportion of the population in each
#'   one-year age group.
#'
#' @export

define_demography <- function(project,
                              demography = dgeom(1:100, prob = 1/20)) {
  
  # check inputs
  assert_custom_class(project, "covfefe_project")
  assert_pos(demography)
  
  # modify project
  project$sim_parameters$demography <- demography/sum(demography)
  
  # return
  invisible(project)
}

#------------------------------------------------
#' @title Define simulation migration parameters
#'
#' @description Define simulation migration parameters
#'
#' @param project covfefe project.
#' @param migration_matrix TODO.
#'
#' @export

define_migration <- function(project,
                             migration_matrix = matrix(1)) {
  
  # check inputs
  assert_custom_class(project, "covfefe_project")
  
  # modify project
  project$sim_parameters$migration_matrix <- migration_matrix
  
  # return
  invisible(project)
}

#------------------------------------------------
#' @title Simulate from simple individual-based model
#'
#' @description Simulate from simple individual-based model
#'
#' @param max_time run simulation for this many days.
#' @param output_counts whether to output daily counts of key quantities, such
#'   as the number of infected hosts and the EIR.
#' @param output_age_times a vector of times at which complete age distributions
#'   are output.
#' @param output_infection_history whether to output complete infection history.
#' @param silent whether to write messages to console.
#'
#' @export

run_sim <- function(project,
                    max_time = 365,
                    output_counts = TRUE,
                    output_age_times = max_time,
                    output_infection_history = FALSE,
                    silent = FALSE) {
  
  # check inputs
  assert_custom_class(project, "covfefe_project")
  assert_single_pos_int(max_time)
  assert_single_logical(output_counts)
  assert_vector(output_age_times)
  assert_pos_int(output_age_times)
  assert_leq(output_age_times, max_time)
  assert_single_logical(output_infection_history)
  assert_single_logical(silent)
  
  # create argument list
  args <- project$sim_parameters
  args$run_parameters <- list(max_time = max_time,
                              output_counts = output_counts,
                              output_age_times = output_age_times,
                              output_infection_history = output_infection_history,
                              silent = silent)
  
  # start timer
  t0 <- Sys.time()
  
  # run efficient C++ function
  output_raw <- indiv_sim_cpp(args)
  
  message(sprintf("completed in %s seconds", round(Sys.time() - t0, 2))) 
  
  return(output_raw)
  
  # create return object
  ret <- list()
  class(ret) <- "covfefe_indiv"
  
  # process counts
  if (output_counts) {
    counts <- list()
    for (k in 1:n_demes) {
      Sh_store <- output_raw$Sh_store[[k]]
      Ih_store <- output_raw$Ih_store[[k]]
      EIR_store <- output_raw$EIR_store[[k]]
      counts[[k]] <- cbind(time = 1:max_time, 
                           Sh = Sh_store,
                           Ih = Ih_store,
                           EIR = EIR_store)
    }
    names(counts) <- paste0("deme", 1:n_demes)
    ret$counts <- counts
  }
  
  # process innoculations
  if (output_innoculations) {
    innoculations <- list()
    for (k in 1:n_demes) {
      raw_innoculations <- rcpp_to_mat(output_raw$innoculations[[k]])
      colnames(raw_innoculations) <- paste0("Ih", 0:max_innoculations)
      innoculations[[k]] <- cbind(time = 1:max_time,
                                  raw_innoculations)
    }
    names(innoculations) <- paste0("deme", 1:n_demes)
    ret$innoculations <- innoculations
  }
  
  # process infection history
  if (output_infection_history) {
    ret$infection_history <- output_raw$infection_history
  }
  
  # end timer
  if (!silent) {
    message(sprintf("completed in %s seconds", round(Sys.time() - t0, 2))) 
  }
  
  # return invisibly
  invisible(ret)
}

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
