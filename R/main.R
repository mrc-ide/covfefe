
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
#' @param prob_infection probability a human becomes infected after being bitten
#'   by an infected mosquito.
#' @param prob_acute probability an infection goes through an acute phase.
#' @param prob_AC probability of acute infection transitioning to chronic before
#'   clearing, as opposed to clearing directly.
#' @param duration_acute vector or list specifying probability distribution of
#'   time (in days) of acute phase of disease. If a list then the first element
#'   specifies the distribution for the first incident of acute disease, the
#'   second element for the second incident of acute disease and so on (the
#'   final distribution is used for all remaining incidents). If a vector then
#'   the same distribution is used for all incidents of acute disease.
#' @param duration_chronic equivalent to \code{duration_acute} but for chronic
#'   phase of disease.
#' @param infectivity_acute probability a mosquito becomes infected after biting
#'   a human host in the acute phase.
#' @param infectivity_chronic probability a mosquito becomes infected after
#'   biting a human host in the chronic phase.
#' @param max_innoculations maximum number of innoculations that an individual
#'   can hold simultaneously.
#'
#' @export

define_epi_parameters <- function(project,
                                  a = 0.3,
                                  p = 0.85,
                                  mu = -log(p),
                                  u = 12,
                                  v = 10,
                                  g = 12,
                                  prob_infection = 0.6,
                                  prob_acute = c(1,0.5,0),
                                  prob_AC = 0.2,
                                  duration_acute = dgeom(1:25, 1/5),
                                  duration_chronic = dgeom(1:1000, 1/200),
                                  infectivity_acute = 0.07,
                                  infectivity_chronic = 0.07,
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
  assert_pos(prob_infection)
  assert_bounded(prob_infection)
  assert_pos(prob_acute)
  assert_bounded(prob_acute)
  assert_single_pos(prob_AC)
  assert_bounded(prob_AC)
  if (!is.list(duration_acute)) {
    duration_acute <- list(duration_acute)
  }
  mapply(assert_pos, duration_acute)
  if (!is.list(duration_chronic)) {
    duration_chronic <- list(duration_chronic)
  }
  mapply(assert_pos, duration_chronic)
  assert_pos(infectivity_acute)
  assert_bounded(infectivity_acute)
  assert_pos(infectivity_chronic)
  assert_bounded(infectivity_chronic)
  assert_single_pos_int(max_innoculations, zero_allowed = FALSE)
  
  # normalise distributions
  for (i in 1:length(duration_acute)) {
    duration_acute[[i]] <- duration_acute[[i]]/sum(duration_acute[[i]])
  }
  for (i in 1:length(duration_chronic)) {
    duration_chronic[[i]] <- duration_chronic[[i]]/sum(duration_chronic[[i]])
  }
  
  # modify project
  project$sim_parameters$epi_parameters <- list(a = a,
                                                p = p,
                                                mu = mu,
                                                u = u,
                                                v = v,
                                                g = g,
                                                prob_infection = prob_infection,
                                                prob_acute = prob_acute,
                                                prob_AC = prob_AC,
                                                duration_acute = duration_acute,
                                                duration_chronic = duration_chronic,
                                                infectivity_acute = infectivity_acute,
                                                infectivity_chronic = infectivity_chronic,
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
#' @description Define deography used in simulation. The raw input is a life
#'   table giving the probability of death in each one-year age group. This in
#'   turn is used to derive the distribution of age of death, and the stable age
#'   distribution.
#'
#' @param project covfefe project.
#' @param life_table vector specifying probability of death in each one-year age
#'   group. Final value must be 1 to ensure a closed population.
#'
#' @export

define_demograpy <- function(project,
                             life_table = NULL) {
  
  # check inputs
  assert_custom_class(project, "covfefe_project")
  if (is.null(life_table)) {
    life_table <- covfefe_file("mali_life_table.rds")
  }
  assert_pos(life_table)
  assert_bounded(life_table)
  assert_eq(life_table[length(life_table)], 1, message = "the final value in the life table must be 1, representing a 100%% chance of dying, to ensure a closed population")
  
  # compute distribution of age of death
  n <- length(life_table)
  age_death <- rep(0,n)
  remaining <- 1
  for (i in 1:n) {
    age_death[i] <- remaining*life_table[i]
    remaining <- remaining*(1 - life_table[i])
  }
  
  # convert life table to transition matrix
  m <- matrix(0,n,n)
  m[col(m) == (row(m)+1)] <- 1 - life_table[1:(n-1)]
  m[,1] <- 1 - rowSums(m)
  
  # convert to rates
  r = m - diag(n)
  
  # compute Eigenvalues of the rate matrix
  E = eigen(t(r))
  
  # there should be one Eigenvalue that is zero (up to limit of computational
  # precision). Find which Eigenvalue this is
  w <- which.min(abs(E$values))
  
  # the stable solution is the corresponding Eigenvector, suitably normalised
  age_stable <- Re(E$vectors[,w]/sum(E$vectors[,w]))
  
  # modify project
  project$sim_parameters$demography <- list(life_table = life_table,
                                            age_death = age_death,
                                            age_stable = age_stable)
  
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
#' @param output_daily_counts whether to output daily counts of key quantities,
#'   such as the number of infected hosts and the EIR.
#' @param output_age_distributions whether to output complete age distributions
#' @param output_age_times a vector of times at which complete age distributions
#'   are output.
#' @param output_infection_history whether to output complete infection history.
#' @param silent whether to write messages to console.
#'
#' @export

run_sim <- function(project,
                    max_time = 365,
                    output_daily_counts = TRUE,
                    output_age_distributions = TRUE,
                    output_age_times = max_time,
                    output_infection_history = FALSE,
                    filepath_migration,
                    silent = FALSE) {
  
  # check inputs
  assert_custom_class(project, "covfefe_project")
  assert_single_pos_int(max_time)
  assert_single_logical(output_daily_counts)
  assert_single_logical(output_age_distributions)
  assert_vector(output_age_times)
  assert_pos_int(output_age_times)
  assert_leq(output_age_times, max_time)
  assert_single_logical(output_infection_history)
  assert_single_string(filepath_migration)
  assert_single_logical(silent)
  
  # get useful quantities
  n_demes <- length(project$sim_parameters$deme_parameters$H)
  n_output_age_times <- length(output_age_times)
  
  # create argument list
  args <- project$sim_parameters
  args$run_parameters <- list(max_time = max_time,
                              output_daily_counts = output_daily_counts,
                              output_age_distributions = output_age_distributions,
                              output_age_times = output_age_times,
                              output_infection_history = output_infection_history,
                              filepath_migration = filepath_migration,
                              silent = silent)
  
  # run efficient C++ function
  output_raw <- indiv_sim_cpp(args)
  
  # process daily counts etc.
  daily_counts <- NULL
  if (output_daily_counts) {
    daily_counts <- list()
    for (k in 1:n_demes) {
      daily_counts[[k]] <- data.frame(H = output_raw$H_store[[k]],
                                      Sh = output_raw$Sh_store[[k]],
                                      Lh = output_raw$Lh_store[[k]],
                                      Ah = output_raw$Ah_store[[k]],
                                      Ch = output_raw$Ch_store[[k]],
                                      Sv = output_raw$Sv_store[[k]],
                                      Lv = output_raw$Lv_store[[k]],
                                      Iv = output_raw$Iv_store[[k]],
                                      EIR = output_raw$EIR_store[[k]]*365)
    }
    names(daily_counts) <- paste0("deme", 1:n_demes)
  }
  
  # process age distribution
  age_distributions <- NULL
  if (output_age_distributions) {
    age_distributions <- list()
    for (k in 1:n_demes) {
      age_distributions[[k]] <- list()
      for (i in 1:n_output_age_times) {
        H_age <- output_raw$H_age_store[[k]][[i]]
        Sh_age <- output_raw$Sh_age_store[[k]][[i]]
        Lh_age <- output_raw$Lh_age_store[[k]][[i]]
        Ah_age <- output_raw$Ah_age_store[[k]][[i]]
        Ch_age <- output_raw$Ch_age_store[[k]][[i]]
        prev_Lh_age <- Lh_age/H_age
        prev_Ah_age <- Ah_age/H_age
        prev_Ch_age <- Ch_age/H_age
        inc_Lh_age <- output_raw$inc_Lh_age_store[[k]][[i]]*365
        inc_Ah_age <- output_raw$inc_Ah_age_store[[k]][[i]]*365
        
        age_distributions[[k]][[i]] <- data.frame(H = H_age,
                                                  Sh = Sh_age,
                                                  Lh = Lh_age,
                                                  Ah = Ah_age,
                                                  Ch = Ch_age,
                                                  prev_Lh = prev_Lh_age,
                                                  prev_Ah = prev_Ah_age,
                                                  prev_Ch = prev_Ch_age,
                                                  inc_Lh = inc_Lh_age,
                                                  inc_Ah = inc_Ah_age)
      }
      names(age_distributions[[k]]) <- paste0("time", 1:n_output_age_times)
    }
    names(age_distributions) <- paste0("deme", 1:n_demes)
  }
  
  ret <- list(daily_counts = daily_counts,
              age_distributions = age_distributions)
  
  return(ret)
  
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
