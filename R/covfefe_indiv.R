
#------------------------------------------------
#' @title Simulate from simple individual-based model
#'
#' @description Simulate from simple individual-based model
#' 
#' @details TODO
#'
#' @param max_time run simulation for this many days.
#' @param a human blood feeding rate. The proportion of mosquitoes that feed on
#'   humans each day.
#' @param p mosquito probability of surviving one day.
#' @param mu mosquito instantaneous death rate. mu = -log(p) unless otherwise
#'   specified.
#' @param u intrinsic incubation period. The number of days from infection to
#'   blood-stage infection in a human host
#' @param v extrinsic incubation period. The number of days from infection to
#'   becoming infectious in a mosquito.
#' @param g lag time between human blood-stage infection and production of
#'   gametocytes.
#' @param r daily recovery rate.
#' @param b probability a human becomes infected after being bitten by an
#'   infected mosquito.
#' @param c probability a mosquito becomes infected after biting an infected
#'   human.
#' @param Eh_init vector or scalar specifying the initial number of infected
#'   humans in each deme. If scalar then the same value is used for all demes.
#' @param H vector or scalar specifying human population size in each deme. If 
#'   scalar then the same value is used for all demes. Alternatively, if
#'   \code{H_auto} is \code{TRUE} then this parameter is ignored and human
#'   population sizes are chosed automatically from the migration matrix.
#' @param M vector or scalar specifying mosquito population size (strictly the
#'   number of adult female mosquitoes) in each deme. If scalar then the same
#'   value is used for all demes.
#' @param max_innoculations maximum number of innoculations that an individual
#'   can hold simultaneously.
#' @param migration_matrix matrix giving the probability of migrating from 
#'   current deme (in rows) to new deme (in columns). Diagonal elements are
#'   ignored.
#' @param H_auto if TRUE then human population sizes are chosen automatically
#'   based on the migration matrix. In this case the total number of humans is
#'   equal to sum(H).
#' @param demog distribution of the number of hosts in each 1-year age bin.
#' @param output_counts whether to output daily counts of key quantities, such
#'   as the number of infected hosts and the EIR.
#' @param output_innoculations whether to output daily counts of the number of 
#'   individuals with each possible number of innoculations (from 0 up to
#'   \code{max_innoculations}).
#' @param output_infection_history whether to output complete infection history.
#' @param silent whether to write messages to console.
#'
#' @export
#' @examples
#' # TODO

sim_indiv <- function(max_time = 365, a = 0.3, p = 0.9, mu = -log(p), u = 10, v = 10, g = 10, r = 1/20, b = 1, c = 1, Eh_init = 10, H = 1000, M = 1000, max_innoculations = 5, migration_matrix = matrix(1), H_auto = FALSE, demog = dgeom(1:100, prob=1/20), output_counts = TRUE, output_innoculations = TRUE, output_infection_history = FALSE, silent = FALSE) {
  
  # check arguments
  assert_pos_int(max_time, zero_allowed = FALSE)
  assert_pos(a)
  assert_bounded(p, inclusive_left = FALSE, inclusive_right = FALSE)
  assert_pos(mu)
  assert_pos_int(u, zero_allowed = FALSE)
  assert_pos_int(v, zero_allowed = FALSE)
  assert_pos_int(g, zero_allowed = FALSE)
  assert_pos(r)
  assert_bounded(b)
  assert_bounded(c)
  assert_pos_int(Eh_init)
  assert_pos_int(H, zero_allowed = FALSE)
  assert_pos_int(M)
  assert_pos_int(max_innoculations, zero_allowed = FALSE)
  assert_square_matrix(migration_matrix)
  assert_bounded(migration_matrix)
  assert_that(all(rowSums(migration_matrix) <= 1), msg = "row sums of migration matrix cannot exceed 1")
  n_demes <- nrow(migration_matrix)
  assert_in(length(Eh_init), c(1, n_demes))
  assert_in(length(H), c(1, n_demes))
  assert_in(length(M), c(1, n_demes))
  assert_logical(H_auto)
  assert_pos(demog)
  assert_logical(output_counts)
  assert_logical(output_innoculations)
  assert_logical(output_infection_history)
  
  # force some scalar inputs to vector
  Eh_init <- force_vector(Eh_init, n_demes)
  H <- force_vector(H, n_demes)
  M <- force_vector(M, n_demes)
  assert_leq(Eh_init, H)
  
  # fill in migration matrix diagonal elements and test if there is any migration
  diag(migration_matrix) <- 0
  diag(migration_matrix) <- 1 - rowSums(migration_matrix)
  any_migration <- (sum(diag(migration_matrix)) < n_demes)
  
  # if any migration then calculate expected human population sizes
  if (any_migration) {
    
    # calculate expected human population sizes from equilibrium solution to
    # given migration matrix
    mig_eigen <- Re(eigen(t(migration_matrix))$vectors[,1])
    mig_eigen <- mig_eigen/sum(mig_eigen)
    H_expected <- round(sum(H)*mig_eigen)
    
    # it is still possible that H_expected is not stable due to rounding errors,
    # therefore step forward until H_expected no longer changes
    H_prev <- rep(0,n_demes)
    while(any(H_expected != H_prev)) {
      H_prev <- H_expected
      H_expected <- round(as.vector(crossprod(H_expected, migration_matrix)))
    }
  } else {  # if no migration
    
    # share expected human hosts equally between demes
    H_expected <- round(rep(sum(H)/n_demes, n_demes))
  }
  
  # set human population sizes automatically
  if (H_auto) {
    H <- H_expected
  }
  
  # check that human population sizes unchanged after migration
  H_stepforward <- round(as.vector(crossprod(H, migration_matrix)))
  if (any(H != H_stepforward)) {
    stop(paste0("migration_matrix and H incompatible. Migration with these rates will lead to H changing over time. Example stable proportions are H <- c(", paste(H_expected, collapse = ", "), ")."))
  }
  
  # calculate number of individuals that move each generation
  delta_mig <- round(outer(H, rep(1,n_demes)) * migration_matrix * (1-diag(1,n_demes)))
  
  # it is still possible that delta_mig does not represent stable migration due 
  # to rounding errors Calculate discrepancies and alter matrix until completely
  # stable
  if (any_migration) {
    discrep <- rowSums(delta_mig) - colSums(delta_mig)
    for (k in 1:(n_demes-1)) {
      discrep_diff <- abs(discrep[k] + discrep[(k+1):n_demes])
      best_diff <- k + which.min(discrep_diff)
      delta_mig[k,best_diff] <- delta_mig[k,best_diff] - discrep[k]
      discrep <- rowSums(delta_mig) - colSums(delta_mig)
    }
  }
  
  # define argument list
  args <- list(max_time = max_time,
               a = a,
               mu = mu,
               u = u,
               v = v,
               g = g,
               r = r,
               b = b,
               c = c,
               Eh_init = Eh_init,
               H = H,
               M = M,
               max_innoculations = max_innoculations,
               delta_mig = mat_to_rcpp(delta_mig),
               demog = demog,
               output_counts = output_counts,
               output_innoculations = output_innoculations,
               output_infection_history = output_infection_history,
               R_draws = R_draws
  )
  
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
# overload print() function for covfefe_indiv
# (not exported)
#' @noRd
print.covfefe_indiv <- function(x, ...) {

  # print without class name
  print(unclass(x))
  
  # return invisibly
  invisible(x)
}

#------------------------------------------------
# overload summary() function for covfefe_indiv
# (not exported)
#' @noRd
summary.covfefe_indiv <- function(x, ...) {

  print("TODO")

}

#------------------------------------------------
# determine if object is of class covfefe_indiv
# (not exported)
#' @noRd
is.covfefe_indiv <- function(x) {
  inherits(x, "covfefe_indiv")
}
