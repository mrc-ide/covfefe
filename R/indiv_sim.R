
# -----------------------------------
#' @title indiv_sim
#'
#' @description Draws from individual-based model
#' 
#' @details TODO
#'
#' @param max_time run simulation for this many days
#' @param a human blood feeding rate. The proportion of mosquitoes that feed on humans each day.
#' @param p mosquito probability of surviving one day.
#' @param mu mosquito instantaneous death rate. mu = -log(p) unless otherwise specified.
#' @param u intrinsic incubation period. The number of days from infection to becoming infectious in a human host.
#' @param v extrinsic incubation period. The number of days from infection to becoming infectious in a mosquito.
#' @param g lag time between human bloodstage infection and production of gametocytes.
#' @param r daily recovery rate.
#' @param b probability a human becomes infected after being bitten by an infected mosquito.
#' @param c probability a mosquito becomes infected after biting an infected human.
#' @param Ih_init initial number of infectious humans in each deme.
#' @param Iv_init initial number of infectious mosquitoes in each deme.
#' @param H human population size in each deme.
#' @param M mosquito population size (number of adult female mosquitoes) in each deme.
#' @param migration_matrix matrix giving the probability of migrating from current deme (rows) to new deme (columns).
#' @param H_auto if TRUE then human population sizes are chosen automatically based on the migration matrix. The total number of humans is equal to sum(H).
#'
#' @export

indiv_sim <- function(max_time = 100, a = 0.3, p = 0.9, mu = -log(p), u = 22, v = 10, g = 10, r = 1/200, b = 1, c = 1, Ih_init = 10, H = 100, M = 100, max_clonality = 5, migration_matrix = matrix(1), H_auto = FALSE) {
  
  # check arguments
  assert_pos_int(max_time, zero_allowed = FALSE)
  assert_pos(a)
  assert_bounded(p, left = 0, inclusive_left = FALSE)
  assert_pos(mu)
  assert_pos_int(u, zero_allowed = FALSE)
  assert_pos_int(v, zero_allowed = FALSE)
  assert_pos_int(g, zero_allowed = FALSE)
  assert_pos(r)
  assert_bounded(b)
  assert_bounded(c)
  assert_pos_int(Ih_init)
  assert_pos_int(H, zero_allowed = FALSE)
  assert_pos_int(M, zero_allowed = FALSE)
  assert_pos_int(max_clonality, zero_allowed = FALSE)
  if (!H_auto) {
    assert_same_length(Ih_init, H, M)
  }
  assert_same_length(Ih_init, M)
  assert_that(all(Ih_init < H))
  assert_bounded(migration_matrix)
  demes <- length(M)
  assert_that(nrow(migration_matrix) == demes, msg = sprintf("migration matrix must have %s rows and columns to match other inputs", demes))
  assert_that(isTRUE(all.equal(rowSums(migration_matrix), rep(1,demes))), msg = "each row of migration matrix must sum to 1")
  
  # calculate expected human population sizes for given migration matrix
  mig_eigen <- Re(eigen(t(migration_matrix))$vectors[,1])
  mig_eigen <- mig_eigen/sum(mig_eigen)
  H_expected <- round(sum(H)*mig_eigen)
  if (H_auto) {
    H <- H_expected
    message(paste0("auto human population size H = c(", paste(H, collapse=", "), ")"))
  }
  
  # check that equilibrium solution of migration matrix is compatible with human
  # population sizes
  if (!isTRUE(all.equal(H, H_expected))) {
    stop(paste0("migration_matrix and H incompatible. Migration with these rates will lead to H deviating over time. To rectify this, the equilibrium solution of the migration_matrix must lead to the same proportions as found in H. For the given migration matrix, the required proportions are H <- c(", paste(H_expected, collapse = ", "), ")."))
  }
  
  # calculate number of individuals that move each generation
  delta_mig <- round(outer(H, rep(1,demes)) * migration_matrix * (1-diag(1,demes)))
  
  # it is still possible that delta_mig does not represent stable migration due
  # to rounding errors Calculate discrepancies and alter matrix until 
  # completely stable
  discrep <- rep(1,demes)
  for (k in 1:demes) {
    # calculate discrepancies between rows and columns
    discrep[k] <- sum(delta_mig[k,]) - sum(delta_mig[,k])
  }
  for (k in 1:(demes-1)) {
    # adjust delta_mig to reduce discrepancies
    discrep_diff <- abs(discrep[k] + discrep[(k+1):demes])
    best_diff <- k + which.min(discrep_diff)
    delta_mig[k,best_diff] <- delta_mig[k,best_diff] - discrep[k]
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
               Ih_init = Ih_init,
               H = H,
               M = M,
               max_clonality = max_clonality,
               delta_mig = mat_to_rcpp(delta_mig),
               demes = demes
               )
  
  t0 <- Sys.time()
  
  # run efficient C++ function
  output_raw <- indiv_sim_cpp(args)
  
  # process output
  clonality <- list()
  for (k in 1:demes) {
    clonality[[k]] <- rcpp_to_mat(output_raw$clonality[[k]])
  }
  line_list <- rcpp_to_mat(output_raw$line_list)
  
  message(sprintf("completed in %s seconds", round(Sys.time() - t0, 2)))
  
  # return as list
  ret <- list(clonality = clonality,
              line_list = line_list)
  return(ret)
}
