
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

indiv_sim <- function(max_time = 100, a = 0.3, p = 0.9, mu = -log(p), u = 22, v = 10, g = 10, r = 1/200, b = 1, c = 1, Ih_init = 10, H = 100, M = 100, max_infections = 5, migration_matrix = matrix(1), H_auto = FALSE) {
  
  # check arguments
  assert_pos_int(max_time, zero_allowed = FALSE)
  assert_pos(a)
  assert_bounded(p, left = 0, inclusive_left = FALSE)
  assert_pos(mu)
  assert_pos_int(u, zero_allowed = FALSE)
  assert_that(max_time>=u)  # TODO - remove?
  assert_pos_int(v, zero_allowed = FALSE)
  assert_pos_int(g, zero_allowed = FALSE)
  assert_pos(r)
  assert_bounded(b)
  assert_bounded(c)
  assert_pos_int(Ih_init)
  assert_pos_int(H, zero_allowed = FALSE)
  assert_pos_int(M, zero_allowed = TRUE)
  assert_pos_int(max_infections, zero_allowed = FALSE)
  if (!H_auto) {
    assert_same_length(Ih_init, H, M)
  }
  assert_same_length(Ih_init, M)
  assert_that(all(Ih_init < H))
  assert_bounded(migration_matrix)
  demes <- length(M)
  assert_that(nrow(migration_matrix) == demes, msg = sprintf("migration matrix must have %s rows and columns to match other inputs", demes))
  assert_that(all(rowSums(migration_matrix)<=1), msg = "row sums of migration matrix cannot exceed 1")
  
  # fill in migration matrix diagonal elements and test if there is any migration
  diag(migration_matrix) <- 0
  diag(migration_matrix) <- 1 - rowSums(migration_matrix)
  any_migration <- sum(diag(migration_matrix)) < demes
  
  # calculate expected human population sizes
  if (any_migration) {
    
    # calculate expected human population sizes from equilibrium solution to
    # given migration matrix
    mig_eigen <- Re(eigen(t(migration_matrix))$vectors[,1])
    mig_eigen <- mig_eigen/sum(mig_eigen)
    H_expected <- round(sum(H)*mig_eigen)
    
    # it is still possible that H_expected is not stable due to rounding errors.
    # Step forward until H_expected no longer changes
    H_prev <- rep(0,demes)
    while(any(H_expected != H_prev)) {
      H_prev <- H_expected
      H_expected <- round(as.vector(crossprod(H_expected, migration_matrix)))
    }
  } else {
    
    # share human hosts equally between demes
    H_expected <- round(rep(sum(H)/demes, demes))
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
  delta_mig <- round(outer(H, rep(1,demes)) * migration_matrix * (1-diag(1,demes)))
  
  # it is still possible that delta_mig does not represent stable migration due
  # to rounding errors Calculate discrepancies and alter matrix until 
  # completely stable
  discrep <- rowSums(delta_mig) - colSums(delta_mig)
  for (k in 1:(demes-1)) {
    # adjust delta_mig to reduce discrepancies
    discrep_diff <- abs(discrep[k] + discrep[(k+1):demes])
    best_diff <- k + which.min(discrep_diff)
    delta_mig[k,best_diff] <- delta_mig[k,best_diff] - discrep[k]
    discrep <- rowSums(delta_mig) - colSums(delta_mig)
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
               max_infections = max_infections,
               migration_matrix = mat_to_rcpp(migration_matrix),
               demes = demes
               )
  
  t0 <- Sys.time()
  
  #print(delta_mig)
  #print(discrep)
  #print(cbind(rowSums(delta_mig), colSums(delta_mig), H))
  #return(args)
  
  # run efficient C++ function
  output_raw <- indiv_sim_cpp(args)
  
  # process output
  n_bloodstage <- list()
  for (k in 1:demes) {
    n_bloodstage[[k]] <- rcpp_to_mat(output_raw$n_bloodstage[[k]])
  }
  #line_list <- rcpp_to_mat(output_raw$line_list)
  line_list <- -9
  
  message(sprintf("completed in %s seconds", round(Sys.time() - t0, 2)))
  
  # return as list
  ret <- list(n_bloodstage = n_bloodstage,
              line_list = line_list)
  return(ret)
}
