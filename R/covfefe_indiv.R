
#------------------------------------------------
#' @title Simulate from simple individual-based model
#'
#' @description Simulate from simple individual-based model
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
#' @param H human population size in each deme.
#' @param M mosquito population size (number of adult female mosquitoes) in each deme.
#' @param max_infections maximum number of infections that an individual can hold simultaneously.
#' @param migration_matrix matrix giving the probability of migrating from current deme (rows) to new deme (columns).
#' @param H_auto if TRUE then human population sizes are chosen automatically based on the migration matrix. In this case the total number of humans is equal to sum(H).
#'
#' @export
#' @examples
#' # TODO

sim_indiv <- function(max_time = 365, a = 0.3, p = 0.9, mu = -log(p), u = 22, v = 10, g = 10, r = 1/20, b = 1, c = 1, Ih_init = 10, H = 1000, M = 500, max_infections = 5, migration_matrix = matrix(1), H_auto = FALSE) {
  
  # check arguments
  assert_pos_int(max_time, zero_allowed = FALSE)
  assert_pos(a)
  assert_bounded(p, left = 0, inclusive_left = FALSE)
  assert_pos(mu)
  assert_pos_int(u, zero_allowed = FALSE)
  assert_that(max_time>=u)  # TODO - remove constraint?
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
  assert_that(all(Ih_init <= H))
  assert_bounded(migration_matrix)
  demes <- length(M)
  assert_that(nrow(migration_matrix) == demes, msg = sprintf("migration matrix must have %s rows and columns to match other inputs", demes))
  assert_that(all(rowSums(migration_matrix)<=1), msg = "row sums of migration matrix cannot exceed 1")
  
  # fill in migration matrix diagonal elements and test if there is any migration
  diag(migration_matrix) <- 0
  diag(migration_matrix) <- 1 - rowSums(migration_matrix)
  any_migration <- sum(diag(migration_matrix)) < demes
  
  # if any migration then calculate expected human population sizes
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
  } else {  # if no migration
    
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
  # to rounding errors Calculate discrepancies and alter matrix until completely
  # stable
  if (any_migration) {
    discrep <- rowSums(delta_mig) - colSums(delta_mig)
    for (k in 1:(demes-1)) {
      discrep_diff <- abs(discrep[k] + discrep[(k+1):demes])
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
               Ih_init = Ih_init,
               H = H,
               M = M,
               max_infections = max_infections,
               delta_mig = mat_to_rcpp(delta_mig),
               demes = demes
  )
  
  # start timer
  t0 <- Sys.time()
  
  # run efficient C++ function
  output_raw <- indiv_sim_cpp(args)
  
  # process output
  daily_counts <- list()
  for (k in 1:demes) {
    daily_counts[[k]] <- cbind(1:max_time, rcpp_to_mat(output_raw$daily_counts[[k]]))
    colnames(daily_counts[[k]]) <- c("time", "Sh", paste0("Ih", 1:max_infections))
  }
  names(daily_counts) <- paste0("deme", 1:demes)
  
  # create custom class
  ret <- list(H = H,
              daily_counts = daily_counts,
              infection_history = output_raw$infection_history)
  class(ret) <- "covfefe_indiv"
  
  # end timer
  message(sprintf("completed in %s seconds", round(Sys.time() - t0, 2)))
  
  # return invisibly
  invisible(ret)
}

#------------------------------------------------
# overload print() function for covfefe_indiv
# (not exported)

print.covfefe_indiv <- function(x, ...) {

  # print without class name
  print(unclass(x))
  
  # return invisibly
  invisible(x)
}

#------------------------------------------------
# overload summary() function for covfefe_indiv
# (not exported)

summary.covfefe_indiv <- function(x, ...) {

  print("TODO")

}

#------------------------------------------------
# determine if object is of class covfefe_indiv
# (not exported)

is.covfefe_indiv <- function(x) {
  inherits(x, "covfefe_indiv")
}
