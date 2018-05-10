
# -----------------------------------
#' @title RM1_stochastic_sync
#'
#' @description Draws from synchronous stochastic Ross-Macdonald model
#' 
#' @details Note that this model is stochastic and uses a daily time step, whereas the classic Ross-Macdonald system is deterministic and assumes continuous time. The stochastic nature of this model means we may see fade-out at low transmission, causing the average behaviour to differ from the deterministic prediction. The daily time step also means results may not agree exactly with classic Ross-Macdonald model.
#'
#' @param max_time run simulation for this many days
#' @param a human blood feeding rate. The proportion of mosquitoes that feed on humans each day.
#' @param p mosquito probability of surviving one day.
#' @param mu mosquito instantaneous death rate. mu = -log(p) unless otherwise specified.
#' @param u intrinsic incubation period. The number of days from infection to becoming infectious in a human host.
#' @param v extrinsic incubation period. The number of days from infection to becoming infectious in a mosquito .
#' @param r daily recovery rate.
#' @param b probability a human becomes infected after being bitten by an infected mosquito.
#' @param c probability a mosquito becomes infected after biting an infected human.
#' @param Ih_init initial number of infectious humans in each deme.
#' @param Iv_init initial number of infectious mosquitoes in each deme.
#' @param H human population size in each deme.
#' @param M mosquito population size (number of adult female mosquitoes) in each deme.
#' @param migration_matrix matrix giving the probability of migrating from current deme (rows) to new deme (columns).
#'
#' @export

ross_macdonald <- function(max_time = 100, a = 0.3, p = 0.9, mu = -log(p), u = 22, v = 10, r = 1/200, b = 1, c = 1, Ih_init = 10, Iv_init = 0, H = 100, M = 100, migration_matrix = matrix(0)) {
  
  # check arguments
  assert_pos_int(max_time, zero_allowed = FALSE)
  assert_pos(a)
  assert_bounded(p, left = 0, inclusive_left = FALSE)
  assert_pos(mu)
  assert_pos_int(u, zero_allowed = FALSE)
  assert_pos_int(v, zero_allowed = FALSE)
  assert_pos(r)
  assert_bounded(b)
  assert_bounded(c)
  assert_pos_int(Ih_init)
  assert_pos_int(Iv_init)
  assert_pos_int(H, zero_allowed = FALSE)
  assert_pos_int(M, zero_allowed = FALSE)
  assert_same_length(Ih_init, Iv_init, H, M)
  assert_bounded(migration_matrix)
  assert_symmetric_matrix(migration_matrix)
  assert_that(nrow(migration_matrix) == length(H))
  
  # get useful values
  demes <- length(H)
  
  # define argument list
  args <- list(max_time = max_time,
               a = a,
               mu = mu,
               u = u,
               v = v,
               r = r,
               b = b,
               c = c,
               Ih_init = Ih_init,
               Iv_init = Iv_init,
               H = H,
               M = M,
               demes = demes
               )
  
  # run efficient C++ function
  output_raw <- ross_macdonald_cpp(args)
  
  # get model states in all demes
  model_states <- list()
  for (k in 1:demes) {
    model_states[[k]] <- data.frame(Sh = output_raw$Sh[[k]],
                           Eh = output_raw$Eh[[k]],
                           Ih = output_raw$Ih[[k]],
                           Sv = output_raw$Sv[[k]],
                           Ev = output_raw$Ev[[k]],
                           Iv = output_raw$Iv[[k]],
                           H = H[k],
                           M = M[k])
    model_states[[k]][model_states[[k]]<0] <- NA
  }
  names(model_states) <- paste0("deme", 1:demes)
  
  # get durations
  durations = -9
  
  # return as list
  ret <- list(model_states = model_states,
              durations = durations)
  return(ret)
}
