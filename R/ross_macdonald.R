
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
#' @param Eh_init initial number of infected but not infectious humans.
#' @param Ih_init initial number of infectious humans.
#' @param Ev_init initial number of infected but not yet infectious mosquitoes.
#' @param Iv_init initial number of infectious mosquitoes.
#' @param H human population size.
#' @param M mosquito population size (number of adult female mosquitoes).
#'
#' @export

ross_macdonald <- function(max_time = 100, a = 0.3, p = 0.9, mu = -log(p), u = 22, v = 10, r = 1/200, b = 1, c = 1, Eh_init = 0, Ih_init = 10, Ev_init = 0, Iv_init = 0, H = 100, M = 100) {
  
  # check arguments
  assert_pos_int(max_time, zero_allowed = FALSE)
  assert_pos(a)
  assert_bounded(p, left = 0, inclusive_left = FALSE)
  assert_pos(mu)
  assert_pos_int(u)
  assert_pos_int(v)
  assert_pos(r)
  assert_bounded(b)
  assert_bounded(c)
  assert_pos_int(Eh_init)
  assert_pos_int(Ih_init)
  assert_pos_int(Ev_init)
  assert_pos_int(Iv_init)
  assert_pos_int(H, zero_allowed = FALSE)
  assert_pos_int(M, zero_allowed = FALSE)
  
  # TODO - check u, v = 0 OK
  
  # define argument list
  args <- list(max_time = max_time,
               a = a,
               mu = mu,
               u = u,
               v = v,
               r = r,
               b = b,
               c = c,
               Eh_init = Eh_init,
               Ih_init = Ih_init,
               Ev_init = Ev_init,
               Iv_init = Iv_init,
               H = H,
               M = M
               )
  
  # run efficient C++ function
  output_raw <- ross_macdonald_cpp(args)
  
  # process output
  output <- as.data.frame(output_raw)
  output$H <- H
  output$M <- M
  output[output<0] <- NA
  
  return(output)
}
