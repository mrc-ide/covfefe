
# -----------------------------------
# default ggplot parameters
# (not exported)
#' @noRd
ggplot_default <- function(x = 1) {
  if (x==1) {
    ret <- ggplot() + theme_bw()
  }
  invisible(ret)
}

#------------------------------------------------
#' @title Plot daily number of innoculated infections
#'
#' @description Plot daily number of hosts with each possible number of innoculated infections.
#'
#' @details TODO
#' 
#' @param x TODO
#' @param deme TODO
#' @param ymax TODO
#'
#' @export
#' @examples
#' # TODO

plot_innoculations <- function(x, deme = 1, ymax = NULL) {
  
  # check inputs
  assert_covfefe_indiv(x)
  assert_non_null(x$innoculations)
  demes <- length(x$innoculations)
  assert_leq(deme, demes)
  
  # split raw innoculation counts
  innoculations_raw <- as.data.frame(x$innoculations[[deme]])
  times <- innoculations_raw[,1]
  innoculations <- innoculations_raw[,-(1:2),drop=FALSE]
  
  # set default ymax
  ymax <- define_default(ymax, max(rowSums(innoculations)))
  
  # get data into ggplot format
  n <- ncol(innoculations)
  df <- data.frame(time = rep(times, n),
                   value = unlist(innoculations),
                   n_innoculations = rep(1:n, each = length(times)))
  df$n_innoculations <- factor(df$n_innoculations, levels = as.character(n:1))
  
  # define colour palette
  cr <- colorRampPalette(c("firebrick4", "indianred1"))
  
  # produce plot
  ret <- ggplot_default()
  ret <- ret + coord_cartesian(ylim = c(0, ymax))
  ret <- ret + geom_area(aes(x = time, y = value, fill = n_innoculations), data = df)
  ret <- ret + scale_fill_manual(values = cr(n))
  
  # return plot object
  return(ret)
}

#------------------------------------------------
#' @title Plot \code{covfefe_indiv} object
#'
#' @description Default plot for \code{covfefe_indiv} object.
#'
#' @details TODO
#' 
#' @param x TODO
#' @param y TODO
#' @param ... TODO
#'
#' @export
#' @examples
#' # TODO

plot.covfefe_indiv <- function(x, y, ...) {
  plot_innoculations(x, 1, ...)
}

#------------------------------------------------
#' @title Simulate and plot EIR vs. innoculations
#'
#' @description Simulate and plot EIR vs. innoculations
#'
#' @details TODO
#' 
#' @param x TODO
#' @param y TODO
#' @param ... TODO
#'
#' @export
#' @examples
#' # TODO

sim_plot_EIR_innoculations <- function(max_time = 365, a = 0.3, p = 0.9, mu = -log(p), u = 22, v = 10, g = 10, r = 1/20, b = 1, c = 1, H = 100, M = c(5e1, 1e2, 2e2, 5e2, 1e3, 2e3, 5e3, 1e4), max_innoculations = 5, reps = 10) {
  
  # object for storing results of simulation
  n <- length(M)
  EIR <- (1:n)/1e6
  innoculations <- matrix(0, n, max_innoculations)
  
  # loop over mosquito densities and reps
  for (i in 1:n) {
    for (rep in 1:reps) {
      
      # simulate from model
      sim <- sim_indiv(max_time = max_time,
                       a = a,
                       p = p,
                       mu = mu,
                       u = u,
                       v = v,
                       g = g,
                       r = r,
                       b = b,
                       c = c,
                       Ih_init = H,
                       H = H,
                       M = M[i],
                       max_innoculations = max_innoculations,
                       migration_matrix = diag(1),
                       H_auto = FALSE,
                       output_counts = TRUE,
                       output_innoculations = TRUE,
                       output_infection_history = FALSE,
                       silent = TRUE)
      
      # get EIR and innoculations at final time point
      EIR[i] <- EIR[i] + sim$counts$deme1[max_time,"EIR"]/reps
      innoculations[i,] <- innoculations[i,] + 100*sim$innoculations$deme1[max_time,-(1:2)]/H/reps
    }
    
  }
  
  # get results into ggplot format
  df <- data.frame(EIR = rep(EIR, times=max_innoculations),
                   value = as.vector(innoculations),
                   n_innoculations = rep(1:max_innoculations, each=n))
  df$n_innoculations <- factor(df$n_innoculations, levels = as.character(n:1))
  
  # define colour palette
  cr <- colorRampPalette(c("firebrick4", "indianred1"))
  
  # produce plot
  plot1 <- ggplot_default()
  plot1 <- plot1 + coord_cartesian(ylim = c(0, 100))
  plot1 <- plot1 + geom_area(aes(x = EIR, y = value, fill = n_innoculations), data = df)
  plot1 <- plot1 + scale_fill_manual(values = cr(max_innoculations+1), name = "number of\ninnoculations")
  plot1 <- plot1 + scale_x_continuous(trans='log10')
  plot1 <- plot1 + xlab("EIR") + ylab("prevalence")
  print(plot1)
  
  # return list
  ret <- list(data = df,
              plot = plot1)
  return(ret)
}



