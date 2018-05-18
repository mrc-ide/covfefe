
#pragma once

#include <Rcpp.h>

//------------------------------------------------
// simulate population and draw blood stage hosts at designated times
Rcpp::List draw_hosts(Rcpp::List &infection_history, Rcpp::List &args);

//------------------------------------------------
// prune infection tree
Rcpp::List prune_cpp(Rcpp::List &infection_history, Rcpp::List &samp_hosts_raw, Rcpp::List &args);

//------------------------------------------------
// simulate genotypes
Rcpp::List sim_genotypes_cpp(Rcpp::List &samp_hosts_raw, Rcpp::List &args);

