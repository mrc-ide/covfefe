
#pragma once

#include <Rcpp.h>

//------------------------------------------------
// simulate population and draw blood stage hosts at designated times
std::vector<std::vector<std::map<int, int>>> draw_hosts(Rcpp::List &line_list, Rcpp::List &args);

//------------------------------------------------
// prune infection tree and simulate genotypes
Rcpp::List sim_genotypes_cpp(Rcpp::List line_list, Rcpp::List args);

