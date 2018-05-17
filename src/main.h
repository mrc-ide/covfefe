
#pragma once

#include <Rcpp.h>

//------------------------------------------------
// simulate population and draw blood stage hosts at designated times
std::vector<std::vector<std::map<int, int>>> draw_hosts(Rcpp::List &line_list, Rcpp::List &args);

//------------------------------------------------
// prune infection tree
void prune(std::vector<int> &pruned_de_novo, std::vector<std::vector<std::pair<int, int>>> &pruned_infection, std::vector<std::vector<int>> &pruned_delete, Rcpp::List &line_list, std::vector<int> &samp_times, std::vector<std::vector<std::map<int, int>>> &samp_hosts, int demes);

//------------------------------------------------
// simulate genotypes
Rcpp::List sim_genotypes_cpp(Rcpp::List line_list, Rcpp::List args);

