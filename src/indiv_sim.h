
#pragma once

#include <Rcpp.h>

//------------------------------------------------
// draw from individual-based model
Rcpp::List indiv_sim_cpp(Rcpp::List args);

void draw_migration2(std::vector<std::vector<int>> delta_mig, std::vector<std::vector<double>> &mig, std::vector<int> &H, std::vector<int> &mig_order);

//------------------------------------------------
// draw migration events
void draw_migration(std::vector<std::vector<int>> delta_mig, std::vector<std::vector<double>> &mig, std::vector<int> &H, std::vector<int> &mig_order, std::vector<int> &rowsum, std::vector<int> &colsum);