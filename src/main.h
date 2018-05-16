
#pragma once

#include <Rcpp.h>

//------------------------------------------------
// simplify line list
Rcpp::List simplify_line_list_cpp(Rcpp::List line_list, Rcpp::List args);

//------------------------------------------------
// simulate genotypes
Rcpp::List sim_genotypes_cpp(Rcpp::List args);
