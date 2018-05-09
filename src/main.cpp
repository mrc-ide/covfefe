
#include "main.h"

using namespace std;

//------------------------------------------------
// Simulate genotypes
// [[Rcpp::export]]
Rcpp::List sim_genotypes_cpp(Rcpp::List proj, Rcpp::List args) {

  // print message to console
  Rcpp::Rcout << "running C++ dummy1_cpp function \n";
  
  // get input arguments
  Rcpp::List proj_parameters = proj("parameters");
  
  // get inputs from Rcpp format to base C++ format
  double r = Rcpp::as<double>(proj_parameters("recom_rate"));
  //vector<double> x = Rcpp::as<vector<double>>(args("x"));

  // do something
  //for (int i=0; i<int(x.size()); i++) {
  //  x[i] *= a;
  //}

  Rcpp::stop("foobar");
  
  // return as Rcpp list
  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("x_modified") = r);
  return ret;
}
