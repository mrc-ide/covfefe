
#include "main.h"

using namespace std;

//------------------------------------------------
// Example function
// [[Rcpp::export]]
Rcpp::List dummy1_cpp(Rcpp::List args) {

  // print message to console
  Rcpp::Rcout << "running C++ dummy1_cpp function \n";

  // get inputs from Rcpp format to base C++ format
  double a = Rcpp::as<double>(args("a"));
  vector<double> x = Rcpp::as<vector<double>>(args("x"));

  // do something
  for (int i=0; i<int(x.size()); i++) {
    x[i] *= a;
  }

  // return as Rcpp list
  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("x_modified") = x);
  return ret;
}
