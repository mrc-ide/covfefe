
#include "host.h"
#include "misc.h"

using namespace std;

//------------------------------------------------
// Simulate genotypes
// [[Rcpp::export]]
Rcpp::List sim_genotypes_cpp(Rcpp::List args) {

  // print message to console
  Rcpp::Rcout << "running C++ dummy1_cpp function \n";
  
  // get input arguments
  Rcpp::List args_parameters = args("parameters");
  
  double r = Rcpp::as<double>(args_parameters("recom_rate"));
  int demes = Rcpp::as<int>(args("demes"));
  int max_time = Rcpp::as<int>(args("max_time"));
  vector<vector<vector<int>>> durations = rcpp_to_array_int(args("durations"));
  vector<vector<vector<int>>> migrations = rcpp_to_array_int(args("migrations"));
  
  // initialise human hosts in each deme
  vector<vector<host>> y(demes);
  for (int k=0; k<demes; k++) {
    int n_inf = durations[k][0].size();
    for (int i=0; i<n_inf; i++) {
      y[k].push_back( host(durations[k][0][i]) );
    }
  }
  
  //-------------------------
  
  // step through time
  for (int t=0; t<max_time; t++) {
    for (int k=0; k<demes; k++) {
      
      // implement migration
      
      
    }
  }
  
  
  //Rcpp::stop("foobar");
  
  
  //-------------------------
  
  // return as Rcpp list
  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("x_modified") = r);
  return ret;
}
