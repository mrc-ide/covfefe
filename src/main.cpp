
#include "host.h"
#include "misc.h"

using namespace std;

//------------------------------------------------
// simplify line list
// [[Rcpp::export]]
Rcpp::List simplify_line_list_cpp(Rcpp::List line_list, Rcpp::List args) {
  
  vector<int> samp_times = rcpp_to_vector_int(args["samp_times"]);
  vector<vector<int>> samp_demes = rcpp_to_mat_int(args["samp_demes"]);
  vector<vector<int>> samp_num = rcpp_to_mat_int(args["samp_num"]);
  int demes = rcpp_to_int(args["demes"]);
  int n_samp = samp_times.size();
  int max_samp_time = samp_times[n_samp-1];
  
  vector<vector<unordered_set<int>>> samp_hosts(n_samp, vector<unordered_set<int>>(demes));
  vector<map<int, int>> hosts(demes);
  
  int i = 0;
  vector<int> this_line;
  while (i<max_samp_time*4) {
    
    // migration
    this_line = rcpp_to_vector_int(line_list[i]);
    for (int j=0; j<int(this_line.size())/3; j++) {
      int this_host = this_line[3*j];
      int this_deme = this_line[3*j+1];
      int new_deme = this_line[3*j+2];
      hosts[new_deme][this_host] = hosts[this_deme][this_host];
      hosts[this_deme].erase(this_host);
    }
    i++;
    
    // infection
    this_line = rcpp_to_vector_int(line_list[i]);
    for (int j=0; j<int(this_line.size())/5; j++) {
      int this_host = this_line[5*j];
      int this_deme = this_line[5*j+1];
      if (hosts[this_deme].count(this_host)==0) { // new infection
        hosts[this_deme][this_host] = 1;
      } else {  // superinfection
        hosts[this_deme][this_host]++;
      }
    }
    i += 2; // (skip bloodstage)
    
    // recovery
    this_line = rcpp_to_vector_int(line_list[i]);
    for (int j=0; j<int(this_line.size())/3; j++) {
      int this_host = this_line[3*j];
      int this_deme = this_line[3*j+1];
      if (hosts[this_deme][this_host]==1) {
        hosts[this_deme].erase(this_host);
      } else {
        hosts[this_deme][this_host]--;
      }
    }
    i++;
    
    // sample individuals
    
    
    //for (auto it = hosts[0].begin(); it!=hosts[0].end(); ++it) {
    //  Rcpp::Rcout << it->second << " ";
    //}
    //Rcpp::Rcout << "\n";
  }
  
  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("foo") = samp_times);
  return ret;
}

//------------------------------------------------
// simulate genotypes
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
