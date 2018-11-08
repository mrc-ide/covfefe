
#pragma once

#include "misc.h"

#ifdef RCPP_ACTIVE
#include <Rcpp.h>
#endif

//------------------------------------------------
// parameters of individual-based simulation model
class Parameters {
  
public:
  
  // PUBLIC OBJECTS
  
  // epidemiological parameters
  static double a;
  static double mu;
  static int u;
  static int v;
  static int g;
  static double r;
  static std::vector<double> b;
  static double c;
  static int max_innoculations;
  
  // deme parameters
  static std::vector<int> H_vec;
  static std::vector<int> seed_infections;
  static std::vector<int> M_vec;
  static int n_demes;
  
  // demography
  static std::vector<double> demography;
  static int n_demography;
  
  // migration
  // TODO
  
  // run parameters
  static int max_time;
  static bool output_counts;
  static std::vector<int> output_age_times;
  static int n_output_age_times;
  static bool output_infection_history;
  
  // misc parameters
  static double prob_v_death;  // daily probability of mosquito death
  static double prob_h_recovery;  // daily probability of human recovery
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Parameters();
#ifdef RCPP_ACTIVE
  Parameters(const Rcpp::List &args);
#endif
  
  // methods
  void print_summary();
  
};
