
#pragma once

#include <Rcpp.h>

//------------------------------------------------
// parameters of individual-based simulation model
class Parameters {
  
public:
  
  // PUBLIC OBJECTS
  
  // model parameters
  static int max_time;
  static double a;
  static double mu;
  static int u;
  static int v;
  static int g;
  static double r;
  static std::vector<double> b;
  static double c;
  static std::vector<int> Ih_vec;
  static std::vector<int> H_vec;
  static std::vector<int> M_vec;
  static int max_innoculations;
  static std::vector<std::vector<double>> delta_mig;
  static std::vector<double> demog;
  static bool output_counts;
  static bool output_innoculations;
  static bool output_infection_history;
  static int n_demes;
  
  // daily probability of mosquito death
  static double prob_v_death;
  static double prob_h_recovery;
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Parameters() {};
  Parameters(const Rcpp::List &args);
  
};