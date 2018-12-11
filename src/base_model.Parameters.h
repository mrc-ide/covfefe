
#pragma once

#include "misc.h"

#ifdef RCPP_ACTIVE
#include <Rcpp.h>
#endif

#ifdef OLD_VERSION


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
  static std::vector<double> b;
  static int n_b;
  static std::vector<double> prob_acute;
  static int n_prob_acute;
  static double prob_AC;
  static std::vector<std::vector<double>> duration_acute;
  static int n_duration_acute;
  static std::vector<std::vector<double>> duration_chronic;
  static int n_duration_chronic;
  static double c;
  static int max_innoculations;
  
  // deme parameters
  static std::vector<int> H_vec;
  static std::vector<int> seed_infections;
  static std::vector<int> M_vec;
  static int n_demes;
  
  // demography
  static std::vector<double> life_table;
  static std::vector<double> age_death;
  static std::vector<double> age_stable;
  static int n_age;
  
  // migration
  // TODO
  
  // run parameters
  static int max_time;
  static bool output_daily_counts;
  static bool output_age_distributions;
  static std::vector<int> output_age_times;
  static int n_output_age_times;
  static bool output_infection_history;
  
  // misc parameters
  static double prob_v_death;  // daily probability of mosquito death
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Parameters();
#ifdef RCPP_ACTIVE
  Parameters(const Rcpp::List &args);
#endif
  
  // methods
  void print_summary();
  
};

#else

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
  static std::vector<double> prob_infection;
  static int n_prob_infection;
  static std::vector<double> prob_acute;
  static int n_prob_acute;
  static double prob_AC;
  static std::vector<std::vector<double>> duration_acute;
  static int n_duration_acute;
  static std::vector<std::vector<double>> duration_chronic;
  static int n_duration_chronic;
  static std::vector<double> infectivity_acute;
  static int n_infectivity_acute;
  static std::vector<double> infectivity_chronic;
  static int n_infectivity_chronic;
  static int max_innoculations;
  
  // deme parameters
  static std::vector<int> H_vec;
  static std::vector<int> seed_infections;
  static std::vector<int> M_vec;
  static int n_demes;
  
  // demography
  static std::vector<double> life_table;
  static std::vector<double> age_death;
  static std::vector<double> age_stable;
  static int n_age;
  
  // migration
  // TODO
  
  // run parameters
  static int max_time;
  static bool output_daily_counts;
  static bool output_age_distributions;
  static std::vector<int> output_age_times;
  static int n_output_age_times;
  static bool output_infection_history;
  static std::string filepath_migration;
  
  
  // misc parameters
  static double prob_v_death;  // daily probability of mosquito death
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Parameters();
#ifdef RCPP_ACTIVE
  Parameters(const Rcpp::List &args);
#endif
  
  // methods
  void print_summary();
  
};

#endif
