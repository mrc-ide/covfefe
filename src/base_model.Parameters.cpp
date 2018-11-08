
#include "base_model.Parameters.h"

#include <cmath>

using namespace std;

//------------------------------------------------
// declare static member variables

// epidemiological parameters
double Parameters::a;
double Parameters::mu;
int Parameters::u;
int Parameters::v;
int Parameters::g;
double Parameters::r;
vector<double> Parameters::b;
double Parameters::c;
int Parameters::max_innoculations;

// deme parameters
vector<int> Parameters::H_vec;
vector<int> Parameters::seed_infections;
vector<int> Parameters::M_vec;
int Parameters::n_demes;

// demography
vector<double> Parameters::demography;
int Parameters::n_demography;

// migration
// TODO

// run parameters
int Parameters::max_time;
bool Parameters::output_counts;
vector<int> Parameters::output_age_times;
int Parameters::n_output_age_times;
bool Parameters::output_infection_history;

// misc parameters
double Parameters::prob_v_death;
double Parameters::prob_h_recovery;

//------------------------------------------------
// constructor
#ifdef RCPP_ACTIVE
Parameters::Parameters() {}
#else
Parameters::Parameters() {
  
  // HARD-CODE PARAMETERS
  
  // epidemiological parameters
  a = 0.3;
  mu = -log(0.9);
  u = 10;
  v = 10;
  g = 10;
  r = 1.0/20;
  b = vector<double>(1, 1.0);
  c = 1.0;
  max_innoculations = 5;
  
  // deme parameters
  int n_demes0 = 1;
  H_vec = vector<int>(n_demes0, 1000);
  seed_infections = vector<int>(n_demes0, 100);
  M_vec = vector<int>(n_demes0, 1000);
  n_demes = n_demes0;
  
  // demography
  demography = {100};
  n_demography = int(demography.size());
  
  // migration
  // (TODO)
  
  // run parameters
  max_time = 365*5;
  output_counts = true;
  output_age_times = {max_time};
  n_output_age_times = int(output_age_times.size());
  output_infection_history = false;
  
  // misc parameters
  prob_v_death = 1 - exp(-mu);  // daily probability of mosquito death
  prob_h_recovery = 1 - exp(-r);  // daily probability of human recovery
  
}
#endif

//------------------------------------------------
// constructor
#ifdef RCPP_ACTIVE
Parameters::Parameters(const Rcpp::List &args) {
  
  // epidemiological parameters
  Rcpp::List args_epi_parameters = args["epi_parameters"];
  a = rcpp_to_double(args_epi_parameters["a"]);
  mu = rcpp_to_double(args_epi_parameters["mu"]);
  u = rcpp_to_int(args_epi_parameters["u"]);
  v = rcpp_to_int(args_epi_parameters["v"]);
  g = rcpp_to_int(args_epi_parameters["g"]);
  r = rcpp_to_double(args_epi_parameters["r"]);
  b = rcpp_to_vector_double(args_epi_parameters["b"]);
  c = rcpp_to_double(args_epi_parameters["c"]);
  max_innoculations = rcpp_to_int(args_epi_parameters["max_innoculations"]);
  
  // deme parameters
  Rcpp::List args_deme_parameters = args["deme_parameters"];
  H_vec = rcpp_to_vector_int(args_deme_parameters["H"]);
  seed_infections = rcpp_to_vector_int(args_deme_parameters["seed_infections"]);
  M_vec = rcpp_to_vector_int(args_deme_parameters["M"]);
  n_demes = int(H_vec.size());
  
  // demography
  demography = rcpp_to_vector_double(args["demography"]);
  n_demography = int(demography.size());
  
  // migration
  // (TODO)
  
  // run parameters
  Rcpp::List args_run_parameters = args["run_parameters"];
  max_time = rcpp_to_int(args_run_parameters["max_time"]);
  output_counts = rcpp_to_bool(args_run_parameters["output_counts"]);
  output_age_times = rcpp_to_vector_int(args_run_parameters["output_age_times"]);
  n_output_age_times = int(output_age_times.size());
  output_infection_history = rcpp_to_bool(args_run_parameters["output_infection_history"]);
  
  // misc parameters
  prob_v_death = 1 - exp(-mu);  // daily probability of mosquito death
  prob_h_recovery = 1 - exp(-r);  // daily probability of human recovery
  
}
#endif

//------------------------------------------------
// print summary
void Parameters::print_summary() {
  
  // epidemiological parameters
  print("-- epidemiological parameters --");
  print("a: ", a);
  print("mu: ", mu);
  print("u: ", u);
  print("v: ", v);
  print("g: ", g);
  print("r: ", r);
  print("b: ");
  print_vector(b);
  print("c: ", c);
  print("max_innoculations: ", max_innoculations);
  print("");
  
  // deme parameters
  print("-- deme parameters --");
  print("H_vec: ");
  print_vector(H_vec);
  print("seed_infections: ");
  print_vector(seed_infections);
  print("M_vec: ");
  print_vector(M_vec);
  print("n_demes: ", n_demes);
  print("");
  
  // demography
  print("-- demography --");
  print("demography: ");
  print_vector(demography);
  print("n_demography: ", n_demography);
  print("");
  
  // migration
  print("-- migration --");
  print("migration: (TODO)");
  print("");
  
  // run parameters
  print("-- run parameters --");
  print("max_time: ", max_time);
  print("output_counts: ", output_counts);
  print("output_age_times: ");
  print_vector(output_age_times);
  print("n_output_age_times: ", n_output_age_times);
  print("output_infection_history: ", output_infection_history);
  print("");
  
  // misc parameters
  print("-- misc parameters --");
  print("prob_v_death: ", prob_v_death);
  print("prob_h_recovery: ", prob_h_recovery);
  print("");
}
