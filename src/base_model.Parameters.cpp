
#include "base_model.Parameters.h"

#include <cmath>

using namespace std;

//------------------------------------------------
// declare static member variables

int Parameters::max_time;
double Parameters::a;
double Parameters::mu;
int Parameters::u;
int Parameters::v;
int Parameters::g;
double Parameters::r;
vector<double> Parameters::b;
double Parameters::c;
vector<int> Parameters::Eh_vec;
vector<int> Parameters::H_vec;
vector<int> Parameters::M_vec;
int Parameters::max_innoculations;
vector<vector<double>> Parameters::delta_mig;
vector<double> Parameters::demog;
bool Parameters::output_counts;
bool Parameters::output_innoculations;
bool Parameters::output_infection_history;
int Parameters::n_demes;

double Parameters::prob_v_death;
double Parameters::prob_h_recovery;

//------------------------------------------------
// constructor
#ifdef RCPP_ACTIVE
Parameters::Parameters() {}
#else
Parameters::Parameters() {
  
  // extract model parameters
  max_time = 365*5;
  a = 0.3;
  mu = -log(0.9);
  u = 10;
  v = 10;
  g = 10;
  r = 1.0/20;
  b = vector<double>(1, 1.0);
  c = 1.0;
  int n_demes0 = 1;
  Eh_vec = vector<int>(n_demes0, 100);
  H_vec = vector<int>(n_demes0, 1000);
  M_vec = vector<int>(n_demes0, 1000);
  max_innoculations = 5;
  //delta_mig = rcpp_to_matrix_double(args["delta_mig"]);
  demog = {100};
  output_counts = true;
  output_innoculations = true;
  output_infection_history = false;
  n_demes = int(H_vec.size());
  
  // daily probability of mosquito death
  prob_v_death = 1 - exp(-mu);
  prob_h_recovery = 1 - exp(-r);
}
#endif

//------------------------------------------------
// constructor
#ifdef RCPP_ACTIVE
Parameters::Parameters(const Rcpp::List &args) {
  
  // extract model parameters
  max_time = rcpp_to_int(args["max_time"]);
  a = rcpp_to_double(args["a"]);
  mu = rcpp_to_double(args["mu"]);
  u = rcpp_to_int(args["u"]);
  v = rcpp_to_int(args["v"]);
  g = rcpp_to_int(args["g"]);
  r = rcpp_to_double(args["r"]);
  b = rcpp_to_vector_double(args["b"]);
  c = rcpp_to_double(args["c"]);
  Eh_vec = rcpp_to_vector_int(args["Eh_init"]);
  H_vec = rcpp_to_vector_int(args["H"]);
  M_vec = rcpp_to_vector_int(args["M"]);
  max_innoculations = rcpp_to_int(args["max_innoculations"]);
  delta_mig = rcpp_to_matrix_double(args["delta_mig"]);
  demog = rcpp_to_vector_double(args["demog"]);
  output_counts = rcpp_to_bool(args["output_counts"]);
  output_innoculations = rcpp_to_bool(args["output_innoculations"]);
  output_infection_history = rcpp_to_bool(args["output_infection_history"]);
  n_demes = H_vec.size();
  
  // daily probability of mosquito death
  prob_v_death = 1 - exp(-mu);
  prob_h_recovery = 1 - exp(-r);
}
#endif

