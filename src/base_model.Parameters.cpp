
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
vector<double> Parameters::prob_infection;
int Parameters::n_prob_infection;
vector<double> Parameters::prob_acute;
int Parameters::n_prob_acute;
double Parameters::prob_AC;
vector<vector<double>> Parameters::duration_acute;
int Parameters::n_duration_acute;
vector<vector<double>> Parameters::duration_chronic;
int Parameters::n_duration_chronic;
vector<double> Parameters::infectivity_acute;
int Parameters::n_infectivity_acute;
vector<double> Parameters::infectivity_chronic;
int Parameters::n_infectivity_chronic;
int Parameters::max_innoculations;

// deme parameters
vector<int> Parameters::H_vec;
vector<int> Parameters::seed_infections;
vector<int> Parameters::M_vec;
int Parameters::n_demes;

// demography
vector<double> Parameters::life_table;
vector<double> Parameters::age_death;
vector<double> Parameters::age_stable;
int Parameters::n_age;

// migration
// TODO

// run parameters
int Parameters::max_time;
bool Parameters::output_daily_counts;
bool Parameters::output_age_distributions;
vector<int> Parameters::output_age_times;
int Parameters::n_output_age_times;
bool Parameters::output_infection_history;
string Parameters::filepath_migration;

// misc parameters
double Parameters::prob_v_death;

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
  prob_infection = vector<double>(1, 1.0);
  n_prob_infection = int(prob_infection.size());
  prob_acute = vector<double>(1, 1.0);
  n_prob_acute = int(prob_acute.size());
  prob_AC = 0.5;
  duration_acute = vector<vector<double>>(1, vector<double>(0.1, 10));
  n_duration_acute = int(duration_clinical.size());
  duration_chronic = vector<vector<double>>(1, vector<double>(0.1, 10));
  n_duration_chronic = int(duration_chronic.size());
  infectivity_acute = vector<double>(1, 1.0);
  n_infectivity_acute = int(infectivity_acute.size());
  infectivity_chronic = vector<double>(1, 1.0);
  n_infectivity_chronic = int(infectivity_chronic.size());
  max_innoculations = 5;
  
  // deme parameters
  int n_demes0 = 1;
  H_vec = vector<int>(n_demes0, 1000);
  seed_infections = vector<int>(n_demes0, 100);
  M_vec = vector<int>(n_demes0, 1000);
  n_demes = n_demes0;
  
  // demography
  life_table = {1};
  age_death = {1};
  age_stable = {1};
  n_age = int(age_stable.size());
  
  // migration
  // (TODO)
  
  // run parameters
  max_time = 365*5;
  output_daily_counts = true;
  output_age_distributions = true;
  output_age_times = {max_time};
  n_output_age_times = int(output_age_times.size());
  output_infection_history = false;
  
  // misc parameters
  prob_v_death = 1 - exp(-mu);  // daily probability of mosquito death
  
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
  prob_infection = rcpp_to_vector_double(args_epi_parameters["prob_infection"]);
  n_prob_infection = int(prob_infection.size());
  prob_acute = rcpp_to_vector_double(args_epi_parameters["prob_acute"]);
  n_prob_acute = int(prob_acute.size());
  prob_AC = rcpp_to_double(args_epi_parameters["prob_AC"]);
  duration_acute = rcpp_to_matrix_double(args_epi_parameters["duration_acute"]);
  n_duration_acute = int(duration_acute.size());
  duration_chronic = rcpp_to_matrix_double(args_epi_parameters["duration_chronic"]);
  n_duration_chronic = int(duration_chronic.size());
  infectivity_acute = rcpp_to_vector_double(args_epi_parameters["infectivity_acute"]);
  n_infectivity_acute = int(infectivity_acute.size());
  infectivity_chronic = rcpp_to_vector_double(args_epi_parameters["infectivity_chronic"]);
  n_infectivity_chronic = int(infectivity_chronic.size());
  max_innoculations = rcpp_to_int(args_epi_parameters["max_innoculations"]);
  
  // deme parameters
  Rcpp::List args_deme_parameters = args["deme_parameters"];
  H_vec = rcpp_to_vector_int(args_deme_parameters["H"]);
  seed_infections = rcpp_to_vector_int(args_deme_parameters["seed_infections"]);
  M_vec = rcpp_to_vector_int(args_deme_parameters["M"]);
  n_demes = int(H_vec.size());
  
  // demography
  Rcpp::List args_demography = args["demography"];
  life_table = rcpp_to_vector_double(args_demography["life_table"]);
  age_death = rcpp_to_vector_double(args_demography["age_death"]);
  age_stable = rcpp_to_vector_double(args_demography["age_stable"]);
  n_age = int(age_stable.size());
  
  // migration
  // (TODO)
  
  // run parameters
  Rcpp::List args_run_parameters = args["run_parameters"];
  max_time = rcpp_to_int(args_run_parameters["max_time"]);
  output_daily_counts = rcpp_to_bool(args_run_parameters["output_daily_counts"]);
  output_age_distributions = rcpp_to_bool(args_run_parameters["output_age_distributions"]);
  output_age_times = rcpp_to_vector_int(args_run_parameters["output_age_times"]);
  n_output_age_times = int(output_age_times.size());
  output_infection_history = rcpp_to_bool(args_run_parameters["output_infection_history"]);
  filepath_migration = rcpp_to_string(args_run_parameters["filepath_migration"]);
  print(filepath_migration);
  
  // misc parameters
  prob_v_death = 1 - exp(-mu);  // daily probability of mosquito death
  
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
  print("prob_infection: ");
  print_vector(prob_infection);
  print("prob_acute: ");
  print_vector(prob_acute);
  print("prob_AC: ", prob_AC);
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
  print("life_table: ");
  print_vector(life_table);
  print("age_death: ");
  print_vector(age_death);
  print("age_stable: ");
  print_vector(age_stable);
  print("n_age: ", n_age);
  print("");
  
  // migration
  print("-- migration --");
  print("migration: (TODO)");
  print("");
  
  // run parameters
  print("-- run parameters --");
  print("max_time: ", max_time);
  print("output_daily_counts: ", output_daily_counts);
  print("output_age_distributions: ", output_age_distributions);
  print("output_age_times: ");
  print_vector(output_age_times);
  print("n_output_age_times: ", n_output_age_times);
  print("output_infection_history: ", output_infection_history);
  print("");
  
  // misc parameters
  print("-- misc parameters --");
  print("prob_v_death: ", prob_v_death);
  print("");
}

