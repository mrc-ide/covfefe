
#include "base_model.Host.h"
#include "probability.h"

using namespace std;


//------------------------------------------------
// constructor
Host::Host(int max_innoculations) {
  
  // initialise innoculation objects
  innoc_active = vector<bool>(max_innoculations);
  innoc_status = vector<Status>(max_innoculations);
  innoc_status_prev_update_time = vector<int>(max_innoculations);
  innoc_status_next_update_time = vector<int>(max_innoculations);
  innoc_infective_start_acute = vector<int>(max_innoculations);
  innoc_infective_start_chronic = vector<int>(max_innoculations);
  innoc_infective_acute_chronic = vector<int>(max_innoculations);
  innoc_infective_stop_acute = vector<int>(max_innoculations);
  innoc_infective_stop_chronic = vector<int>(max_innoculations);
  
}

//------------------------------------------------
// reset host parameters
void Host::reset(int ID, int deme, int birth_day, int death_day) {
  
  // unique ID and record of current deme
  this->ID = ID;
  this->deme = deme;
  
  // indices relating to global distributions
  prob_infection_index = 0;
  prob_acute_index = 0;
  duration_acute_index = 0;
  duration_chronic_index = 0;
  infectivity_acute_index = 0;
  infectivity_chronic_index = 0;
  
  // dates of birth and death
  this->birth_day = birth_day;
  this->death_day = death_day;
  
  // innoculation objects
  fill(innoc_active.begin(), innoc_active.end(), false);
  fill(innoc_status.begin(), innoc_status.end(), Inactive);
  fill(innoc_status_prev_update_time.begin(), innoc_status_prev_update_time.end(), 0);
  fill(innoc_status_next_update_time.begin(), innoc_status_next_update_time.end(), 0);
  fill(innoc_infective_start_acute.begin(), innoc_infective_start_acute.end(), 0);
  fill(innoc_infective_start_chronic.begin(), innoc_infective_start_chronic.end(), 0);
  fill(innoc_infective_acute_chronic.begin(), innoc_infective_acute_chronic.end(), 0);
  fill(innoc_infective_stop_acute.begin(), innoc_infective_stop_acute.end(), 0);
  fill(innoc_infective_stop_chronic.begin(), innoc_infective_stop_chronic.end(), 0);
  
  // innoculation counts
  cumulative_n_innoculations = 0;
  n_latent = 0;
  n_acute = 0;
  n_chronic = 0;
  n_infective_acute = 0;
  n_infective_chronic = 0;
}

//------------------------------------------------
// get total number of innoculations
int Host::get_n_innoculations() {
  return n_latent + n_acute + n_chronic + n_infective_acute + n_infective_chronic;
}

//------------------------------------------------
// get total number of asexual innoculations
int Host::get_n_asexual() {
  return n_latent + n_acute + n_chronic;
}

//------------------------------------------------
// get total number of infective innoculations
int Host::get_n_infective() {
  return n_infective_acute + n_infective_chronic;
}

//------------------------------------------------
// print summary
void Host::summary() {
  print("ID: ", ID);
  print("deme: ", deme);
  print("prob_infection_index: ", prob_infection_index);
  print("prob_acute_index: ", prob_acute_index);
  print("duration_acute_index: ", duration_acute_index);
  print("birth_day: ", birth_day);
  print("death_day: ", death_day);
  print("");
  
  print("cumulative_n_innoculations: ", cumulative_n_innoculations);
  print("n_latent: ", n_latent);
  print("n_acute: ", n_acute);
  print("n_chronic: ", n_chronic);
  print("n_infective_acute: ", n_infective_acute);
  print("n_infective_chronic: ", n_infective_chronic);
  print("");
  
  print_vector(innoc_active);
  print_vector(innoc_status);
  print_vector(innoc_status_prev_update_time);
  print_vector(innoc_status_next_update_time);
  print_vector(innoc_infective_start_acute);
  print_vector(innoc_infective_start_chronic);
  print_vector(innoc_infective_acute_chronic);
  print_vector(innoc_infective_stop_acute);
  print_vector(innoc_infective_stop_chronic);
  print("");
}
