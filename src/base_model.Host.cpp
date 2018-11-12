
#include "base_model.Host.h"
#include "probability.h"

using namespace std;

//------------------------------------------------
// constructor
Host::Host(int max_innoculations) {
  
  // initialise innoculation objects
  innoc_active = vector<bool>(max_innoculations);
  innoc_status = vector<Status>(max_innoculations);
  innoc_status_update_time = vector<int>(max_innoculations);
  innoc_infective = vector<bool>(max_innoculations);
  innoc_infective_start_time = vector<int>(max_innoculations);
  innoc_infective_stop_time = vector<int>(max_innoculations);
  
}

//------------------------------------------------
// reset host parameters
void Host::reset(int ID, int deme, int birth_day, int death_day) {
  
  // unique ID and record of current deme
  this->ID = ID;
  this->deme = deme;
  
  // host properties
  b_index = 0;
  prob_acute_index = 0;
  duration_acute_index = 0;
  duration_chronic_index = 0;
  this->birth_day = birth_day;
  this->death_day = death_day;
  
  // innoculation objects
  fill(innoc_active.begin(), innoc_active.end(), false);
  fill(innoc_status.begin(), innoc_status.end(), Inactive);
  fill(innoc_status_update_time.begin(), innoc_status_update_time.end(), 0);
  fill(innoc_infective.begin(), innoc_infective.end(), false);
  fill(innoc_infective_start_time.begin(), innoc_infective_start_time.end(), 0);
  fill(innoc_infective_stop_time.begin(), innoc_infective_stop_time.end(), 0);
  
  // innoculation counts
  cumulative_n_innoculations = 0;
  n_latent = 0;
  n_acute = 0;
  n_chronic = 0;
  n_infective = 0;
}

//------------------------------------------------
// get next free innoculation slot. Return -1 if no free slot
int Host::get_innoculation_slot() {
  
  for (int i=0; i<int(innoc_active.size()); ++i) {
    if (!innoc_active[i]) {
      return i;
    }
  }
  return -1;
}

//------------------------------------------------
// new innoculation, scheduled to transition to blood-stage at time t
void Host::new_innoculation(int slot, int t) {
  
  // update latent count
  n_latent++;
  
  // add new innoculation
  innoc_active[slot] = true;
  innoc_status[slot] = Latent;
  innoc_status_update_time[slot] = t;
  innoc_infective[slot] = false;
  
}

//------------------------------------------------
// print summary
void Host::summary() {
  print("ID: ", ID);
  print("deme: ", deme);
  print("b_index: ", b_index);
  print("prob_acute_index: ", prob_acute_index);
  print("duration_acute_index: ", duration_acute_index);
  print("birth_day: ", birth_day);
  print("death_day: ", death_day);
  print("");
  
  print("cumulative_n_innoculations: ", cumulative_n_innoculations);
  print("n_latent: ", n_latent);
  print("n_acute: ", n_acute);
  print("n_chronic: ", n_chronic);
  print("n_infective: ", n_infective);
  print("");
  
  print_vector(innoc_active);
  print_vector(innoc_status);
  print_vector(innoc_status_update_time);
  print_vector(innoc_infective);
  print_vector(innoc_infective_start_time);
  print_vector(innoc_infective_stop_time);
  print("");
}

//------------------------------------------------
// get total number of innoculations
int Host::get_n_innoculations() {
  return n_latent + n_acute + n_chronic + n_infective;
}

//------------------------------------------------
// get total number of asexual innoculations
int Host::get_n_asexual() {
  return n_latent + n_acute + n_chronic;
}
