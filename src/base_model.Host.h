
#pragma once

#include "misc.h"

#ifdef RCPP_ACTIVE
#include <Rcpp.h>
#endif

#include <list>

#ifdef OLD_VERSION


//------------------------------------------------
// enumerate possible host status
enum Status {Inactive, Latent, Acute, Chronic, Treated};

//------------------------------------------------
// class defining host
class Host {
  
public:
  
  // PUBLIC OBJECTS
  
  // unique ID and record of current deme
  int ID;
  int deme;
  
  // host properties
  int b_index;
  int prob_acute_index;
  int duration_acute_index;
  int duration_chronic_index;
  int birth_day;
  int death_day;
  
  // innoculation objects
  std::vector<bool> innoc_active;
  std::vector<Status> innoc_status;
  std::vector<int> innoc_status_update_time;
  std::vector<bool> innoc_infective;
  std::vector<int> innoc_infective_start_time;
  std::vector<int> innoc_infective_stop_time;
  
  // innoculation counts
  int cumulative_n_innoculations;
  int n_latent;
  int n_acute;
  int n_chronic;
  int n_infective;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Host() {};
  Host(int max_innoculations);
  
  // methods
  void reset(int ID, int deme, int birth_day, int death_day);
  int get_innoculation_slot();
  void new_innoculation(int slot, int t);
  void summary();
  int get_n_innoculations();
  int get_n_asexual();
  
};

#else

//------------------------------------------------
// enumerate possible host status
enum Status {Null, Inactive, Latent, Acute, Chronic, Treated};

//------------------------------------------------
// class defining host
class Host {
  
public:
  
  // PUBLIC OBJECTS
  
  // unique ID and record of current deme
  int ID;
  int deme;
  
  // indices relating to global distributions
  int b_index;
  int prob_acute_index;
  int duration_acute_index;
  int duration_chronic_index;
  int infectivity_acute_index;
  int infectivity_chronic_index;
  
  // dates of birth and death
  int birth_day;
  int death_day;
  
  // innoculation objects
  std::vector<Status> innoc_status;
  std::vector<Status> innoc_next_status;
  std::vector<int> innoc_status_update_time;
  std::vector<Status> innoc_infective_status;
  std::vector<Status> innoc_next_infective_status;
  std::vector<int> innoc_infective_status_update_time;
  
  // innoculation counts
  int cumulative_n_innoculations;
  int n_latent;
  int n_acute;
  int n_chronic;
  int n_infective_acute;
  int n_infective_chronic;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Host() {};
  Host(int max_innoculations);
  
  // methods
  void reset(int ID, int deme, int birth_day, int death_day);
  int get_n_innoculations();
  int get_n_asexual();
  int get_n_infective();
  void summary();
  
};

#endif
