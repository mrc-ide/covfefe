
#pragma once

#include "misc.h"
#include "base_model.Parameters.h"
#include "base_model.Host.h"
#include "base_model.Mosquito.h"
#include "base_model.Sampler.h"

#ifdef RCPP_ACTIVE
#include <Rcpp.h>
#endif

#include <set>
#include <map>


//------------------------------------------------
// class defining individual-based simulation model
class Dispatcher : public Parameters {
  
public:
  
  // PUBLIC OBJECTS
  
  // objects for sampling from probability distributions
  Sampler sampler_age_stable;
  Sampler sampler_age_death;
  std::vector<Sampler> sampler_duration_acute;
  std::vector<Sampler> sampler_duration_chronic;
  
  // scheduler objects
  std::vector<std::set<int>> schedule_death;
  std::vector<std::vector<std::pair<int, int>>> schedule_status_update;
  std::vector<std::vector<std::pair<int, int>>> schedule_infective_start_acute;
  std::vector<std::vector<std::pair<int, int>>> schedule_infective_start_chronic;
  std::vector<std::vector<std::pair<int, int>>> schedule_infective_acute_chronic;
  std::vector<std::vector<std::pair<int, int>>> schedule_infective_stop_acute;
  std::vector<std::vector<std::pair<int, int>>> schedule_infective_stop_chronic;
  
  // counts of host types
  int H_total;
  std::vector<int> H;
  std::vector<int> Sh;
  std::vector<int> Lh;
  std::vector<int> Ah;
  std::vector<int> Ch;
  
  // other quantities to keep track of
  std::vector<double> EIR;
  
  // population of human hosts over all demes
  std::vector<Host> hosts;
  int next_host_ID;
  
  // store the integer index of hosts in each deme
  std::vector<std::vector<int>> host_vec;
  std::vector<std::vector<int>> host_infective_vec;
  
  // counts of mosquito types
  std::vector<int> Sv;
  std::vector<int> Lv;
  std::vector<int> Iv;
  
  // population of mosquitoes (see Indiv_deme.cpp for details)
  std::vector<std::vector<std::vector<Mosquito>>> Ev_mosq;
  std::vector<std::vector<Mosquito>> Iv_mosq;
  std::vector<std::vector<int>> Ev_death;
  int ringtime;
  
  // objects for storing daily counts etc.
  std::vector<std::vector<int>> H_store;
  std::vector<std::vector<int>> Sh_store;
  std::vector<std::vector<int>> Lh_store;
  std::vector<std::vector<int>> Ah_store;
  std::vector<std::vector<int>> Ch_store;
  std::vector<std::vector<int>> Sv_store;
  std::vector<std::vector<int>> Lv_store;
  std::vector<std::vector<int>> Iv_store;
  std::vector<std::vector<double>> EIR_store;
  
  // objects for storing age distributions
  std::vector<std::vector<std::vector<int>>> H_age_store;
  std::vector<std::vector<std::vector<int>>> Sh_age_store;
  std::vector<std::vector<std::vector<int>>> Lh_age_store;
  std::vector<std::vector<std::vector<int>>> Ah_age_store;
  std::vector<std::vector<std::vector<int>>> Ch_age_store;
  std::vector<std::vector<std::vector<double>>> inc_Lh_age_store;
  std::vector<std::vector<std::vector<double>>> inc_Ah_age_store;
  
  // misc
  double prob_h_infectious_bite;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Dispatcher();
  
  // methods
  void new_infection(int this_host, Mosquito &mosq, int t);
  void latent_acute(int this_host, int this_slot, int t);
  void latent_chronic(int this_host, int this_slot, int t);
  void acute_chronic(int this_host, int this_slot, int t);
  void acute_recover(int this_host, int this_slot, int t);
  void chronic_recover(int this_host, int this_slot, int t);
  void update_age_incidence(int this_host, int this_deme, int output_age_time_index, int t);
  void simulate();
  
};

