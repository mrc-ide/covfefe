
#pragma once

#include "base_model.Parameters.h"
#include "base_model.Population.h"
#include "base_model.Scheduler.h"
#include "base_model.Host.h"
#include "base_model.Mosquito.h"
#include <Rcpp.h>

//------------------------------------------------
// single deme of individual-based simulation model
class Deme : public Population, public Scheduler {
  
public:
  
  // PUBLIC OBJECTS
  
  //std::vector<int> demog_draws;
  
  // index of this deme
  int this_deme;
  
  // counts of host types
  int H;
  int Sh;
  int Ih;
  
  // indices of hosts in this deme
  std::vector<int> hosts_uninfective;
  std::vector<int> hosts_infective;
  
  // total beta values (host infectiousness)
  double beta_uninfective;
  double beta_infective;
  double beta_total;
  
  // counts of mosquito types
  int M;
  int Sv;
  
  // population of mosquitoes (see Indiv_deme.cpp for details)
  std::vector<std::vector<Mosquito>> Ev_mosq;
  std::vector<Mosquito> Iv_mosq;
  std::vector<int> Ev_death;
  int ringtime;
  
  // objects for storing daily counts
  std::vector<int> Sh_store;
  std::vector<int> Ih_store;
  std::vector<double> EIR_store;
  
  // object for storing the number of hosts that carry each possible 
  // number of innoculations
  std::vector<std::vector<int>> innoculations_store;
  
  // infection history list for storing complete history of events
  /*
  std::vector<int> history_migration;
  std::vector<int> history_infection;
  std::vector<int> history_bloodstage;
  std::vector<int> history_recover;
  std::vector<std::vector<std::vector<int>>> history_full;
  */
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Deme() {};
  Deme(int this_deme);
  
  // methods
  void get_beta();
  void init(std::vector<int> &hosts_uninfective, int Ih);
  void human_infection(int host_index, Mosquito &mosq, int t);
  void step_forward(int t);
  
};