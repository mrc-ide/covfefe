
#pragma once

#include "misc.h"
#include "base_model.Parameters.h"
#include "base_model.Population.h"
#include "base_model.Host.h"
#include "base_model.Mosquito.h"

#ifdef RCPP_ACTIVE
#include <Rcpp.h>
#endif

//------------------------------------------------
// single deme of individual-based simulation model
class Deme : public Population, public Parameters {
  
public:
  
  // PUBLIC OBJECTS
  
  //std::vector<int> demog_draws;
  
  // index of this deme
  int this_deme;
  
  // number of hosts
  int H;
  
  // store the integer index of hosts in this deme
  std::vector<int> host_vec;
  std::vector<int> host_infective_vec;
  
  // counts of mosquito types
  int M;
  int Sv;
  
  // population of mosquitoes (see Indiv_deme.cpp for details)
  std::vector<std::vector<Mosquito>> Ev_mosq;
  std::vector<Mosquito> Iv_mosq;
  std::vector<int> Ev_death;
  int ringtime;
  
  // infection history list for storing complete history of events
  /*
  std::vector<int> history_migration;
  std::vector<int> history_infection;
  std::vector<int> history_bloodstage;
  std::vector<int> history_recover;
  std::vector<std::vector<std::vector<int>>> history_full;
  */
  
  // misc
  double EIR;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Deme() {};
  Deme(int this_deme);
  
  // methods
  void init(std::vector<int> &host_vec0, int Eh);
  void human_infection(int this_host, Mosquito &mosq, int t);
  void step_forward(int t);
  
};
