
#pragma once

#include "misc.h"
#include "base_model.Parameters.h"
#include "base_model.Population.h"
#include "base_model.Deme.h"

#ifdef RCPP_ACTIVE
#include <Rcpp.h>
#endif

//------------------------------------------------
// class defining individual-based simulation model
class Dispatcher : public Parameters {
  
public:
  
  // PUBLIC OBJECTS
  
  // population of human hosts over all demes
  int H;
  Population population;
  
  // vector of demes
  std::vector<Deme> demes;
  
  // initialise objects for implementing migration
  //std::vector<std::vector<std::vector<int>>> mig_noninf_hosts;
  //std::vector<std::vector<std::vector<int>>> mig_inf_hosts;
  
  // objects for storing daily counts
  std::vector<std::vector<int>> Sh_store;
  std::vector<std::vector<int>> Eh_store;
  std::vector<std::vector<int>> Ih_store;
  //std::vector<std::vector<double>> EIR_store;
  
  // object for storing the number of hosts that carry each possible
  // number of innoculations
  //std::vector<std::vector<int>> innoculations_store;
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Dispatcher();
  
  // methods
  void simulate();
  //void migrate();
  
};
