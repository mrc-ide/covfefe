
#pragma once

#include "base_model.Parameters.h"
#include "base_model.Population.h"
#include "base_model.Scheduler.h"
#include "base_model.Deme.h"
#include <Rcpp.h>

//------------------------------------------------
// class defining individual-based simulation model
class Dispatcher : public Parameters {
  
public:
  
  // PUBLIC OBJECTS
  
  // population of human hosts over all demes
  Population population;
  
  // event scheduler
  Scheduler scheduler;
  
  // vector of demes
  std::vector<Deme> demes;
  
  // initialise objects for implementing migration
  //std::vector<std::vector<std::vector<int>>> mig_noninf_hosts;
  //std::vector<std::vector<std::vector<int>>> mig_inf_hosts;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Dispatcher();
  
  // methods
  void simulate();
  //void migrate();
  
};