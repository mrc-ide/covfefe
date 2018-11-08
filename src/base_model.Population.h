
#pragma once

#include "misc.h"
#include "base_model.Host.h"
#include "base_model.Parameters.h"

#ifdef RCPP_ACTIVE
#include <Rcpp.h>
#endif

//------------------------------------------------
// class defining host population
class Population : public Parameters {
  
public:
  
  // PUBLIC OBJECTS
  
  // static counter of next host ID
  static int next_host_ID;
  
  // static vector of hosts
  static std::vector<Host> hosts;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Population() {};
  Population(int H);
  
  // methods
  void enact_death(int this_host, int t);
  
};
