
#pragma once

#include "misc.h"
#include "base_model.Host.h"

#ifdef RCPP_ACTIVE
#include <Rcpp.h>
#endif

//------------------------------------------------
// class defining host population
class Population {
  
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
  
};
