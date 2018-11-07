
#pragma once

#include <Rcpp.h>

//------------------------------------------------
// class defining mosquito
class Mosquito {
  
public:
  
  // PUBLIC OBJECTS
  
  // properties
  int host_ID;
  int host_infections;
  int infection_time;
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Mosquito();
  Mosquito(int host_ID, int host_infections, int infection_time);
  
};
