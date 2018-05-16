
#pragma once

#include <Rcpp.h>

//------------------------------------------------
// class defining host
class indiv_host {
  
public:
  
  // PUBLIC OBJECTS
  
  // individual properties
  int deme;
  int total_infections;
  int n_latent;
  int n_bloodstage;
  int n_infective;
  
  // PUBLIC FUNCTIONS
  
  // constructors
  indiv_host();
  
};
