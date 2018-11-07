
#pragma once

#include "base_model.Parameters.h"
#include <Rcpp.h>

//------------------------------------------------
// class defining host
class Host : public Parameters {
  
public:
  
  // PUBLIC OBJECTS
  
  // unique ID and record of current deme
  int ID;
  int deme;
  
  // infection properties
  int total_infections;
  int n_latent;
  int n_bloodstage;
  int n_infective;
  double beta;
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Host();
  
  // methods
  void new_infection();
  void transition_bloodstage();
  void transition_infective();
  int get_n_innoculations();
  
};
