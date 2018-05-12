
#pragma once

#include <Rcpp.h>

//------------------------------------------------
// class defining host
class indiv_host {
  
public:
  
  // PUBLIC OBJECTS
  
  // individual properties
  int ID;
  int n_infections;
  int n_status3;
  
  // vector of infections. Each infection has elements:
  //  1. status (1=liver stage, 2=blood stage, 3=blood stage + gametocytes)
  //  2. infection time
  //  3. recovery time
  std::vector<std::vector<int>> infections;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  indiv_host();
  
  // other functions
  void new_infection(int infection_time, int recovery_time);
  void step_forward(int t, int u, int g, std::unordered_set<int> &h_infectious);
};
