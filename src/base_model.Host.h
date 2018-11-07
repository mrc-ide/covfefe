
#pragma once

#include "misc.h"
#include "base_model.Parameters.h"

#ifdef RCPP_ACTIVE
#include <Rcpp.h>
#endif

#include <list>

//------------------------------------------------
// class defining host
class Host : public Parameters {
  
public:
  
  // PUBLIC OBJECTS
  
  // unique ID and record of current deme
  int ID;
  int deme;
  
  // host properties
  double beta;
  int next_event_time;
  
  // list of innoculations
  //std::list<Innoculation> innoculations;
  
  // innoculation objects
  std::vector<bool> innoc_active;
  std::vector<int> innoc_ID;
  std::vector<int> innoc_status;
  std::vector<int> innoc_status_update_time;
  std::vector<bool> innoc_infective;
  std::vector<int> innoc_infective_start_time;
  std::vector<int> innoc_infective_stop_time;
  
  // innoculation counts
  int cumulative_n_innoculations;
  int n_innoculations;
  int n_latent;
  int n_bloodstage;
  int n_infective;
  int n_asexual;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Host();
  
  // methods
  void new_innoculation(int t);
  void enact_events(int t, std::vector<int> &host_infective_vec, int this_host);
  void summary();
  
};
