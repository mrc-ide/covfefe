
#pragma once

#include "base_model.Parameters.h"
#include <Rcpp.h>

//------------------------------------------------
// class defining scheduler
class Scheduler : public Parameters {
  
public:
  
  // PUBLIC OBJECTS
  
  // objects for scheduling future events
  static std::vector<std::vector<std::pair<int, int>>> schedule_bloodstage;
  static std::vector<std::vector<std::pair<int, int>>> schedule_infective;
  static std::vector<std::vector<std::tuple<int, int, bool>>> schedule_recover;
  //std::vector<std::vector<std::vector<int>>> schedule_recover;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  Scheduler();
  
  // methods
  void init();
  void schedule_future_events(int host_index, int infection_ID, int t);
  
};
