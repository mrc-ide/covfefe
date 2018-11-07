
#include "base_model.Scheduler.h"
#include "probability.h"
#include "misc.h"

using namespace std;

//------------------------------------------------
// declare static member variables

vector<vector<pair<int, int>>> Scheduler::schedule_bloodstage;
vector<vector<pair<int, int>>> Scheduler::schedule_infective;
vector<vector<tuple<int, int, bool>>> Scheduler::schedule_recover;

//------------------------------------------------
// constructor
Scheduler::Scheduler() {
  
  // schedule move to blood-stage using pair of values:
  // first = index of host, second = infection ID
  schedule_bloodstage = vector<vector<pair<int, int>>>(max_time);
  
  // schedule move to infective stage using pair of values:
  // first = index of host, second = infection ID
  schedule_infective = vector<vector<pair<int, int>>>(max_time);
  
  // schedule move to infective stage using tuple of three values:
  // first = index of host, second = infection ID, third = whether infection is in infective stage
  schedule_recover = vector<vector<tuple<int, int, bool>>>(max_time);
  
}

//------------------------------------------------
// schedule future events
void Scheduler::schedule_future_events(int host_index, int infection_ID, int t) {
  
  // schedule move to blood-stage
  if ((t+u) < max_time) {
    schedule_bloodstage[t+u].emplace_back(host_index, infection_ID);
  }
  
  // schedule becoming infective and/or recovery
  int dur = rgeom1(prob_h_recovery);
  if (dur > g) {  // if become infective prior to recovery
    
    if ((t+u+dur) < max_time) { // if both becoming infective and recovery occur within max_time
      schedule_infective[t+u+g].emplace_back(host_index, infection_ID);
      schedule_recover[t+u+dur].emplace_back(host_index, infection_ID, true);
    } else if ((t+u+g) < max_time) {  // if becoming infective but not recovery occurs within max_time
      schedule_infective[t+u+g].emplace_back(host_index, infection_ID);
    }
    
  } else {  // if recover prior to becoming infective
    if ((t+u+dur)<max_time) { // if recovery occurs within max_time
      schedule_recover[t+u+dur].emplace_back(host_index, infection_ID, false);
    }
  }
}
