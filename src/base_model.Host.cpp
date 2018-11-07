
#include "base_model.Host.h"
#include "probability.h"

using namespace std;

//------------------------------------------------
// default constructor for host class
Host::Host() {
  
  // unique ID and record of current deme
  ID = 0;
  deme = 0;
  
  // host properties
  beta = b[0];
  next_event_time = max_time+1;
  
  // innoculation objects
  innoc_active = vector<bool>(max_innoculations, false);
  innoc_ID = vector<int>(max_innoculations);
  innoc_status = vector<int>(max_innoculations);
  innoc_status_update_time = vector<int>(max_innoculations);
  innoc_infective = vector<bool>(max_innoculations, false);
  innoc_infective_start_time = vector<int>(max_innoculations);
  innoc_infective_stop_time = vector<int>(max_innoculations);
  
  // innoculation counts
  cumulative_n_innoculations = 0;
  n_innoculations = 0;
  n_latent = 0;
  n_bloodstage = 0;
  n_infective = 0;
  n_asexual = 0;
}

//------------------------------------------------
// new innoculation
void Host::new_innoculation(int t) {
  
  // abort if reached max innoculations
  if (n_innoculations == max_innoculations) {
    return;
  }
  
  // find next free innoculation slot
  int next_innoculation = 0;
  for (int i=0; i<max_innoculations; ++i) {
    if (!innoc_active[i]) {
      next_innoculation = i;
      break;
    }
  }
  
  // add new innoculation
  innoc_active[next_innoculation] = true;
  innoc_ID[next_innoculation] = cumulative_n_innoculations;
  innoc_status[next_innoculation] = 0;
  innoc_status_update_time[next_innoculation] = t+u;
  innoc_infective[next_innoculation] = false;
  innoc_infective_start_time[next_innoculation] = t+u+g;
  innoc_infective_stop_time[next_innoculation] = max_time+1;
  
  // update next event time
  if (t+u < next_event_time) {
    next_event_time = t+u;
  }
  
  // update counts
  cumulative_n_innoculations++;
  n_innoculations++;
  n_latent++;
  n_asexual++;
  
  // update infection probability
  if (cumulative_n_innoculations < b.size()) {
    beta = b[cumulative_n_innoculations];
  } else {
    beta = b[b.size()-1];
  }
}

//------------------------------------------------
// enact scheduled events
void Host::enact_events(int t, vector<int> &host_infective_vec, int this_host) {
  
  // enact events scheduled for time t
  if (t != next_event_time) {
    return;
  }
  
  //if (ID == 0) {
  //  print(t);
  //  summary();
  //}
  
  // store whether infective before updating
  bool infective_before = (n_infective > 0);
  
  // update all active innoculations
  for (int i=0; i<max_innoculations; ++i) {
    if (innoc_active[i]) {
      
      // if due status update
      if (innoc_status_update_time[i] == t) {
        switch (innoc_status[i]) {
          case 0:  // latent become blood-stage
            innoc_status_update_time[i] += rgeom1(r)+1;
            n_latent--;
            n_bloodstage++;
            break;
            
          case 1:  // blood-stage recover
            innoc_status_update_time[i] = max_time+1;
            innoc_infective_stop_time[i] = t+g;
            n_bloodstage--;
            n_asexual--;
            break;
            
          default:
            break;
        }
        innoc_status[i]++;
      }
      
      // if due infective status update
      if (innoc_infective_start_time[i] == t) {
        innoc_infective[i] = true;
        n_infective++;
      }
      if (innoc_infective_stop_time[i] == t) {
        innoc_active[i] = false;
        n_infective--;
        n_innoculations--;
      }
      
    }
  }
  
  // recalculate next_event_time
  next_event_time = max_time+1;
  for (int i=0; i<max_innoculations; ++i) {
    if (innoc_active[i]) {
      if (innoc_status_update_time[i] < next_event_time) {
        next_event_time = innoc_status_update_time[i];
      }
      if (innoc_infective[i]) {
        if (innoc_infective_stop_time[i] < next_event_time) {
          next_event_time = innoc_infective_stop_time[i];
        }
      } else {
        if (innoc_infective_start_time[i] < next_event_time) {
          next_event_time = innoc_infective_start_time[i];
        }
      }
    }
  }
  
  //if (ID == 0) {
  //  summary();
  //}
  
  // apply change in infective status
  bool infective_after = (n_infective > 0);
  if (!infective_before && infective_after) {  // newly infective
    host_infective_vec.push_back(this_host);
  }
  if (infective_before && !infective_after) {  // newly uninfective
    host_infective_vec.erase(remove(host_infective_vec.begin(), host_infective_vec.end(), this_host));
  }
  
}

//------------------------------------------------
// print summary
void Host::summary() {
  print("ID: ", ID);
  print("deme: ", deme);
  print("beta: ", beta);
  print("next_event_time: ", next_event_time);
  print("cumulative_n_innoculations: ", cumulative_n_innoculations);
  print("n_innoculations: ", n_innoculations);
  print("n_latent: ", n_latent);
  print("n_bloodstage: ", n_bloodstage);
  print("n_infective: ", n_infective);
  print("n_asexual: ", n_asexual);
  
  print_vector(innoc_active);
  print_vector(innoc_status);
  print_vector(innoc_status_update_time);
  print_vector(innoc_infective);
  print_vector(innoc_infective_start_time);
  print_vector(innoc_infective_stop_time);
  print("");
}
