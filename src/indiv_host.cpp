
#include "indiv_host.h"
#include "misc.h"

using namespace std;

//------------------------------------------------
// default constructor for host class
indiv_host::indiv_host() {
  ID = 0;
  n_infections = 0;
  n_status3 = 0;
}

//------------------------------------------------
// new infection
void indiv_host::new_infection(int infection_time, int recovery_time) {
  
  //  Each infection has elements:
  //  1. status (1=liver stage, 2=blood stage, 3=infectious stage)
  //  2. infection time
  //  3. recovery time
  n_infections++;
  infections.push_back(vector<int>{1, infection_time, recovery_time});
  
}

//------------------------------------------------
// step time forward
void indiv_host::step_forward(int t, int u, int g, unordered_set<int> &h_infectious) {
  
  // return if no infection
  if (n_infections == 0) {
    return;
  }
  
  // loop through infections. Establish whether event needs to happen at this
  // time point
  int i = 0;
  while (i<n_infections) {
    if (t == infections[i][2]) {  // recovery
      if (infections[i][0]==3) {  // decrease n_status3 if this infection in status 3
        n_status3--;
        if (n_status3 == 0) { // change contribution to infectiousness if necessary
          h_infectious.erase(ID);
        }
      }
      infections.erase(infections.begin() + i); // drop this infection
      n_infections--;
      continue;
    } else if (t == (infections[i][1] + u)) { // move to blood stage
      infections[i][0] = 2;
    } else if (t == (infections[i][1] + u + g)) { // move to infectious stage
      infections[i][0] = 3;
      if (n_status3 == 0) { // change contribution to infectiousness if necessary
        h_infectious.insert(ID);
      }
      n_status3++;
    }
    i++;
  }
  
}

