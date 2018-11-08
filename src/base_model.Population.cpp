
#include "base_model.Population.h"
#include "probability.h"

using namespace std;

//------------------------------------------------
// declare static member variables

int Population::next_host_ID;

vector<Host> Population::hosts;

//------------------------------------------------
// constructor
Population::Population(int H) {
  
  // create vector of hosts
  hosts = vector<Host>(H);
  
  // initialise hosts
  next_host_ID = 0;
  for (int i=0; i<H; i++) {
    
    // draw life duration from demography distribution
    int life_years = sample1(demography);
    int life_days = life_years*365 + sample2(1, 365);
    
    // draw age at time 0
    int age_days = sample2(0, life_days);
    
    // initialise host
    hosts[i].reset(next_host_ID++, -age_days, life_days-age_days);
  }
  
}

//------------------------------------------------
// check for host death
void Population::enact_death(int this_host, int t) {
  
  // reset host on death
  if (hosts[this_host].death_day == t) {
    
    // draw life duration from demography distribution
    int life_years = sample1(demography);
    int life_days = life_years*365 + sample2(1, 365);
    hosts[this_host].reset(next_host_ID++, t, t + life_days);
  }
  
}
