
#include "base_model.Population.h"

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
  
  // assign unique IDs
  next_host_ID = 0;
  for (int i=0; i<H; i++) {
    hosts[i].ID = next_host_ID++;
  }
  
}
