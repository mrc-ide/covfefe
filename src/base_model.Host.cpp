
#include "base_model.Host.h"
#include "misc.h"

using namespace std;

//------------------------------------------------
// default constructor for host class
Host::Host() {
  
  // unique ID and record of current deme
  ID = 0;
  deme = 0;
  
  // infection properties
  total_infections = 0;
  n_latent = 0;
  n_bloodstage = 0;
  n_infective = 0;
  beta = b[0];
}

//------------------------------------------------
// new infection
void Host::new_infection() {
  total_infections++;
  n_latent++;
  if (total_infections < b.size()) {
    beta = b[total_infections];
  } else {
    beta = b[b.size()-1];
  }
}

//------------------------------------------------
// transition to blood-stage
void Host::transition_bloodstage() {
  n_latent--;
  n_bloodstage++;
}

//------------------------------------------------
// transition to infective stage
void Host::transition_infective() {
  n_bloodstage--;
  n_infective++;
}

//------------------------------------------------
// get current number of innoculations, including all latent and blood-stage
int Host::get_n_innoculations() {
  return n_latent + n_bloodstage + n_infective;
}
