
#include "base_model.Deme.h"
#include "probability.h"

using namespace std;

//------------------------------------------------
// constructor for Deme class
Deme::Deme(int this_deme) {
  
  // draw from demography distribution
  //Rcpp::Function R_draws = args["R_draws"];
  //demog_draws = rcpp_to_vector_int(R_draws(demog, sum(H)));
  
  // set index of this deme
  this->this_deme = this_deme;
  
  // initialise population of mosquitoes. An infective mosquito consists of
  // three values: 1) the ID of the host that infected it, 2) the number of
  // observable infections in that host at the time of biting, 3) the time the
  // mosquito became infected by that host. If the infective mosquito ever bites
  // a new host and passes on the infection, these values will be written to the
  // infection history record, providing a link between the original host and
  // the newly infected host.
  // 
  // Typically a large number of mosquitoes become infected each time step, but
  // most of them die before making it through the lag phase (the extrinsic
  // incubation period). To avoid creating a large number of infective
  // mosquitoes that die before they could pass on infection, we instead draw
  // the number of infective mosquitoes that die at each day of the lag phase,
  // and only mosquitoes that make it all the way through without dying are
  // instantiated. The mosquitoes that die in the lag phase must re-enter the
  // population of susceptible mosquitoes on their day of death (i.e. we assume
  // constant mosquito population size). Ev_death stores the number of deaths of
  // infective mosquitoes at each day of the lag phase going forward until day
  // v. These deaths will be added back into the susceptible pool on that day.
  // Ev_mosq stores the full details of infective mosquitoes that are due to
  // emerge from the lag phase at each day going forward. These infective
  // mosquitoes will be added to Iv_mosq on that day. Finally, both Ev_death and
  // Ev_mosq are implemented as ring-buffers, meaning they wrap around.
  // ringtime stores the current index of the ring-buffer, which
  // wraps back round to 0 once it hits v days.
  Ev_death = vector<int>(v);
  Ev_mosq = vector<vector<Mosquito>>(v);
  ringtime = 0;
  
  // initialise infection history list for storing complete history of events
  //if (output_infection_history) {
  //  history_full = array_int(max_time);
  //}
  
}

//------------------------------------------------
// initialise deme. Assign hosts and seed infections
void Deme::init(vector<int> &host_vec0, int Eh) {
  
  // get number of hosts
  H = int(host_vec0.size());
  
  // initialise counts of mosquito types
  M = M_vec[this_deme];
  Sv = M;
  
  // assign hosts to this deme
  host_vec = host_vec0;
  for (int i=0; i<H; i++) {
    int this_host = host_vec[i];
    hosts[this_host].deme = this_deme;
  }
  
  // seeding infections are de novo, and do not stem from a mosquito. Hence we
  // need a "null" mosquito to carry out this infection. This null mosquito will
  // be recorded in the infection history using -1 values
  Mosquito null_mosquito;
  
  // seed initial infections
  reshuffle(host_vec);
  for (int i=0; i<Eh; i++) {
    int this_host = host_vec[i];
    human_infection(this_host, null_mosquito, 0);
  }
  
  /*
  // save initial infection history
  if (output_infection_history) {
    history_full[0].push_back(history_migration);
    history_full[0].push_back(history_infection);
    history_full[0].push_back(history_bloodstage);
    history_full[0].push_back(history_recover);
  }
  */
  
  // misc
  EIR = 0;
}

//------------------------------------------------
// human new infection
void Deme::human_infection(int this_host, Mosquito &mosq, int t) {
  
  // get ID of this infection
  //int infection_ID = hosts[host_index].cumulative_n_innoculations;
  
  // update host object
  hosts[this_host].new_innoculation(t);
  
  // add to infection history
  //if (output_infection_history) {
    //vector<int> tmp = {host.ID, host.deme, infection_ID, mosq.host_ID, mosq.host_infections, mosq.infection_time};
    //push_back_multiple(history_infection, tmp);
  //}
  
}

//------------------------------------------------
// move time forward one step
void Deme::step_forward(int t) {
  
  /*
  // clear history lists
  if (output_infection_history) {
    history_migration.clear();
    history_infection.clear();
    history_bloodstage.clear();
    history_recover.clear();
  }
  */
  
  
  //-------- MOSQUITO EVENTS --------
  
  // update ring buffer index
  ringtime = (ringtime == v-1) ? 0 : ringtime+1;
  
  // move Ev into Iv
  if (Ev_mosq[ringtime].size() > 0) {
    push_back_multiple(Iv_mosq, Ev_mosq[ringtime]);
    Ev_mosq[ringtime].clear();
  }
  
  // deaths in Iv
  int death_Iv = rbinom1(int(Iv_mosq.size()), prob_v_death);
  for (int i=0; i<death_Iv; i++) {
    int rnd1 = sample2(0, int(Iv_mosq.size())-1);
    quick_erase(Iv_mosq, rnd1);
  }
  
  // draw number of new infections
  double rate_v_infection = a*c*host_infective_vec.size()/double(H); // rate of mosquito infection
  double prob_v_infection_or_death = 1 - exp(-(rate_v_infection + mu)); // mosquito infection or death in Sv state (competing hazards)
  double prob_v_infection = rate_v_infection/(rate_v_infection + mu); // relative rate of mosquito infection vs. death
  int v_infection_or_death = rbinom1(Sv, prob_v_infection_or_death); // number mosquito infection or death in Sv state
  int v_infection = rbinom1(v_infection_or_death, prob_v_infection); // number mosquito infection
  
  // update Sv
  Sv += -v_infection + Ev_death[ringtime] + death_Iv;
  Ev_death[ringtime] = 0; // once these lag stage deaths have been moved into Sv, reset this entry
  
  // the majority of new mosquito infections will die in lag phase. Schedule
  // these deaths to move back into Sv in future steps
  if (v_infection > 0) {
    
    // schedule deaths in Ev state
    int j2 = ringtime;
    for (int j=0; j<v; j++) {
      
      // wrap around ring buffer
      j2 = (j2 == (v-1)) ? 0 : j2++;
      
      // draw number of newly infected mosquitoes that die j steps forward from
      // current time
      int Ev_j_deaths = rbinom1(v_infection, prob_v_death);
      
      // update objects for storing mosquito deaths
      Ev_death[j2] += Ev_j_deaths;
      v_infection -= Ev_j_deaths;
      if (v_infection == 0) {
        break;
      }
    }
    
    // add infections that make it through lag phase to Ev
    for (int i=0; i<v_infection; i++) {
      
      // choose which human host caused the infection
      int rnd1 = sample2(0, int(host_infective_vec.size())-1);
      int this_host = host_infective_vec[rnd1];
      int this_n = 1; // TODO //hosts_infective[rnd1].n_bloodstage + hosts_infective[rnd1].n_infective;
      Ev_mosq[ringtime].emplace_back(this_host, this_n, t);
    }
    
  }
  
  
  //-------- HUMAN EVENTS --------
  
  // get number of new infectious bites on humans
  EIR = a*Iv_mosq.size()/double(H);
  double prob_h_infectious_bite = 1 - exp(-EIR);  // probability of new infectious bite per host
  int h_infectious_bite = rbinom1(H, prob_h_infectious_bite);  // total number of new infectious bites
  
  // apply new infectios bites
  for (int i=0; i<h_infectious_bite; i++) {
    
    // choose host at random
    int rnd1 = sample2(0, H-1);
    int this_host = host_vec[rnd1];
    
    // determine whether infectious bite is successful
    if (rbernoulli1(hosts[this_host].beta)) {
      
      // choose mosquito at random and carry out infection
      int rnd2 = sample2(0, int(Iv_mosq.size())-1);
      human_infection(this_host, Iv_mosq[rnd2], t);
    }
  }
  
  
  //-------- SCHEDULED EVENTS --------
  
  // carry out scheduled events
  for (int i=0; i<int(host_vec.size()); i++) {
    int this_host = host_vec[i];
    
    // check for death
    enact_death(this_host, t);
    
    // enact epidemiological events
    hosts[this_host].enact_events(t, host_infective_vec, i);
  }
  
  
  /*
  // add to infection history
  if (output_infection_history) {
    history_full[t].push_back(history_migration);
    history_full[t].push_back(history_infection);
    history_full[t].push_back(history_bloodstage);
    history_full[t].push_back(history_recover);
  }
  */
}
