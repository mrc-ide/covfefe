
#include "base_model.Deme.h"
#include "probability.h"
#include "misc.h"

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
  // u. These deaths will be added back into the susceptible pool on that day.
  // Ev_mosq stores the full details of infective mosquitoes that are due to
  // emerge from the lag phase at each day going forward. These infective
  // mosquitoes will be added to Iv_mosq on that day. Finally, both Ev_death and
  // Ev_mosq are implemented as ring-buffers, meaning they wrap around.
  // ringtime stores the current index of the ring-buffer, which
  // wraps back round to 0 once it hits u days.
  Ev_death = vector<int>(u);
  Ev_mosq = vector<vector<Mosquito>>(u+1);
  ringtime = 0;
  
  // initialise objects for storing daily counts of key quantities
  if (output_counts) {
    Sh_store = vector<int>(max_time);
    Ih_store = vector<int>(max_time);
    EIR_store = vector<double>(max_time);
  }
  
  // initialise object for storing the number of hosts that carry each possible
  // number of innoculations (from 0 up to max_innoculations) at each time step
  if (output_innoculations) {
    innoculations_store = matrix_int(max_time, max_innoculations + 1);
  }
  
  // initialise infection history list for storing complete history of events
  //if (output_infection_history) {
  //  history_full = array_int(max_time);
  //}
  
}

//------------------------------------------------
// get total beta value (infectiousness) of all hosts
void Deme::get_beta() {
  
  // get beta_uninfective
  beta_uninfective = 0;
  for (int i=0; i<int(hosts_uninfective.size()); i++) {
    int this_host = hosts_uninfective[i];
    beta_uninfective += hosts[this_host].beta;
  }
  
  // get beta_infective
  beta_infective = 0;
  for (int i=0; i<int(hosts_infective.size()); i++) {
    int this_host = hosts_infective[i];
    beta_infective += hosts[this_host].beta;
  }
  
  // get beta_total
  beta_total = beta_uninfective + beta_infective;
}

//------------------------------------------------
// initialise deme. Assign hosts and seed infections
void Deme::init(vector<int> &hosts_uninfective, int Ih) {
  
  // initialise counts of host types
  H = hosts_uninfective.size();
  this->Ih = Ih;
  Sh = H - Ih;
  
  // initialise counts of mosquito types
  M = M_vec[this_deme];
  Sv = M;
  
  // assign hosts to this deme
  this->hosts_uninfective = hosts_uninfective;
  for (int i=0; i<H; i++) {
    int this_host = hosts_uninfective[i];
    hosts[this_host].deme = this_deme;
  }
  
  // seeding infections are de novo, and do not stem from a mosquito. Hence we
  // need a "null" mosquito to carry out this infection. This null mosquito will
  // be recorded in the infection history
  Mosquito null_mosquito;
  
  // seed initial infections
  reshuffle(hosts_uninfective);
  for (int i=0; i<Ih; i++) {
    int this_host = hosts_uninfective[i];
    human_infection(this_host, null_mosquito, 0);
  }
  
  // save initial counts
  if (output_counts) {
    Sh_store[0] = Sh;
    Ih_store[0] = Ih;
    EIR_store[0] = 0;
  }
  
  // save initial innoculations
  if (output_innoculations) {
    innoculations_store[0][0] = Sh;
    innoculations_store[0][1] = Ih;
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
}

//------------------------------------------------
// human new infection
void Deme::human_infection(int host_index, Mosquito &mosq, int t) {
  
  // skip if reached max innoculations
  int n_innoculations = hosts[host_index].get_n_innoculations();
  if (n_innoculations == max_innoculations) {
    return;
  }
  
  // get ID of this infection
  int infection_ID = hosts[host_index].total_infections;
  
  // add to infection history
  //if (output_infection_history) {
    //vector<int> tmp = {host.ID, host.deme, infection_ID, mosq.host_ID, mosq.host_infections, mosq.infection_time};
    //push_back_multiple(history_infection, tmp);
  //}
  
  // schedule future events
  schedule_future_events(host_index, infection_ID, t);
  
  // update host object
  hosts[host_index].new_infection();
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
  ringtime = (ringtime == v) ? 0 : ringtime+1;
  
  // move Ev into Iv
  if (Ev_mosq[ringtime].size() > 0) {
    push_back_multiple(Iv_mosq, Ev_mosq[ringtime]);
    Ev_mosq[ringtime].clear();
  }
  
  // deaths in Iv
  int death_Iv = rbinom1(Iv_mosq.size(), prob_v_death);
  for (int i=0; i<death_Iv; i++) {
    int rnd1 = sample2(0, Iv_mosq.size()-1);
    quick_erase(Iv_mosq, rnd1);
  }
  
  // draw number of new infections
  double rate_v_infection = a*c*hosts_infective.size()/double(H); // rate of mosquito infection
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
      int rnd1 = sample2(0, hosts_infective.size()-1);
      int this_host = hosts_infective[rnd1];
      int this_n = 1; // TODO //hosts_infective[rnd1].n_bloodstage + hosts_infective[rnd1].n_infective;
      Ev_mosq[ringtime].emplace_back(this_host, this_n, t);
    }
    
  }
  
  //print(t, rate_v_infection, v_infection, Iv_mosq.size());
  
  
  //-------- HUMAN EVENTS --------
  /*
  // get beta values
  get_beta();
  
  // get number of new infections
  double EIR = a*Iv_mosq.size()/double(H);
  double rate_h_infection = beta_total*EIR;
  double prob_h_infection = 1 - exp(-rate_h_infection);
  int h_infection = rbinom1(H, prob_h_infection);  // total number new infections
  //double prob_h_infection_noninf = hosts_noninfective.size()/double(H);  // proportion that occur in non-infective hosts
  */
  //print(t, h_infection);
  print(t);
  
  /*
  // apply new infections
  for (int i=0; i<h_infection; i++) {
    
    // choose host at random
    int host_ID;
    if (rbernoulli1(prob_h_infection_noninf)) {
      int rnd1 = sample2(0, hosts_noninfective[k].size()-1);
      host_ID = hosts_noninfective[k][rnd1];
    } else {
      int rnd1 = sample2(0, hosts_infective[k].size()-1);
      host_ID = hosts_infective[k][rnd1];
    }
    
    // choose mosquito at random
    int rnd2 = sample2(0, Iv[k].size()-1);
    
    // carry out infection
    human_infection(host_ID, k, Iv[k][rnd2], t);
    
  } // end loop through new infections
  
  // store counts
  if (output_counts) {
    EIR_store[k][t] = EIR;
  }
  */
  
  /*
  //-------- IMPLEMENT SCHEDULED EVENTS --------
  
  // recovery
  for (int i=0; i<int(schedule_recover[t].size()); i++) {
    int host_ID = schedule_recover[t][i][0];
    int inf_ID = schedule_recover[t][i][1];
    int recovery_type = schedule_recover[t][i][2];
    int host_deme = hosts[host_ID].deme;
    
    // add to recover list
    if (output_infection_history) {
      vector<int> tmp = {host_ID, host_deme, inf_ID};
      push_back_multiple(history_recover, tmp);
    }
    
    // type 0 = blood stage recovery, type 1 = infective recovery
    if (recovery_type==0) {
      hosts[host_ID].n_bloodstage--;
    } else {
      hosts[host_ID].n_infective--;
      
      // move from infective to noninfective list
      if (hosts[host_ID].n_infective==0) {
        hosts_infective[host_deme].erase(remove(hosts_infective[host_deme].begin(), hosts_infective[host_deme].end(), host_ID));
        hosts_noninfective[host_deme].push_back(host_ID);
      }
    }
  }
  */
  //-------- STORE RESULTS --------
  
  // store daily counts of key quantities
  if (output_counts) {
    Sh_store[t] = Sh;
    Ih_store[t] = Ih;
    EIR_store[t] = 0;//EIR;
  }
  
  /*
  // TODO - increment counts as events happen, rather than looping through all
  // hosts?
  if (output_innoculations) {
    for (int k=0; k<demes; k++) {
      int s = hosts_infective[k].size();
      for (int i=0; i<H[k]; i++) {
        int host_ID;
        if (i<s) {
          host_ID = hosts_infective[k][i];
        } else {
          host_ID = hosts_noninfective[k][i-s];
        }
        int nb = hosts[host_ID].n_bloodstage + hosts[host_ID].n_infective;
        if (nb<=max_innoculations) {
          innoculations_full[k][t][nb]++;
        }
      }
    }
  }
  
  // add to infection history
  if (output_infection_history) {
    history_full[t].push_back(history_migration);
    history_full[t].push_back(history_infection);
    history_full[t].push_back(history_bloodstage);
    history_full[t].push_back(history_recover);
  }
  */
}
