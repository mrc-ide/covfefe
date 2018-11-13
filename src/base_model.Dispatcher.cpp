
#include "base_model.Dispatcher.h"
#include "probability.h"

using namespace std;

#ifdef OLD_VERSION

//------------------------------------------------
// default constructor for Dispatcher class
Dispatcher::Dispatcher() {
  
  // objects for sampling from probability distributions
  sampler_age_stable = Sampler(age_stable, 1000);
  sampler_age_death = Sampler(age_death, 1000);
  sampler_duration_acute = vector<Sampler>(n_duration_acute);
  for (int i=0; i<n_duration_acute; ++i) {
    sampler_duration_acute[i] = Sampler(duration_acute[i], 1000);
  }
  sampler_duration_chronic = vector<Sampler>(n_duration_chronic);
  for (int i=0; i<n_duration_chronic; ++i) {
    sampler_duration_chronic[i] = Sampler(duration_chronic[i], 1000);
  }
  
  // events are enacted using scheduler objects. New events (e.g. infection) are
  // generated in the current time step, and future events (e.g. transition to
  // blood-stage) are scheduled for future time steps using these objects. This
  // avoids the need to loop through every host in every time step, as we only
  // need to modify the hosts for which we have scheduled events.
  //
  // some scheduled events may need to be modified - for example, if a host
  // randomly seeks treatment on second innoculation then this will shorten the
  // duration of all innoculations, meaning the scheduled recovery time of the
  // first innoculation must be changed. For this reason, hosts also contain a
  // record of their own event timings, which can be used to locate the
  // corresponding records in the scheduler objects.
  schedule_death = vector<set<int>>(max_time+1);
  schedule_status_update = vector<vector<pair<int, int>>>(max_time+1);
  schedule_infective_start = vector<vector<pair<int, int>>>(max_time+1);
  schedule_infective_stop = vector<vector<pair<int, int>>>(max_time+1);
  
  // counts of host types
  H_total = sum(H_vec);
  H = H_vec;
  Sh = H;
  Lh = vector<int>(n_demes);
  Ah = vector<int>(n_demes);
  Ch = vector<int>(n_demes);
  
  // other quantities to keep track of
  EIR = vector<double>(n_demes);
  
  // initialise single population of human hosts over all demes. This is
  // preferable to using separate vectors of hosts for each deme, as this would
  // mean moving hosts around due to migration. With a single static population
  // we can simply change the "deme" attribute of a host to represent migration
  hosts = vector<Host>(H_total);
  
  // for each deme, store the integer index of all hosts in that deme
  host_vec = vector<vector<int>>(n_demes);
  host_infective_vec = vector<vector<int>>(n_demes);
  int tmp1 = 0;
  for (int k=0; k<n_demes; ++k) {
    host_vec[k] = seq_int(tmp1, tmp1+H[k]-1);
    tmp1 += H[k];
  }
  
  // initialise hosts
  next_host_ID = 0;
  for (int k=0; k<n_demes; ++k) {
    for (int i=0; i<H[k]; i++) {
      int this_host = host_vec[k][i];
      
      // draw age from demography distribution
      int age_years = sampler_age_stable.draw() - 1;
      int extra_days = sample2(1, 365);
      int age_days = age_years*365 + extra_days;
      
      // draw duration of life from demography distribution looking forward from
      // current age. This is tricky, as we must account for the fact that if we
      // are already part way into an age group we have a reduced probability of
      // dying within that age group.
      int life_days = 0;
      double prop_year_remaining = 1 - extra_days/365.0;
      double prob_die_this_year = life_table[age_years]*prop_year_remaining;
      if (rbernoulli1(prob_die_this_year) || age_years == (n_age-1)) {
        life_days = age_years*365 + sample2(extra_days, 365);
      } else {
        for (int i=(age_years+1); i<n_age; ++i) {
          if (rbernoulli1(life_table[i])) {
            life_days = i*365 + sample2(1, 365);
            break;
          }
        }
      }
      
      // convert to final birth and death days
      int birth_day = -age_days;
      int death_day = life_days - age_days;
      if (death_day == 0) {  // in the unlikely even that due to die on day zero, delay death by one day
        death_day++;
      }
      
      // initialise host
      hosts[this_host] = Host(max_innoculations);
      hosts[this_host].reset(next_host_ID++, k, birth_day, death_day);
      
      // schedule death
      if (death_day <= max_time) {
        schedule_death[death_day].insert(this_host);
      }
      
    }
  }
  
  // seeding infections are de novo, and do not stem from a mosquito. Hence we
  // need a "null" mosquito to carry out this infection. This null mosquito will
  // be recorded in the infection history using -1 values
  Mosquito null_mosquito;
  
  // seed initial infections
  for (int k=0; k<n_demes; ++k) {
    reshuffle(host_vec[k]);
    for (int i=0; i<seed_infections[k]; i++) {
      int this_host = host_vec[k][i];
      human_infection(this_host, null_mosquito, 0);
    }
  }
  
  // initialise population of mosquitoes in all demes. An infective mosquito
  // consists of three values: 1) the ID of the host that infected it, 2) the
  // number of observable infections in that host at the time of biting, 3) the
  // time the mosquito became infected by that host. If the infective mosquito
  // ever bites a new host and passes on the infection, these values will be
  // written to the infection history record, thereby providing a link between
  // the original host and the newly infected host.
  // 
  // typically a large number of mosquitoes become infected each time step, but
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
  // mosquitoes will be added to Iv_mosq on that day. Both Ev_death and Ev_mosq
  // are implemented as ring-buffers, meaning they wrap around. ringtime stores
  // the current index of the ring-buffer, which wraps back round to 0 once it
  // hits v days. Ev_death, Ev_mosq and Iv_mosq are duplicated over all demes.
  // The same ringtime is used for all demes.
  Ev_death = vector<vector<int>>(n_demes, vector<int>(v));
  Ev_mosq = vector<vector<vector<Mosquito>>>(n_demes, vector<vector<Mosquito>>(v));
  Iv_mosq = vector<vector<Mosquito>>(n_demes);
  ringtime = 0;
  
  // initialise counts of mosquito types
  Sv = M_vec;
  
  // objects for storing daily counts etc.
  if (output_daily_counts) {
    H_store = vector<vector<int>>(n_demes, vector<int>(max_time+1));
    Sh_store = vector<vector<int>>(n_demes, vector<int>(max_time+1));
    Lh_store = vector<vector<int>>(n_demes, vector<int>(max_time+1));
    Ah_store = vector<vector<int>>(n_demes, vector<int>(max_time+1));
    Ch_store = vector<vector<int>>(n_demes, vector<int>(max_time+1));
    EIR_store = vector<vector<double>>(n_demes, vector<double>(max_time+1));
    
    // counts at time 0
    for (int k=0; k<n_demes; ++k) {
      H_store[k][0] = H[k];
      Sh_store[k][0] = Sh[k];
      Lh_store[k][0] = Lh[k];
    }
  }
  
  // initialise objects for storing age distributions
  if (output_age_distributions) {
    vector<vector<vector<int>>> tmp_mat_int = vector<vector<vector<int>>>(n_demes, vector<vector<int>>(n_output_age_times, vector<int>(n_age)));
    vector<vector<vector<double>>> tmp_mat_double = vector<vector<vector<double>>>(n_demes, vector<vector<double>>(n_output_age_times, vector<double>(n_age)));
    H_age_store = tmp_mat_int;
    prev_Sh_age_store = tmp_mat_double;
    prev_Lh_age_store = tmp_mat_double;
    prev_Ah_age_store = tmp_mat_double;
    prev_Ch_age_store = tmp_mat_double;
    inc_Lh_age_store = tmp_mat_double;
    inc_Ah_age_store = tmp_mat_double;
  }
  
}

//------------------------------------------------
// human new infection
void Dispatcher::human_infection(int this_host, Mosquito &mosq, int t) {
  
  // register infection
  hosts[this_host].cumulative_n_innoculations++;
  
  // get next free innoculation slot. If none free then return without creating
  // new innoculation
  int slot = hosts[this_host].get_innoculation_slot();
  if (slot == -1) {
    return;
  }
  
  // update deme counts
  int this_deme = hosts[this_host].deme;
  if (hosts[this_host].get_n_asexual() == 0) {
    Sh[this_deme]--;
  }
  if (hosts[this_host].n_latent == 0) {
    Lh[this_deme]++;
  }
  
  // update host with new innoculation
  hosts[this_host].new_innoculation(slot, t+u);
  
  // schedule move to blood-stage
  if (t+u <= max_time) {
    schedule_status_update[t+u].emplace_back(this_host, slot);
  }
  
  // get ID of this infection
  //int infection_ID = hosts[host_index].cumulative_n_innoculations;
  
  // add to infection history
  //if (output_infection_history) {
  //vector<int> tmp = {host.ID, host.deme, infection_ID, mosq.host_ID, mosq.host_infections, mosq.infection_time};
  //push_back_multiple(history_infection, tmp);
  //}
  
}

//------------------------------------------------
// run simulation
void Dispatcher::simulate() {
  
  
  //-------- TODO - MIGRATION --------
  
  
  // loop through daily time steps
  int output_age_time_index = 0;
  for (int t=1; t<=max_time; t++) {
    
    // loop through demes
    for (int k=0; k<n_demes; ++k) {
      
      
      //-------- MOSQUITO EVENTS --------
      
      // update ring buffer index
      ringtime = (ringtime == v-1) ? 0 : ringtime+1;
      
      // move Ev into Iv
      if (Ev_mosq[k][ringtime].size() > 0) {
        push_back_multiple(Iv_mosq[k], Ev_mosq[k][ringtime]);
        Ev_mosq[k][ringtime].clear();
      }
      
      // deaths in Iv
      int death_Iv = rbinom1(int(Iv_mosq[k].size()), prob_v_death);
      for (int i=0; i<death_Iv; i++) {
        int rnd1 = sample2(0, int(Iv_mosq[k].size())-1);
        quick_erase(Iv_mosq[k], rnd1);
      }
      
      // draw number of new infections
      double rate_v_infection = a*c*host_infective_vec[k].size()/double(H[k]); // rate of mosquito infection
      double prob_v_infection_or_death = 1 - exp(-(rate_v_infection + mu)); // probability of mosquito infection or death in Sv state (competing hazards)
      double prob_v_infection = rate_v_infection/(rate_v_infection + mu); // relative probability of mosquito infection vs. death
      int v_infection_or_death = rbinom1(Sv[k], prob_v_infection_or_death); // number of mosquito infection or death in Sv state
      int v_infection = rbinom1(v_infection_or_death, prob_v_infection); // number of mosquito infection in Sv state
      
      // update Sv
      Sv[k] += -v_infection + Ev_death[k][ringtime] + death_Iv;
      Ev_death[k][ringtime] = 0; // once these lag stage deaths have been moved into Sv, reset this entry
      
      // the majority of new mosquito infections will die in lag phase. Schedule
      // these deaths to move back into Sv in future steps
      if (v_infection > 0) {
        
        // schedule deaths in Ev state
        int j2 = ringtime;
        for (int j=0; j<v; j++) {
          
          // wrap around ring buffer
          j2 = (j2 == (v-1)) ? 0 : j2+1;
          
          // draw number of newly infected mosquitoes that die j steps forward from
          // current time
          int Ev_j_deaths = rbinom1(v_infection, prob_v_death);
          
          // update objects for storing mosquito deaths
          Ev_death[k][j2] += Ev_j_deaths;
          v_infection -= Ev_j_deaths;
          if (v_infection == 0) {
            break;
          }
        }
        
        // add infections that make it through lag phase to Ev
        for (int i=0; i<v_infection; i++) {
          
          // choose which human host caused the infection
          int rnd1 = sample2(0, int(host_infective_vec[k].size())-1);
          int this_host = host_infective_vec[k][rnd1];
          int this_n = 1; // TODO //hosts_infective[rnd1].n_bloodstage + hosts_infective[rnd1].n_infective;
          Ev_mosq[k][ringtime].emplace_back(this_host, this_n, t);
        }
        
      }
      
      
      //-------- NEW HUMAN EVENTS --------
      
      // get number of new infectious bites on humans
      EIR[k] = a*Iv_mosq[k].size()/double(H[k]);
      double prob_h_infectious_bite = 1 - exp(-EIR[k]);  // probability of new infectious bite per host
      int h_infectious_bite = rbinom1(H[k], prob_h_infectious_bite);  // total number of new infectious bites
      
      // update complete age distributions of incidence
      if (output_age_distributions) {
        if (output_age_times[output_age_time_index] == t) {
          
          // loop through all hosts
          for (int i=0; i<H[k]; ++i) {
            int this_host = host_vec[k][i];
            int age_days = t - hosts[this_host].birth_day;
            int age_years = floor(age_days/365.0);
            double this_b = b[hosts[this_host].b_index];
            double this_prob_acute = prob_acute[hosts[this_host].prob_acute_index];
            
            // update raw number of hosts
            H_age_store[k][output_age_time_index][age_years]++;
            
            // update incidence of infection
            inc_Lh_age_store[k][output_age_time_index][age_years] += prob_h_infectious_bite*this_b;
            
            // update incidence of acute disease
            inc_Ah_age_store[k][output_age_time_index][age_years] += prob_h_infectious_bite*this_b*this_prob_acute;
            
          }
        }
      }
      
      // apply new infectios bites
      for (int i=0; i<h_infectious_bite; i++) {
        
        // choose host at random
        int rnd1 = sample2(0, H[k]-1);
        int this_host = host_vec[k][rnd1];
        
        // determine whether infectious bite is successful
        double host_b = b[hosts[this_host].b_index];
        if (rbernoulli1(host_b)) {
          
          // choose mosquito at random and carry out infection
          int rnd2 = sample2(0, int(Iv_mosq[k].size())-1);
          human_infection(this_host, Iv_mosq[k][rnd2], t);
        }
        
        // update b_index irrespective of whether infection successful
        if (hosts[this_host].b_index < (n_b-1)) {
          hosts[this_host].b_index++;
        }
        
      }
      
    }  // end loop through demes
    
    
    //-------- SCHEDULED HUMAN EVENTS --------
    
    // scheduled deaths
    for (auto it = schedule_death[t].begin(); it != schedule_death[t].end(); ++it) {
      int this_host = *it;
      int this_home_deme = hosts[this_host].deme;
      int this_deme = hosts[this_host].deme;
      
      // update deme counts
      if (hosts[this_host].get_n_asexual() > 0) {
        Sh[this_deme]++;
      }
      if (hosts[this_host].n_latent > 0) {
        Lh[this_deme]--;
      }
      if (hosts[this_host].n_acute > 0) {
        Ah[this_deme]--;
      }
      if (hosts[this_host].n_chronic > 0) {
        Ch[this_deme]--;
      }
      
      // drop from infective list if necessary
      if (hosts[this_host].n_infective > 0) {
        host_infective_vec[this_deme].erase(remove(host_infective_vec[this_deme].begin(), host_infective_vec[this_deme].end(), this_host));
      }
      
      // draw life duration from demography distribution
      int life_years = sampler_age_death.draw() - 1;
      int life_days = life_years*365 + sample2(1, 365);
      int death_day = t + life_days;
      hosts[this_host].reset(next_host_ID++, this_home_deme, t, death_day);
      
      // schedule new death
      if (death_day <= max_time) {
        schedule_death[death_day].insert(this_host);
      }
    }
    
    // scheduled status updates
    for (auto it = schedule_status_update[t].begin(); it != schedule_status_update[t].end(); ++it) {
      int this_host = it->first;
      int this_slot = it->second;
      int this_deme = hosts[this_host].deme;
      int this_death_day = hosts[this_host].death_day;
      
      // switch on current status
      switch (hosts[this_host].innoc_status[this_slot]) {
      
      // latent become clinical or asymptomatic
      case Latent: {
        
        // draw whether becomes acute or chronic
        double host_p_acute = prob_acute[hosts[this_host].prob_acute_index];
        bool become_acute = rbernoulli1(host_p_acute);
        
        // update p_symptomatic_index
        if (hosts[this_host].prob_acute_index < (n_prob_acute-1)) {
          hosts[this_host].prob_acute_index++;
        }
        
        // split method between clinical and asymptomatic
        int duration = 0;
        if (become_acute) {
          
          // update status
          hosts[this_host].innoc_status[this_slot] = Acute;
          
          // update status time
          duration = sampler_duration_acute[hosts[this_host].duration_acute_index].draw() + 1;
          hosts[this_host].innoc_status_update_time[this_slot] = t + duration;
          
          // update duration_acute_index
          if (hosts[this_host].duration_acute_index < (n_duration_acute-1)) {
            hosts[this_host].duration_acute_index++;
          }
          
          // update host counts
          hosts[this_host].n_latent--;
          hosts[this_host].n_acute++;
          
          // update deme counts
          if (hosts[this_host].n_latent == 0) {
            Lh[this_deme]--;
          }
          if (hosts[this_host].n_acute == 1) {
            Ah[this_deme]++;
          }
          
        } else {
          
          // update status
          hosts[this_host].innoc_status[this_slot] = Chronic;
          
          // update status time
          duration = sampler_duration_chronic[hosts[this_host].duration_chronic_index].draw() + 1;
          hosts[this_host].innoc_status_update_time[this_slot] = t + duration;
          
          // update duration_acute_index
          if (hosts[this_host].duration_chronic_index < (n_duration_chronic-1)) {
            hosts[this_host].duration_chronic_index++;
          }
          
          // update host counts
          hosts[this_host].n_latent--;
          hosts[this_host].n_chronic++;
          
          // update deme counts
          if (hosts[this_host].n_latent == 0) {
            Lh[this_deme]--;
          }
          if (hosts[this_host].n_chronic == 1) {
            Ch[this_deme]++;
          }
        }
        
        // update infective start time
        hosts[this_host].innoc_infective_start_time[this_slot] = t+g;
        
        // schedule new events
        if (t+duration <= max_time && t+duration < this_death_day) {
          schedule_status_update[t+duration].emplace_back(this_host, this_slot);
        }
        if (t+g <= max_time && t+g < this_death_day) {
          schedule_infective_start[t+g].emplace_back(this_host, this_slot);
        }
        
        break;
      } 
        
        // acute phase recover
      case Acute: {
        
        // draw whether recovers or becomes chronic
        bool become_chronic = rbernoulli1(prob_AC);
        
        // split method between recovery and chronic
        int duration = 0;
        if (become_chronic) {
          
          // update status
          hosts[this_host].innoc_status[this_slot] = Chronic;
          
          // update status time
          duration = sampler_duration_chronic[hosts[this_host].duration_chronic_index].draw() + 1;
          hosts[this_host].innoc_status_update_time[this_slot] = t + duration;
          
          // update duration_acute_index
          if (hosts[this_host].duration_chronic_index < (n_duration_chronic-1)) {
            hosts[this_host].duration_chronic_index++;
          }
          
          // update host counts
          hosts[this_host].n_acute--;
          hosts[this_host].n_chronic++;
          
          // update deme counts
          if (hosts[this_host].n_acute == 0) {
            Ah[this_deme]--;
          }
          if (hosts[this_host].n_chronic == 1) {
            Ch[this_deme]++;
          }
          
          // schedule new events
          if (t+duration <= max_time && t+duration < this_death_day) {
            schedule_status_update[t+duration].emplace_back(this_host, this_slot);
          }
          
        } else {
          
          // update status
          hosts[this_host].innoc_status[this_slot] = Inactive;
          
          // update host counts
          hosts[this_host].n_acute--;
          
          // update deme counts
          if (hosts[this_host].get_n_asexual() == 0) {
            Sh[this_deme]++;
          }
          if (hosts[this_host].n_acute == 0) {
            Ah[this_deme]--;
          }
          
          // update infective stop time
          hosts[this_host].innoc_infective_stop_time[this_slot] = t+g;
          
          // schedule new events
          if (t+g <= max_time && t+g < this_death_day) {
            schedule_infective_stop[t+g].emplace_back(this_host, this_slot);
          }
        }
        
        break;
      }
        
        // asymptomatic recover
      case Chronic: {
        
        // update status
        hosts[this_host].innoc_status[this_slot] = Inactive;
        
        // update host counts
        hosts[this_host].n_chronic--;
        
        // update deme counts
        if (hosts[this_host].get_n_asexual() == 0) {
          Sh[this_deme]++;
        }
        if (hosts[this_host].n_chronic == 0) {
          Ch[this_deme]--;
        }
        
        // update infective stop time
        hosts[this_host].innoc_infective_stop_time[this_slot] = t+g;
        
        // schedule new events
        if (t+g <= max_time && t+g < this_death_day) {
          schedule_infective_stop[t+g].emplace_back(this_host, this_slot);
        }
        
        break;
      }
        
      default:
        break;
      }
      
    }
    
    // scheduled infection start events
    for (auto it = schedule_infective_start[t].begin(); it != schedule_infective_start[t].end(); ++it) {
      int this_host = it->first;
      int slot = it->second;
      
      // if newly infective then add to infectives list
      if (hosts[this_host].n_infective == 0) {
        int this_deme = hosts[this_host].deme;
        host_infective_vec[this_deme].push_back(this_host);
      }
      
      // update infective status and counts
      hosts[this_host].innoc_infective[slot] = true;
      hosts[this_host].n_infective++;
    }
    
    // scheduled infection stop events
    for (auto it = schedule_infective_stop[t].begin(); it != schedule_infective_stop[t].end(); ++it) {
      int this_host = it->first;
      int slot = it->second;
      
      // update infective status and counts
      hosts[this_host].innoc_active[slot] = false;
      hosts[this_host].n_infective--;
      
      // if no longer infective then drop from infectives list
      if (hosts[this_host].n_infective == 0) {
        int this_deme = hosts[this_host].deme;
        host_infective_vec[this_deme].erase(remove(host_infective_vec[this_deme].begin(), host_infective_vec[this_deme].end(), this_host));
      }
    }
    
    
    //-------- STORE RESULTS --------
    
    // store daily counts etc.
    if (output_daily_counts) {
      for (int k=0; k<n_demes; ++k) {
        H_store[k][t] = H[k];
        Sh_store[k][t] = Sh[k];
        Lh_store[k][t] = Lh[k];
        Ah_store[k][t] = Ah[k];
        Ch_store[k][t] = Ch[k];
        EIR_store[k][t] = EIR[k];
      }
    }
    
    // update complete age distributions
    if (output_age_distributions) {
      if (output_age_times[output_age_time_index] == t) {
        for (int k=0; k<n_demes; ++k) {
          
          // loop through all hosts
          for (int i=0; i<H[k]; ++i) {
            int this_host = host_vec[k][i];
            int age_days = t - hosts[this_host].birth_day;
            int age_years = floor(age_days/365.0);
            
            // update susceptible prevalence
            if (hosts[this_host].get_n_innoculations() == 0) {
              prev_Sh_age_store[k][output_age_time_index][age_years] += 1.0;
            }
            
            // update latent prevalence
            if (hosts[this_host].n_latent > 0) {
              prev_Lh_age_store[k][output_age_time_index][age_years] += 1.0;
            }
            
            // update acute prevalence
            if (hosts[this_host].n_acute > 0) {
              prev_Ah_age_store[k][output_age_time_index][age_years] += 1.0;
            }
            
            // update chronic prevalence
            if (hosts[this_host].n_chronic > 0) {
              prev_Ch_age_store[k][output_age_time_index][age_years] += 1.0;
            }
            
          }
          
          // divide by number of hosts per age bin
          for (int i=0; i<n_age; ++i) {
            prev_Sh_age_store[k][output_age_time_index][i] /= H_age_store[k][output_age_time_index][i];
            prev_Lh_age_store[k][output_age_time_index][i] /= H_age_store[k][output_age_time_index][i];
            prev_Ah_age_store[k][output_age_time_index][i] /= H_age_store[k][output_age_time_index][i];
            prev_Ch_age_store[k][output_age_time_index][i] /= H_age_store[k][output_age_time_index][i];
            
            inc_Lh_age_store[k][output_age_time_index][i] /= H_age_store[k][output_age_time_index][i];
            inc_Ah_age_store[k][output_age_time_index][i] /= H_age_store[k][output_age_time_index][i];
          }
        }
        
        // update output_age_time_index
        if (output_age_time_index < (output_age_times.size()-1)) {
          output_age_time_index++;
        }
      }
    }
    
  } // end time loop
  
}

#else

//------------------------------------------------
// default constructor for Dispatcher class
Dispatcher::Dispatcher() {
  
  // objects for sampling from probability distributions
  sampler_age_stable = Sampler(age_stable, 1000);
  sampler_age_death = Sampler(age_death, 1000);
  sampler_duration_acute = vector<Sampler>(n_duration_acute);
  for (int i=0; i<n_duration_acute; ++i) {
    sampler_duration_acute[i] = Sampler(duration_acute[i], 1000);
  }
  sampler_duration_chronic = vector<Sampler>(n_duration_chronic);
  for (int i=0; i<n_duration_chronic; ++i) {
    sampler_duration_chronic[i] = Sampler(duration_chronic[i], 1000);
  }
  
  // events are enacted using scheduler objects. New events (e.g. infection) are
  // generated in the current time step, and future events (e.g. transition to
  // blood-stage) are scheduled for future time steps using these objects. This
  // avoids the need to loop through every host in every time step, as we only
  // need to modify the hosts for which we have scheduled events.
  //
  // some scheduled events may need to be modified - for example, if a host
  // randomly seeks treatment on second innoculation then this will shorten the
  // duration of all innoculations, meaning the scheduled recovery time of the
  // first innoculation must be changed. For this reason, hosts also contain a
  // record of their own event timings, which can be used to locate the
  // corresponding records in the scheduler objects.
  schedule_death = vector<set<int>>(max_time+1);
  schedule_status_update = vector<vector<pair<int, int>>>(max_time+1);
  schedule_infective_status_update = vector<vector<pair<int, int>>>(max_time+1);
  
  // counts of host types
  H_total = sum(H_vec);
  H = H_vec;
  Sh = H;
  Lh = vector<int>(n_demes);
  Ah = vector<int>(n_demes);
  Ch = vector<int>(n_demes);
  
  // other quantities to keep track of
  EIR = vector<double>(n_demes);
  
  // initialise single population of human hosts over all demes. This is
  // preferable to using separate vectors of hosts for each deme, as this would
  // mean moving hosts around due to migration. With a single static population
  // we can simply change the "deme" attribute of a host to represent migration
  hosts = vector<Host>(H_total);
  
  // for each deme, store the integer index of all hosts in that deme
  host_vec = vector<vector<int>>(n_demes);
  host_infective_vec = vector<vector<int>>(n_demes);
  int tmp1 = 0;
  for (int k=0; k<n_demes; ++k) {
    host_vec[k] = seq_int(tmp1, tmp1+H[k]-1);
    tmp1 += H[k];
  }
  
  // initialise hosts
  next_host_ID = 0;
  for (int k=0; k<n_demes; ++k) {
    for (int i=0; i<H[k]; i++) {
      int this_host = host_vec[k][i];
      
      // draw age from demography distribution
      int age_years = sampler_age_stable.draw() - 1;
      int extra_days = sample2(1, 365);
      int age_days = age_years*365 + extra_days;
      
      // draw duration of life from demography distribution looking forward from
      // current age. This is tricky, as we must account for the fact that if we
      // are already part way into an age group we have a reduced probability of
      // dying within that age group.
      int life_days = 0;
      double prop_year_remaining = 1 - extra_days/365.0;
      double prob_die_this_year = life_table[age_years]*prop_year_remaining;
      if (rbernoulli1(prob_die_this_year) || age_years == (n_age-1)) {
        life_days = age_years*365 + sample2(extra_days, 365);
      } else {
        for (int i=(age_years+1); i<n_age; ++i) {
          if (rbernoulli1(life_table[i])) {
            life_days = i*365 + sample2(1, 365);
            break;
          }
        }
      }
      
      // convert to final birth and death days
      int birth_day = -age_days;
      int death_day = life_days - age_days;
      if (death_day == 0) {  // in the unlikely even that due to die on day zero, delay death by one day
        death_day++;
      }
      
      // initialise host
      hosts[this_host] = Host(max_innoculations);
      hosts[this_host].reset(next_host_ID++, k, birth_day, death_day);
      
      // schedule death
      if (death_day <= max_time) {
        schedule_death[death_day].insert(this_host);
      }
      
    }
  }
  
  // seeding infections are de novo, and do not stem from a mosquito. Hence we
  // need a "null" mosquito to carry out this infection. This null mosquito will
  // be recorded in the infection history using -1 values
  Mosquito null_mosquito;
  
  // seed initial infections
  for (int k=0; k<n_demes; ++k) {
    reshuffle(host_vec[k]);
    for (int i=0; i<seed_infections[k]; i++) {
      int this_host = host_vec[k][i];
      human_infection(this_host, null_mosquito, 0);
    }
  }
  
  // initialise population of mosquitoes in all demes. An infective mosquito
  // consists of three values: 1) the ID of the host that infected it, 2) the
  // number of observable infections in that host at the time of biting, 3) the
  // time the mosquito became infected by that host. If the infective mosquito
  // ever bites a new host and passes on the infection, these values will be
  // written to the infection history record, thereby providing a link between
  // the original host and the newly infected host.
  // 
  // typically a large number of mosquitoes become infected each time step, but
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
  // mosquitoes will be added to Iv_mosq on that day. Both Ev_death and Ev_mosq
  // are implemented as ring-buffers, meaning they wrap around. ringtime stores
  // the current index of the ring-buffer, which wraps back round to 0 once it
  // hits v days. Ev_death, Ev_mosq and Iv_mosq are duplicated over all demes.
  // The same ringtime is used for all demes.
  Ev_death = vector<vector<int>>(n_demes, vector<int>(v));
  Ev_mosq = vector<vector<vector<Mosquito>>>(n_demes, vector<vector<Mosquito>>(v));
  Iv_mosq = vector<vector<Mosquito>>(n_demes);
  ringtime = 0;
  
  // initialise counts of mosquito types
  Sv = M_vec;
  
  // objects for storing daily counts etc.
  if (output_daily_counts) {
    H_store = vector<vector<int>>(n_demes, vector<int>(max_time+1));
    Sh_store = vector<vector<int>>(n_demes, vector<int>(max_time+1));
    Lh_store = vector<vector<int>>(n_demes, vector<int>(max_time+1));
    Ah_store = vector<vector<int>>(n_demes, vector<int>(max_time+1));
    Ch_store = vector<vector<int>>(n_demes, vector<int>(max_time+1));
    EIR_store = vector<vector<double>>(n_demes, vector<double>(max_time+1));
    
    // counts at time 0
    for (int k=0; k<n_demes; ++k) {
      H_store[k][0] = H[k];
      Sh_store[k][0] = Sh[k];
      Lh_store[k][0] = Lh[k];
    }
  }
  
  // initialise objects for storing age distributions
  if (output_age_distributions) {
    vector<vector<vector<int>>> tmp_mat_int = vector<vector<vector<int>>>(n_demes, vector<vector<int>>(n_output_age_times, vector<int>(n_age)));
    vector<vector<vector<double>>> tmp_mat_double = vector<vector<vector<double>>>(n_demes, vector<vector<double>>(n_output_age_times, vector<double>(n_age)));
    H_age_store = tmp_mat_int;
    prev_Sh_age_store = tmp_mat_double;
    prev_Lh_age_store = tmp_mat_double;
    prev_Ah_age_store = tmp_mat_double;
    prev_Ch_age_store = tmp_mat_double;
    inc_Lh_age_store = tmp_mat_double;
    inc_Ah_age_store = tmp_mat_double;
  }
  
}

//------------------------------------------------
// human new infection
void Dispatcher::human_infection(int this_host, Mosquito &mosq, int t) {
  
  // update cumulative innoculations, irrespective of whether this innoculation
  // passes the next series of checks
  hosts[this_host].cumulative_n_innoculations++;
  
  // get next free innoculation slot. If none free then return without creating
  // new innoculation
  int this_slot = 0;
  for (int i=0; i<max_innoculations; ++i) {
    if (hosts[this_host].innoc_status[i] == Inactive && hosts[this_host].innoc_infective_status[i] == Inactive) {
      break;
    }
    this_slot++;
  }
  if (this_slot == max_innoculations) {
    return;
  }
  
  // update deme counts
  int this_deme = hosts[this_host].deme;
  if (hosts[this_host].get_n_asexual() == 0) {
    Sh[this_deme]--;
  }
  if (hosts[this_host].n_latent == 0) {
    Lh[this_deme]++;
  }
  
  // update host counts
  hosts[this_host].n_latent++;
  
  // add new innoculation and reset other categories
  hosts[this_host].innoc_status[this_slot] = Latent;
  hosts[this_host].innoc_next_status[this_slot] = Null;
  hosts[this_host].innoc_status_update_time[this_slot] = t+u;
  hosts[this_host].innoc_infective_status[this_slot] = Inactive;
  hosts[this_host].innoc_next_infective_status[this_slot] = Null;
  hosts[this_host].innoc_infective_status_update_time[this_slot] = 0;
  
  // schedule status update
  if (t+u <= max_time) {
    schedule_status_update[t+u].emplace_back(this_host, this_slot);
  }
  
  // get ID of this infection
  //int infection_ID = hosts[host_index].cumulative_n_innoculations;
  
  // add to infection history
  //if (output_infection_history) {
  //vector<int> tmp = {host.ID, host.deme, infection_ID, mosq.host_ID, mosq.host_infections, mosq.infection_time};
  //push_back_multiple(history_infection, tmp);
  //}
  
}

//------------------------------------------------
// run simulation
void Dispatcher::simulate() {
  
  
  //-------- TODO - MIGRATION --------
  
  
  // loop through daily time steps
  int output_age_time_index = 0;
  for (int t=1; t<=max_time; t++) {
    
    // loop through demes
    for (int k=0; k<n_demes; ++k) {
      
      
      //-------- MOSQUITO EVENTS --------
      
      // update ring buffer index
      ringtime = (ringtime == v-1) ? 0 : ringtime+1;
      
      // move Ev into Iv
      if (Ev_mosq[k][ringtime].size() > 0) {
        push_back_multiple(Iv_mosq[k], Ev_mosq[k][ringtime]);
        Ev_mosq[k][ringtime].clear();
      }
      
      // deaths in Iv
      int death_Iv = rbinom1(int(Iv_mosq[k].size()), prob_v_death);
      for (int i=0; i<death_Iv; i++) {
        int rnd1 = sample2(0, int(Iv_mosq[k].size())-1);
        quick_erase(Iv_mosq[k], rnd1);
      }
      
      // draw number of new infections
      double rate_v_bite_infective = a*host_infective_vec[k].size()/double(H[k]); // rate of mosquito biting infective hosts
      double prob_v_bite_infective_or_death = 1 - exp(-(rate_v_bite_infective + mu)); // probability of mosquito biting infective host or dying in Sv state (competing hazards)
      double prob_v_bite_infective = rate_v_bite_infective/(rate_v_bite_infective + mu); // relative probability of mosquito biting infective host vs. dying
      int v_bite_infective_or_death = rbinom1(Sv[k], prob_v_bite_infective_or_death); // number of mosquito biting infective hosts or dying in Sv state
      int v_bite_infective = rbinom1(v_bite_infective_or_death, prob_v_bite_infective); // number of mosquito biting infective hosts
      
      // update Sv
      Sv[k] += -v_bite_infective + Ev_death[k][ringtime] + death_Iv;
      Ev_death[k][ringtime] = 0; // once these lag stage deaths have been moved into Sv, reset this entry
      
      // loop through bites on infective hosts
      for (int i=0; i<v_bite_infective; ++i) {
        
        // choose host at random
        int rnd1 = sample2(0, host_infective_vec[k].size()-1);
        int this_host = host_infective_vec[k][rnd1];
        
        // determine whether infectious bite is successful
        bool host_acute = (hosts[this_host].n_acute > 0);
        double prob_infective = 0; 
        if (host_acute) {
          prob_infective = infectivity_acute[hosts[this_host].infectivity_acute_index];
        } else {
          prob_infective = infectivity_chronic[hosts[this_host].infectivity_chronic_index];
        }
        prob_infective = 1.0;
        if (rbernoulli1(prob_infective)) {
          
          // the majority of new mosquito infections will die in lag phase.
          // Schedule these deaths to move back into Sv in future steps.
          // Otherwise add to Ev_mosq
          int v_time_death = rgeom1(prob_v_death) + 1;
          if (v_time_death <= v) {
            Ev_death[k][(ringtime+v_time_death) % v]++;
          } else {
            int this_n = 1; // TODO //hosts_infective[rnd1].n_bloodstage + hosts_infective[rnd1].n_infective;
            Ev_mosq[k][ringtime].emplace_back(this_host, this_n, t);
          }
        }
      }
      
      
      //-------- NEW HUMAN EVENTS AND STORE INCIDENCE --------
      
      // get number of new infectious bites on humans
      EIR[k] = a*Iv_mosq[k].size()/double(H[k]);
      double prob_h_infectious_bite = 1 - exp(-EIR[k]);  // probability of new infectious bite per host
      int h_infectious_bite = rbinom1(H[k], prob_h_infectious_bite);  // total number of new infectious bites
      
      // update complete age distributions of incidence
      if (output_age_distributions) {
        if (output_age_times[output_age_time_index] == t) {
          
          // loop through all hosts
          for (int i=0; i<H[k]; ++i) {
            int this_host = host_vec[k][i];
            int age_days = t - hosts[this_host].birth_day;
            int age_years = floor(age_days/365.0);
            double this_b = b[hosts[this_host].b_index];
            double this_prob_acute = prob_acute[hosts[this_host].prob_acute_index];
            
            // update raw number of hosts
            H_age_store[k][output_age_time_index][age_years]++;
            
            // update incidence of infection
            inc_Lh_age_store[k][output_age_time_index][age_years] += prob_h_infectious_bite*this_b;
            
            // update incidence of acute disease
            inc_Ah_age_store[k][output_age_time_index][age_years] += prob_h_infectious_bite*this_b*this_prob_acute;
            
          }
        }
      }
      
      // apply new infectious bites
      for (int i=0; i<h_infectious_bite; i++) {
        
        // choose host at random
        int rnd1 = sample2(0, H[k]-1);
        int this_host = host_vec[k][rnd1];
        
        // determine whether infectious bite is successful
        double host_b = b[hosts[this_host].b_index];
        if (rbernoulli1(host_b)) {
          
          // choose mosquito at random and carry out infection
          int rnd2 = sample2(0, int(Iv_mosq[k].size())-1);
          human_infection(this_host, Iv_mosq[k][rnd2], t);
        }
        
        // update b_index irrespective of whether infection successful
        if (hosts[this_host].b_index < (n_b-1)) {
          hosts[this_host].b_index++;
        }
        
      }
      
    }  // end loop through demes
    
    
    //-------- SCHEDULED HUMAN EVENTS --------
    
    // scheduled deaths
    for (auto it = schedule_death[t].begin(); it != schedule_death[t].end(); ++it) {
      int this_host = *it;
      int this_home_deme = hosts[this_host].deme;
      int this_deme = hosts[this_host].deme;
      
      // update deme counts
      if (hosts[this_host].get_n_asexual() > 0) {
        Sh[this_deme]++;
      }
      if (hosts[this_host].n_latent > 0) {
        Lh[this_deme]--;
      }
      if (hosts[this_host].n_acute > 0) {
        Ah[this_deme]--;
      }
      if (hosts[this_host].n_chronic > 0) {
        Ch[this_deme]--;
      }
      
      // drop from infective list if necessary
      if (hosts[this_host].get_n_infective() > 0) {
        host_infective_vec[this_deme].erase(remove(host_infective_vec[this_deme].begin(), host_infective_vec[this_deme].end(), this_host));
      }
      
      // draw life duration from demography distribution
      int life_years = sampler_age_death.draw() - 1;
      int life_days = life_years*365 + sample2(1, 365);
      int death_day = t + life_days;
      hosts[this_host].reset(next_host_ID++, this_home_deme, t, death_day);
      
      // schedule new death
      if (death_day <= max_time) {
        schedule_death[death_day].insert(this_host);
      }
    }
    
    // scheduled status updates
    for (auto it = schedule_status_update[t].begin(); it != schedule_status_update[t].end(); ++it) {
      int this_host = it->first;
      int this_slot = it->second;
      int this_deme = hosts[this_host].deme;
      int this_death_day = hosts[this_host].death_day;
      
      // switch depending on current status
      switch (hosts[this_host].innoc_status[this_slot]) {
        
        // latent become acute or chronic
        case Latent: {
          
          // draw whether becomes acute or chronic
          double host_prob_acute = prob_acute[hosts[this_host].prob_acute_index];
          bool become_acute = rbernoulli1(host_prob_acute);
          
          // update prob_acute_index
          if (hosts[this_host].prob_acute_index < (n_prob_acute-1)) {
            hosts[this_host].prob_acute_index++;
          }
          
          // split method between acute and chronic
          if (become_acute) {
            
            // update status
            hosts[this_host].innoc_status[this_slot] = Acute;
            
            // update status update time
            int duration = sampler_duration_acute[hosts[this_host].duration_acute_index].draw() + 1;
            hosts[this_host].innoc_status_update_time[this_slot] = t + duration;
            
            // update duration_acute_index
            if (hosts[this_host].duration_acute_index < (n_duration_acute-1)) {
              hosts[this_host].duration_acute_index++;
            }
            
            // update infectivity_acute_index
            if (hosts[this_host].infectivity_acute_index < (n_infectivity_acute-1)) {
              hosts[this_host].infectivity_acute_index++;
            }
            
            // update next infective status
            hosts[this_host].innoc_next_infective_status[this_slot] = Acute;
            
            // update infective status update time
            hosts[this_host].innoc_infective_status_update_time[this_slot] = t+g;
            
            // update host counts
            hosts[this_host].n_latent--;
            hosts[this_host].n_acute++;
            
            // update deme counts
            if (hosts[this_host].n_latent == 0) {
              Lh[this_deme]--;
            }
            if (hosts[this_host].n_acute == 1) {
              Ah[this_deme]++;
            }
            
            // schedule new events
            if (t+duration <= max_time && t+duration < this_death_day) {
              schedule_status_update[t+duration].emplace_back(this_host, this_slot);
            }
            if (t+g <= max_time && t+g < this_death_day) {
              schedule_infective_status_update[t+g].emplace_back(this_host, this_slot);
            }
            
          } else {
            
            // update status
            hosts[this_host].innoc_status[this_slot] = Chronic;
            
            // update status update time
            int duration = sampler_duration_chronic[hosts[this_host].duration_chronic_index].draw() + 1;
            hosts[this_host].innoc_status_update_time[this_slot] = t + duration;
            
            // update duration_chronic_index
            if (hosts[this_host].duration_chronic_index < (n_duration_chronic-1)) {
              hosts[this_host].duration_chronic_index++;
            }
            
            // update infectivity_chronic_index
            if (hosts[this_host].infectivity_chronic_index < (n_infectivity_chronic-1)) {
              hosts[this_host].infectivity_chronic_index++;
            }
            
            // update next infective status
            hosts[this_host].innoc_next_infective_status[this_slot] = Chronic;
            
            // update infective status update time
            hosts[this_host].innoc_infective_status_update_time[this_slot] = t+g;
            
            // update host counts
            hosts[this_host].n_latent--;
            hosts[this_host].n_chronic++;
            
            // update deme counts
            if (hosts[this_host].n_latent == 0) {
              Lh[this_deme]--;
            }
            if (hosts[this_host].n_chronic == 1) {
              Ch[this_deme]++;
            }
            
            // schedule new events
            if (t+duration <= max_time && t+duration < this_death_day) {
              schedule_status_update[t+duration].emplace_back(this_host, this_slot);
            }
            if (t+g <= max_time && t+g < this_death_day) {
              schedule_infective_status_update[t+g].emplace_back(this_host, this_slot);
            }
            
          }  // end if
          
          break;
        } 
      
        // acute phase become chronic or recover
        case Acute: {
          
          // draw whether recovers or becomes chronic
          bool become_chronic = rbernoulli1(prob_AC);
          
          // split method between chronic disease and recovery
          if (become_chronic) {
            
            // update status
            hosts[this_host].innoc_status[this_slot] = Chronic;
            
            // update status time
            int duration = sampler_duration_chronic[hosts[this_host].duration_chronic_index].draw() + 1;
            hosts[this_host].innoc_status_update_time[this_slot] = t + duration;
            
            // update duration_chronic_index
            if (hosts[this_host].duration_chronic_index < (n_duration_chronic-1)) {
              hosts[this_host].duration_chronic_index++;
            }
            
            // update next infective status
            hosts[this_host].innoc_next_infective_status[this_slot] = Chronic;
            
            // update infective status update time
            hosts[this_host].innoc_infective_status_update_time[this_slot] = t+g;
            
            // update host counts
            hosts[this_host].n_acute--;
            hosts[this_host].n_chronic++;
            
            // update deme counts
            if (hosts[this_host].n_acute == 0) {
              Ah[this_deme]--;
            }
            if (hosts[this_host].n_chronic == 1) {
              Ch[this_deme]++;
            }
            
            // schedule new events
            if (t+duration <= max_time && t+duration < this_death_day) {
              schedule_status_update[t+duration].emplace_back(this_host, this_slot);
            }
            if (t+g <= max_time && t+g < this_death_day) {
              schedule_infective_status_update[t+g].emplace_back(this_host, this_slot);
            }
            
          } else {
            
            // update status
            hosts[this_host].innoc_status[this_slot] = Inactive;
            
            // update next infective status
            hosts[this_host].innoc_next_infective_status[this_slot] = Inactive;
            
            // update infective status update time
            hosts[this_host].innoc_infective_status_update_time[this_slot] = t+g;
            
            // update host counts
            hosts[this_host].n_acute--;
            
            // update deme counts
            if (hosts[this_host].get_n_asexual() == 0) {
              Sh[this_deme]++;
            }
            if (hosts[this_host].n_acute == 0) {
              Ah[this_deme]--;
            }
            
            // schedule new events
            if (t+g <= max_time && t+g < this_death_day) {
              schedule_infective_status_update[t+g].emplace_back(this_host, this_slot);
            }
            
          }  // end if
          
          break;
        }
        
        // chronic recover
        case Chronic: {
          
          // update status
          hosts[this_host].innoc_status[this_slot] = Inactive;
          
          // update next infective status
          hosts[this_host].innoc_next_infective_status[this_slot] = Inactive;
          
          // update infective status update time
          hosts[this_host].innoc_infective_status_update_time[this_slot] = t+g;
          
          // update host counts
          hosts[this_host].n_chronic--;
          
          // update deme counts
          if (hosts[this_host].get_n_asexual() == 0) {
            Sh[this_deme]++;
          }
          if (hosts[this_host].n_chronic == 0) {
            Ch[this_deme]--;
          }
          
          // schedule new events
          if (t+g <= max_time && t+g < this_death_day) {
            schedule_infective_status_update[t+g].emplace_back(this_host, this_slot);
          }
          
          break;
        }
        
        default:
          break;
      }  // end switch
      
    }
    
    // scheduled infectivity events
    for (auto it = schedule_infective_status_update[t].begin(); it != schedule_infective_status_update[t].end(); ++it) {
      int this_host = it->first;
      int this_slot = it->second;
      int this_deme = hosts[this_host].deme;
      bool infective_before = (hosts[this_host].get_n_infective() > 0);
      
      // subtract from current counts
      if (hosts[this_host].innoc_infective_status[this_slot] == Acute) {
        hosts[this_host].n_infective_acute--;
      }
      if (hosts[this_host].innoc_infective_status[this_slot] == Chronic) {
        hosts[this_host].n_infective_chronic--;
      }
      
      // add to new counts
      if (hosts[this_host].innoc_next_infective_status[this_slot] == Acute) {
        hosts[this_host].n_infective_acute++;
      }
      if (hosts[this_host].innoc_next_infective_status[this_slot] == Chronic) {
        hosts[this_host].n_infective_chronic++;
      }
      
      // step forward infective status
      hosts[this_host].innoc_infective_status[this_slot] = hosts[this_host].innoc_next_infective_status[this_slot];
      
      // if newly infective then add to infectives list
      bool infective_after = (hosts[this_host].get_n_infective() > 0);
      if (!infective_before && infective_after) {
        host_infective_vec[this_deme].push_back(this_host);
      }
      
      // if no longer infective then drop from infectives list
      if (infective_before && !infective_after) {
        host_infective_vec[this_deme].erase(remove(host_infective_vec[this_deme].begin(), host_infective_vec[this_deme].end(), this_host));
      }
      
    }
    
    
    //-------- STORE RESULTS --------
    
    // store daily counts etc.
    if (output_daily_counts) {
      for (int k=0; k<n_demes; ++k) {
        H_store[k][t] = H[k];
        Sh_store[k][t] = Sh[k];
        Lh_store[k][t] = Lh[k];
        Ah_store[k][t] = Ah[k];
        Ch_store[k][t] = Ch[k];
        EIR_store[k][t] = EIR[k];
      }
    }
    
    // update complete age distributions
    if (output_age_distributions) {
      if (output_age_times[output_age_time_index] == t) {
        for (int k=0; k<n_demes; ++k) {
          
          // loop through all hosts
          for (int i=0; i<H[k]; ++i) {
            int this_host = host_vec[k][i];
            int age_days = t - hosts[this_host].birth_day;
            int age_years = floor(age_days/365.0);
            
            // update susceptible prevalence
            if (hosts[this_host].get_n_innoculations() == 0) {
              prev_Sh_age_store[k][output_age_time_index][age_years] += 1.0;
            }
            
            // update latent prevalence
            if (hosts[this_host].n_latent > 0) {
              prev_Lh_age_store[k][output_age_time_index][age_years] += 1.0;
            }
            
            // update acute prevalence
            if (hosts[this_host].n_acute > 0) {
              prev_Ah_age_store[k][output_age_time_index][age_years] += 1.0;
            }
            
            // update chronic prevalence
            if (hosts[this_host].n_chronic > 0) {
              prev_Ch_age_store[k][output_age_time_index][age_years] += 1.0;
            }
            
          }
          
          // divide by number of hosts per age bin
          for (int i=0; i<n_age; ++i) {
            prev_Sh_age_store[k][output_age_time_index][i] /= H_age_store[k][output_age_time_index][i];
            prev_Lh_age_store[k][output_age_time_index][i] /= H_age_store[k][output_age_time_index][i];
            prev_Ah_age_store[k][output_age_time_index][i] /= H_age_store[k][output_age_time_index][i];
            prev_Ch_age_store[k][output_age_time_index][i] /= H_age_store[k][output_age_time_index][i];
            
            inc_Lh_age_store[k][output_age_time_index][i] /= H_age_store[k][output_age_time_index][i];
            inc_Ah_age_store[k][output_age_time_index][i] /= H_age_store[k][output_age_time_index][i];
          }
        }
        
        // update output_age_time_index
        if (output_age_time_index < (output_age_times.size()-1)) {
          output_age_time_index++;
        }
      }
    }
    
  } // end time loop
  
}

#endif

/*
//------------------------------------------------
// carry out migration
void Dispatcher::migrate() {
  
  // schedule hosts to move from deme k1 to deme k2
  for (int k1=0; k1<demes; k1++) {
    for (int k2=0; k2<demes; k2++) {
      if (k1==k2 || delta_mig[k1][k2]==0) {
        continue;
      }
      
      // loop through all migration events
      for (int i=0; i<delta_mig[k1][k2]; i++) {
        
        // calculate probability migration event is in non-infective vs. infective host
        int n_noninf = hosts_noninfective[k1].size();
        int n_inf = hosts_infective[k1].size();
        double prob_h_migration_noninf = n_noninf/double(n_noninf+n_inf);  // proportion of migrations in non-infetive hosts
        
        // migration in either non-infective or infective
        int host_ID;
        if (rbernoulli1(prob_h_migration_noninf)) { // schedule non-infective to move
          
          int rnd1 = sample2(0, hosts_noninfective[k1].size()-1);
          host_ID = hosts_noninfective[k1][rnd1];
          mig_noninf_hosts[k1][k2].push_back(host_ID);
          hosts_noninfective[k1].erase(hosts_noninfective[k1].begin()+rnd1);
          
        } else {  // schedule infective to move
          
          int rnd1 = sample2(0, hosts_infective[k1].size()-1);
          host_ID = hosts_infective[k1][rnd1];
          mig_inf_hosts[k1][k2].push_back(host_ID);
          hosts_infective[k1].erase(hosts_infective[k1].begin()+rnd1);
        }
        
        // if host infected then add migration event to infection history
        if (output_infection_history) {
          int nb = hosts[host_ID].n_latent + hosts[host_ID].n_bloodstage + hosts[host_ID].n_infective;
          if (nb>0) {
            vector<int> tmp = {host_ID, hosts[host_ID].deme, k2};
            push_back_multiple(history_migration, tmp);
          }
        }
        
        // change deme of migrant
        hosts[host_ID].deme = k2;
      }
    }
  }
  // update non-infective and infective lists with new hosts
  for (int k1=0; k1<demes; k1++) {
    for (int k2=0; k2<demes; k2++) {
      if (k1==k2 || delta_mig[k1][k2]==0) {
        continue;
      }
      push_back_multiple(hosts_noninfective[k2], mig_noninf_hosts[k1][k2]);
      push_back_multiple(hosts_infective[k2], mig_inf_hosts[k1][k2]);
      
      mig_noninf_hosts[k1][k2].clear();
      mig_inf_hosts[k1][k2].clear();
    }
  }
  
}
*/
