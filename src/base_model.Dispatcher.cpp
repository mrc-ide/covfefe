
#include "base_model.Dispatcher.h"
#include "probability.h"

using namespace std;


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
  schedule_infective_start_acute = vector<vector<pair<int, int>>>(max_time+1);
  schedule_infective_start_chronic = vector<vector<pair<int, int>>>(max_time+1);
  schedule_infective_acute_chronic = vector<vector<pair<int, int>>>(max_time+1);
  schedule_infective_stop_acute = vector<vector<pair<int, int>>>(max_time+1);
  schedule_infective_stop_chronic = vector<vector<pair<int, int>>>(max_time+1);
  
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
      new_infection(this_host, null_mosquito, 0);
      
      // update prob_infection_index
      if (hosts[this_host].prob_infection_index < (n_prob_infection-1)) {
        hosts[this_host].prob_infection_index++;
      }
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
  Lv = vector<int>(n_demes);
  Iv = vector<int>(n_demes);
  
  // objects for storing daily counts etc.
  if (output_daily_counts) {
    H_store = vector<vector<int>>(n_demes, vector<int>(max_time+1));
    Sh_store = vector<vector<int>>(n_demes, vector<int>(max_time+1));
    Lh_store = vector<vector<int>>(n_demes, vector<int>(max_time+1));
    Ah_store = vector<vector<int>>(n_demes, vector<int>(max_time+1));
    Ch_store = vector<vector<int>>(n_demes, vector<int>(max_time+1));
    Sv_store = vector<vector<int>>(n_demes, vector<int>(max_time+1));
    Lv_store = vector<vector<int>>(n_demes, vector<int>(max_time+1));
    Iv_store = vector<vector<int>>(n_demes, vector<int>(max_time+1));
    EIR_store = vector<vector<double>>(n_demes, vector<double>(max_time+1));
    
    // counts at time 0
    for (int k=0; k<n_demes; ++k) {
      H_store[k][0] = H[k];
      Sh_store[k][0] = Sh[k];
      Lh_store[k][0] = Lh[k];
      Sv_store[k][0] = Sv[k];
    }
  }
  
  // initialise objects for storing age distributions
  if (output_age_distributions) {
    vector<vector<vector<int>>> tmp_mat_int = vector<vector<vector<int>>>(n_demes, vector<vector<int>>(n_output_age_times, vector<int>(n_age)));
    vector<vector<vector<double>>> tmp_mat_double = vector<vector<vector<double>>>(n_demes, vector<vector<double>>(n_output_age_times, vector<double>(n_age)));
    H_age_store = tmp_mat_int;
    Sh_age_store = tmp_mat_int;
    Lh_age_store = tmp_mat_int;
    Ah_age_store = tmp_mat_int;
    Ch_age_store = tmp_mat_int;
    inc_Lh_age_store = tmp_mat_double;
    inc_Ah_age_store = tmp_mat_double;
  }
  
  // misc
  prob_h_infectious_bite = 0;
  
}

//------------------------------------------------
// human new infection
void Dispatcher::new_infection(int this_host, Mosquito &mosq, int t) {
  
  // update cumulative innoculations, irrespective of whether this innoculation
  // passes the next series of checks
  hosts[this_host].cumulative_n_innoculations++;
  
  // get next free innoculation slot. If none free then return without creating
  // new innoculation
  int this_slot = 0;
  for (int i=0; i<max_innoculations; ++i) {
    if (!hosts[this_host].innoc_active[i]) {
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
  hosts[this_host].innoc_active[this_slot] = true;
  hosts[this_host].innoc_status[this_slot] = Latent;
  hosts[this_host].innoc_status_prev_update_time[this_slot] = 0;
  hosts[this_host].innoc_status_next_update_time[this_slot] = t+u;
  hosts[this_host].innoc_infective_start_acute[this_slot] = 0;
  hosts[this_host].innoc_infective_start_chronic[this_slot] = 0;
  hosts[this_host].innoc_infective_acute_chronic[this_slot] = 0;
  hosts[this_host].innoc_infective_stop_acute[this_slot] = 0;
  hosts[this_host].innoc_infective_stop_chronic[this_slot] = 0;
  
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
// latent transition to acute infection
void Dispatcher::latent_acute(int this_host, int this_slot, int t) {
  
  // get host properties
  int this_deme = hosts[this_host].deme;
  int this_death_day = hosts[this_host].death_day;
  
  // update status
  hosts[this_host].innoc_status[this_slot] = Acute;
  
  // update status update time
  int duration = sampler_duration_acute[hosts[this_host].duration_acute_index].draw() + 1;
  hosts[this_host].innoc_status_prev_update_time[this_slot] = hosts[this_host].innoc_status_next_update_time[this_slot];
  hosts[this_host].innoc_status_next_update_time[this_slot] = t + duration;
  
  // update duration_acute_index
  if (hosts[this_host].duration_acute_index < (n_duration_acute-1)) {
    hosts[this_host].duration_acute_index++;
  }
  
  // update infectivity_acute_index
  if (hosts[this_host].infectivity_acute_index < (n_infectivity_acute-1)) {
    hosts[this_host].infectivity_acute_index++;
  }
  
  // schedule status update
  if (t+duration <= max_time && t+duration < this_death_day) {
    schedule_status_update[t+duration].emplace_back(this_host, this_slot);
  }
  
  // update infective status
  hosts[this_host].innoc_infective_start_acute[this_slot] = t+g;
  
  // schedule infective update
  if (t+g <= max_time && t+g < this_death_day) {
    schedule_infective_start_acute[t+g].emplace_back(this_host, this_slot);
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
  
}

//------------------------------------------------
// latent transition to chronic infection
void Dispatcher::latent_chronic(int this_host, int this_slot, int t) {
  
  // get host properties
  int this_deme = hosts[this_host].deme;
  int this_death_day = hosts[this_host].death_day;
  
  // update status
  hosts[this_host].innoc_status[this_slot] = Chronic;
  
  // update status update time
  int duration = sampler_duration_chronic[hosts[this_host].duration_chronic_index].draw() + 1;
  hosts[this_host].innoc_status_prev_update_time[this_slot] = hosts[this_host].innoc_status_next_update_time[this_slot];
  hosts[this_host].innoc_status_next_update_time[this_slot] = t + duration;
  
  // update duration_chronic_index
  if (hosts[this_host].duration_chronic_index < (n_duration_chronic-1)) {
    hosts[this_host].duration_chronic_index++;
  }
  
  // update infectivity_chronic_index
  if (hosts[this_host].infectivity_chronic_index < (n_infectivity_chronic-1)) {
    hosts[this_host].infectivity_chronic_index++;
  }
  
  // schedule status update
  if (t+duration <= max_time && t+duration < this_death_day) {
    schedule_status_update[t+duration].emplace_back(this_host, this_slot);
  }
  
  // update infective status
  hosts[this_host].innoc_infective_start_chronic[this_slot] = t+g;
  
  // schedule infective update
  if (t+g <= max_time && t+g < this_death_day) {
    schedule_infective_start_chronic[t+g].emplace_back(this_host, this_slot);
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

//------------------------------------------------
// acute transition to chronic infection
void Dispatcher::acute_chronic(int this_host, int this_slot, int t) {
  
  // get host properties
  int this_deme = hosts[this_host].deme;
  int this_death_day = hosts[this_host].death_day;
  
  // update status
  hosts[this_host].innoc_status[this_slot] = Chronic;
  
  // update status update time
  int duration = sampler_duration_chronic[hosts[this_host].duration_chronic_index].draw() + 1;
  hosts[this_host].innoc_status_prev_update_time[this_slot] = hosts[this_host].innoc_status_next_update_time[this_slot];
  hosts[this_host].innoc_status_next_update_time[this_slot] = t + duration;
  
  // update duration_chronic_index
  if (hosts[this_host].duration_chronic_index < (n_duration_chronic-1)) {
    hosts[this_host].duration_chronic_index++;
  }
  
  // schedule status update
  if (t+duration <= max_time && t+duration < this_death_day) {
    schedule_status_update[t+duration].emplace_back(this_host, this_slot);
  }
  
  // update infective status
  hosts[this_host].innoc_infective_acute_chronic[this_slot] = t+g;
  
  // schedule infective update
  if (t+g <= max_time && t+g < this_death_day) {
    schedule_infective_acute_chronic[t+g].emplace_back(this_host, this_slot);
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
  
}

//------------------------------------------------
// acute transition to recovery
void Dispatcher::acute_recover(int this_host, int this_slot, int t) {
  
  // get host properties
  int this_deme = hosts[this_host].deme;
  int this_death_day = hosts[this_host].death_day;
  
  // update status
  hosts[this_host].innoc_status[this_slot] = Inactive;
  
  // update status update time (no more updates)
  hosts[this_host].innoc_status_prev_update_time[this_slot] = hosts[this_host].innoc_status_next_update_time[this_slot];
  hosts[this_host].innoc_status_next_update_time[this_slot] = 0;
  
  // update infective status
  hosts[this_host].innoc_infective_stop_acute[this_slot] = t+g;
  
  // schedule infective update
  if (t+g <= max_time && t+g < this_death_day) {
    schedule_infective_stop_acute[t+g].emplace_back(this_host, this_slot);
  }
  
  // update host counts
  hosts[this_host].n_acute--;
  
  // update deme counts
  if (hosts[this_host].get_n_asexual() == 0) {
    Sh[this_deme]++;
  }
  if (hosts[this_host].n_acute == 0) {
    Ah[this_deme]--;
  }
  
}

//------------------------------------------------
// chronic transition to recovery
void Dispatcher::chronic_recover(int this_host, int this_slot, int t) {
  
  // get host properties
  int this_deme = hosts[this_host].deme;
  int this_death_day = hosts[this_host].death_day;
  
  // update status
  hosts[this_host].innoc_status[this_slot] = Inactive;
  
  // update status update time (no more updates)
  hosts[this_host].innoc_status_prev_update_time[this_slot] = hosts[this_host].innoc_status_next_update_time[this_slot];
  hosts[this_host].innoc_status_next_update_time[this_slot] = 0;
  
  // update infective status
  hosts[this_host].innoc_infective_stop_chronic[this_slot] = t+g;
  
  // schedule infective update
  if (t+g <= max_time && t+g < this_death_day) {
    schedule_infective_stop_chronic[t+g].emplace_back(this_host, this_slot);
  }
  
  // update host counts
  hosts[this_host].n_chronic--;
  
  // update deme counts
  if (hosts[this_host].get_n_asexual() == 0) {
    Sh[this_deme]++;
  }
  if (hosts[this_host].n_chronic == 0) {
    Ch[this_deme]--;
  }
  
}

//------------------------------------------------
// update age incidence vectors
void Dispatcher::update_age_incidence(int this_host, int this_deme, int output_age_time_index, int t) {
  
  // get host properties
  int age_days = t - hosts[this_host].birth_day;
  int age_years = floor(age_days/365.0);
  double this_prob_infection = prob_infection[hosts[this_host].prob_infection_index];
  double this_prob_acute = prob_acute[hosts[this_host].prob_acute_index];
  
  // update raw number of hosts
  H_age_store[this_deme][output_age_time_index][age_years]++;
  
  // if host susceptible
  if (hosts[this_host].get_n_innoculations() == 0) {
    
    // update incidence of infection
    inc_Lh_age_store[this_deme][output_age_time_index][age_years] += prob_h_infectious_bite*this_prob_infection;
    
    // update incidence of acute disease
    inc_Ah_age_store[this_deme][output_age_time_index][age_years] += prob_h_infectious_bite*this_prob_infection*this_prob_acute;
     
  }
  
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
        int delta_Ev = int(Ev_mosq[k][ringtime].size());
        Lv[k] -= delta_Ev;
        Iv[k] += delta_Ev;
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
      
      // update counts with Lv and Iv deaths. Once lag stage deaths have been
      // moved into Sv, reset this entry
      Sv[k] += Ev_death[k][ringtime] + death_Iv;
      Lv[k] -= Ev_death[k][ringtime];
      Iv[k] -= death_Iv;
      Ev_death[k][ringtime] = 0;
      
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
        if (rbernoulli1(prob_infective)) {
          
          // update Sv and Lv
          Sv[k]--;
          Lv[k]++;
          
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
      prob_h_infectious_bite = 1 - exp(-EIR[k]);  // probability of new infectious bite per host
      int h_infectious_bite = rbinom1(H[k], prob_h_infectious_bite);  // total number of new infectious bites
      
      // update complete age distributions of incidence
      if (output_age_distributions) {
        if (output_age_times[output_age_time_index] == t) {
          for (int i=0; i<H[k]; ++i) {
            int this_host = host_vec[k][i];
            update_age_incidence(this_host, k, output_age_time_index, t);
          }
        }
      }
      
      // apply new infectious bites
      for (int i=0; i<h_infectious_bite; i++) {
        
        // choose host at random
        int rnd1 = sample2(0, H[k]-1);
        int this_host = host_vec[k][rnd1];
        
        // determine whether infectious bite is successful
        double host_prob_infection = prob_infection[hosts[this_host].prob_infection_index];
        if (rbernoulli1(host_prob_infection)) {
          
          // choose mosquito at random and carry out infection
          int rnd2 = sample2(0, int(Iv_mosq[k].size())-1);
          new_infection(this_host, Iv_mosq[k][rnd2], t);
        }
        
        // update b_index irrespective of whether infection successful
        if (hosts[this_host].prob_infection_index < (n_prob_infection-1)) {
          hosts[this_host].prob_infection_index++;
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
          
          // apply transition
          if (become_acute) {
            latent_acute(this_host, this_slot, t);
          } else {
            latent_chronic(this_host, this_slot, t);
          }
          
          break;
        } 
      
        // acute phase become chronic or recover
        case Acute: {
          
          // draw whether recovers or becomes chronic
          bool become_chronic = rbernoulli1(prob_AC);
          
          // split method between chronic disease and recovery
          if (become_chronic) {
            acute_chronic(this_host, this_slot, t);
          } else {
            acute_recover(this_host, this_slot, t);
          }
          
          break;
        }
        
        // chronic recover
        case Chronic: {
          chronic_recover(this_host, this_slot, t);
          break;
        }
        
        default:
          break;
      }  // end switch
    }  // end status updates
    
    // evaluate new acute infective
    for (auto it = schedule_infective_start_acute[t].begin(); it != schedule_infective_start_acute[t].end(); ++it) {
      int this_host = it->first;
      int this_deme = hosts[this_host].deme;
      
      // update counts
      hosts[this_host].n_infective_acute++;
      
      // if newly infective then add to infectives list
      if (hosts[this_host].get_n_infective() == 1) {
        host_infective_vec[this_deme].push_back(this_host);
      }
    }
    
    // evaluate new chronic infective
    for (auto it = schedule_infective_start_chronic[t].begin(); it != schedule_infective_start_chronic[t].end(); ++it) {
      int this_host = it->first;
      int this_deme = hosts[this_host].deme;
      
      // update counts
      hosts[this_host].n_infective_chronic++;
      
      // if newly infective then add to infectives list
      if (hosts[this_host].get_n_infective() == 1) {
        host_infective_vec[this_deme].push_back(this_host);
      }
    }
    
    // evaluate acute infective become chronic
    for (auto it = schedule_infective_acute_chronic[t].begin(); it != schedule_infective_acute_chronic[t].end(); ++it) {
      int this_host = it->first;
      
      // update counts
      hosts[this_host].n_infective_acute--;
      hosts[this_host].n_infective_chronic++;
    }
    
    // evaluate recovery of acute infective
    for (auto it = schedule_infective_stop_acute[t].begin(); it != schedule_infective_stop_acute[t].end(); ++it) {
      int this_host = it->first;
      int this_slot = it->second;
      int this_deme = hosts[this_host].deme;
      
      // update counts
      hosts[this_host].n_infective_acute--;
      
      // re-activate this slot
      hosts[this_host].innoc_active[this_slot] = false;
      
      // if no longer infective then drop from infectives list
      if (hosts[this_host].get_n_infective() == 0) {
        host_infective_vec[this_deme].erase(remove(host_infective_vec[this_deme].begin(), host_infective_vec[this_deme].end(), this_host));
      }
    }
    
    // evaluate recovery of chronic infective
    for (auto it = schedule_infective_stop_chronic[t].begin(); it != schedule_infective_stop_chronic[t].end(); ++it) {
      int this_host = it->first;
      int this_slot = it->second;
      int this_deme = hosts[this_host].deme;
      
      // update counts
      hosts[this_host].n_infective_chronic--;
      
      // re-activate this slot
      hosts[this_host].innoc_active[this_slot] = false;
      
      // if no longer infective then drop from infectives list
      if (hosts[this_host].get_n_infective() == 0) {
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
        Sv_store[k][t] = Sv[k];
        Lv_store[k][t] = Lv[k];
        Iv_store[k][t] = Iv[k];
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
              Sh_age_store[k][output_age_time_index][age_years]++;
            }
            
            // update latent prevalence
            if (hosts[this_host].n_latent > 0) {
              Lh_age_store[k][output_age_time_index][age_years]++;
            }
            
            // update acute prevalence
            if (hosts[this_host].n_acute > 0) {
              Ah_age_store[k][output_age_time_index][age_years]++;
            }
            
            // update chronic prevalence
            if (hosts[this_host].n_chronic > 0) {
              Ch_age_store[k][output_age_time_index][age_years]++;
            }
            
          }
          
          // divide incidence by number of hosts per age bin
          for (int i=0; i<n_age; ++i) {
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
