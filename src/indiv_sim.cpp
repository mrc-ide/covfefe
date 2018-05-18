
#include "indiv_sim.h"
#include "indiv_host.h"
#include "probability.h"
#include "misc.h"

using namespace std;

//------------------------------------------------
// Draws from individual-based model
// [[Rcpp::export]]
Rcpp::List indiv_sim_cpp(Rcpp::List args) {
  
  // extract arguments
  int max_time = rcpp_to_int(args["max_time"]);
  double a = rcpp_to_double(args["a"]);
  double mu = rcpp_to_double(args["mu"]);
  int u = rcpp_to_int(args["u"]);
  int v = rcpp_to_int(args["v"]);
  int g = rcpp_to_int(args["g"]);
  double r = rcpp_to_double(args["r"]);
  double b = rcpp_to_double(args["b"]);
  double c = rcpp_to_double(args["c"]);
  vector<int> Ih_init = rcpp_to_vector_int(args["Ih_init"]);
  vector<int> H = rcpp_to_vector_int(args["H"]);
  vector<int> M = rcpp_to_vector_int(args["M"]);
  int max_infections = rcpp_to_int(args["max_infections"]);
  vector<vector<double>> delta_mig = rcpp_to_matrix_double(args["delta_mig"]);
  int demes = rcpp_to_int(args["demes"]);
  
  // define probability of human recovery and mosquito death
  double prob_h_recovery = 1 - exp(-r);
  double prob_v_death = 1 - exp(-mu);
  
  // initialise infection history list for storing complete history of events
  vector<int> history_migration;
  vector<int> history_infection;
  vector<int> history_bloodstage;
  vector<int> history_recover;
  vector<vector<vector<int>>> history_full(max_time);
  
  // initialise daily_counts for storing the number of hosts with each possible
  // number of blood stage infections (from 0 up to max_infections) at each time
  // step and in each deme. At time 0 all hosts have 0 blood stage infections.
  vector<vector<vector<int>>> daily_counts(demes, vector<vector<int>>(max_time, vector<int>(max_infections+1)));
  for (int k=0; k<demes; k++) {
    daily_counts[k][0][0] = H[k];
  }
  
  // initialise objects for scheduling future events
  vector<vector<pair<int, int>>> schedule_bloodstage(max_time); // first=host ID, second=infection ID
  vector<vector<pair<int, int>>> schedule_infective(max_time); // first=host ID, second=infection ID
  vector<vector<tuple<int, int, int>>> schedule_recover(max_time);  // first=host ID, second=infection ID, third=infection is in infective stage? (0=false, 1=true)
  
  // create population of human hosts. Create objects for indicating which hosts
  // are infective vs. non-infective in each deme. This is useful because only 
  // infective hosts are relevant to the mosquito population. All hosts start
  // out non-infective.
  vector<indiv_host> host(sum(H));
  vector<vector<int>> host_noninfective(demes);
  vector<vector<int>> host_infective(demes);
  int j = 0;
  for (int k=0; k<demes; k++) {
    host_noninfective[k] = vector<int>(H[k]);
    for (int i=0; i<H[k]; i++) {
      host[j].deme = k;
      host_noninfective[k][i] = j++;
    }
  }
  
  // seed initial infections
  for (int k=0; k<demes; k++) {
    for (int i=0; i<Ih_init[k]; i++) {
      int host_ID = host_noninfective[k][i];
      int inf_ID = 0; // start with infection ID 0
      
      // add to infection list
      vector<int> tmp = {host_ID, k, inf_ID, -1, -1, -1}; // the -1 values indicate that this is a de-novo infection, i.e. there is no source host
      push_back_multiple(history_infection, tmp);
      
      // schedule blood stage for u steps ahead
      schedule_bloodstage[u].push_back({host_ID, inf_ID});
      
      // record infection in host
      host[host_ID].total_infections++;
      host[host_ID].n_latent++;
    }
  }
  
  // add first entries to infection history
  history_full[0].push_back(history_migration);
  history_full[0].push_back(history_infection);
  history_full[0].push_back(history_bloodstage);
  history_full[0].push_back(history_recover);
  
  // initialise objects representing mosquitoes. An infective mosquito is 
  // represented as a tuple of three values: 1) the ID of the host that infected
  // it, 2) the number of observable infections in that host (i.e. blood stage
  // or infective stage), 3) the time the mosquito became infected by that host.
  // If the infective mosquito ever bites a new host and passes on the
  // infection, these values will make it into the infection history.
  //
  // Typically a large number of mosquitoes become infected each time step, but 
  // most of them die before making it through the lag phase (the extrinsic 
  // incubation period). To avoid creating a large number of infective 
  // mosquitoes that die before they could pass on infection, we instead draw
  // the number of infective mosquitoes that die at each day of the lag phase,
  // and only those mosquitoes that make it all the way through are instantiated
  // as tuples. n_Ev_death stores the number of scheduled deaths at each day
  // going forward until day u. These deaths will be added back into the
  // susceptible pool on that day. Likewise, Ev stores the full tuples of
  // infective mosquitoes that are due to emerge from the lag phase at each day
  // going forward. These infective mosquitoes will be added to Iv on that day.
  // Finally, both n_Ev_death and Ev are implemented as ring-buffers, meaning
  // they wrap around. v_ringbuffer_thistime stores the current day, which wraps
  // back round to 0 once it hits u days.
  vector<int> n_Sv = M; // number of susceptible mosquitoes in each deme (all mosquitoes initially susceptible)
  vector<vector<int>> n_Ev_death(demes, vector<int>(u));
  vector<vector<vector<tuple<int, int, int>>>> Ev(demes, vector<vector<tuple<int, int, int>>>(u+1));
  vector<vector<tuple<int, int, int>>> Iv(demes);
  int v_ringbuffer_thistime = 0;
  
  // initialise objects for implementing migration
  vector<vector<vector<int>>> mig_noninf_hosts(demes, vector<vector<int>>(demes));  // to store IDs of non-infective hosts that will move demes in a given time step
  vector<vector<vector<int>>> mig_inf_hosts(demes, vector<vector<int>>(demes));  // to store IDs of infective hosts that will move demes in a given time step
  
  
  
  // SIMULATION ############################################################
  
  for (int t=1; t<max_time; t++) {
    
    // clear history lists
    history_migration.clear();
    history_infection.clear();
    history_bloodstage.clear();
    history_recover.clear();
    
    // update ring buffer indices
    v_ringbuffer_thistime = (v_ringbuffer_thistime==v) ? 0 : v_ringbuffer_thistime+1;
    
    
    //#### MIGRATION
    // schedule hosts to move from deme k1 to deme k2
    for (int k1=0; k1<demes; k1++) {
      for (int k2=0; k2<demes; k2++) {
        if (k1==k2 || delta_mig[k1][k2]==0) {
          continue;
        }
        
        // loop through all migration events
        for (int i=0; i<delta_mig[k1][k2]; i++) {
          
          // calculate probability migration event is in non-infective vs. infective host
          int n_noninf = host_noninfective[k1].size();
          int n_inf = host_infective[k1].size();
          double prob_h_migration_noninf = n_noninf/double(n_noninf+n_inf);  // proportion of migrations in non-infetive hosts
          
          // migration in either non-infective or infective
          int host_ID;
          if (rbernoulli1(prob_h_migration_noninf)) { // schedule non-infective to move
            
            int rnd1 = sample2(0, host_noninfective[k1].size()-1);
            host_ID = host_noninfective[k1][rnd1];
            mig_noninf_hosts[k1][k2].push_back(host_ID);
            host_noninfective[k1].erase(host_noninfective[k1].begin()+rnd1);
            
          } else {  // schedule infective to move
            
            int rnd1 = sample2(0, host_infective[k1].size()-1);
            host_ID = host_infective[k1][rnd1];
            mig_inf_hosts[k1][k2].push_back(host_ID);
            host_infective[k1].erase(host_infective[k1].begin()+rnd1);
          }
          
          // if host infected then add migration event to infection history
          int nb = host[host_ID].n_latent + host[host_ID].n_bloodstage + host[host_ID].n_infective;
          if (nb>0) {
            vector<int> tmp = {host_ID, host[host_ID].deme, k2};
            push_back_multiple(history_migration, tmp);
          }
          
          // change deme of migrant
          host[host_ID].deme = k2;
        }
      }
    }
    // update non-infective and infective lists with new hosts
    for (int k1=0; k1<demes; k1++) {
      for (int k2=0; k2<demes; k2++) {
        if (k1==k2 || delta_mig[k1][k2]==0) {
          continue;
        }
        push_back_multiple(host_noninfective[k2], mig_noninf_hosts[k1][k2]);
        push_back_multiple(host_infective[k2], mig_inf_hosts[k1][k2]);
        
        mig_noninf_hosts[k1][k2].clear();
        mig_inf_hosts[k1][k2].clear();
      }
    }
    
    
    // loop through all demes
    for (int k=0; k<demes; k++) {
      
      //#### MOSQUITO EVENTS
      // move Ev into Iv
      push_back_multiple(Iv[k], Ev[k][v_ringbuffer_thistime]);
      Ev[k][v_ringbuffer_thistime].clear();
      
      // deaths in Iv
      int v_death_Iv = rbinom1(Iv[k].size(), prob_v_death);
      for (int i=0; i<v_death_Iv; i++) {
        int rnd1 = sample2(0,Iv[k].size()-1);
        Iv[k].erase(Iv[k].begin()+rnd1);
      }
      
      // draw number of new infections
      double rate_v_infection = a*c*host_infective[k].size()/double(H[k]); // rate of mosquito infection
      double prob_v_infection_or_death = 1 - exp(-(rate_v_infection + mu)); // mosquito infection or death in Sv state (competing hazards)
      double prob_v_infection = rate_v_infection/(rate_v_infection + mu); // relative rate of mosquito infection vs. death
      int v_infection_or_death = rbinom1(n_Sv[k], prob_v_infection_or_death); // number mosquito infection or death in Sv state
      int v_infection = rbinom1(v_infection_or_death, prob_v_infection); // number mosquito infection
      
      // update n_Sv with events so far
      n_Sv[k] += -v_infection + n_Ev_death[k][v_ringbuffer_thistime] + v_death_Iv; // update Sv
      n_Ev_death[k][v_ringbuffer_thistime] = 0; // once these lag stage deaths have been moved into Sv, reset this entry
      
      // the majority of new mosquito infections will die in lag phase. Schedule
      // these deaths to move back into Sv in future steps
      if (v_infection>0) {
        
        // schedule deaths in Ev state
        int j2 = v_ringbuffer_thistime;
        for (int j=0; j<v; j++) {
          j2 = (j2==(v-1)) ? 0 : j2++;  // wrap around ring buffer
          int Ev_j_deaths = rbinom1(v_infection, prob_v_death);
          n_Ev_death[k][j2] += Ev_j_deaths;
          v_infection -= Ev_j_deaths;
          if (v_infection==0) {
            break;
          }
        }
        
        // add remaining infections to Ev
        for (int i=0; i<v_infection; i++) {
          int rnd1 = sample2(0,host_infective[k].size()-1);
          int this_host = host_infective[k][rnd1];
          int this_n = host[this_host].n_bloodstage + host[this_host].n_infective;
          Ev[k][v_ringbuffer_thistime].push_back({this_host, this_n, t});
        }
      }
      
      
      //#### HUMAN EVENTS
      // get number of new infections
      double rate_h_infection = a*b*Iv[k].size()/double(H[k]);
      double prob_h_infection = 1 - exp(-rate_h_infection);
      int h_infection = rbinom1(H[k], prob_h_infection);  // total number new infections
      double prob_h_infection_noninf = host_noninfective[k].size()/double(H[k]);  // proportion that occur in non-infetive hosts
      
      // apply new infections
      for (int i=0; i<h_infection; i++) {
        
        // choose host at random
        int host_ID;
        if (rbernoulli1(prob_h_infection_noninf)) {
          int rnd1 = sample2(0, host_noninfective[k].size()-1);
          host_ID = host_noninfective[k][rnd1];
        } else {
          int rnd1 = sample2(0, host_infective[k].size()-1);
          host_ID = host_infective[k][rnd1];
        }
        
        // skip if reached max infections
        int nb = host[host_ID].n_latent + host[host_ID].n_bloodstage + host[host_ID].n_infective;
        if (nb==max_infections) {
          continue;
        }
        
        // choose mosquito at random
        int rnd2 = sample2(0, Iv[k].size()-1);
        
        // add to infection list
        vector<int> tmp = {host_ID, k, host[host_ID].total_infections, get<0>(Iv[k][rnd2]), get<1>(Iv[k][rnd2]), get<2>(Iv[k][rnd2])};
        push_back_multiple(history_infection, tmp);
        
        // schedule move to blood stage
        if ((t+u)<max_time) {
          schedule_bloodstage[t+u].push_back({host_ID, host[host_ID].total_infections});
        }
        
        // infection
        host[host_ID].total_infections++;
        host[host_ID].n_latent++;
        
      }
      
    } // end loop through demes
    
    
    //#### IMPLEMENT SCHEDULED EVENTS
    // move to blood stage
    for (int i=0; i<int(schedule_bloodstage[t].size()); i++) {
      int host_ID = schedule_bloodstage[t][i].first;
      int inf_ID = schedule_bloodstage[t][i].second;
      int host_deme = host[host_ID].deme;
      
      // add to bloodstage list
      vector<int> tmp = {host_ID, host_deme, inf_ID};
      push_back_multiple(history_bloodstage, tmp);
      
      // new blood stage infection
      host[host_ID].n_latent--;
      host[host_ID].n_bloodstage++;
      
      // schedule infective and/or recovery
      int dur = rgeom1(prob_h_recovery);
      if (dur>g) {  // if become infective prior to recovery
        
        if ((t+dur)<max_time) { // both infective and recovery occur within max_time
          schedule_infective[t+g].push_back({host_ID, inf_ID});
          schedule_recover[t+dur].push_back({host_ID, inf_ID, 1});
        } else if ((t+g)<max_time) {  // infective but not recovery occur within max_time
          schedule_infective[t+g].push_back({host_ID, inf_ID});
        }
        
      } else {  // if recover prior to become infective
        if ((t+dur)<max_time) { // recovery occurs within max_time
          schedule_recover[t+dur].push_back({host_ID, inf_ID, 0});
        }
      }
    }
    
    // blood stage become infective
    for (int i=0; i<int(schedule_infective[t].size()); i++) {
      int host_ID = schedule_infective[t][i].first;
      
      // move to infective
      host[host_ID].n_bloodstage--;
      host[host_ID].n_infective++;
      
      // move from noninfective to infective list
      if (host[host_ID].n_infective==1) {
        int host_deme = host[host_ID].deme;
        host_noninfective[host_deme].erase(remove(host_noninfective[host_deme].begin(), host_noninfective[host_deme].end(), host_ID));
        host_infective[host_deme].push_back(host_ID);
      }
    }
    
    // recovery
    for (int i=0; i<int(schedule_recover[t].size()); i++) {
      int host_ID = get<0>(schedule_recover[t][i]);
      int inf_ID = get<1>(schedule_recover[t][i]);
      int recovery_type = get<2>(schedule_recover[t][i]);
      int host_deme = host[host_ID].deme;
      
      // add to recover list
      vector<int> tmp = {host_ID, host_deme, inf_ID};
      push_back_multiple(history_recover, tmp);
      
      // type 0 = blood stage recovery, type 1 = infective recovery
      if (recovery_type==0) {
        host[host_ID].n_bloodstage--;
      } else {
        host[host_ID].n_infective--;
        
        // move from infective to noninfective list
        if (host[host_ID].n_infective==0) {
          host_infective[host_deme].erase(remove(host_infective[host_deme].begin(), host_infective[host_deme].end(), host_ID));
          host_noninfective[host_deme].push_back(host_ID);
        }
      }
    }
    
    //#### STORE RESULTS
    for (int k=0; k<demes; k++) {
      int s = host_infective[k].size();
      for (int i=0; i<H[k]; i++) {
        int host_ID;
        if (i<s) {
          host_ID = host_infective[k][i];
        } else {
          host_ID = host_noninfective[k][i-s];
        }
        int nb = host[host_ID].n_bloodstage + host[host_ID].n_infective;
        if (nb<=max_infections) {
          daily_counts[k][t][nb]++;
        }
      }
    }
    
    // add to infection history
    history_full[t].push_back(history_migration);
    history_full[t].push_back(history_infection);
    history_full[t].push_back(history_bloodstage);
    history_full[t].push_back(history_recover);
    
  } // end simulation loop
  
  // return as list
  return Rcpp::List::create(Rcpp::Named("daily_counts") = daily_counts,
                            Rcpp::Named("infection_history") = history_full);
}

