
#include "ross_macdonald.h"
#include "probability.h"
#include "misc.h"

using namespace std;

//------------------------------------------------
// Draws from synchronous stochastic Ross-Macdonald model
// [[Rcpp::export]]
Rcpp::List ross_macdonald_cpp(Rcpp::List args) {
  
  // define parameters
  int max_time = Rcpp::as<int>(args["max_time"]);
  double a = Rcpp::as<double>(args["a"]);
  double mu = Rcpp::as<double>(args["mu"]);
  int u = Rcpp::as<int>(args["u"]);
  int v = Rcpp::as<int>(args["v"]);
  double r = Rcpp::as<double>(args["r"]);
  double b = Rcpp::as<double>(args["b"]);
  double c = Rcpp::as<double>(args["c"]);
  vector<int> Ih_init = Rcpp::as<vector<int>>(args["Ih_init"]);
  vector<int> Iv_init = Rcpp::as<vector<int>>(args["Iv_init"]);
  vector<int> H = Rcpp::as<vector<int>>(args["H"]);
  vector<int> M = Rcpp::as<vector<int>>(args["M"]);
  vector<vector<double>> mig = Rcpp::as<vector<vector<double>>>(args["migration_matrix"]);
  int demes = Rcpp::as<int>(args["demes"]);
  
  // store duration of all infections. infection_time is a vector over demes, 
  // then over times. First element holds time of infection, second element 
  // holds deme of infection. These values can eventually be used to calculate
  // the duration of infection
  vector<vector<pair<int, int>>> infection_time(demes);
  vector<vector<vector<int>>> durations(demes, vector<vector<int>>(max_time));
  
  // define and initialise model states
  vector<int> Sh(demes);
  vector<int> Eh(demes);
  vector<int> Ih = Ih_init;
  vector<int> Sv(demes);
  vector<int> Ev(demes);
  vector<int> Iv = Iv_init;
  for (int k=0; k<demes; k++) {
    Sh[k] = H[k] - Ih_init[k];
    Sv[k] = M[k] - Iv_init[k];
    for (int j=0; j<Ih_init[k]; j++) {
      infection_time[k].push_back(make_pair(0,k));
    }
  }
  
  // vectors for storing model states in all demes, at all points in time. Note
  // that -1 acts as an indicator that these values should be replaced with NA
  // once we are back in R
  vector<vector<double>> Sh_store(demes, vector<double>(max_time, -1));
  vector<vector<double>> Eh_store(demes, vector<double>(max_time, -1));
  vector<vector<double>> Ih_store(demes, vector<double>(max_time, -1));
  vector<vector<double>> Sv_store(demes, vector<double>(max_time, -1));
  vector<vector<double>> Ev_store(demes, vector<double>(max_time, -1));
  vector<vector<double>> Iv_store(demes, vector<double>(max_time, -1));
  for (int k=0; k<demes; k++) {
    Sh_store[k][0] = Sh[k];
    Eh_store[k][0] = Eh[k];
    Ih_store[k][0] = Ih[k];
    Sv_store[k][0] = Sv[k];
    Ev_store[k][0] = Ev[k];
    Iv_store[k][0] = Iv[k];
  }
  
  // create vectors to keep track of lag states. Eh_list and Ev_list are 
  // ring-buffers that store the number of agents in the lag state for each 
  // number of days back in time. u_buffer_index is the current time, and 
  // u_buffer_index_delay is u steps back in time. Note that
  // u_buffer_index_delay is actually u steps BEHIND u_buffer_index, but due to
  // the ring-buffer looping round this actually places it one step in front of
  // u_buffer_index at all times.
  vector<vector<int>> Eh_list(demes, vector<int>(u+1));
  vector<vector<int>> Ev_list(demes, vector<int>(v+1));
  int u_buffer_index = 0;
  int u_buffer_index_delay = 1;
  int v_buffer_index = 0;
  
  // initialise rates and probabilities
  double rate_h_infection, rate_v_infection;
  double prob_h_infection, prob_v_infection_or_death, prob_v_infection;
  double prob_h_recovery = 1 - exp(-r);
  double prob_v_death = 1 - exp(-mu);
  
  // initialise objects for implementing migration
  vector<vector<int>> mig_events(demes, vector<int>(demes));
  vector<vector<vector<int>>> mig_events_store(max_time, mig_events);
  vector<vector<pair<int, int>>> mig_infection_time(demes);
  
  // carry out simulation loop
  for (int i=1; i<max_time; i++) {
    
    // update ring buffer indices
    u_buffer_index = (u_buffer_index==u) ? 0 : u_buffer_index+1;
    u_buffer_index_delay = (u_buffer_index_delay==u) ? 0 : u_buffer_index_delay+1;
    v_buffer_index = (v_buffer_index==v) ? 0 : v_buffer_index+1;
    
    // implement migration
    // clear temporary object for storing infection times of migrating hosts
    for (int k=0; k<demes; k++) {
      mig_infection_time[k].clear();
    }
    // schedule human hosts to be moved
    for (int k1=0; k1<demes; k1++) {
      mig_events[k1] = rmultinom1(Ih[k1], mig[k1]);
      mig_events[k1][k1] = 0;
      
      // move this many infection times to the temporary object
      for (int k2=0; k2<demes; k2++) {
        if (k2==k1) {
          continue;
        }
        for (int j=0; j<mig_events[k1][k2]; j++) {
          int rnd1 = sample2(0, infection_time[k1].size()-1);
          mig_infection_time[k2].push_back(infection_time[k1][rnd1]);
          infection_time[k1].erase(infection_time[k1].begin() + rnd1);
        }
      }
    }
    mig_events_store[i] = mig_events;
    // move human hosts
    for (int k1=0; k1<demes; k1++) {
      for (int k2=0; k2<demes; k2++) {
        Ih[k1] += mig_events[k2][k1] - mig_events[k1][k2];
      }
      infection_time[k1].insert(infection_time[k1].end(), mig_infection_time[k1].begin(), mig_infection_time[k1].end());
    }
    
    // loop through all demes
    for (int k=0; k<demes; k++) {
      
      // calculate rates of events
      rate_h_infection = a*b*Iv[k]/double(H[k]); // human infection (move to Eh state)
      rate_v_infection = a*c*Ih[k]/double(H[k]); // mosquito infection (move to Ev state)
      
      // convert to probabilities, allowing for competing hazards
      prob_h_infection = 1 - exp(-rate_h_infection); // human infection (move to Eh state)
      prob_v_infection_or_death = 1 - exp(-(rate_v_infection + mu)); // mosquito infection or death in Sv state (competing hazards)
      prob_v_infection = rate_v_infection/(rate_v_infection + mu); // relative rate of mosquito infection
      
      // human events
      int h_infection = rbinom1(Sh[k], prob_h_infection); // human infection (move to Eh state)
      int h_recovery = rbinom1(Ih[k], prob_h_recovery); // human recovery
      Sh[k] += -h_infection + h_recovery; // update Sh
      
      Eh_list[k][u_buffer_index] = h_infection; // add infecteds to Eh list
      Eh[k] += Eh_list[k][u_buffer_index] - Eh_list[k][u_buffer_index_delay]; // update Eh
      
      Ih[k] += Eh_list[k][u_buffer_index_delay] - h_recovery; // update Ih
      
      // update infection times. Add new infections and drop those that recover
      for (int j=0; j<h_infection; j++) {
        infection_time[k].push_back(make_pair(i,k));
      }
      for (int j=0; j<h_recovery; j++) {
        int rnd1 = sample2(0, infection_time[k].size()-1);
        int time_interval = i - infection_time[k][rnd1].first;
        durations[infection_time[k][rnd1].second][infection_time[k][rnd1].first].push_back(time_interval);
        infection_time[k].erase(infection_time[k].begin() + rnd1);
      }
      
      // mosquito events
      int v_infection_or_death = rbinom1(Sv[k], prob_v_infection_or_death); // mosquito infection or death in Sv state (competing hazards)
      int v_infection = rbinom1(v_infection_or_death, prob_v_infection); // mosquito infection
      int v_death_Iv = rbinom1(Iv[k], prob_v_death); // mosquito death in Iv state
      Sv[k] += -v_infection + v_death_Iv; // update Sv
      
      Ev_list[k][v_buffer_index] = v_infection; // add infecteds to Ev list
      Ev[k] += v_infection; // add new infecteds to Ev
      
      // deaths in Ev state
      int j2 = v_buffer_index;
      for (int j=0; j<v; j++) { // loop through all previous entries in Ev list
        j2 = (j2==0) ? v : j2-1;
        if (Ev_list[k][j2]>0) {
          int v_death_Ev = rbinom1(Ev_list[k][j2], prob_v_death); // mosquito death in this Ev state
          Ev_list[k][j2] -= v_death_Ev; // subtract mosquito deaths from this element
          Ev[k] -= v_death_Ev; // subtract mosquito deaths from Ev counter
          Sv[k] += v_death_Ev; // respawn mosquitoes in Sv state
        }
      }
      
      // at this stage j2 is behind v_buffer_index by v steps. Move final Ev
      // list entries from Ev to Iv
      Ev[k] -= Ev_list[k][j2];
      Iv[k] += Ev_list[k][j2] - v_death_Iv;
      
      // store values
      Sh_store[k][i] = Sh[k];
      Eh_store[k][i] = Eh[k];
      Ih_store[k][i] = Ih[k];
      Sv_store[k][i] = Sv[k];
      Ev_store[k][i] = Ev[k];
      Iv_store[k][i] = Iv[k];
      
    } // end loop over demes
  } // end loop over time
  
  // finalise durations
  for (int k=0; k<demes; k++) {
    for (int j=0; j<int(infection_time[k].size()); j++) {
      durations[infection_time[k][j].second][infection_time[k][j].first].push_back(-1);
    }
  }
  
  // return values
  return Rcpp::List::create(Rcpp::Named("Sh")=Sh_store,
                            Rcpp::Named("Eh")=Eh_store,
                            Rcpp::Named("Ih")=Ih_store,
                            Rcpp::Named("Sv")=Sv_store,
                            Rcpp::Named("Ev")=Ev_store,
                            Rcpp::Named("Iv")=Iv_store,
                            Rcpp::Named("durations")=durations,
                            Rcpp::Named("migrations")=mig_events_store);
}
