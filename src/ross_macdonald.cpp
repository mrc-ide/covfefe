
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
  int demes = Rcpp::as<int>(args["demes"]);
  
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
  
  // carry out simulation loop
  for (int i=1; i<max_time; i++) {
    
    // update ring buffer indices
    u_buffer_index = (u_buffer_index==u) ? 0 : u_buffer_index+1;
    u_buffer_index_delay = (u_buffer_index_delay==u) ? 0 : u_buffer_index_delay+1;
    v_buffer_index = (v_buffer_index==v) ? 0 : v_buffer_index+1;
    
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
  
  // return values
  return Rcpp::List::create(Rcpp::Named("Sh")=Sh_store,
                            Rcpp::Named("Eh")=Eh_store,
                            Rcpp::Named("Ih")=Ih_store,
                            Rcpp::Named("Sv")=Sv_store,
                            Rcpp::Named("Ev")=Ev_store,
                            Rcpp::Named("Iv")=Iv_store);
}
