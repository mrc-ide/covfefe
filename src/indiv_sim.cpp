
#include "indiv_sim.h"
#include "indiv_host.h"
#include "probability.h"
#include "misc.h"

using namespace std;

//------------------------------------------------
// Draws from individual-based model
// [[Rcpp::export]]
Rcpp::List indiv_sim_cpp(Rcpp::List args) {
  
  // define parameters
  int max_time = Rcpp::as<int>(args["max_time"]);
  double a = Rcpp::as<double>(args["a"]);
  double mu = Rcpp::as<double>(args["mu"]);
  int u = Rcpp::as<int>(args["u"]);
  int v = Rcpp::as<int>(args["v"]);
  int g = Rcpp::as<int>(args["g"]);
  double r = Rcpp::as<double>(args["r"]);
  double b = Rcpp::as<double>(args["b"]);
  double c = Rcpp::as<double>(args["c"]);
  vector<int> Ih_init = Rcpp::as<vector<int>>(args["Ih_init"]);
  vector<int> H = Rcpp::as<vector<int>>(args["H"]);
  vector<int> M = Rcpp::as<vector<int>>(args["M"]);
  int max_clonality = Rcpp::as<int>(args["max_clonality"]);
  vector<vector<double>> delta_mig = Rcpp::as<vector<vector<double>>>(args["delta_mig"]);
  int demes = Rcpp::as<int>(args["demes"]);
  
  double prob_h_recovery = 1 - exp(-r);
  double prob_v_death = 1 - exp(-mu);
  
  // store results
  vector<vector<int>> line_list;
  vector<vector<vector<int>>> clonality(demes, vector<vector<int>>(max_time, vector<int>(max_clonality+1)));
  
  // create population of human hosts
  vector<vector<indiv_host>> h(demes);
  int j = 0;
  for (int k=0; k<demes; k++) {
    h[k] = vector<indiv_host>(H[k]);
    for (int i=0; i<H[k]; i++) {
      h[k][i].ID = j++;
    }
  }
  
  // seed initial infections
  for (int k=0; k<demes; k++) {
    for (int i=0; i<Ih_init[k]; i++) {
      h[k][i].new_infection(0, u + rgeom1(prob_h_recovery));
      line_list.push_back(vector<int>{0, h[k][i].ID, 1, k, -1});
    }
  }
  vector<unordered_set<int>> h_infectious(demes); // set of infectious hosts in each deme
  
  // create objects representing mosquitoes
  vector<int> n_Sv = M;
  vector<int> n_Iv(demes);
  vector<vector<int>> n_Ev(demes, vector<int>(u+1));
  vector<vector<vector<pair<int, int>>>> Ev(demes, vector<vector<pair<int, int>>>(u+1));
  vector<vector<pair<int, int>>> Iv(demes);
  int v_buffer_index = 0;
  
  // initialise objects for implementing migration
  vector<vector<vector<indiv_host>>> mig_hosts(demes, vector<vector<indiv_host>>(demes));
  
  // carry out simulation loop
  for (int t=0; t<max_time; t++) {
    
    // update ring buffer indices
    v_buffer_index = (v_buffer_index==v) ? 0 : v_buffer_index+1;
    
    // implement migration
    // schedule hosts to move
    for (int k1=0; k1<demes; k1++) {
      for (int k2=0; k2<demes; k2++) {
        if (k1==k2) {
          continue;
        }
        mig_hosts[k1][k2].clear();
        for (int i=0; i<delta_mig[k1][k2]; i++) {
          int rnd1 = sample2(0, h[k1].size()-1);
          mig_hosts[k1][k2].push_back(h[k1][rnd1]);
          h[k1].erase(h[k1].begin()+rnd1);
        }
      }
    }
    // move hosts
    for (int k1=0; k1<demes; k1++) {
      for (int k2=0; k2<demes; k2++) {
        push_back_multiple(h[k2], mig_hosts[k1][k2]);
      }
    }
    
    // loop through all demes
    for (int k=0; k<demes; k++) {
      
      // human events
      double rate_h_infection = a*b*n_Iv[k]/double(H[k]);
      double prob_h_infection = 1 - exp(-rate_h_infection);
      int h_infection = rbinom1(H[k], prob_h_infection);
      for (int i=0; i<h_infection; i++) {
        int rnd1 = sample2(0, H[k]-1);
        if (h[k][rnd1].n_infections<max_clonality) {
          int rnd2 = sample2(0, n_Iv[k]-1);
          h[k][rnd1].new_infection(t, t + u + rgeom1(prob_h_recovery));
          line_list.push_back(vector<int>{t, h[k][rnd1].ID, 2, Iv[k][rnd2].first, Iv[k][rnd2].second});
        }
      }
      
      // mosquito events
      int v_death_Iv = rbinom1(n_Iv[k], prob_v_death); // mosquito death in Iv state
      for (int i=0; i<v_death_Iv; i++) {
        int rnd1 = sample2(0,n_Iv[k]-1-i);
        Iv[k].erase(Iv[k].begin()+rnd1);
      }
      n_Iv[k] -= v_death_Iv;
      
      int n_infectious = h_infectious[k].size();
      double rate_v_infection = a*c*n_infectious/double(H[k]);
      double prob_v_infection_or_death = 1 - exp(-(rate_v_infection + mu)); // mosquito infection or death in Sv state (competing hazards)
      double prob_v_infection = rate_v_infection/(rate_v_infection + mu); // relative rate of mosquito infection
      int v_infection_or_death = rbinom1(n_Sv[k], prob_v_infection_or_death); // mosquito infection or death in Sv state (competing hazards)
      int v_infection = rbinom1(v_infection_or_death, prob_v_infection); // mosquito infection
      n_Sv[k] += -v_infection + v_death_Iv; // update Sv
      n_Ev[k][v_buffer_index] = v_infection; // add infecteds to Ev list
      
      // fill in Ev[k][v_buffer_index] with source hosts
      if (v_infection>0) {
        vector<int> v_infection_draws(v_infection);
        for (int i=0; i<v_infection; i++) {
          v_infection_draws[i] = sample2(0,n_infectious-1);
        }
        sort(v_infection_draws.begin(), v_infection_draws.end());
        int i = 0;
        int j = 0;
        for (auto it=h_infectious[k].begin(); it!=h_infectious[k].end(); ++it) {
          while (v_infection_draws[j]==i) {
            Ev[k][v_buffer_index].push_back(pair<int, int>{*it, t});
            j++;
          }
          if (j==v_infection) {
            break;
          }
          i++;
        }
      }
      
      // deaths in Ev state
      int j2 = v_buffer_index;
      for (int j=0; j<v; j++) { // loop through all previous entries in Ev list
        j2 = (j2==0) ? v : j2-1;
        if (n_Ev[k][j2]>0) {
          int v_death_Ev = rbinom1(n_Ev[k][j2], prob_v_death); // mosquito death in this Ev state
          for (int i=0; i<v_death_Ev; i++) {
            int rnd1 = sample2(0,n_Ev[k][j2]-1-i);
            Ev[k][j2].erase(Ev[k][j2].begin()+rnd1);
          }
          n_Ev[k][j2] -= v_death_Ev; // subtract mosquito deaths from this element
          n_Sv[k] += v_death_Ev; // respawn mosquitoes in Sv state
        }
      }
      
      // at this stage j2 is behind v_buffer_index by v steps. Move final Ev
      // list entries from Ev to Iv
      for (int i=0; i<n_Ev[k][j2]; i++) {
        Iv[k].push_back(Ev[k][j2][i]);
      }
      n_Iv[k] += n_Ev[k][j2];
      n_Ev[k][j2] = 0;
      Ev[k][j2].clear();
      
      // step time forward and store counts
      for (int i=0; i<int(h[k].size()); i++) {
        // step time forwards
        h[k][i].step_forward(t, u, g, h_infectious[k]);
        
        // record clonality
        if (h[k][i].n_infections<max_clonality) {
          clonality[k][t][h[k][i].n_infections]++;
        } else {
          clonality[k][t][max_clonality]++;
        }
      }
      
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("clonality") = clonality,
                            Rcpp::Named("line_list") = line_list);
}
