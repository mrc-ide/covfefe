
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
  int max_infections = Rcpp::as<int>(args["max_infections"]);
  vector<vector<double>> delta_mig = Rcpp::as<vector<vector<double>>>(args["delta_mig"]);
  int demes = Rcpp::as<int>(args["demes"]);
  
  double prob_h_recovery = 1 - exp(-r);
  double prob_v_death = 1 - exp(-mu);
  
  // store results
  vector<vector<int>> line_list;
  vector<int> line_list_migration;
  vector<int> line_list_infection;
  vector<int> line_list_bloodstage;
  vector<int> line_list_recover;
  vector<vector<vector<int>>> store_n_bloodstage(demes, vector<vector<int>>(max_time, vector<int>(max_infections+1)));
  
  // create schedules
  vector<vector<pair<int, int>>> schedule_bloodstage(max_time);
  vector<vector<pair<int, int>>> schedule_infective(max_time);
  vector<vector<tuple<int, int, int>>> schedule_recover(max_time);  // first=host ID, second=infection ID, third=is infective?(0=false, 1=true)
  
  // create population of human hosts
  vector<indiv_host> host(sum(H));
  vector<vector<int>> host_noninfective(demes);
  vector<vector<int>> host_infective(demes);
  int j = 0;
  for (int k=0; k<demes; k++) {
    host_noninfective[k] = vector<int>(H[k]);
    for (int i=0; i<H[k]; i++) {
      host[j].deme = k;
      host_noninfective[k][i] = j;
      j++;
    }
  }
  
  // seed initial infections
  for (int k=0; k<demes; k++) {
    for (int i=0; i<Ih_init[k]; i++) {
      
      // choose host
      int host_ID = host_noninfective[k][i];
      
      // add to infection list
      vector<int> tmp = {host_ID, k, host[host_ID].total_infections, -1, -1};
      push_back_multiple(line_list_infection, tmp);
      
      // schedule blood stage
      schedule_bloodstage[u].push_back({host_ID, 0});
      
      // infect
      host[host_ID].total_infections++;
      host[host_ID].n_latent++;
    }
  }
  
  // add first entries to line list
  line_list.push_back(line_list_migration);
  line_list.push_back(line_list_infection);
  line_list.push_back(line_list_bloodstage);
  line_list.push_back(line_list_recover);
  
  // store initial results
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
        store_n_bloodstage[k][0][nb]++;
      }
    }
  }
  
  // create objects representing mosquitoes
  vector<int> n_Sv = M;
  vector<vector<int>> n_Ev_death(demes, vector<int>(u));
  vector<vector<vector<pair<int, int>>>> Ev(demes, vector<vector<pair<int, int>>>(u+1));
  vector<vector<pair<int, int>>> Iv(demes);
  int v_ringbuffer_thistime = 0;
  
  // initialise objects for implementing migration
  vector<vector<vector<int>>> mig_noninf_hosts(demes, vector<vector<int>>(demes));
  vector<vector<vector<int>>> mig_inf_hosts(demes, vector<vector<int>>(demes));
  
  // carry out simulation loop
  for (int t=1; t<max_time; t++) {
    
    // clear lists
    line_list_migration.clear();
    line_list_infection.clear();
    line_list_bloodstage.clear();
    line_list_recover.clear();
    
    // update ring buffer indices
    v_ringbuffer_thistime = (v_ringbuffer_thistime==v) ? 0 : v_ringbuffer_thistime+1;
    
    
    //#### MIGRATION
    // schedule hosts to move
    for (int k1=0; k1<demes; k1++) {
      for (int k2=0; k2<demes; k2++) {
        if (k1==k2 || delta_mig[k1][k2]==0) {
          continue;
        }
        
        // split migrations between non-infectives and infectives
        for (int i=0; i<delta_mig[k1][k2]; i++) {
          int n_noninf = host_noninfective[k1].size();
          int n_inf = host_infective[k1].size();
          double prob_h_migration_noninf = n_noninf/double(n_noninf+n_inf);  // proportion of migrations in non-infetive hosts
          
          int host_ID;
          if (rbernoulli1(prob_h_migration_noninf)) { // schedule non-infectives to move
            int rnd1 = sample2(0, host_noninfective[k1].size()-1);
            host_ID = host_noninfective[k1][rnd1];
            mig_noninf_hosts[k1][k2].push_back(host_ID);
            host_noninfective[k1].erase(host_noninfective[k1].begin()+rnd1);
          } else {  // schedule infectives to move
            int rnd1 = sample2(0, host_infective[k1].size()-1);
            host_ID = host_infective[k1][rnd1];
            mig_inf_hosts[k1][k2].push_back(host_ID);
            host_infective[k1].erase(host_infective[k1].begin()+rnd1);
          }
          
          // add migration event to line list
          int nb = host[host_ID].n_latent + host[host_ID].n_bloodstage + host[host_ID].n_infective;
          if (nb>0) {
            vector<int> tmp = {host_ID, host[host_ID].deme, k2};
            push_back_multiple(line_list_migration, tmp);
          }
          
          host[host_ID].deme = k2;
        }
      }
    }
    // move hosts
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
          Ev[k][v_ringbuffer_thistime].push_back({host_infective[k][rnd1], t});
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
        vector<int> tmp = {host_ID, k, host[host_ID].total_infections, Iv[k][rnd2].first, Iv[k][rnd2].second};
        push_back_multiple(line_list_infection, tmp);
        
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
      int host_deme = host[host_ID].deme;
      int inf_ID = schedule_bloodstage[t][i].second;
      
      // add to bloodstage list
      vector<int> tmp = {host_ID, host_deme, inf_ID};
      push_back_multiple(line_list_bloodstage, tmp);
      
      // new blood stage infection
      host[host_ID].n_latent--;
      host[host_ID].n_bloodstage++;
      
      // schedule infective and/or recovery
      int dur = rgeom1(prob_h_recovery);
      if (dur>g) {  // if become infective prior to recovery
        if ((t+dur)<max_time) {
          schedule_infective[t+g].push_back({host_ID, inf_ID});
          schedule_recover[t+dur].push_back({host_ID, inf_ID, 1});
        } else if ((t+g)<max_time) {
          schedule_infective[t+g].push_back({host_ID, inf_ID});
        }
      } else {  // if recover prior to become infective
        if ((t+dur)<max_time) {
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
      int host_deme = host[host_ID].deme;
      int inf_ID = get<1>(schedule_recover[t][i]);
      int inf_type = get<2>(schedule_recover[t][i]);
      
      // add to recover list
      vector<int> tmp = {host_ID, host_deme, inf_ID};
      push_back_multiple(line_list_recover, tmp);
      
      // type 0 = blood stage recovery, type 1 = infective recovery
      if (inf_type==1) {
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
          store_n_bloodstage[k][t][nb]++;
        }
      }
    }
    
    // to line list
    line_list.push_back(line_list_migration);
    line_list.push_back(line_list_infection);
    line_list.push_back(line_list_bloodstage);
    line_list.push_back(line_list_recover);
    
  } // end loop over time
  
  return Rcpp::List::create(Rcpp::Named("n_bloodstage") = store_n_bloodstage,
                            Rcpp::Named("line_list") = line_list);
}

