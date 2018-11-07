
#include "base_model.Dispatcher.h"
#include "probability.h"

using namespace std;

//------------------------------------------------
// default constructor for Dispatcher class
Dispatcher::Dispatcher() {
  
  // initialise single static population of human hosts over all demes. This is
  // preferable to using separate vectors of hosts for each deme, as this would
  // mean moving hosts around due to migration. With a single static population
  // we can simply change the "deme" attribute of a host to represent migration.
  H = sum(H_vec);
  population = Population(H);
  
  // create vector of demes
  demes = vector<Deme>(n_demes);
  for (int k=0; k<n_demes; ++k) {
    demes[k] = Deme(k);
  }
  
  // initialise all demes (assign hosts and seed infections)
  int tmp = 0;
  for (int k=0; k<n_demes; ++k) {
    vector<int> hosts_k = seq_int(tmp, tmp+H_vec[k]-1);
    tmp += H_vec[k];
    demes[k].init(hosts_k, Eh_vec[k]);
  }
  
  // initialise objects for implementing migration. Store IDs of infective and
  // non-infective hosts that will move demes in a given time step
  //mig_inf_hosts = array_int(demes, demes);
  //mig_noninf_hosts = array_int(demes, demes);
  
  // initialise objects for storing daily counts of key quantities
  if (output_counts) {
    Sh_store = vector<vector<int>>(n_demes, vector<int>(max_time));
    Eh_store = vector<vector<int>>(n_demes, vector<int>(max_time));
    Ih_store = vector<vector<int>>(n_demes, vector<int>(max_time));
    //EIR_store = vector<vector<double>>(n_demes, vector<double>(max_time));
    
    for (int k=0; k<n_demes; ++k) {
      Sh_store[k][0] = H_vec[k] - Eh_vec[k];
      Eh_store[k][0] = Eh_vec[k];
    }
  }
  
}

//------------------------------------------------
// run simulation
void Dispatcher::simulate() {
  
  // loop through daily time steps
  for (int t=1; t<max_time; t++) {
    
    // step forward in all demes
    for (int k=0; k<n_demes; ++k) {
      demes[k].step_forward(t);
    }
    
    // store daily counts of key quantities
    if (output_counts) {
      for (int i=0; i<H; ++i) {
        int this_deme = population.hosts[i].deme;
        
        // store susceptible count
        if (population.hosts[i].n_asexual == 0) {
          Sh_store[this_deme][t]++;
        }
        
        // store latent count
        if (population.hosts[i].n_latent > 0) {
          Eh_store[this_deme][t]++;
        }
        
        // store bloodstage count
        if (population.hosts[i].n_bloodstage > 0) {
          Ih_store[this_deme][t]++;
        }
      }
      //print(t, Sh_store[0][t]);
    }
    
    if (t == 1000) {
      //population.hosts[0].summary();
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
