
#include "base_model.Dispatcher.h"
#include "probability.h"
#include "misc.h"

using namespace std;

//------------------------------------------------
// default constructor for Dispatcher class
Dispatcher::Dispatcher() {
  
  // initialise single static population of human hosts over all demes. This is
  // preferable to using separate vectors of hosts for each deme, as this would
  // mean moving hosts around due to migration. With a single static population
  // we can simply change the "deme" attribute of a host to represent migration.
  population = Population(sum(H_vec));
  
  // initialise a single static scheduler object. This will be used to schedule
  // future events for all individuals in all demes. The reason for using a
  // single scheduler over all demes is that an event could be scheduled in one
  // deme, but the host could have migrated before the scheduled time, meaning
  // it should be enacted in a different deme.
  scheduler = Scheduler();
  
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
    demes[k].init(hosts_k, Ih_vec[k]);
  }
  
  
  // initialise objects for implementing migration. Store IDs of infective and
  // non-infective hosts that will move demes in a given time step
  //mig_inf_hosts = array_int(demes, demes);
  //mig_noninf_hosts = array_int(demes, demes);
  
}

//------------------------------------------------
// run simulation
void Dispatcher::simulate() {
  
  // loop through daily time steps
  for (int t=1; t<max_time; t++) {
    
    //print_stars();
    //print(t);
    
    // step forward in all demes
    for (int k=0; k<n_demes; ++k) {
      demes[k].step_forward(t);
    }
    
    // implement scheduled transitions to blood-stage
    for (int i=0; i<int(scheduler.schedule_bloodstage[t].size()); i++) {
      int host_index = scheduler.schedule_bloodstage[t][i].first;
      //int infection_ID = scheduler.schedule_bloodstage[t][i].second;
      //int host_deme = population.hosts[host_index].deme;
      
      // transition to blood stage
      population.hosts[host_index].transition_bloodstage();
      
      /*
      // add to bloodstage list
      if (output_infection_history) {
        vector<int> tmp = {host_ID, host_deme, inf_ID};
        push_back_multiple(history_bloodstage, tmp);
      }
      */
    }
    
    // blood stage become infective
    for (int i=0; i<int(scheduler.schedule_infective[t].size()); i++) {
      int host_index = scheduler.schedule_infective[t][i].first;
      
      // move from noninfective to infective list within the host's deme
      if (population.hosts[host_index].n_infective == 0) {
        int host_deme = population.hosts[host_index].deme;
        demes[host_deme].hosts_uninfective.erase(remove(demes[host_deme].hosts_uninfective.begin(), demes[host_deme].hosts_uninfective.end(), host_index));
        demes[host_deme].hosts_infective.push_back(host_index);
      }
      
      // transition to infective stage
      population.hosts[host_index].transition_infective();
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