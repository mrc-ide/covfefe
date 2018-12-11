
#include "main.h"
#include "base_model.Parameters.h"
#include "base_model.Dispatcher.h"
#include "probability.h"
#include "genotype.h"

#include <chrono>

using namespace std;

//------------------------------------------------
// main function (when not run using Rcpp)
#ifndef RCPP_ACTIVE
int main(int argc, const char * argv[]) {
  
  // run simulation
  indiv_sim_cpp();
  
}
#endif

//------------------------------------------------
// draw from simple individual-based model
#ifdef RCPP_ACTIVE
// [[Rcpp::export]]
Rcpp::List indiv_sim_cpp(Rcpp::List args) {
  
  // start timer
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  
  // extract model parameters into separate class
  Parameters parameters(args);
  //parameters.print_summary();
  
  // create simulation dispatcher object
  Dispatcher dispatcher;
  
  // carry out simulation
  dispatcher.simulate();
  
  // end timer
  chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
  chrono::duration<double> time_span = chrono::duration_cast< chrono::duration<double> >(t2-t1);
  print("completed in", time_span.count(), "seconds\n");
  
  return Rcpp::List::create(Rcpp::Named("H_store") = dispatcher.H_store,
                            Rcpp::Named("Sh_store") = dispatcher.Sh_store,
                            Rcpp::Named("Lh_store") = dispatcher.Lh_store,
                            Rcpp::Named("Ah_store") = dispatcher.Ah_store,
                            Rcpp::Named("Ch_store") = dispatcher.Ch_store,
                            Rcpp::Named("Sv_store") = dispatcher.Sv_store,
                            Rcpp::Named("Lv_store") = dispatcher.Lv_store,
                            Rcpp::Named("Iv_store") = dispatcher.Iv_store,
                            Rcpp::Named("EIR_store") = dispatcher.EIR_store,
                            Rcpp::Named("H_age_store") = dispatcher.H_age_store,
                            Rcpp::Named("Sh_age_store") = dispatcher.Sh_age_store,
                            Rcpp::Named("Lh_age_store") = dispatcher.Lh_age_store,
                            Rcpp::Named("Ah_age_store") = dispatcher.Ah_age_store,
                            Rcpp::Named("Ch_age_store") = dispatcher.Ch_age_store,
                            Rcpp::Named("inc_Lh_age_store") = dispatcher.inc_Lh_age_store,
                            Rcpp::Named("inc_Ah_age_store") = dispatcher.inc_Ah_age_store);
}
#else
int indiv_sim_cpp() {
  
  // start timer
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
  
  // extract model parameters into separate class
  Parameters parameters;
  
  // create simulation dispatcher object
  Dispatcher dispatcher;
  
  // carry out simulation
  dispatcher.simulate();
  
  // end timer
  chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
  chrono::duration<double> time_span = chrono::duration_cast< chrono::duration<double> >(t2-t1);
  print("completed in", time_span.count(), "seconds\n");
  
  return 0;
}
#endif

/*
//------------------------------------------------
// simulate population and draw blood stage hosts at designated times
// [[Rcpp::export]]
Rcpp::List draw_hosts_cpp(Rcpp::List &infection_history, Rcpp::List &args) {
  
  // extract arguments
  vector<int> samp_times = rcpp_to_vector_int(args["samp_times"]);
  vector<vector<int>> samp_demes = rcpp_to_matrix_int(args["samp_demes"]);
  vector<vector<int>> samp_num = rcpp_to_matrix_int(args["samp_num"]);
  int demes = rcpp_to_int(args["demes"]);
  bool recover_pop_counts = rcpp_to_bool(args["recover_pop_counts"]);
  int max_infections = rcpp_to_int(args["max_infections"]);
  int n_samp = samp_times.size();
  int max_samp_time = samp_times[n_samp-1];
  int max_time = infection_history.size();
  
  // create output object. For each deme and each designated time point, create 
  // a map of hosts. Map key is host ID, map value is the number of observable
  // (i.e. bloodstage or infective stage) infections at that time.
  vector<vector<map<int, int>>> samp_hosts(demes, vector<map<int, int>>(n_samp));
  
  // for each deme, create running map of hosts. Map key is host ID, map value
  // is a second map of infections. Second map key is the infection ID, second
  // map value is whether the infection is observable (i.e. bloodstage or
  // infective stage).
  vector<map<int, map<int, bool>>> hosts(demes);
  
  // optionally store counts of number of blood stage infections at each time
  // point
  vector<vector<vector<int>>> pop_counts(demes, vector<vector<int>>(max_time, vector<int>(max_infections)));
  
  // run simulation until max_time, or until max_samp_time if not recovering
  // population counts
  int max_t = recover_pop_counts ? max_time : max_samp_time;
  
  
  // SIMULATION ############################################################
  
  int samp_i = 0; // step through sample times
  Rcpp::List infection_history_t; // stores each time point of infection history
  vector<int> this_line;  // stores each sub-list of infection history
  
  for (int t=0; t<max_t; t++) {
    infection_history_t = infection_history[t];
    
    //#### MIGRATION
    this_line = rcpp_to_vector_int(infection_history_t[0]);
    for (int j=0; j<int(this_line.size())/3; j++) {
      int host_ID = this_line[3*j];
      int host_deme = this_line[3*j+1];
      int new_deme = this_line[3*j+2];
      
      // add host to new deme and erase from old
      hosts[new_deme][host_ID] = hosts[host_deme][host_ID];
      hosts[host_deme].erase(host_ID);
    }
    
    //#### INFECTION
    this_line = rcpp_to_vector_int(infection_history_t[1]);
    for (int j=0; j<int(this_line.size())/6; j++) {
      int host_ID = this_line[6*j];
      int host_deme = this_line[6*j+1];
      int inf_ID = this_line[6*j+2];
      
      // add new non-observable infection to host
      hosts[host_deme][host_ID][inf_ID] = false;
    }
    
    //#### BLOODSTAGE
    this_line = rcpp_to_vector_int(infection_history_t[2]);
    for (int j=0; j<int(this_line.size())/3; j++) {
      int host_ID = this_line[3*j];
      int host_deme = this_line[3*j+1];
      int inf_ID = this_line[3*j+2];
      
      // infection becomes observable
      hosts[host_deme][host_ID][inf_ID] = true;
    }
    
    //#### RECOVERY
    this_line = rcpp_to_vector_int(infection_history_t[3]);
    for (int j=0; j<int(this_line.size())/3; j++) {
      int host_ID = this_line[3*j];
      int host_deme = this_line[3*j+1];
      int inf_ID = this_line[3*j+2];
      
      // if host has a single infection then this represents complete recovery
      // (i.e. drop host). Otherwise drop just this one infection
      if (hosts[host_deme][host_ID].size()==1) {
        hosts[host_deme].erase(host_ID);
      } else {
        hosts[host_deme][host_ID].erase(inf_ID);
      }
    }
    
    //#### ADD TO OUTPUT
    // recover counts of number of blood stage infections
    if (recover_pop_counts) {
      for (int k=0; k<demes; k++) {
        for (auto it = hosts[k].begin(); it!=hosts[k].end(); ++it) {
          int host_ID = it->first;
          
          // count observable infections
          int n_observable = 0;
          for (auto it2 = hosts[k][host_ID].begin(); it2!=hosts[k][host_ID].end(); ++it2) {
            int inf_ID = it2->first;
            if (hosts[k][host_ID][inf_ID] == true) {
              n_observable++;
            }
          }
          
          // add to population counts
          if (n_observable>0 && n_observable<=max_infections) {
            pop_counts[k][t][n_observable-1]++;
          }
        }
      }
    }
    // if designated sampling time then sample hosts
    if (t==(samp_times[samp_i]-1)) {
      for (int ki=0; ki<int(samp_demes[samp_i].size()); ki++) { // loop through demes
        int k = samp_demes[samp_i][ki]-1;
        int n = samp_num[samp_i][ki];
        
        // loop through all hosts in this deme
        vector<int> v;
        for (auto it = hosts[k].begin(); it!=hosts[k].end(); ++it) {
          int host_ID = it->first;
          
          // loop through infections to see if any are observable
          bool is_observable = false;
          for (auto it2 = hosts[k][host_ID].begin(); it2!=hosts[k][host_ID].end(); ++it2) {
            int inf_ID = it2->first;
            if (hosts[k][host_ID][inf_ID] == true) {
              is_observable = true;
              break;
            }
          }
          
          // if host is observable
          if (is_observable) {
            // if n<0 then sampling all individuals, so add directly to sample. 
            // Otherwise add to vector v which will be randomly subsampled
            if (n<0) {
              samp_hosts[k][samp_i][host_ID] = hosts[k][host_ID].size();
            } else {
              v.push_back(host_ID);
            }
          }
        }
        
        // randomly subsample v into final output
        if (n>0) {
          n = (n < v.size()) ? n : v.size();  // limit n at total number observable infections
          for (int i=0; i<n; i++) {
            int rnd1 = sample2(0, v.size()-1);
            int host_ID = v[rnd1];
            samp_hosts[k][samp_i][host_ID] = hosts[k][host_ID].size();
            v.erase(v.begin()+rnd1);
          }
        }
      } // close ki loop
      samp_i++;
      
    } // close if (t == samp time) condition
    
  } // end simulation
  
  // return
  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("samp_hosts") = samp_hosts,
                                      Rcpp::Named("pop_counts") = pop_counts);
  return ret;
}

//------------------------------------------------
// prune infection tree
// [[Rcpp::export]]
Rcpp::List prune_cpp(Rcpp::List &infection_history, Rcpp::List &samp_hosts_raw, Rcpp::List &args) {
  
  // extract arguments
  vector<int> samp_times = rcpp_to_vector_int(args["samp_times"]);
  int demes = samp_hosts_raw.size();
  int n_samp = samp_times.size();
  int max_samp_time = samp_times[n_samp-1];
  
  // recover samp_hosts object from samp_hosts_raw. samp_hosts is a map of hosts
  // for each deme and each designated time point. The map key is the host ID,
  // the map value is the number of observable (i.e. bloodstage or infective
  // stage) infections at that time.
  vector<vector<map<int, int>>> samp_hosts(demes, vector<map<int, int>>(n_samp));
  Rcpp::List samp_hosts_raw_k;
  vector<int> v;
  for (int k=0; k<demes; k++) {
    samp_hosts_raw_k = samp_hosts_raw[k];
    for (int i=0; i<n_samp; i++) {
      v = rcpp_to_vector_int(samp_hosts_raw_k[i]);
      int n = v.size()/2;
      for (int j=0; j<n; j++) {
        samp_hosts[k][i][v[j]] = v[n+j];
      }
    }
  }
  
  // object for storing results. At each time step, records hosts that require 
  // de-novo creation, hosts that require infection and from whom (two values), 
  // and hosts that can be safely deleted as they are never needed again.
  vector<vector<vector<int>>> pruned(max_samp_time, vector<vector<int>>(3));
  
  // create a schedule of hosts that become relevant at each time point. 
  // Relevant hosts are those that are in the history of the sampled hosts, i.e.
  // hosts for which we will ultimately need to simulate genotypes. Start by
  // adding samp_hosts to the schedule.
  vector<map<int, int>> schedule_hosts(max_samp_time);
  for (int k=0; k<demes; k++) {
    for (int i=0; i<n_samp; i++) {
      int this_time = samp_times[i]-1;
      for (auto it = samp_hosts[k][i].begin(); it!=samp_hosts[k][i].end(); ++it) {
        int host_ID = it->first;
        schedule_hosts[this_time][host_ID] = samp_hosts[k][i][host_ID];
      }
    }
  }
  
  // running list of relevant hosts
  map<int, int> hosts;
  
  // store the latest occurring entry for each relevant host
  map<int, int> latest_entry;
  
  
  // SIMULATION ############################################################
  
  // loop backwards through infection history
  Rcpp::List infection_history_t; // stores each time point of infection history
  vector<int> this_line;
  for (int t=max_samp_time-1; t>=0; t--) {
    infection_history_t = infection_history[t];
    
    //#### ADD SCHEDULED HOSTS TO RUNNING LIST
    for (auto it = schedule_hosts[t].begin(); it!=schedule_hosts[t].end(); ++it) {
      int host_ID = it->first;
      hosts[host_ID] = schedule_hosts[t][host_ID];
    }
    
    //#### INFECTION
    this_line = rcpp_to_vector_int(infection_history_t[1]);
    for (int j=0; j<int(this_line.size())/6; j++) {
      int host_ID = this_line[6*j];
      
      // if this host is in list of relevant hosts
      if (hosts.count(host_ID)>0) {
        int source_ID = this_line[6*j+3];
        
        // if this is a de-novo infection then add to pruned_de_novo list and skip over
        if (source_ID==-1) {
          pruned[t][0].push_back(host_ID);
          hosts.erase(host_ID);
          continue;
        }
        
        // add source to schedule
        int source_n = this_line[6*j+4];
        int source_time = this_line[6*j+5];
        schedule_hosts[source_time][source_ID] = source_n;
        
        // add to list of infections
        pruned[source_time][1].push_back(source_ID);
        pruned[source_time][1].push_back(host_ID);
        
        // update latest entry for this source as needed
        if (latest_entry.count(source_ID)==0) {
          latest_entry[source_ID] = source_time;
        } else if (latest_entry[source_ID]<source_time) {
          latest_entry[source_ID] = source_time;
        }
        
        // reduce infections of host_ID by 1, and drop if none remaining
        hosts[host_ID]--;
        if (hosts[host_ID]==0) {
          hosts.erase(host_ID);
        }
      }
    }
    
  } // end loop through infection history
  
  // remove samp_hosts from latest entries
  for (int k=0; k<demes; k++) {
    for (int i=0; i<n_samp; i++) {
      for (auto it = samp_hosts[k][i].begin(); it!=samp_hosts[k][i].end(); ++it) {
        int host_ID = it->first;
        latest_entry.erase(host_ID);
      }
    }
  }
  
  // fill in deletion list
  for (auto it = latest_entry.begin(); it!=latest_entry.end(); ++it) {
    int host_ID = it->first;
    int delete_time = it->second;
    pruned[delete_time][2].push_back(host_ID);
  }
  
  // return
  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("pruned") = pruned);
  return ret;
}

//------------------------------------------------
// simulate genotypes
// [[Rcpp::export]]
Rcpp::List sim_genotypes_cpp(Rcpp::List &samp_hosts_raw, Rcpp::List &args) {
  
  // extract arguments
  vector<int> samp_times = rcpp_to_vector_int(args["samp_times"]);
  vector<vector<vector<int>>> pruned = rcpp_to_array_int(args["pruned"]);
  vector<double> dist_oocysts = rcpp_to_vector_double(args["dist_oocysts"]);
  vector<double> dist_hepatocytes = rcpp_to_vector_double(args["dist_hepatocytes"]);
  vector<vector<int>> loci = rcpp_to_matrix_int(args["loci"]);
  vector<bool> loci_null = rcpp_to_vector_bool(args["loci_null"]);
  double recom_rate = rcpp_to_double(args["recom_rate"]);
  int n_samp = samp_times.size();
  int max_samp_time = pruned.size();
  int demes = samp_hosts_raw.size();
  
  // get sample host IDs in each deme and at each sample time
  vector<vector<vector<int>>> samp_hosts(demes, vector<vector<int>>(n_samp));
  Rcpp::List samp_hosts_raw_k;
  vector<int> v;
  for (int k=0; k<demes; k++) {
    samp_hosts_raw_k = samp_hosts_raw[k];
    for (int i=0; i<n_samp; i++) {
      v = rcpp_to_vector_int(samp_hosts_raw_k[i]);
      int n = v.size()/2;
      samp_hosts[k][i] = vector<int>(n);
      for (int j=0; j<n; j++) {
        samp_hosts[k][i][j] = v[j];
      }
    }
  }
  
  // drop null chromosomes
  for (int i=0; i<int(loci_null.size()); i++) {
    if (loci_null[i] == true) {
      loci[i].clear();
    }
  }
  
  // map for storing genotypes
  map<int, genotype> gen_map;
  
  
  // SIMULATION  ############################################################
  
  // loop through time
  int next_denovo = 0;
  for (int t=0; t<max_samp_time; t++) {
    
    //#### DE-NOVO GENOTYPES
    for (int i=0; i<int(pruned[t][0].size()); i++) {
      int host_ID = pruned[t][0][i];
      gen_map[host_ID].de_novo(loci,next_denovo++);
    }
    
    //#### INFECTION
    for (int i=0; i<int(pruned[t][1].size())/2; i++) {
      int host_ID = pruned[t][1][2*i];
      int targ_ID = pruned[t][1][2*i+1];
      genotype_infect(gen_map[host_ID], gen_map[targ_ID], dist_oocysts, dist_hepatocytes, loci, recom_rate);
    }
    
    //#### DELETE REDUNDANT
    for (int i=0; i<int(pruned[t][2].size()); i++) {
      int host_ID = pruned[t][2][i];
      gen_map.erase(host_ID);
    }
    
  } // end loop through time
  
  
  // save genotypes in list
  vector<vector<vector<vector<vector<vector<int>>>>>> genotypes(demes);
  for (int k=0; k<demes; k++) {
    genotypes[k] = vector<vector<vector<vector<vector<int>>>>>(n_samp);
    for (int i=0; i<n_samp; i++) {
      int n_hosts = samp_hosts[k][i].size();
      genotypes[k][i] = vector<vector<vector<vector<int>>>>(n_hosts);
      for (int j=0; j<n_hosts; j++) {
        int host_ID = samp_hosts[k][i][j];
        genotypes[k][i][j] = gen_map[host_ID].x;
      }
    }
  }
  
  // return
  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("genotypes") = genotypes);
  return ret;
}
*/

