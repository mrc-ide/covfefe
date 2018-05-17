
#include "misc.h"
#include "probability.h"
#include "genotype.h"

using namespace std;

//------------------------------------------------
// simulate population and draw blood stage hosts at designated times
vector<vector<map<int, int>>> draw_hosts(Rcpp::List &line_list, Rcpp::List &args) {
  
  // extract arguments
  vector<int> samp_times = rcpp_to_vector_int(args["samp_times"]);
  vector<vector<int>> samp_demes = rcpp_to_matrix_int(args["samp_demes"]);
  vector<vector<int>> samp_num = rcpp_to_matrix_int(args["samp_num"]);
  int demes = rcpp_to_int(args["demes"]);
  int n_samp = samp_times.size();
  int max_samp_time = samp_times[n_samp-1];
  
  // create output object. For each deme and each designated time point, stores
  // a map of hosts. First element is host ID, second element is the number of
  // observable (i.e. bloodstage or infective stage) infections at the time of
  // sampling.
  vector<vector<map<int, int>>> store_hosts(demes, vector<map<int, int>>(n_samp));
  
  // create running population of hosts as a map. First element is the host ID, 
  // second element is a map of infections. First element is the infection ID, 
  // second element is whether the infection is observable
  vector<map<int, map<int, bool>>> hosts(demes);
  
  
  // SIMULATION ############################################################
  
  int i = 0;  // step through line-list
  int samp_i = 0; // step through sample times
  vector<int> this_line;  // stores each line of line-list
  while (i<max_samp_time*4) {
    int t = i/4;  // for entries per time step
    
    //#### MIGRATION
    this_line = rcpp_to_vector_int(line_list[i]);
    for (int j=0; j<int(this_line.size())/3; j++) {
      int host_ID = this_line[3*j];
      int host_deme = this_line[3*j+1];
      int new_deme = this_line[3*j+2];
      // add host to new deme and erase from old
      hosts[new_deme][host_ID] = hosts[host_deme][host_ID];
      hosts[host_deme].erase(host_ID);
    }
    i++;
    
    //#### INFECTION
    this_line = rcpp_to_vector_int(line_list[i]);
    for (int j=0; j<int(this_line.size())/6; j++) {
      int host_ID = this_line[6*j];
      int host_deme = this_line[6*j+1];
      int inf_ID = this_line[6*j+2];
      // add new non-observable infection to host
      hosts[host_deme][host_ID][inf_ID] = false;
    }
    i++;
    
    //#### BLOODSTAGE
    this_line = rcpp_to_vector_int(line_list[i]);
    for (int j=0; j<int(this_line.size())/3; j++) {
      int host_ID = this_line[3*j];
      int host_deme = this_line[3*j+1];
      int inf_ID = this_line[3*j+2];
      // infection becomes observable
      hosts[host_deme][host_ID][inf_ID] = true;
    }
    i++;
    
    //#### RECOVERY
    this_line = rcpp_to_vector_int(line_list[i]);
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
    i++;
    
    //#### ADD TO OUTPUT
    if (t==(samp_times[samp_i]-1)) {  // if designated sampling time
      for (int ki=0; ki<int(samp_demes[samp_i].size()); ki++) {
        int k = samp_demes[samp_i][ki]-1;
        int n = samp_num[samp_i][ki];
        
        // loop through all hosts
        vector<int> v;
        for (auto it = hosts[k].begin(); it!=hosts[k].end(); ++it) {
          int host_ID = it->first;
          // loop through infections to see if any are observable
          bool is_observable = false;
          for (auto it2 = hosts[k][host_ID].begin(); it2!=hosts[k][host_ID].end(); ++it2) {
            if (hosts[k][host_ID][it2->first]) {
              is_observable = true;
              break;
            }
          }
          // if observable
          if (is_observable) {
            // if n<0 then sampling all individuals, so add directly to sample. 
            // Otherwise add to vector v which will be randomly subsampled
            if (n<0) {
              store_hosts[k][samp_i][host_ID] = hosts[k][host_ID].size();
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
            store_hosts[k][samp_i][host_ID] = hosts[k][host_ID].size();
            v.erase(v.begin()+rnd1);
          }
        }
      } // close ki loop
      samp_i++;
      
    } // close if (t == samp time) condition
    
  } // end simulation
  
  // return
  return store_hosts;
}

//------------------------------------------------
// prune infection tree and simulate genotypes
// [[Rcpp::export]]
Rcpp::List sim_genotypes_cpp(Rcpp::List line_list, Rcpp::List args) {
  
  // extract arguments
  vector<int> samp_times = rcpp_to_vector_int(args["samp_times"]);
  int demes = rcpp_to_int(args["demes"]);
  int n_samp = samp_times.size();
  int max_samp_time = samp_times[n_samp-1];
  
  // draw hosts from given sampling times
  vector<vector<map<int, int>>> samp_hosts = draw_hosts(line_list, args);
  /*
  for (int i=0; i<n_samp; i++) {
    int this_time = samp_times[i]-1;
    Rcpp::Rcout << this_time << ":   ";
    for (auto it = samp_hosts[0][i].begin(); it!=samp_hosts[0][i].end(); ++it) {
      int host_ID = it->first;
      Rcpp::Rcout << host_ID << " ";
    }
    Rcpp::Rcout << "\n";
  }
  */
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
  
  // save events that make up pruned tree
  vector<int> de_novo;
  vector<vector<pair<int, int>>> list_infection(max_samp_time);
  map<int, int> latest_entry;
  
  
  // SIMULATION 1 ############################################################
  // prune infection tree
  
  // loop backwards through line list
  vector<int> this_line;
  int i = max_samp_time*4;
  while (i>0) {
    int t = i/4-1;
    
    //#### ADD SCHEDULED HOSTS
    for (auto it = schedule_hosts[t].begin(); it!=schedule_hosts[t].end(); ++it) {
      int host_ID = it->first;
      hosts[host_ID] = schedule_hosts[t][host_ID];
    }
    
    //#### INFECTION
    i -= 3;
    this_line = rcpp_to_vector_int(line_list[i]);
    for (int j=0; j<int(this_line.size())/6; j++) {
      int host_ID = this_line[6*j];
      
      // if this host is in list of relevant hosts
      if (hosts.count(host_ID)>0) {
        int source_ID = this_line[6*j+3];
        
        // if this is a de-novo infection then add to de_novo list and skip over
        if (source_ID==-1) {
          de_novo.push_back(host_ID);
          hosts.erase(host_ID);
          continue;
        }
        
        // add source to schedule
        int source_n = this_line[6*j+4];
        int source_time = this_line[6*j+5];
        schedule_hosts[source_time][source_ID] = source_n;
        
        // add to list of infections
        list_infection[source_time].push_back({source_ID, host_ID});
        
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
    i--;
    
    /*
    Rcpp::Rcout << t+1 << ":  ";
    for (auto it = hosts.begin(); it!=hosts.end(); ++it) {
      int host_ID = it->first;
      int inf_n = hosts[host_ID];
      Rcpp::Rcout << host_ID << " ";
    }
    Rcpp::Rcout << "\n";
    */
  }
  
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
  vector<vector<int>> list_delete(max_samp_time);
  for (auto it = latest_entry.begin(); it!=latest_entry.end(); ++it) {
    int host_ID = it->first;
    int delete_time = it->second;
    list_delete[delete_time].push_back(host_ID);
  }
  /*
  for (int i=0; i<max_samp_time; i++) {
    Rcpp::Rcout << i+1 << ":  ";
    for (auto it = list_infection[i].begin(); it!=list_infection[i].end(); ++it) {
      int host_ID = it->first;
      int targ_ID = it->second;
      Rcpp::Rcout << host_ID << "-" << targ_ID << " ";
    }
    Rcpp::Rcout << "  |  ";
    for (int j=0; j<int(list_delete[i].size()); j++) {
      Rcpp::Rcout << list_delete[i][j] << " ";
    }
    Rcpp::Rcout << "\n";
  }
  */
  
  
  // SIMULATION 2 ############################################################
  // simulate genetic data
  
  //vector<int> de_novo;
  //vector<vector<pair<int, int>>> list_infection(max_samp_time);
  //vector<vector<int>> list_delete(max_samp_time);
  
  // extract arguments
  vector<vector<int>> loci = rcpp_to_matrix_int(args["loci"]);
  int n_chrom = loci.size();
  
  map<int, genotype> gen_map;
  
  // initialise de-novo infections
  for (int i=0; i<int(de_novo.size()); i++) {
    int host_ID = de_novo[i];
    gen_map[host_ID] = genotype();
    
  }
  
  
  Rcpp::List ret = Rcpp::List::create(Rcpp::Named("foo") = samp_times);
  return ret;
}
