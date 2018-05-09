
#include "host.h"
#include "misc.h"

using namespace std;

//------------------------------------------------
// default constructor for host class
host::host() {}

//------------------------------------------------
// constructor for host class
host::host(int duration_) {
  
  // parameters
  duration = duration_;
  //print(duration);
  
}

//------------------------------------------------
/*
// generate scaffold groupings
void MCMC::scaffold_mcmc(Rcpp::List &args) {
  
  // print header
  if (print_console) {
    print("Generating", scaf_n, "scaffolds");
  }
  
  // extract R functions
  Rcpp::Function test_convergence = args["test_convergence"];
  Rcpp::Function update_progress = args["update_progress"];
  
  // initialise objects needed for generating scaffolds
  int batch_size = 10;
  int max_batches = 10;
  vector<double> batch_vec(batch_size);
  int dummy_accept = 0;
  
  // generate multiple scaffolds
  for (int scaf_rep=0; scaf_rep<scaf_n; scaf_rep++) {
    
    // reset particles
    for (int r=0; r<rungs; r++) {
      particle_vec[r].reset();
      particle_vec[r].beta = beta_vec[r];
    }
    rung_order = seq_int(0,rungs-1);
    
    // store log-likelihoods
    vector<vector<double>> scaf_log_like(rungs, vector<double>(batch_size));
    
    // loop through iterations in batches
    for (int batch=0; batch<max_batches; batch++) {
      
      // iterations of this batch
      for (int rep=0; rep<batch_size; rep++) {
        
        // update particles
        for (int r=0; r<rungs; r++) {
          int rung = rung_order[r];
          
          // update parameters
          particle_vec[rung].update_mu();
          particle_vec[rung].update_group();
          
          // calculate log-likelihood
          particle_vec[rung].calc_log_like();
          
          // split-merge step
          if (splitmerge_on) {
            particle_vec[rung].splitmerge_propose(dummy_accept);
          }
        }
        
        // Metropolis-coupling
        if (coupling_on) {
          metropolis_coupling();
        }
        
        // store log-likelihoods
        for (int r=0; r<rungs; r++) {
          int rung = rung_order[r];
          scaf_log_like[rung][batch*batch_size + rep] = particle_vec[rung].loglike;
        }
        
      } // end iterations of this batch
      
      // break if converged
      bool all_converged = true;
      for (int r=0; r<rungs; r++) {
        int rung = rung_order[r];
        bool rung_converged = test_convergence(scaf_log_like[rung]);
        if (!rung_converged) {
          all_converged = false;
          break;
        }
      }
      if (all_converged) {
        break;
      }
      
      // if not converged expand scaf_log_like
      for (int r=0; r<rungs; r++) {
        push_back_multiple(scaf_log_like[r], batch_vec);
      }
      
      // throw warning if still not converged at end of batches
      if (batch == (max_batches-1)) {
        print("warning: scaffold did not converge");
      }
      
    } // loop over batches
    
    // loop over rungs
    for (int r=0; r<rungs; r++) {
      int rung = rung_order[r];
      
      // re-order group allocation to be always-increasing
      particle_vec[rung].group_increasing();
      
      // store this particle
      scaf_group[rung][scaf_rep] = particle_vec[rung].group;
      for (int i=0; i<n; i++) {
        scaf_counts[rung][scaf_rep][scaf_group[rung][scaf_rep][i]] ++;
        scaf_x_sum[rung][scaf_rep][scaf_group[rung][scaf_rep][i]] += (*x_ptr)[i];
      }
    }
    
    // update scaffold progress bar
    if (print_console) {
      update_progress(args, 1, scaf_rep+1, scaf_n);
    }
    
  } // loop over scaf_n
  
  // load these scaffolds back into particles
  for (int r=0; r<rungs; r++) {
    particle_vec[r].scaf_group = scaf_group[r];
    particle_vec[r].scaf_counts = scaf_counts[r];
    particle_vec[r].scaf_x_sum = scaf_x_sum[r];
  }
  
}
*/

