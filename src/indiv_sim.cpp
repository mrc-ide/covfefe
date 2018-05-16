
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
  vector<vector<double>> mig = Rcpp::as<vector<vector<double>>>(args["migration_matrix"]);
  int demes = Rcpp::as<int>(args["demes"]);
  
  double prob_h_recovery = 1 - exp(-r);
  double prob_v_death = 1 - exp(-mu);
  
  // store results
  vector<vector<int>> line_list;
  vector<vector<vector<int>>> store_n_bloodstage(demes, vector<vector<int>>(max_time+1, vector<int>(max_infections+1)));
  
  // create schedules
  vector<vector<int>> schedule_bloodstage(max_time+1);
  vector<vector<int>> schedule_infective(max_time+1);
  vector<vector<pair<int, int>>> schedule_recover(max_time+1);
  
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
      int rnd1 = sample2(0,H[k]-1);
      int host_ID = host_noninfective[k][rnd1];
      host[host_ID].total_infections++;
      schedule_bloodstage[u].push_back(host_ID);
    }
  }
  
  // create objects representing mosquitoes
  vector<int> n_Sv = M;
  vector<vector<int>> n_Ev_death(demes, vector<int>(u));
  vector<vector<vector<pair<int, int>>>> Ev(demes, vector<vector<pair<int, int>>>(u+1));
  vector<vector<pair<int, int>>> Iv(demes);
  int v_ringbuffer_thistime = 0;
  
  // initialise objects for implementing migration
  vector<vector<int>> delta_mig(demes, vector<int>(demes));
  vector<vector<vector<int>>> mig_noninf_hosts(demes, vector<vector<int>>(demes));
  vector<vector<vector<int>>> mig_inf_hosts(demes, vector<vector<int>>(demes));
  vector<int> mig_order = seq_int(0,demes-1);
  vector<int> rowsum(demes);
  vector<int> colsum(demes);
  
  // carry out simulation loop
  for (int t=0; t<=max_time; t++) {
    
    // update ring buffer indices
    v_ringbuffer_thistime = (v_ringbuffer_thistime==v) ? 0 : v_ringbuffer_thistime+1;
    
    
    //#### MIGRATION
    // draw number of migration events
    //draw_migration(delta_mig, mig, H, mig_order, rowsum, colsum);
    draw_migration2(delta_mig, mig, H, mig_order);
    
    
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
          
          if (rbernoulli1(prob_h_migration_noninf)) { // schedule non-infectives to move
            int rnd1 = sample2(0, host_noninfective[k1].size()-1);
            int this_ID = host_noninfective[k1][rnd1];
            host[this_ID].deme = k2;
            mig_noninf_hosts[k1][k2].push_back(this_ID);
            host_noninfective[k1].erase(host_noninfective[k1].begin()+rnd1);
          } else {  // schedule infectives to move
            int rnd1 = sample2(0, host_infective[k1].size()-1);
            int this_ID = host_infective[k1][rnd1];
            host[this_ID].deme = k2;
            mig_inf_hosts[k1][k2].push_back(this_ID);
            host_infective[k1].erase(host_infective[k1].begin()+rnd1);
          }
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
        int this_ID;
        if (rbernoulli1(prob_h_infection_noninf)) {
          int rnd1 = sample2(0, host_noninfective[k].size()-1);
          this_ID = host_noninfective[k][rnd1];
        } else {
          int rnd1 = sample2(0, host_infective[k].size()-1);
          this_ID = host_infective[k][rnd1];
        }
        
        // skip if reached max infections
        int nb = host[this_ID].n_bloodstage + host[this_ID].n_infective;
        if (nb==max_infections) {
          continue;
        }
        
        // increase total infections
        host[this_ID].total_infections++;
        
        // schedule move to blood stage
        if ((t+u)<=max_time) {
          schedule_bloodstage[t+u].push_back(this_ID);
        }
      }
      
    } // end loop through demes
    
    
    //#### IMPLEMENT SCHEDULED EVENTS
    // move to blood stage
    for (int i=0; i<int(schedule_bloodstage[t].size()); i++) {
      int this_ID = schedule_bloodstage[t][i];
      
      // new blood stage infection
      host[this_ID].n_bloodstage++;
      
      // schedule infective and/or recovery
      int dur = rgeom1(prob_h_recovery);
      if (dur>g) {  // if become infective prior to recovery
        if ((t+dur)<=max_time) {
          schedule_infective[t+g].push_back(this_ID);
          schedule_recover[t+dur].push_back({this_ID,2});
        } else if ((t+g)<=max_time) {
          schedule_infective[t+g].push_back(this_ID);
        }
      } else {  // if recover prior to become infective
        if ((t+dur)<=max_time) {
          schedule_recover[t+dur].push_back({this_ID,1});
        }
      }
    }
    // blood stage become infective
    for (int i=0; i<int(schedule_infective[t].size()); i++) {
      int this_ID = schedule_infective[t][i];
      
      // move to infective
      host[this_ID].n_bloodstage--;
      host[this_ID].n_infective++;
      
      // move from noninfective to infective list
      if (host[this_ID].n_infective==1) {
        int this_deme = host[this_ID].deme;
        host_noninfective[this_deme].erase(remove(host_noninfective[this_deme].begin(), host_noninfective[this_deme].end(), this_ID));
        host_infective[this_deme].push_back(this_ID);
      }
    }
    // recovery
    for (int i=0; i<int(schedule_recover[t].size()); i++) {
      int this_ID = schedule_recover[t][i].first;
      
      // type 1 = blood stage recovery, type 2 = infective recovery
      int this_type = schedule_recover[t][i].second;
      if (this_type==1) {
        host[this_ID].n_bloodstage--;
      } else {
        host[this_ID].n_infective--;
        
        // move from infective to noninfective list
        if (host[this_ID].n_infective==0) {
          int this_deme = host[this_ID].deme;
          host_infective[this_deme].erase(remove(host_infective[this_deme].begin(), host_infective[this_deme].end(), this_ID));
          host_noninfective[this_deme].push_back(this_ID);
        }
      }
    }
    
    
    //#### STORE RESULTS
    for (int k=0; k<demes; k++) {
      int s = host_infective[k].size();
      for (int i=0; i<H[k]; i++) {
        int this_ID;
        if (i<s) {
          this_ID = host_infective[k][i];
        } else {
          this_ID = host_noninfective[k][i-s];
        }
        int nb = host[this_ID].n_bloodstage + host[this_ID].n_infective;
        if (nb<=max_infections) {
          store_n_bloodstage[k][t][nb]++;
        }
      }
    }
    
  } // end loop over time
  
  return Rcpp::List::create(Rcpp::Named("n_bloodstage") = store_n_bloodstage,
                            Rcpp::Named("line_list") = -2);
}

void draw_migration2(vector<vector<int>> delta_mig, vector<vector<double>> &mig, vector<int> &H, vector<int> &mig_order) {
  
  // 
  print_stars("", 50);
  int demes = mig.size();
  vector<int> rowsum(demes);
  vector<int> colsum(demes);
  vector<double> p(demes,1.0);
  vector<double> q(demes);
  for (int i=0; i<3; i++) {
    
    // draw values for row i and column i of delta_mig from binomial
    // distribution
    for (int j=i; j<demes; j++) {
      int rnd1 = rbinom1(H[i]-rowsum[i], mig[i][j]/p[i]);
      delta_mig[i][j] = rnd1;
      rowsum[i] += rnd1;
      colsum[j] += rnd1;
      p[i] -= mig[i][j];
    }
    for (int j=i+1; j<demes; j++) {
      int rnd1 = rbinom1(H[j]-rowsum[j], mig[j][i]/p[j]);
      delta_mig[j][i] = rnd1;
      rowsum[j] += rnd1;
      colsum[i] += rnd1;
      p[j] -= mig[j][i];
    }
    
    print_matrix(delta_mig);
    print_vector(rowsum);
    print_vector(colsum);
    print("");
    if (i==2) {return;}
    
    // modify delta_mig
    if (colsum[i] < rowsum[i]) {  // if need more flow in
      int diff = rowsum[i] - colsum[i];
      fill(q.begin(), q.end(), 0);
      for (int j=i+1; j<demes; j++) {
        q[j] = delta_mig[i][j];
      }
      double qsum = sum(q);
      for (int z=0; z<diff; z++) {
        int rnd1 = sample1(q, qsum) - 1;
        delta_mig[i][rnd1]--;
        rowsum[i]--;
        colsum[rnd1]--;
        q[rnd1] -= 1;
        qsum -= 1;
      }
      /*
      int diff = rowsum[i] - colsum[i];
      fill(q.begin(), q.end(), 0);
      for (int j=i+1; j<demes; j++) {
        q[j] = mig[j][i]*(H[j]-rowsum[j]);
      }
      double qsum = sum(q);
      for (int z=0; z<diff; z++) {
        int rnd1 = sample1(q, qsum) - 1;
        delta_mig[rnd1][i]++;
        rowsum[rnd1]++;
        colsum[i]++;
        q[rnd1] -= mig[rnd1][i];
        qsum -= mig[rnd1][i];
      }
      */
    } else {  // if need less flow in
      int diff = colsum[i] - rowsum[i];
      fill(q.begin(), q.end(), 0);
      for (int j=i+1; j<demes; j++) {
        q[j] = delta_mig[j][i];
      }
      double qsum = sum(q);
      for (int z=0; z<diff; z++) {
        int rnd1 = sample1(q, qsum) - 1;
        delta_mig[rnd1][i]--;
        rowsum[rnd1]--;
        colsum[i]--;
        q[rnd1] -= 1;
        qsum -= 1;
      }
    }
    print_matrix(delta_mig);
    print_vector(rowsum);
    print_vector(colsum);
    print("");
    
  }
  
}

//------------------------------------------------
// draw migration events
void draw_migration(vector<vector<int>> delta_mig, vector<vector<double>> &mig, vector<int> &H, vector<int> &mig_order, vector<int> &rowsum, vector<int> &colsum) {
  
  // start by drawing number of events from multinomial distribution
  int demes = mig.size();
  for (int k1=0; k1<demes; k1++) {
    delta_mig[k1] = rmultinom1(H[k1], mig[k1]);
    delta_mig[k1][k1] = 0;
  }
  
  // calculate excess flow out of demes
  vector<int> excess_out(demes);
  row_sums(rowsum, delta_mig);
  col_sums(colsum, delta_mig);
  int total_loss = 0;
  for (int i=0; i<demes; i++) {
    excess_out[i] = rowsum[i] - colsum[i];
    total_loss += excess_out[i]*excess_out[i];
  }
  
  
  // re-shuffle order in which demes are modified
  //reshuffle(mig_order);
  
  // modify delta_mig to preserve constant human population size
  int cutout = 0;
  while (total_loss>0) {
    cutout++;
    
    // search all cells of delta_mig to find best option to modify
    int best_i = 0;
    int best_j = 0;
    int best_loss_improve = 0;
    int best_sub = 0;
    for (int i=0; i<demes; i++) {
      for (int j=0; j<demes; j++) {
        if (i==j) {
          continue;
        }
        if (excess_out[i]>0) {  // if excess flow out
          // get optimal value to subtract from delta_mig[i][j]
          int to_sub = round(0.5*(excess_out[i]-excess_out[j]));
          to_sub = (to_sub<=delta_mig[i][j]) ? to_sub : delta_mig[i][j];
          
          // calculate difference in loss function
          int new_loss = pow(excess_out[i]-to_sub, 2) + pow(excess_out[j]+to_sub, 2);
          int old_loss = pow(excess_out[i], 2) + pow(excess_out[j], 2);
          int loss_improve = old_loss - new_loss;
          
          // if best so far then keep
          if (loss_improve >= best_loss_improve) {
            best_loss_improve = loss_improve;
            best_sub = to_sub;
            best_i = i;
            best_j = j;
          }
        } else if (excess_out[i]<0) {
          // get optimal value to add to delta_mig[i][j]
          int to_add = round(0.5*(-excess_out[i]+excess_out[j]));
          to_add = (to_add<=(H[i]-rowsum[i])) ? to_add : (H[i]-rowsum[i]);
          
          // calculate difference in loss function
          int new_loss = pow(excess_out[i]+to_add, 2) + pow(excess_out[j]-to_add, 2);
          int old_loss = pow(excess_out[i], 2) + pow(excess_out[j], 2);
          int loss_improve = old_loss - new_loss;
          
          // if best so far then keep
          if (loss_improve >= best_loss_improve) {
            best_loss_improve = loss_improve;
            best_sub = -to_add;
            best_i = i;
            best_j = j;
          }
        }
      }
    }
    
    if (best_loss_improve==0 && false) {
      print("oh dear");
      
      print_matrix(delta_mig);
      print_vector(excess_out);
      print_vector(H);
      
      Rcpp::stop("foobar");
      
    }
    
    // modify migration matrix
    delta_mig[best_i][best_j] -= best_sub;
    
    // recalculate excess flow out of demes
    row_sums(rowsum, delta_mig);
    col_sums(colsum, delta_mig);
    total_loss = 0;
    for (int i=0; i<demes; i++) {
      excess_out[i] = rowsum[i] - colsum[i];
      total_loss += excess_out[i]*excess_out[i];
    }
    
    if (cutout>100) {
      
      print_matrix(delta_mig);
      print_vector(excess_out);
      print_vector(H);
      
      Rcpp::stop("foobar");
    }
    
  }
  
}


