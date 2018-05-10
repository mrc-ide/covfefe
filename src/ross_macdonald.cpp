
#include "ross_macdonald.h"
#include "probability.h"
#include "misc.h"

using namespace std;

//------------------------------------------------
// Draws from synchronous stochastic Ross-Macdonald model
// [[Rcpp::export]]
Rcpp::List ross_macdonald_cpp(Rcpp::List args) {
  
  // get input arguments
  int max_time = Rcpp::as<int>(args["max_time"]);
  double a = Rcpp::as<double>(args["a"]);
  double mu = Rcpp::as<double>(args["mu"]);
  int u = Rcpp::as<int>(args["u"]);
  int v = Rcpp::as<int>(args["v"]);
  double r = Rcpp::as<double>(args["r"]);
  double b = Rcpp::as<double>(args["b"]);
  double c = Rcpp::as<double>(args["c"]);
  int Eh_init = Rcpp::as<int>(args["Eh_init"]);
  int Ih_init = Rcpp::as<int>(args["Ih_init"]);
  int Ev_init = Rcpp::as<int>(args["Ev_init"]);
  int Iv_init = Rcpp::as<int>(args["Iv_init"]);
  int H = Rcpp::as<int>(args["H"]);
  int M = Rcpp::as<int>(args["M"]);
  
  // setup some initial parameters
  int Sh = H-Eh_init-Ih_init;
  int Eh = Eh_init;
  int Ih = Ih_init;
  int Sv = M-Ev_init-Iv_init;
  int Ev = Ev_init;
  int Iv = Iv_init;
  double H_inv = 1/double(H);
  vector<double> Sh_vec(max_time, -1); // -1 acts as an indicator that these values should be replaced with NA in the R function
  vector<double> Eh_vec(max_time, -1);
  vector<double> Ih_vec(max_time, -1);
  vector<double> Sv_vec(max_time, -1);
  vector<double> Ev_vec(max_time, -1);
  vector<double> Iv_vec(max_time, -1);
  Sh_vec[0] = Sh;
  Eh_vec[0] = Eh;
  Ih_vec[0] = Ih;
  Sv_vec[0] = Sv;
  Ev_vec[0] = Ev;
  Iv_vec[0] = Iv;
  
  // create vectors to store lag states
  vector<double> Eh_list(u+1);
  vector<double> Ev_list(v+1);
  Eh_list[0] = Eh;
  Ev_list[0] = Ev;
  int uBuffer_index = 0;
  int uBuffer_index_delay = 1; // note that uBuffer_index_delay is actually u_step steps BEHIND uBuffer_index, but due to the ring-buffer looping round this actually places it one step in front of uBuffer_index at all times.
  int vBuffer_index = 0;
  
  // initialise variables
  double rate_h_infection, rate_v_infection;
  double prob_h_infection, prob_v_infection_or_death, prob_v_infection;
  double prob_h_recovery = 1 - exp(-r);
  double prob_v_death = 1 - exp(-mu);
  int h_infection, h_recovery, v_infection_or_death, v_infection, v_death_Ev, v_death_Iv;
  int j2;
  
  // carry out simulation loop
  for (int i=1; i<max_time; i++) {
    
    // calculate rates of events
    rate_h_infection = a*b*Iv*H_inv; // human infection (move to Eh state)
    rate_v_infection = a*c*Ih*H_inv; // mosquito infection (move to Ev state)
    
    // convert to probabilities, allowing for competing hazards
    prob_h_infection = 1 - exp(-rate_h_infection); // human infection (move to Eh state)
    prob_v_infection_or_death = 1 - exp(-(rate_v_infection + mu)); // mosquito infection or death in Sv state (competing hazards)
    prob_v_infection = rate_v_infection/(rate_v_infection + mu); // relative rate of mosquito infection
    
    // update ring buffer indices
    uBuffer_index = (uBuffer_index==u) ? 0 : uBuffer_index+1;
    uBuffer_index_delay = (uBuffer_index_delay==u) ? 0 : uBuffer_index_delay+1;
    vBuffer_index = (vBuffer_index==v) ? 0 : vBuffer_index+1;
    
    // human events
    h_infection = rbinom1(Sh, prob_h_infection); // human infection (move to Eh state)
    h_recovery = rbinom1(Ih, prob_h_recovery); // human recovery
    Sh += -h_infection + h_recovery; // update Sh
    
    Eh_list[uBuffer_index] = h_infection; // add infecteds to Eh list
    Eh += Eh_list[uBuffer_index] - Eh_list[uBuffer_index_delay]; // update Eh
    
    Ih += Eh_list[uBuffer_index_delay] - h_recovery; // update Ih
    
    // mosquito events
    v_infection_or_death = rbinom1(Sv, prob_v_infection_or_death); // mosquito infection or death in Sv state (competing hazards)
    v_infection = rbinom1(v_infection_or_death, prob_v_infection); // mosquito infection
    v_death_Iv = rbinom1(Iv, prob_v_death); // mosquito death in Iv state
    Sv += -v_infection + v_death_Iv; // update Sv
    
    Ev_list[vBuffer_index] = v_infection; // add infecteds to Ev list
    Ev += v_infection; // add new infecteds to Ev
    
    j2 = vBuffer_index;
    for (int j=0; j<v; j++) { // loop through all previous entries in Ev list
      j2 = (j2==0) ? v : j2-1;
      if (Ev_list[j2]>0) {
        v_death_Ev = rbinom1(Ev_list[j2], prob_v_death); // mosquito death in this Ev state
        Ev_list[j2] -= v_death_Ev; // subtract mosquito deaths from this element
        Ev -= v_death_Ev; // subtract mosquito deaths from Ev counter
        Sv += v_death_Ev; // respawn mosquitoes in Sv state
      }
    }
    Ev -= Ev_list[j2]; // at this stage j2 is behind vBuffer_index by v_step steps. Subtract final Ev list entries from Ev counter as they transition from lag state
    
    Iv += Ev_list[j2] - v_death_Iv; // update Iv
    
    // store values
    Sh_vec[i] = Sh;
    Eh_vec[i] = Eh;
    Ih_vec[i] = Ih;
    Sv_vec[i] = Sv;
    Ev_vec[i] = Ev;
    Iv_vec[i] = Iv;
    
  }
  
  // return values
  return Rcpp::List::create(Rcpp::Named("Sh")=Sh_vec,
                            Rcpp::Named("Eh")=Eh_vec,
                            Rcpp::Named("Ih")=Ih_vec,
                            Rcpp::Named("Sv")=Sv_vec,
                            Rcpp::Named("Ev")=Ev_vec,
                            Rcpp::Named("Iv")=Iv_vec
  );
}
