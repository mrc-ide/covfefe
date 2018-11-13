
#include "probability.h"

#include <math.h>

using namespace std;

#ifndef RCPP_ACTIVE
//-- set random seed --
random_device rd;
default_random_engine generator(rd());

uniform_real_distribution<double> uniform_0_1(0.0,1.0);
#endif

//------------------------------------------------
// draw from continuous uniform distribution on interval [0,1)
double runif_0_1() {
#ifdef RCPP_ACTIVE
  return R::runif(0,1);
#else
  return uniform_0_1(generator);
#endif
}

//------------------------------------------------
// draw from continuous uniform distribution on interval [a,b)
double runif1(const double a, const double b) {
#ifdef RCPP_ACTIVE
  return R::runif(a,b);
#else
  uniform_real_distribution<double> uniform_a_b(a,b);
  return uniform_a_b(generator);
#endif
}

//------------------------------------------------
// draw from Bernoulli(p) distribution
bool rbernoulli1(const double p) {
#ifdef RCPP_ACTIVE
  return R::rbinom(1, p);
#else
  bernoulli_distribution dist_bernoulli(p);
  return dist_bernoulli(generator);
#endif
}

//------------------------------------------------
// draw from Binomial(n, p) distribution
int rbinom1(const int n, const double p) {
  if (n == 0 || p >= 1) {
    return n;
  }
#ifdef RCPP_ACTIVE
  return R::rbinom(n, p);
#else
  binomial_distribution<int> dist_binomial(n,p);
  return dist_binomial(generator);
#endif
}

//------------------------------------------------
// draw from Multinomial(n, p) distribution
vector<int> rmultinom1(const int n, vector<double> &p) {
  int k = int(p.size());
  vector<int> ret(k);
#ifdef RCPP_ACTIVE
  R::rmultinom(n, &p[0], k, &(ret)[0]);
#else
  int n_remaining = n;
  double p_remaining = sum(p);
  for (int i=0; i<(k-1); ++i) {
    ret[i] = rbinom1(n_remaining, p[i]/p_remaining);
    n_remaining -= ret[i];
    p_remaining -= p[i];
  }
  ret[k-1] = n_remaining;
#endif
  return ret;
}

//------------------------------------------------
// draw from Geometric(p) distribution, with mean (1-p)/p
int rgeom1(const double p) {
#ifdef RCPP_ACTIVE
  return R::rgeom(p);
#else
  geometric_distribution<int> dist_geom(p);
  return dist_geom(generator);
#endif
}

//------------------------------------------------
// draw from exponential(r) distribution
double rexp1(const double r) {
#ifdef RCPP_ACTIVE
  return R::rexp(1/r);
#else
  exponential_distribution<double> dist_exponential(r);
  return dist_exponential(generator);
#endif
}

//------------------------------------------------
// sample single value from given probability vector that sums to p_sum
int sample1(vector<double> &p, const double p_sum) {
  double rand = p_sum*runif_0_1();
  double z = 0;
  for (int i=0; i<int(p.size()); i++) {
    z += p[i];
    if (rand<z) {
      return i+1;
    }
  }
  return 0;
}

//------------------------------------------------
// equivalent to sample1, but ignores first values of p until min_val is reached
int sample1b(vector<double> &p, double p_sum, int min_val) {
  if (min_val == 1) {
    return sample1(p, p_sum);
  }
  for (int i=0; i<(min_val-1); ++i) {
    p_sum -= p[i];
  }
  double rand = p_sum*runif_0_1();
  double z = 0;
  for (int i=(min_val-1); i<int(p.size()); i++) {
    z += p[i];
    if (rand<z) {
      return i+1;
    }
  }
  return 0;
}

//------------------------------------------------
// sample single value x that lies between a and b (inclusive) with equal 
// probability. Works on positive or negative values of a or b, and works 
// irrespective of which of a or b is larger.
int sample2(const int a, const int b) {
  if (a<b) {
    return floor(runif1(a, b+1));
  } else {
    return floor(runif1(b, a+1));
  }
}

//------------------------------------------------
// re-shuffle the order of a vector of integers
void reshuffle(vector<int> &x) {
  int rnd1, tmp1;
  int n = int(x.size());
  for (int i=0; i<n; i++) {
    
    // draw random index from i to end of vector. Note that although runif
    // returns a double, by forcing to int we essentially round this value down
    // to nearest int.
    rnd1 = runif1(i,n);
    
    // temporarily store current value of vector at this position
    tmp1 = x[rnd1];
    
    // swap for value at position i
    x[rnd1] = x[i];
    x[i] = tmp1;
  }
}
