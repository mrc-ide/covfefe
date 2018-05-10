
#include <Rcpp.h>
#include <math.h>
#include "probability.h"
#include "misc.h"

using namespace std;

//------------------------------------------------
// draw from continuous uniform distribution on interval [0,1)
double runif_0_1() {
  return R::runif(0,1);
}

//------------------------------------------------
// draw from continuous uniform distribution on interval [a,b)
double runif1(const double a, const double b) {
  return R::runif(a,b);
}

//------------------------------------------------
// draw from Bernoulli(p) distribution
bool rbernoulli1(const double p) {
  return R::rbinom(1, p);
}

//------------------------------------------------
// draw from Binomial(n, p) distribution
int rbinom1(const int n, const double p) {
  return R::rbinom(n, p);
}

//------------------------------------------------
// draw from Multinomial(n, p) distribution
vector<int> rmultinom1(const int n, vector<double> &p) {
  int k = p.size();
  vector<int> ret(k);
  R::rmultinom(n, &p[0], k, &(ret)[0]);
  return ret;
}

//------------------------------------------------
// sample single value from given probability vector that sums to p_sum
int sample1(const vector<double> &p, const double p_sum) {
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
