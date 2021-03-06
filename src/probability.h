
#pragma once

#include "misc.h"

#include <random>

//------------------------------------------------
// draw from continuous uniform distribution on interval [0,1)
double runif_0_1();

//------------------------------------------------
// draw from continuous uniform distribution on interval [a,b)
double runif1(const double a=0, const double b=1.0);

//------------------------------------------------
// draw from Bernoulli(p) distribution
bool rbernoulli1(const double p);

//------------------------------------------------
// draw from Binomial(n, p) distribution
int rbinom1(const int n, const double p);

//------------------------------------------------
// draw from Multinomial(n, p) distribution
std::vector<int> rmultinom1(const int n, std::vector<double> &p);

//------------------------------------------------
// draw from Geometric(p) distribution
int rgeom1(const double p);

//------------------------------------------------
// draw from exponential(r) distribution
double rexp1(const double r);

//------------------------------------------------
// sample single value from given probability vector that sums to p_sum
int sample1(std::vector<double> &p, const double pSum=1.0);

//------------------------------------------------
// equivalent to sample1, but ignores first values of p until min_val is reached
int sample1b(std::vector<double> &p, double p_sum, int min_val);

//------------------------------------------------
// sample single value x that lies between a and b (inclusive) with equal 
// probability. Works on positive or negative values of a or b, and works 
// irrespective of which of a or b is larger.
int sample2(const int a, const int b);

//------------------------------------------------
// re-shuffle the order of a vector of integers
void reshuffle(std::vector<int> &x);

