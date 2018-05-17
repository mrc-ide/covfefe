
#pragma once

#include <Rcpp.h>

//------------------------------------------------
// class defining genotype
class genotype {
  
public:
  
  // PUBLIC OBJECTS
  
  // parameters
  std::vector<std::vector<std::vector<int>>> x;
  int n_x;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  genotype();
  
  // other functions
  void de_novo(std::vector<std::vector<int>> &loci, int z);
  
};

//------------------------------------------------
// pass infection from source host to dest host via a mosquito
void genotype_infect(genotype &source, genotype &dest, std::vector<std::vector<int>> &loci, double recom_rate);

//------------------------------------------------
// recombine two genotypes
std::vector<std::vector<int>> recombine(std::vector<std::vector<int>> &gen1, std::vector<std::vector<int>> &gen2, std::vector<std::vector<int>> &loci, double recom_rate);

