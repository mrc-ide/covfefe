
#pragma once

#include <Rcpp.h>

//------------------------------------------------
// class defining genotype
class genotype {
  
public:
  
  // PUBLIC OBJECTS
  
  // parameters
  int duration;
  
  
  // PUBLIC FUNCTIONS
  
  // constructors
  genotype();
  genotype(std::vector<std::vector<int>> &loci);
  
};
