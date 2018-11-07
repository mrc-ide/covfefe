
#include "genotype.h"
#include "probability.h"

using namespace std;
/*
//------------------------------------------------
// default constructor for genotype class
genotype::genotype() {
  n_x = 0;
}

//------------------------------------------------
// de novo infection
void genotype::de_novo(vector<vector<int>> &loci, int z) {
  
  // create new genotype with all loci equal to z
  int n_c = loci.size();
  x.clear();
  x.push_back(vector<vector<int>>(n_c));
  for (int c=0; c<n_c; c++) {
    x[0][c] = vector<int>(loci[c].size(), z);
  }
  n_x = 1;
  
}

//------------------------------------------------
// pass infection from source host to dest host via a mosquito
void genotype_infect(genotype &source, genotype &dest, vector<double> &dist_oocysts, vector<double> &dist_hepatocytes, vector<vector<int>> &loci, double recom_rate) {
  
  // if single infection then pass on this infection unchanged
  if (source.n_x==1) {
    dest.x.push_back(source.x[0]);
    dest.n_x++;
    return;
  }
  
  // draw number of oocysts in mosquito
  int n_oocysts = sample1(dist_oocysts);
  
  // draw of a pair of source genotypes for each oocyst
  vector<pair<int, int>> gametocytes(n_oocysts);
  for (int i=0; i<n_oocysts; i++) {
    gametocytes[i].first = sample2(0, source.n_x-1);
    gametocytes[i].second = sample2(0, source.n_x-1);
  }
  
  // draw number of infected hepatocytes
  int n_hepatocytes = sample1(dist_hepatocytes);
  
  // for each hepatocyte, draw an oocyst at random (call this draw r1). Let the
  // two gametocytes that make up this oocyst be g1 and g2. If g1==g2 then there
  // is no need to simulate recombination. Instead, increment
  // original_count[g1], which records how many times each of the original
  // source genotypes make it unchanged into hepatocytes.
  //
  // If g1!=g2 then draw a second value, r2, equally from 1 to 4. The value of 
  // recombinant_index[r1][r2] gives the index of a recombinant genotype (call 
  // this index i1). The value of recombinants[i1] then gives the full 
  // recombinant genotype. Finally, recombinant_count[i1] records how many times
  // this recombinant genotype makes it into hepatocytes.
  //
  // This method of indexing is slightly complex, but ensures that recombinant
  // genotypes are generated only once and at the last possible point they are
  // needed.
  vector<vector<vector<int>>> recombinants;
  map<int, map<int, int>> recombinant_index;
  vector<int> original_count(source.n_x);
  vector<int> recombinant_count;
  for (int i=0; i<n_hepatocytes; i++) {
    int r1 = sample2(0, n_oocysts-1);
    int g1 = gametocytes[r1].first;
    int g2 = gametocytes[r1].second;
    if (g1==g2) {
      original_count[g1]++;
    } else {
      int r2 = sample2(0, 3);
      if (recombinant_index[r1].count(r2)==0) {
        recombinants.push_back(recombine(source.x[g1], source.x[g2], loci, recom_rate));
        recombinant_count.push_back(1);
        recombinant_index[r1][r2] = recombinant_count.size()-1;
      } else {
        recombinant_count[recombinant_index[r1][r2]]++;
      }
    }
  }
  
  // add genotypes to target
  for (int i=0; i<source.n_x; i++) {
    if (original_count[i]>0) {
      dest.x.push_back(source.x[i]);
      dest.n_x++;
    }
  }
  for (int i=0; i<int(recombinant_count.size()); i++) {
    dest.x.push_back(recombinants[i]);
    dest.n_x++;
  }
  
}

//------------------------------------------------
// recombine two genotypes
vector<vector<int>> recombine(vector<vector<int>> &gen1, vector<vector<int>> &gen2, vector<vector<int>> &loci, double recom_rate) {
  
  // initialise return object
  int nc = gen1.size();
  vector<vector<int>> ret(nc);
  
  // loop through chromosomes
  for (int c=0; c<nc; c++) {
    int n_loci = gen1[c].size();
    
    // draw genotype to start copying from at random with equal probability
    int copy_from = sample2(0,1);
    
    // pos0 and pos1 represent the left and right positions of a chromosomal
    // chunk. index0 and index1 represent the corresponding index of the
    // elements immediately to the right of these positions.
    double pos0 = 0;
    double pos1 = 0;
    int index0 = 0;
    int index1 = 0;
    
    // generate chunks until index0 exceeds the number of loci
    while (index0<n_loci) {
      
      // swap which genotype we are copying from
      copy_from = 1 - copy_from;
      
      // draw position of new break point
      pos1 = pos0 + rexp1(recom_rate);
      
      // find index of locus to the right of this position
      vector<int>::iterator it = lower_bound(loci[c].begin(), loci[c].end(), pos1);
      index1 = int(it-loci[c].begin());
      
      // if no change in index then this chunk spans no loci. Update and
      // continue
      if (index1==index0) {
        pos0 = pos1;
        continue;
      }
      
      // cap index1 at maximum number of loci
      index1 = (index1<n_loci) ? index1 : n_loci;
      
      // copy elements index0 to index1-1 of selected genotype
      if (copy_from==0) {
        ret[c].insert(ret[c].end(), gen1[c].begin()+index0, gen1[c].begin()+index1);
      } else {
        ret[c].insert(ret[c].end(), gen2[c].begin()+index0, gen2[c].begin()+index1);
      }
      
      // update position and index
      pos0 = pos1;
      index0 = index1;
    }
    
  } // end loop over chromosomes
  
  // return
  return ret;
}

*/
