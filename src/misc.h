
#pragma once

#include <Rcpp.h>

//------------------------------------------------
// basic sum over elements in a vector (templated for different data types).
template<class TYPE>
TYPE sum(std::vector<TYPE> &x) {
  TYPE output = 0;
  for (int i=0; i<int(x.size()); i++) {
    output += x[i];
  }
  return output;
}

//------------------------------------------------
// sum over rows of a matrix (templated for different data types)
template<class TYPE>
void row_sums(std::vector<TYPE> &ret, std::vector<std::vector<TYPE>> &x) {
  std::fill(ret.begin(), ret.end(), 0);
  for (int i=0; i<int(x.size()); i++) {
    ret[i] = sum(x[i]);
  }
}

//------------------------------------------------
// sum over cols of a matrix (templated for different data types)
template<class TYPE>
void col_sums(std::vector<TYPE> &ret, std::vector<std::vector<TYPE>> &x) {
  std::fill(ret.begin(), ret.end(), 0);
  for (int i=0; i<int(x.size()); i++) {
    for (int j=0; j<int(x[i].size()); j++) {
      ret[j] += x[i][j];
    }
  }
}

//------------------------------------------------
// helper function for printing a single value or series of values (templated for different data types)
template<class TYPE>
void print(TYPE x) {
    Rcpp::Rcout << x << "\n";
    R_FlushConsole();
}
template<class TYPE1, class TYPE2>
void print(TYPE1 x1, TYPE2 x2) {
    Rcpp::Rcout << x1 << " " << x2 << "\n";
    R_FlushConsole();
}
template<class TYPE1, class TYPE2, class TYPE3>
void print(TYPE1 x1, TYPE2 x2, TYPE3 x3) {
    Rcpp::Rcout << x1 << " " << x2 << " " << x3 << "\n";
    R_FlushConsole();
}
template<class TYPE1, class TYPE2, class TYPE3, class TYPE4>
void print(TYPE1 x1, TYPE2 x2, TYPE3 x3, TYPE4 x4) {
    Rcpp::Rcout << x1 << " " << x2 << " " << x3 << " " << x4 << "\n";
    R_FlushConsole();
}
template<class TYPE1, class TYPE2, class TYPE3, class TYPE4, class TYPE5>
void print(TYPE1 x1, TYPE2 x2, TYPE3 x3, TYPE4 x4, TYPE5 x5) {
    Rcpp::Rcout << x1 << " " << x2 << " " << x3 << " " << x4 << " " << x5 << "\n";
    R_FlushConsole();
}
template<class TYPE1, class TYPE2, class TYPE3, class TYPE4, class TYPE5, class TYPE6>
void print(TYPE1 x1, TYPE2 x2, TYPE3 x3, TYPE4 x4, TYPE5 x5, TYPE6 x6) {
    Rcpp::Rcout << x1 << " " << x2 << " " << x3 << " " << x4 << " " << x5 << " " << x6 << "\n";
    R_FlushConsole();
}

//------------------------------------------------
// helper function for printing contents of a vector (templated for different data types)
template<class TYPE>
void print_vector(std::vector<TYPE> &x) {
    for (int i=0; i<x.size(); i++) {
        Rcpp::Rcout << x[i] << " ";
    }
    Rcpp::Rcout << "\n";
    R_FlushConsole();
}

//------------------------------------------------
// helper function for printing contents of a matrix (templated for different data types)
template<class TYPE>
void print_matrix(std::vector< std::vector<TYPE> > &x) {
    for (int i=0; i<x.size(); i++) {
        for (int j=0; j<x[i].size(); j++) {
            Rcpp::Rcout << x[i][j] << " ";
        }
        Rcpp::Rcout << "\n";
    }
    Rcpp::Rcout << "\n";
    R_FlushConsole();
}

//------------------------------------------------
// helper function for printing contents of a 3D array (templated for different data types)
template<class TYPE>
void print_array(std::vector< std::vector< std::vector<TYPE> > > &x) {
    for (int i=0; i<x.size(); i++) {
        Rcpp::Rcout << "--- slice " << i+1 << " ---\n";
        for (int j=0; j<x[i].size(); j++) {
            for (int k=0; k<x[i][j].size(); k++) {
                Rcpp::Rcout << x[i][j][k] << " ";
            }
            Rcpp::Rcout << "\n";
        }
        Rcpp::Rcout << "\n";
    }
    Rcpp::Rcout << "\n";
    R_FlushConsole();
}

//------------------------------------------------
// print simple bar-graph composed of title followed by n stars
void print_stars(std::string title="", int n=10);

//------------------------------------------------
// print "foo", with option number e.g. "foo2"
void foo(int n=0);

//------------------------------------------------
// print "bar", with option number e.g. "bar2"
void bar(int n=0);

//------------------------------------------------
// print "foobar", with option number e.g. "foobar2"
void foobar(int n=0);

//------------------------------------------------
// analogue of R function seq() for integers
std::vector<int> seq_int(int from, int to, int by=1);

//------------------------------------------------
// push back multiple values to vector
template<class TYPE>
void push_back_multiple(std::vector<TYPE> &lhs, std::vector<TYPE> &rhs) {
  lhs.insert(lhs.end(), rhs.begin(), rhs.end());
}

//------------------------------------------------
// converts input from Rcpp::SEXP format to bool format.
int rcpp_to_bool(SEXP x);

//------------------------------------------------
// converts input from Rcpp::SEXP format to int format.
int rcpp_to_int(SEXP x);

//------------------------------------------------
// converts input from Rcpp::SEXP format to double format.
double rcpp_to_double(SEXP x);

//------------------------------------------------
// converts input from Rcpp::List format to vector<int> format.
std::vector<int> rcpp_to_vector_int(SEXP x);

//------------------------------------------------
// converts input from Rcpp::List format to vector<double> format.
std::vector<double> rcpp_to_vector_double(SEXP x);

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<int>> format.
std::vector<std::vector<int>> rcpp_to_matrix_int(Rcpp::List x);

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<double>> format.
std::vector<std::vector<double>> rcpp_to_matrix_double(Rcpp::List x);

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<vector<int>>> format.
std::vector<std::vector<std::vector<int>>> rcpp_to_array_int(Rcpp::List x);

//------------------------------------------------
// converts input from Rcpp::List format to vector<vector<vector<double>>> format.
std::vector<std::vector<std::vector<double>>> rcpp_to_array_double(Rcpp::List x);



