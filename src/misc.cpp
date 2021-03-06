
#include "misc.h"

#include <math.h>

using namespace std;

//------------------------------------------------
// define matrix of integers using nested vectors
vector<vector<int>> matrix_int(int d1) {
  return vector<vector<int>>(d1);
}
vector<vector<int>> matrix_int(int d1, int d2, int x) {
  return vector<vector<int>>(d1, vector<int>(d2, x));
}

//------------------------------------------------
// define matrix of doubles using nested vectors
vector<vector<double>> matrix_double(int d1) {
  return vector<vector<double>>(d1);
}
vector<vector<double>> matrix_double(int d1, int d2, double x) {
  return vector<vector<double>>(d1, vector<double>(d2, x));
}

//------------------------------------------------
// define 3d array of integers using nested vectors
vector<vector<vector<int>>> array_int(int d1) {
  return vector<vector<vector<int>>>(d1);
}
vector<vector<vector<int>>> array_int(int d1, int d2) {
  return vector<vector<vector<int>>>(d1, vector<vector<int>>(d2));
}
vector<vector<vector<int>>> array_int(int d1, int d2, int d3, int x) {
  return vector<vector<vector<int>>>(d1, vector<vector<int>>(d2, vector<int>(d3, x)));
}

//------------------------------------------------
// define 4d array of integers using nested vectors
vector<vector<vector<vector<int>>>> array_int_4d(int d1) {
  return vector<vector<vector<vector<int>>>>(d1);
}
vector<vector<vector<vector<int>>>> array_int_4d(int d1, int d2) {
  return vector<vector<vector<vector<int>>>>(d1, vector<vector<vector<int>>>(d2));
}
vector<vector<vector<vector<int>>>> array_int_4d(int d1, int d2, int d3) {
  return vector<vector<vector<vector<int>>>>(d1, vector<vector<vector<int>>>(d2, vector<vector<int>>(d3)));
}
vector<vector<vector<vector<int>>>> array_int_4d(int d1, int d2, int d3, int d4, int x) {
  return vector<vector<vector<vector<int>>>>(d1, vector<vector<vector<int>>>(d2, vector<vector<int>>(d3, vector<int>(d4, x))));
}

//------------------------------------------------
// basic sum over elements in a vector (templated for different data types)
// sum
// DEFINED IN HEADER

//------------------------------------------------
// sum over rows of a matrix (templated for different data types)
// row_sums
// DEFINED IN HEADER

//------------------------------------------------
// sum over columns of a matrix (templated for different data types)
// col_sums
// DEFINED IN HEADER

//------------------------------------------------
// erases particular element of a vector using efficient method
// quick_erase
// DEFINED IN HEADER

//------------------------------------------------
// helper function for printing a single value (templated for different data types)
// print
// DEFINED IN HEADER

//------------------------------------------------
// helper function for printing contents of a vector (templated for different data types)
// print_vector
// DEFINED IN HEADER

//------------------------------------------------
// helper function for printing contents of a matrix (templated for different data types)
// print_matrix
// DEFINED IN HEADER

//------------------------------------------------
// helper function for printing contents of a 3D array (templated for different data types)
// print_array
// DEFINED IN HEADER

//------------------------------------------------
// print simple bar-graph composed of title followed by n stars
void print_stars(string title, int n) {
#ifdef RCPP_ACTIVE
  Rcpp::Rcout << title << " ";
  for (int i=0; i<n; i++) {
      Rcpp::Rcout << "*";
  }
  Rcpp::Rcout << "\n";
  R_FlushConsole();
#else
  std::cout << title << " ";
  for (int i=0; i<n; i++) {
    std::cout << "*";
  }
  std::cout << "\n";
#endif
}

//------------------------------------------------
// print "foo", with optional number e.g. "foo2"
void foo(int n) {
#ifdef RCPP_ACTIVE
  if (n==0) {
      Rcpp::Rcout << "foo\n";
  } else {
      Rcpp::Rcout << "foo" << n << "\n";
  }
  R_FlushConsole();
#else
  if (n==0) {
    std::cout << "foo\n";
  } else {
    std::cout << "foo" << n << "\n";
  }
#endif
}

//------------------------------------------------
// print "bar", with optional number e.g. "bar2"
void bar(int n) {
#ifdef RCPP_ACTIVE
  if (n==0) {
      Rcpp::Rcout << "bar\n";
  } else {
      Rcpp::Rcout << "bar" << n << "\n";
  }
  R_FlushConsole();
#else
  if (n==0) {
    std::cout << "bar\n";
  } else {
    std::cout << "bar" << n << "\n";
  }
#endif
}

//------------------------------------------------
// print "foobar", with optional number e.g. "foobar2"
void foobar(int n) {
#ifdef RCPP_ACTIVE
  if (n==0) {
      Rcpp::Rcout << "foobar\n";
  } else {
      Rcpp::Rcout << "foobar" << n << "\n";
  }
  R_FlushConsole();
#else
  if (n==0) {
    std::cout << "foobar\n";
  } else {
    std::cout << "foobar" << n << "\n";
  }
#endif
}

//------------------------------------------------
// analogue of R function seq() for integers
vector<int> seq_int(int from, int to, int by) {
  int n = floor((to-from)/double(by)) + 1;
  vector<int> ret(n,from);
  for (int i=1; i<n; i++) {
      from += by;
      ret[i] = from;
  }
  return ret;
}

//------------------------------------------------
// push back multiple values to vector
// push_back_multiple
// DEFINED IN HEADER

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::SEXP format to bool format.
int rcpp_to_bool(SEXP x) {
  return Rcpp::as<bool>(x);
}
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::SEXP format to int format.
int rcpp_to_int(SEXP x) {
  return Rcpp::as<int>(x);
}
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::SEXP format to double format.
double rcpp_to_double(SEXP x) {
  return Rcpp::as<double>(x);
}
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::SEXP format to string format.
string rcpp_to_string(SEXP x) {
  return Rcpp::as<string>(x);
}
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::List format to vector<bool> format.
vector<bool> rcpp_to_vector_bool(SEXP x) {
  return Rcpp::as<vector<bool> >(x);
}
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::List format to vector<int> format.
vector<int> rcpp_to_vector_int(SEXP x) {
    return Rcpp::as<vector<int> >(x);
}
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::List format to vector<double> format.
vector<double> rcpp_to_vector_double(SEXP x) {
    return Rcpp::as<vector<double> >(x);
}
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::List format to vector<vector<int>> format.
vector<vector<int>> rcpp_to_matrix_int(Rcpp::List x) {
    int nrow = int(x.size());
    vector< vector<int> > x_mat(nrow);
    for (int i=0; i<nrow; i++) {
        x_mat[i] = Rcpp::as<vector<int> >(x[i]);
    }
    return x_mat;
}
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::List format to vector<vector<double>> format.
vector<vector<double>> rcpp_to_matrix_double(Rcpp::List x) {
    int nrow = int(x.size());
    vector< vector<double> > x_mat(nrow);
    for (int i=0; i<nrow; i++) {
        x_mat[i] = Rcpp::as<vector<double> >(x[i]);
    }
    return x_mat;
}
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::List format to vector<vector<vector<int>>> format.
vector<vector<vector<int>>> rcpp_to_array_int(Rcpp::List x) {
    int n1 = int(x.size());
    vector< vector< vector<int> > > ret(n1);
    for (int i=0; i<n1; i++) {
        Rcpp::List x_i = x[i];
        int n2 = int(x_i.size());
        ret[i] = vector< vector<int> >(n2);
        for (int j=0; j<n2; j++) {
            ret[i][j] = Rcpp::as<vector<int> >(x_i[j]);
        }
    }
    return ret;
}
#endif

//------------------------------------------------
#ifdef RCPP_ACTIVE
// converts input from Rcpp::List format to vector<vector<vector<double>>> format.
vector<vector<vector<double>>> rcpp_to_array_double(Rcpp::List x) {
    int n1 = int(x.size());
    vector< vector< vector<double> > > ret(n1);
    for (int i=0; i<n1; i++) {
        Rcpp::List x_i = x[i];
        int n2 = int(x_i.size());
        ret[i] = vector< vector<double> >(n2);
        for (int j=0; j<n2; j++) {
            ret[i][j] = Rcpp::as<vector<double> >(x_i[j]);
        }
    }
    return ret;
}
#endif

