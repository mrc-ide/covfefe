// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// indiv_sim_cpp
Rcpp::List indiv_sim_cpp(Rcpp::List args);
RcppExport SEXP _covfefe_indiv_sim_cpp(SEXP argsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type args(argsSEXP);
    rcpp_result_gen = Rcpp::wrap(indiv_sim_cpp(args));
    return rcpp_result_gen;
END_RCPP
}
// sim_genotypes_cpp
Rcpp::List sim_genotypes_cpp(Rcpp::List line_list, Rcpp::List args);
RcppExport SEXP _covfefe_sim_genotypes_cpp(SEXP line_listSEXP, SEXP argsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type line_list(line_listSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type args(argsSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_genotypes_cpp(line_list, args));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_covfefe_indiv_sim_cpp", (DL_FUNC) &_covfefe_indiv_sim_cpp, 1},
    {"_covfefe_sim_genotypes_cpp", (DL_FUNC) &_covfefe_sim_genotypes_cpp, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_covfefe(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
