// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// graph_coloring_dsatur
IntegerVector graph_coloring_dsatur(ListOf<IntegerVector> adj_list);
RcppExport SEXP _graphcoloring_graph_coloring_dsatur(SEXP adj_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< ListOf<IntegerVector> >::type adj_list(adj_listSEXP);
    rcpp_result_gen = Rcpp::wrap(graph_coloring_dsatur(adj_list));
    return rcpp_result_gen;
END_RCPP
}
// graph_coloring_msc
IntegerVector graph_coloring_msc(ListOf<IntegerVector> adj_list);
RcppExport SEXP _graphcoloring_graph_coloring_msc(SEXP adj_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< ListOf<IntegerVector> >::type adj_list(adj_listSEXP);
    rcpp_result_gen = Rcpp::wrap(graph_coloring_msc(adj_list));
    return rcpp_result_gen;
END_RCPP
}
// graph_coloring_lmxrlf
IntegerVector graph_coloring_lmxrlf(ListOf<IntegerVector> adj_list);
RcppExport SEXP _graphcoloring_graph_coloring_lmxrlf(SEXP adj_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< ListOf<IntegerVector> >::type adj_list(adj_listSEXP);
    rcpp_result_gen = Rcpp::wrap(graph_coloring_lmxrlf(adj_list));
    return rcpp_result_gen;
END_RCPP
}
// graph_coloring_hybrid_dsatur_tabucol
IntegerVector graph_coloring_hybrid_dsatur_tabucol(ListOf<IntegerVector> adj_list);
RcppExport SEXP _graphcoloring_graph_coloring_hybrid_dsatur_tabucol(SEXP adj_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< ListOf<IntegerVector> >::type adj_list(adj_listSEXP);
    rcpp_result_gen = Rcpp::wrap(graph_coloring_hybrid_dsatur_tabucol(adj_list));
    return rcpp_result_gen;
END_RCPP
}
// graph_coloring_hybrid_lmxrlf_tabucol
IntegerVector graph_coloring_hybrid_lmxrlf_tabucol(ListOf<IntegerVector> adj_list);
RcppExport SEXP _graphcoloring_graph_coloring_hybrid_lmxrlf_tabucol(SEXP adj_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< ListOf<IntegerVector> >::type adj_list(adj_listSEXP);
    rcpp_result_gen = Rcpp::wrap(graph_coloring_hybrid_lmxrlf_tabucol(adj_list));
    return rcpp_result_gen;
END_RCPP
}
// graph_coloring_tabucol
IntegerVector graph_coloring_tabucol(ListOf<IntegerVector> adj_list, int k, int tabu_size, int rep, int nbmax);
RcppExport SEXP _graphcoloring_graph_coloring_tabucol(SEXP adj_listSEXP, SEXP kSEXP, SEXP tabu_sizeSEXP, SEXP repSEXP, SEXP nbmaxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< ListOf<IntegerVector> >::type adj_list(adj_listSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type tabu_size(tabu_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type rep(repSEXP);
    Rcpp::traits::input_parameter< int >::type nbmax(nbmaxSEXP);
    rcpp_result_gen = Rcpp::wrap(graph_coloring_tabucol(adj_list, k, tabu_size, rep, nbmax));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_graphcoloring_graph_coloring_dsatur", (DL_FUNC) &_graphcoloring_graph_coloring_dsatur, 1},
    {"_graphcoloring_graph_coloring_msc", (DL_FUNC) &_graphcoloring_graph_coloring_msc, 1},
    {"_graphcoloring_graph_coloring_lmxrlf", (DL_FUNC) &_graphcoloring_graph_coloring_lmxrlf, 1},
    {"_graphcoloring_graph_coloring_hybrid_dsatur_tabucol", (DL_FUNC) &_graphcoloring_graph_coloring_hybrid_dsatur_tabucol, 1},
    {"_graphcoloring_graph_coloring_hybrid_lmxrlf_tabucol", (DL_FUNC) &_graphcoloring_graph_coloring_hybrid_lmxrlf_tabucol, 1},
    {"_graphcoloring_graph_coloring_tabucol", (DL_FUNC) &_graphcoloring_graph_coloring_tabucol, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_graphcoloring(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
