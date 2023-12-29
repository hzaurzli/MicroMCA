// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// fastPDist
NumericMatrix fastPDist(NumericMatrix Ar, NumericMatrix Br);
RcppExport SEXP _MicroMCA_fastPDist(SEXP ArSEXP, SEXP BrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Ar(ArSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Br(BrSEXP);
    rcpp_result_gen = Rcpp::wrap(fastPDist(Ar, Br));
    return rcpp_result_gen;
END_RCPP
}
// MCAStep1
List MCAStep1(NumericMatrix X);
RcppExport SEXP _MicroMCA_MCAStep1(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(MCAStep1(X));
    return rcpp_result_gen;
END_RCPP
}
// MCAStep2
List MCAStep2(NumericMatrix Z, NumericMatrix V, NumericVector Dc);
RcppExport SEXP _MicroMCA_MCAStep2(SEXP ZSEXP, SEXP VSEXP, SEXP DcSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type V(VSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Dc(DcSEXP);
    rcpp_result_gen = Rcpp::wrap(MCAStep2(Z, V, Dc));
    return rcpp_result_gen;
END_RCPP
}
// fastOrder
NumericMatrix fastOrder(NumericMatrix Ar, NumericMatrix Br);
RcppExport SEXP _MicroMCA_fastOrder(SEXP ArSEXP, SEXP BrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Ar(ArSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Br(BrSEXP);
    rcpp_result_gen = Rcpp::wrap(fastOrder(Ar, Br));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MicroMCA_fastPDist", (DL_FUNC) &_MicroMCA_fastPDist, 2},
    {"_MicroMCA_MCAStep1", (DL_FUNC) &_MicroMCA_MCAStep1, 1},
    {"_MicroMCA_MCAStep2", (DL_FUNC) &_MicroMCA_MCAStep2, 3},
    {"_MicroMCA_fastOrder", (DL_FUNC) &_MicroMCA_fastOrder, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_MicroMCA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
