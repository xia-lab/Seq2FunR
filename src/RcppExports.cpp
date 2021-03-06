// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// seq2funr
int seq2funr(std::string sampletable, std::string genemap, std::string tfmi, bool outputMappedCleanReads, bool profiling, std::string in1, std::string in2, std::string prefix, std::string mode, int mismatch, int minscore, int minlength, int maxtranslength, int nThreads, bool verbose);
RcppExport SEXP _seq2funR_seq2funr(SEXP sampletableSEXP, SEXP genemapSEXP, SEXP tfmiSEXP, SEXP outputMappedCleanReadsSEXP, SEXP profilingSEXP, SEXP in1SEXP, SEXP in2SEXP, SEXP prefixSEXP, SEXP modeSEXP, SEXP mismatchSEXP, SEXP minscoreSEXP, SEXP minlengthSEXP, SEXP maxtranslengthSEXP, SEXP nThreadsSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type sampletable(sampletableSEXP);
    Rcpp::traits::input_parameter< std::string >::type genemap(genemapSEXP);
    Rcpp::traits::input_parameter< std::string >::type tfmi(tfmiSEXP);
    Rcpp::traits::input_parameter< bool >::type outputMappedCleanReads(outputMappedCleanReadsSEXP);
    Rcpp::traits::input_parameter< bool >::type profiling(profilingSEXP);
    Rcpp::traits::input_parameter< std::string >::type in1(in1SEXP);
    Rcpp::traits::input_parameter< std::string >::type in2(in2SEXP);
    Rcpp::traits::input_parameter< std::string >::type prefix(prefixSEXP);
    Rcpp::traits::input_parameter< std::string >::type mode(modeSEXP);
    Rcpp::traits::input_parameter< int >::type mismatch(mismatchSEXP);
    Rcpp::traits::input_parameter< int >::type minscore(minscoreSEXP);
    Rcpp::traits::input_parameter< int >::type minlength(minlengthSEXP);
    Rcpp::traits::input_parameter< int >::type maxtranslength(maxtranslengthSEXP);
    Rcpp::traits::input_parameter< int >::type nThreads(nThreadsSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(seq2funr(sampletable, genemap, tfmi, outputMappedCleanReads, profiling, in1, in2, prefix, mode, mismatch, minscore, minlength, maxtranslength, nThreads, verbose));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_seq2funR_seq2funr", (DL_FUNC) &_seq2funR_seq2funr, 15},
    {NULL, NULL, 0}
};

RcppExport void R_init_seq2funR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
