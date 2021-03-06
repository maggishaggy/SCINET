// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <../../inst/include/scinet.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// constructNet_summary
List constructNet_summary(mat A, mat net, uvec samples, int rand_sample_no, int rand_sample_size, int thread_no);
RcppExport SEXP _SCINET_constructNet_summary(SEXP ASEXP, SEXP netSEXP, SEXP samplesSEXP, SEXP rand_sample_noSEXP, SEXP rand_sample_sizeSEXP, SEXP thread_noSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< mat >::type net(netSEXP);
    Rcpp::traits::input_parameter< uvec >::type samples(samplesSEXP);
    Rcpp::traits::input_parameter< int >::type rand_sample_no(rand_sample_noSEXP);
    Rcpp::traits::input_parameter< int >::type rand_sample_size(rand_sample_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type thread_no(thread_noSEXP);
    rcpp_result_gen = Rcpp::wrap(constructNet_summary(A, net, samples, rand_sample_no, rand_sample_size, thread_no));
    return rcpp_result_gen;
END_RCPP
}
// constructNet
List constructNet(mat A, mat net, uvec samples, int thread_no);
RcppExport SEXP _SCINET_constructNet(SEXP ASEXP, SEXP netSEXP, SEXP samplesSEXP, SEXP thread_noSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< mat >::type net(netSEXP);
    Rcpp::traits::input_parameter< uvec >::type samples(samplesSEXP);
    Rcpp::traits::input_parameter< int >::type thread_no(thread_noSEXP);
    rcpp_result_gen = Rcpp::wrap(constructNet(A, net, samples, thread_no));
    return rcpp_result_gen;
END_RCPP
}
// constructNet_crossing
List constructNet_crossing(mat A, mat net, uvec rows1, uvec samples1, uvec rows2, uvec samples2, int rand_sample_no, int rand_sample_size, int thread_no);
RcppExport SEXP _SCINET_constructNet_crossing(SEXP ASEXP, SEXP netSEXP, SEXP rows1SEXP, SEXP samples1SEXP, SEXP rows2SEXP, SEXP samples2SEXP, SEXP rand_sample_noSEXP, SEXP rand_sample_sizeSEXP, SEXP thread_noSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< mat >::type net(netSEXP);
    Rcpp::traits::input_parameter< uvec >::type rows1(rows1SEXP);
    Rcpp::traits::input_parameter< uvec >::type samples1(samples1SEXP);
    Rcpp::traits::input_parameter< uvec >::type rows2(rows2SEXP);
    Rcpp::traits::input_parameter< uvec >::type samples2(samples2SEXP);
    Rcpp::traits::input_parameter< int >::type rand_sample_no(rand_sample_noSEXP);
    Rcpp::traits::input_parameter< int >::type rand_sample_size(rand_sample_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type thread_no(thread_noSEXP);
    rcpp_result_gen = Rcpp::wrap(constructNet_crossing(A, net, rows1, samples1, rows2, samples2, rand_sample_no, rand_sample_size, thread_no));
    return rcpp_result_gen;
END_RCPP
}
// identifyDiffEdges
mat identifyDiffEdges(mat mu1, mat sigma1_sq, int n1, mat mu2, mat sigma2_sq, int n2);
RcppExport SEXP _SCINET_identifyDiffEdges(SEXP mu1SEXP, SEXP sigma1_sqSEXP, SEXP n1SEXP, SEXP mu2SEXP, SEXP sigma2_sqSEXP, SEXP n2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< mat >::type mu1(mu1SEXP);
    Rcpp::traits::input_parameter< mat >::type sigma1_sq(sigma1_sqSEXP);
    Rcpp::traits::input_parameter< int >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< mat >::type mu2(mu2SEXP);
    Rcpp::traits::input_parameter< mat >::type sigma2_sq(sigma2_sqSEXP);
    Rcpp::traits::input_parameter< int >::type n2(n2SEXP);
    rcpp_result_gen = Rcpp::wrap(identifyDiffEdges(mu1, sigma1_sq, n1, mu2, sigma2_sq, n2));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SCINET_constructNet_summary", (DL_FUNC) &_SCINET_constructNet_summary, 6},
    {"_SCINET_constructNet", (DL_FUNC) &_SCINET_constructNet, 4},
    {"_SCINET_constructNet_crossing", (DL_FUNC) &_SCINET_constructNet_crossing, 9},
    {"_SCINET_identifyDiffEdges", (DL_FUNC) &_SCINET_identifyDiffEdges, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_SCINET(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
