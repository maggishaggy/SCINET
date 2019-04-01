#include <scinet.h>
#include <RcppArmadillo.h>

#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List constructNet_summary(mat A, mat net, uvec samples, int rand_sample_no, int rand_sample_size, int thread_no) {
	samples --;
	
	field<mat> stats = SCINET::constructNet_summary(A, net, samples, rand_sample_no, rand_sample_size, thread_no);
	
	List res;	
	res["mu"] = stats(0);		
	res["sigma_sq"] = stats(1);
		
	return res;	
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List constructNet(mat A, mat net, uvec samples, int thread_no) {
	samples--;
	
	field<mat> stats = SCINET::constructNet(A, net, samples, thread_no);
	
	List res;	
	res["subs"] = stats(0);		
	res["weights"] = stats(1);
		
	return res;	
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List constructNet_crossing(mat A, mat net, uvec rows1, uvec samples1, uvec rows2, uvec samples2, int rand_sample_no, int rand_sample_size, int thread_no) {

	rows1--; rows2--;
	samples1--; samples2--;
	
	field<mat> stats = SCINET::constructNet_crossing(A, net, rows1, samples1, rows2, samples2, rand_sample_no, rand_sample_size, thread_no);
	
	List res;	
	res["mu"] = stats(0);		
	res["sigma"] = stats(1);
		
	return res;	
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat identifyDiffEdges(mat mu1, mat sigma1_sq, int n1, mat mu2, mat sigma2_sq, int n2) {
  
	mat DiffNet = SCINET::identifyDiffEdges(mu1, sigma1_sq, n1, mu2, sigma2_sq, n2);
		
	return DiffNet;	
}

