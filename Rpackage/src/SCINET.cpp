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
List constructNet(sp_mat A, mat net, IntegerVector samples, int rand_sample_no, int rand_sample_size, int thread_no) {
	uvec samples_uvec(samples.size());
	for (int i = 0; i < samples_uvec.n_elem; i++) {
		samples_uvec(i) = (uword)(samples(i)-1);
	}	  
	
	field<mat> stats = constructNet(A, net, samples_uvec, rand_sample_no, rand_sample_size, thread_no);
	
	List res;	
	res["mu"] = stats(0);		
	res["sigma_sq"] = stats(1);
		
	return res;	
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat constructDiffNet(mat mu1, mat sigma1_sq, int n1, mat mu2, mat sigma2_sq, int n2) {
  
	
	mat DiffNet = identifyDiffEdges(mu1, sigma1_sq, n1, mu2, sigma2_sq, n2);
		
	return DiffNet;	
}

