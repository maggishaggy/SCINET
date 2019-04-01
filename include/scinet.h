#ifndef DECODE_H
#define DECODE_H

//#define ARMA_DONT_USE_WRAPPER
#undef ARMA_BLAS_CAPITALS
#define ARMA_BLAS_UNDERSCORE
#define ARMA_64BIT_WORD
#define ARMA_BLAS_LONG_LONG

#include "armadillo"
using namespace arma;
using namespace std;

#include <omp.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <cmath>
#include <unordered_map>
#include <random>

void cumchi ( double *x, double *df, double *cum, double *ccum );
double r8_normal_01_cdf_inverse ( double p );
double nchoosek(int n, int k);
double r8_choose(int n, int k);

sp_mat read_from_mm(char *path);
sp_mat read_from_table(char *path);
sp_mat read_from_csv(char *path);
uvec read_uvec(char *path);

mat PRImpute(sp_mat profile, sp_mat net, uvec cols, int method, double alpha);

namespace SCINET {
	// Returns a network per condition
	field<mat> constructNet(mat A, mat net, uvec samples, int thread_no);
	
	// Returns mean and std stats for each edge within the group
	field<mat> constructNet_summary(mat A, mat net, uvec samples, int rand_sample_no, int rand_sample_size, int thread_no);
	
	// Performs differential network analysis
	mat identifyDiffEdges(mat mu1, mat sigma1_sq, int n1, mat mu2, mat sigma2_sq, int n2);
	
	// Tissue/cell-type cross-talk
	field<mat> constructNet_crossing(mat A, mat net, uvec rows1, uvec samples1, uvec rows2, uvec samples2, int rand_sample_no, int rand_sample_size, int thread_no);	
}

#endif
