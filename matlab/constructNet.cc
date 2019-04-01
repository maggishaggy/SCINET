#include "armaMex.hpp"
#include <scinet.h>

void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int i;

	mat A;
	int m = (int ) mxGetM(prhs[0]);
	int n = (int ) mxGetN(prhs[0]);
	if(mxIsSparse(prhs[0])) {
		A = mat(armaGetSparseMatrix(prhs[0]));  
	}
	else {
		double* ptr = (double *) reinterpret_cast<double*>(mxGetPr(prhs[0])); 
		
		A = mat(ptr, m, n, true, true);		
	}
	
	mat net;
	if((int ) mxGetM(prhs[1]) != m || (int ) mxGetN(prhs[1]) != m) {
		mexErrMsgTxt("network nodes should match genes in the expression profile.");
	}
	
	if(mxIsSparse(prhs[1])) {
		net = mat(armaGetSparseMatrix(prhs[1]));  
	}
	else {
		double* ptr = (double *) reinterpret_cast<double*>(mxGetPr(prhs[1])); 		
		net = mat(ptr, m, m, true, true);		
	}


	int m2 = (int ) mxGetM(prhs[2]);
	int n2 = (int ) mxGetN(prhs[2]);
	if(min(m2, n2) != 1) {
		mexErrMsgTxt("samples should be a one dimensional vector.");
	}	
	int vector_size = max(m2, n2);
	uvec samples(vector_size);
	if(mxIsSparse(prhs[2])) {
		mat temp = mat(armaGetSparseMatrix(prhs[2]));  
		for (i = 0; i < vector_size; i++) {
			samples(i) = (uword)(temp[i]-1);
		}		
	}
	else {
		double* ptr = (double *) reinterpret_cast<double*>(mxGetPr(prhs[2])); 
		for (i = 0; i < vector_size; i++) {
			samples(i) = (uword)(ptr[i]-1);
		}		
	}

	int thread_no = 8;
					
	if(nrhs > 3) {
		thread_no = mxGetScalar(prhs[3]);
	}


	field<mat> res = SCINET::constructNet(A, net, samples, thread_no);

	mat subs = res(0);
	plhs[0] = mxCreateDoubleMatrix(subs.n_rows, subs.n_cols, mxREAL);
	memcpy(mxGetPr(plhs[0]), subs.memptr(), subs.n_elem * sizeof(double)); 	

	mat weights = res(1);
	plhs[1] = mxCreateDoubleMatrix(weights.n_rows, weights.n_cols, mxREAL);
	memcpy(mxGetPr(plhs[1]), weights.memptr(), weights.n_elem * sizeof(double)); 		
}
