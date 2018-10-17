#include "armaMex.hpp"
#include "scDiffNet.h"

void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int i;

	// Check type of input.
	if ( (mxGetClassID(prhs[0]) != mxDOUBLE_CLASS) )
	mexErrMsgTxt("Input must me of type double.");

	// Check if input is real.
	if ( (mxIsComplex(prhs[0])) )
	mexErrMsgTxt("Input must be real.");

	
	sp_mat A;
	int m = (int ) mxGetM(prhs[0]);
	int n = (int ) mxGetN(prhs[0]);
	if(mxIsSparse(prhs[0])) {
		A = armaGetSparseMatrix(prhs[0]);  
	}
	else {
		double* ptr = (double *) reinterpret_cast<double*>(mxGetPr(prhs[0])); 
		
		A = sp_mat(mat(ptr, m, n, true, true));		
	}
	
		
	sp_mat net;
	if((int ) mxGetM(prhs[1]) != m || (int ) mxGetN(prhs[1]) != m) {
		mexErrMsgTxt("network nodes should match genes in the expression profile.");
	}
	
	if(mxIsSparse(prhs[1])) {
		net = armaGetSparseMatrix(prhs[1]);  
	}
	else {
		double* ptr = (double *) reinterpret_cast<double*>(mxGetPr(prhs[1])); 		
		net = sp_mat(mat(ptr, m, m, true, true));		
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

	// 1: bin/sym, 2: bin/asym, 3:sigmoid/sym, 4:sigmoid/asym
	int method = 1;				
	double alpha = 0.85;				

	if(nrhs > 3) {
		method = mxGetScalar(prhs[3]);
	}

	if(nrhs > 4) {
		alpha = mxGetScalar(prhs[4]);
	}

	mat updated_profile = PRImpute(A, net, samples, method, alpha);

	plhs[0] = mxCreateDoubleMatrix(updated_profile.n_rows, updated_profile.n_cols, mxREAL);
	memcpy(mxGetPr(plhs[0]), updated_profile.memptr(), updated_profile.n_elem * sizeof(double)); 		
}
