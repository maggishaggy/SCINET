#include "armaMex.hpp"
#include "scinet.h"

void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int i;

	mat mu1;
	int m = (int ) mxGetM(prhs[0]);
	int n = (int ) mxGetN(prhs[0]);
	if(mxIsSparse(prhs[0])) {
		mu1 = mat(armaGetSparseMatrix(prhs[0]));  
	}
	else {
		double* ptr = (double *) reinterpret_cast<double*>(mxGetPr(prhs[0])); 
		
		mu1 = (mat(ptr, m, n, true, true));		
	}


	mat sigma1;
	if(mxIsSparse(prhs[1])) {
		sigma1 = mat(armaGetSparseMatrix(prhs[1]));  
	}
	else {
		double* ptr = (double *) reinterpret_cast<double*>(mxGetPr(prhs[1])); 		
		sigma1 = (mat(ptr, m, n, true, true));		
	}
	

	int n1 = mxGetScalar(prhs[2]);
	

	mat mu2;
	if(mxIsSparse(prhs[3])) {
		mu2 = mat(armaGetSparseMatrix(prhs[3]));  
	}
	else {
		double* ptr = (double *) reinterpret_cast<double*>(mxGetPr(prhs[3])); 		
		mu2 = (mat(ptr, m, n, true, true));		
	}



	mat sigma2;
	if(mxIsSparse(prhs[4])) {
		sigma2 = mat(armaGetSparseMatrix(prhs[4]));  
	}
	else {
		double* ptr = (double *) reinterpret_cast<double*>(mxGetPr(prhs[4])); 		
		sigma2 = (mat(ptr, m, n, true, true));		
	}


	int n2 = mxGetScalar(prhs[5]);



	mat diffNet = SCINET::identifyDiffEdges(mu1, sigma1, n1, mu2, sigma2, n2);

	plhs[0] = mxCreateDoubleMatrix(diffNet.n_rows, diffNet.n_cols, mxREAL);
	memcpy(mxGetPr(plhs[0]), diffNet.memptr(), diffNet.n_elem * sizeof(double)); 	
}
