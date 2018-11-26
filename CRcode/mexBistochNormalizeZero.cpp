#include "mex.h"
#include "mexBistochNormalizeZero.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    enum{U, tol, maxIter, Step, OneMat, Xmat};
    enum{V};
       
    int N;
	int M;
    
    double* pU = mxGetPr(prhs[U]);
    double* pTol = mxGetPr(prhs[tol]);
    int* pMaxIter = (int*)mxGetPr(prhs[maxIter]);
    double* pStep = mxGetPr(prhs[Step]);
	double* pOneMat = mxGetPr(prhs[OneMat]);
	double* pXmat = mxGetPr(prhs[Xmat]);

    M = mxGetM(prhs[U]);
	N = mxGetN(prhs[U]);
        
    plhs[V] = mxCreateDoubleMatrix(N, N, mxREAL);
    double* pV = mxGetPr(plhs[V]);
            
    bistocNormZero(pU, N, pTol, pMaxIter, pStep, pOneMat, pXmat, pV);
}