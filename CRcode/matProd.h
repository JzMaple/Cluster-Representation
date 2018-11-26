#include "string.h"

void matProd(double* output, double* A, double* B, int N){
	// memset(output,0,N*N);
	for (int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			for(int inner=0;inner<N;inner++){
				output[i*N+j] += A[i*N+inner]*B[inner*N+j];
			}
		}
	}
}