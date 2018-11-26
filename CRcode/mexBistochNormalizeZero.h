#include "stdlib.h"
#include "string.h"
#include "stdio.h"
#include "matProd.h"

inline double UtruncLower(double U, double LowerBound){
	if (U < LowerBound){
		return LowerBound;
	}
	else{
		return U;
	}
}

inline double UtruncUpper(double U, double UpperBound){
	if (U > UpperBound){
		return UpperBound;
	}
	else{
		return U;
	}
}

inline double max(double a, double b){
	if (a>b){
		return a;
	}
	else{
		return b;
	}
}

inline double min(double a, double b){
	if (a<b){
		return a;
	}
	else{
		return b;
	}
}

void bistocNormZero(double* pU, int N,
		double* pTol,
        int* pMaxIter,
		double* pStep,
		double* pOneMat,
		double* pXmat,
        double* pV)
{
    double tol = (*pTol);
    tol = tol*tol;
    double delta = tol;
    
    int i, nG1, nG2, iter;
    // double *pX2 = new double[N*N];
    double *pTemp1 = new double[N*N];
	double *pTemp2 = new double[N*N];
	double *pTemp3 = new double[N*N];
    
	double *BoundLower = new double[N*N];
	double *BoundUpper = new double[N*N];
	for (int i=0; i<N; i++){
		for (int j=0; j<N; j++){
			BoundLower[i*N+j] = max(-pXmat[i*N+j]/(*pStep), -pXmat[j*N+i]/(*pStep));
			BoundUpper[i*N+j] = min((1-pXmat[i*N+j])/(*pStep), (1-pXmat[j*N+i])/(*pStep));
		}
	}

	for (int i=0; i<N; i++){
		for (int j=0; j<N; j++){
			pV[i*N+j] = pU[i*N+j];
		}
	}

    // double tempSum = 0;
        
    iter = 0;

    while(delta >= tol && iter < (*pMaxIter))
    {
        iter++;
        
		memset(pTemp1,0,N*N*sizeof(double));
		memset(pTemp2,0,N*N*sizeof(double));
		memset(pTemp3,0,N*N*sizeof(double));
		/**
		for (int i=0; i<N; i++){
			for (int j=0; j<N; j++){
				for (int inner=0; inner<N; inner++){
					pTemp1[i*N+j] = pV[i*N+inner]*pOneMat[inner*N+j];
				}
				pTemp1[i*N+j] = pTemp1[i*N+j]/double(N);
			}
		}**/
		matProd(pTemp1,pV,pOneMat,N);
		for (int i=0; i<N; i++){
			for (int j=0; j<N; j++){
				pTemp1[i*N+j] = pTemp1[i*N+j]/double(N);
			}
		}
		/**
		for (int i=0; i<N; i++){
			for (int j=0; j<N; j++){
				for (int inner=0; inner<N; inner++){
					pTemp2[i*N+j] += pOneMat[i*N+inner]*pV[inner*N+j];
				}
				pTemp2[i*N+j] = pTemp2[i*N+j]/double(N);
			}
		}**/
		matProd(pTemp2,pOneMat,pV,N);
		for (int i=0; i<N; i++){
			for (int j=0; j<N; j++){
				pTemp2[i*N+j] = pTemp2[i*N+j]/double(N);
			}
		}
		/**
		for (int i=0; i<N; i++){
			for (int j=0; j<N; j++){
				for (int inner=0; inner<N; inner++){
					pTemp3[i*N+j] += pTemp2[i*N+inner]*pOneMat[inner*N+j];
				}
				pTemp3[i*N+j] = pTemp3[i*N+j]/double(N);
			}
		}**/
		matProd(pTemp3,pOneMat,pTemp1,N);
		for (int i=0; i<N; i++){
			for (int j=0; j<N; j++){
				pTemp3[i*N+j] = pTemp3[i*N+j]/double(N);
			}
		}

		

		for (int i=0; i<N; i++){ // Bregman projection
			for (int j=0; j<N; j++){
				pV[i*N+j] = pV[i*N+j] - pTemp1[i*N+j] - pTemp2[i*N+j] + pTemp3[i*N+j];
			}
		}

		double temp = 0;
		for (int i=0; i<30; i++){
			temp+=pV[i];
		}
		int pp = 0;

		for (int i=0; i<N; i++){ // set diagonal elements 0
			pV[i*N+i] = 0;
		}
		for (int i=0; i<N; i++){
			for (int j=0; j<N; j++){
				if (pV[i*N+j]<BoundLower[i*N+j]){
					pV[i*N+j]=BoundLower[i*N+j];
				}
				else if (pV[i*N+j]>BoundUpper[i*N+j]){
					pV[i*N+j]=BoundUpper[i*N+j];
				}
				// pV[i*N+j] = UtruncLower(pV[i*N+j], BoundLower[i*N+j]);
				// pV[i*N+j] = UtruncUpper(pV[i*N+j], BoundUpper[i*N+j]);
			}
		}

		



        // check the difference for termination criterion

    }

    delete [] pTemp1;
	delete [] pTemp2;
	delete [] pTemp3;
	delete [] BoundLower;
	delete [] BoundUpper;
}

