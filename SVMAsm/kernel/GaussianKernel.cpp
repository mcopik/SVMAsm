/*
 * GaussianKernel.cpp
 *
 *  Created on: Oct 12, 2012
 *      Author: mcopik
 */
#include <cmath>
#include "GaussianKernel.h"

Matrix * GaussianKernel::cacheKernel(Matrix & X,Matrix &, double sigma) {
	return nullptr;
}
double GaussianKernel::kernelFunction(int * x,int * y,int size,double sigma) {
	double sum = 0.0;
	for(int i = 0;i < size;++i) {
		sum -= pow(x[i]-y[i],2.0);
	}
	sum /= 2*pow(sigma,2.0);
	return exp(sum);
}



