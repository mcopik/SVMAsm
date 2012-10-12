/*
 * GaussianKernel.h
 *
 *  Created on: Oct 12, 2012
 *      Author: mcopik
 */

#ifndef GAUSSIANKERNEL_H_
#define GAUSSIANKERNEL_H_
#include "AbstractKernel.h"

class Matrix;

class GaussianKernel: public AbstractKernel {

public:
	Matrix * cacheKernel(Matrix & X,Matrix & Y,double sigma);
	double kernelFunction(int * x,int * y,int size,double sigma);
};


#endif /* GAUSSIANKERNEL_H_ */
