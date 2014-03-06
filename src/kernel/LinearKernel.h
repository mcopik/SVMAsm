/*
 * LinearKernel.h
 *
 *  Created on: Nov 12, 2012
 *      Author: mcopik
 */

#ifndef LINEARKERNEL_H_
#define LINEARKERNEL_H_
#include "../data/Matrix.h"
#include "AbstractKernel.h"

template<class T>
class LinearKernel: public AbstractKernel<T> {
public:
	Matrix<T> cacheKernel(Matrix<T> & X) {
		Matrix<T> retval(X.rows(),X.rows());
		retval = X.multiplyByTranspose(X);
		return std::move(retval);
	}
	double kernelFunction(T * x,T * y,int size) {
		double sum = 0.0;
		for(int i = 0;i < size;++i) {
			sum += x[i]*y[i];
		}
		return sum;
	}
	double kernelFunction(T x,T y) {
		return x*y;
	}
};


#endif /* LINEARKERNEL_H_ */
