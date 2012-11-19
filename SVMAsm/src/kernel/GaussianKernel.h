/*
 * GaussianKernel.h
 *
 *  Created on: Oct 12, 2012
 *      Author: mcopik
 */

#ifndef GAUSSIANKERNEL_H_
#define GAUSSIANKERNEL_H_
#include <cmath>
#include "AbstractKernel.h"

template<class T>
class GaussianKernel: public AbstractKernel<T> {
protected:
	double sigma;
public:
	GaussianKernel(double _sigma = 0):sigma(_sigma) {}
	void setSigma(double _sigma) {
		sigma = _sigma;
	}
	Matrix<T> cacheKernel(Matrix<T> & X) {
		Matrix<T> retval(X.rows,X.rows);
		Vector<T> temp(X.rows);
		/*
		 * temp = sum(X.^2,2);
		 */
		for(unsigned int i = 0;i < X.rows;++i) {
			temp.data[i] = 0.0;
			for(unsigned int j = 0;j < X.cols;++j)
				temp.data[i] += pow(X(i,j),2.0);
		}
		/*
		 * retval = temp + temp' - 2*X*X'
		 * retval = retval .^ kernelFunction(1,0)
		 */
		T a = 1,b = 0;
		double kernelValue = kernelFunction(&a,&b,1);
		Matrix<T> multiplied = X.multiplyByTranspose();
		for(unsigned int i = 0;i < X.rows;++i) {
			for(unsigned int j = 0;j < X.rows;++j) {
				retval(i,j) = temp(i) + temp(j) - 2*multiplied(i,j);
				retval(i,j) = pow(kernelValue,retval(i,j));
			}
		}
		return std::move(retval);
	}
	double kernelFunction(T * x,T * y,int size) {
		double sum = 0.0;
		for(int i = 0;i < size;++i) {
			sum -= pow(x[i]-y[i],2.0);
		}
		sum /= 2*pow(sigma,2.0);
		return exp(sum);
	}
	double kernelFunction(T x,T y) {
		double sum = -pow(x-y,2.0);
		sum /= 2*pow(sigma,2.0);
		return exp(sum);
	}
};


#endif /* GAUSSIANKERNEL_H_ */
