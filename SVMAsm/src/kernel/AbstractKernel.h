/*
 * AbstractKernel.h
 *
 *  Created on: Oct 12, 2012
 *      Author: mcopik
 */

#ifndef ABSTRACTKERNEL_H_
#define ABSTRACTKERNEL_H_
#include "../data/Matrix.h"
#include "../data/Vector.h"

template<class T>
class AbstractKernel {
public:
	virtual Matrix<T> cacheKernel(Matrix<T> & X) = 0;
	virtual double kernelFunction(T * x,T * y,int size) = 0;
	virtual double kernelFunction(T x,T y) = 0;
};


#endif /* ABSTRACTKERNEL_H_ */
