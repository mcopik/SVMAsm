/*
 * AbstractKernel.h
 *
 *  Created on: Oct 12, 2012
 *      Author: mcopik
 */

#ifndef ABSTRACTKERNEL_H_
#define ABSTRACTKERNEL_H_

class Matrix;

class AbstractKernel {
public:
	virtual Matrix * cacheKernel(Matrix & X,Matrix & Y,double sigma) = 0;
	virtual double kernelFunction(int * x,int * y,int size,double sigma) = 0;
};


#endif /* ABSTRACTKERNEL_H_ */
