/*
 * TrainedModel.h
 *
 *  Created on: Oct 12, 2012
 *      Author: mcopik
 */

#ifndef TRAINEDMODEL_H_
#define TRAINEDMODEL_H_
#include "Vector.h"
#include "Matrix.h"

/**
 * Contains model.
 * T - type of training data.
 * U - float/double
 */
template<class T,class U>
class TrainedModel {
public:
	Matrix<T> & X;
	Vector<T> & Y;
	Matrix<U> & cachedKernel;
	double b = 0.0;
	Vector<U> * alphas = nullptr;
	Vector<U> * w = nullptr;
	TrainedModel(Matrix<T> & _x,Vector<T> & _y,Matrix<U> & kernel):X(_x),Y(_y), \
			cachedKernel(kernel){}
	~TrainedModel() {
		delete[] alphas;
		delete[] w;
	}
};


#endif /* TRAINEDMODEL_H_ */
