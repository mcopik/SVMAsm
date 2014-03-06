/*
 * TrainedModel.h
 *
 *  Created on: Oct 12, 2012
 *      Author: mcopik
 */

#ifndef TRAINEDMODEL_H_
#define TRAINEDMODEL_H_
#include <cstdlib>
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
	Matrix<U> * cachedKernel;
	U b = 0.0;
	Vector<U> alphas;
	//Vector<U> * w = nullptr;
	/**
	 * Constructor.
	 */
	TrainedModel(Matrix<T> & _x,Vector<T> & _y)
		:X(_x),Y(_y),cachedKernel(nullptr),alphas(_x.rows()) {
		//just to shut up Valgrind
		for(unsigned int i = 0;i < X.rows();++i)
			alphas(i) = 0;
	}
	/**
	 * Constructor.
	 * Matrix with cached kernel values is copied.
	 */
	TrainedModel(Matrix<T> & _x,Vector<T> & _y,Matrix<U> & kernel)
		:X(_x),Y(_y),cachedKernel(new Matrix<U>(kernel)),alphas(_x.rows()) {
		//just to shut up Valgrind
		for(unsigned int i = 0;i < X.rows();++i)
			alphas(i) = 0;
	}
	/**
	 * Constructor.
	 * Matrix with cached kernel values is moved.
	 */
	TrainedModel(Matrix<T> & _x,Vector<T> & _y,Matrix<U> && kernel)
		:X(_x),Y(_y),cachedKernel(new Matrix<U>(std::move(kernel)))
			,alphas(_x.rows()) {
		//just to shut up Valgrind
		for(unsigned int i = 0;i < X.rows();++i)
			alphas(i) = 0;
	}
	~TrainedModel() {
		//delete w;
		delete cachedKernel;
	}
};


#endif /* TRAINEDMODEL_H_ */
