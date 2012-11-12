/*
 * TrainData.h
 *
 *  Created on: Oct 12, 2012
 *      Author: mcopik
 */

#ifndef TRAINDATA_H_
#define TRAINDATA_H_
#include "Matrix.h"
#include "Vector.h"

/**
 * Contains training data.
 */
template<class T>
class TrainData {
public:
	/**
	 * Training examples MxN
	 * M - number of examples
	 * N - number of features
	 */
	Matrix<T> & X;
	/**
	 * Examples classification Mx1
	 * M - number of examples
	 */
	Vector<T> & Y;
	TrainData(Matrix<T> & _X,Vector<T> & _Y): X(_X),Y(_Y) {}
};


#endif /* TRAINDATA_H_ */
