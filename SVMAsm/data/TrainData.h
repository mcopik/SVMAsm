/*
 * TrainData.h
 *
 *  Created on: Oct 12, 2012
 *      Author: mcopik
 */

#ifndef TRAINDATA_H_
#define TRAINDATA_H_
#include "Matrix.h"

/**
 * Contains training data.
 */
class TrainData {
public:
	/**
	 * Training examples MxN
	 * M - number of examples
	 * N - number of features
	 */
	Matrix & X;
	/**
	 * Examples classification Mx1
	 * M - number of examples
	 */
	Matrix & Y;
	TrainData(Matrix & _X,Matrix & _Y): X(_X),Y(_Y) {}
};


#endif /* TRAINDATA_H_ */
