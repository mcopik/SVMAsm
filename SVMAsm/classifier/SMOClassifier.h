/*
 * SMOClassifier.h
 *
 *  Created on: Oct 12, 2012
 *      Author: mcopik
 */

#ifndef SMOCLASSIFIER_H_
#define SMOCLASSIFIER_H_
#include "AbstractClassifier.h"

template<class T>
class SMOClassifier: public AbstractClassifier<T> {
	using AbstractClassifier<T>::model;
public:
	/**
	 * Trains model on given data.
	 */
	void train(TrainData<T> & train,AbstractKernel<T> & kernel,double sigma) {
		model = new TrainedModel<T>
			(train.X,train.Y,*(kernel.cacheKernel(train.X,train.Y,sigma)));
	}
	/**
	 * Predicts values on test data.
	 */
	void predict(TrainData<T> & test) {

	}
	/**
	 * Destructor.
	 */
	~SMOClassifier() {
		delete model;
	}
};



#endif /* SMOCLASSIFIER_H_ */
