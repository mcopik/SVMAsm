/*
 * AbstractClassifier.h
 *
 *  Created on: Oct 12, 2012
 *      Author: mcopik
 */

#ifndef ABSTRACTCLASSIFIER_H_
#define ABSTRACTCLASSIFIER_H_
#include "../data/TrainData.h"
#include "../data/TrainedModel.h"
#include "../kernel/AbstractKernel.h"


template<class T,class U>
class AbstractClassifier {
protected:
	typedef T inputType;
	typedef U dataType;
	/**
	 * Cost parameter.
	 */
	dataType C;
	/**
	 * Converge condition parameter.
	 */
	dataType error;
	/**
	 * Standard computation/approximation error.
	 */
	dataType epsilon;
	/**
	 * Contains cached kernel values.
	 */
	Matrix<dataType> * cache;
public:
	/**
	 * Base constructor.
	 */
	AbstractClassifier() {
		C = 0;
		error = 0;
		epsilon = 0;
		model = nullptr;
		cache = nullptr;
	}
	/**
	 * Set value of approximation error.
	 * @param epsilon
	 */
	void setEpsilon(dataType epsilon) {
		this->epsilon = epsilon;
	}
	/**
	 * Set cost value.
	 * @param C
	 */
	void setC(dataType C) {
		this->C = C;
	}
	/**
	 * Set converge condition parameter.
	 * @param error
	 */
	void setError(dataType error) {
		this->error = error;
	}
	/**
	 * Set pointer to matrix with cached kernel values.
	 * @param cache
	 */
	void setCachedKernel(Matrix<U> & cache) {
		this->cache = &cache;
	}
	virtual void train(TrainData<inputType> &,AbstractKernel<inputType> &,bool=false) = 0;
	virtual Vector<T> predict(Matrix<inputType> &) = 0;
	virtual ~AbstractClassifier(){}
	TrainedModel<inputType,dataType> * model;
};


#endif /* ABSTRACTCLASSIFIER_H_ */
