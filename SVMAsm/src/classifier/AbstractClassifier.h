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
	U C;
	U error;
public:
	AbstractClassifier() {
		C = 0;
		error = 0;
		model = nullptr;
	}
	void setMaxPasses(int maxPasses) {
		this->maxPasses = maxPasses;
	}
	void setC(U C) {
		this->C = C;
	}
	void setError(U error) {
		this->error = error;
	}
	virtual void train(TrainData<T> &,AbstractKernel<T> &,bool=false) = 0;
	virtual Vector<T> predict(Matrix<T> &) = 0;
	virtual void setCachedKernel(Matrix<U> &) = 0;
	virtual ~AbstractClassifier(){}
	TrainedModel<T,U> * model;
};


#endif /* ABSTRACTCLASSIFIER_H_ */
