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

template<class T>
class AbstractClassifier {
public:
	virtual void train(TrainData<T> &,AbstractKernel<T> &,double) = 0;
	virtual void predict(TrainData<T> &) = 0;
	virtual ~AbstractClassifier(){}
protected:
	TrainedModel<T> * model = nullptr;
};


#endif /* ABSTRACTCLASSIFIER_H_ */
