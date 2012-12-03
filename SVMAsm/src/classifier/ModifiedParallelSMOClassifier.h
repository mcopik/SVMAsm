/*
 * ModifiedParallelSMOClassifier.h
 *
 *  Created on: Dec 3, 2012
 *      Author: mcopik
 */

#ifndef MODIFIEDPARALLELSMOCLASSIFIER_H_
#define MODIFIEDPARALLELSMOCLASSIFIER_H_

#include <iostream>
#include <algorithm>
#include <pthread.h>
#include "AbstractClassifier.h"
#include "../data/TrainedModel.h"
#include "../data/TrainData.h"
#include "../data/ParallelTrainData.h"
#include "../data/Vector.h"
#include "../data/Matrix.h"
#include <random>
#include <cmath>
#include <fstream>
#include <cstdio>

typedef void *(*threadFunction)(void*);//ParallelTrainData<float,float>*);
template<class T,class U>
class ModifiedParallelSMOClassifier: public AbstractClassifier<T,U> {
private:
	typedef T inputType;
	typedef U dataType;
	/**
	 * Pointer to kernel class.
	 */
	AbstractKernel<inputType> * kernel;
	/**
	 * Cached kernel values.
	 */
	Matrix<dataType> * cache = nullptr;
	/**
	 * Cached error function.
	 */
	Vector<dataType> * errorCache = nullptr;
	using AbstractClassifier<inputType,dataType>::C;
	using AbstractClassifier<inputType,dataType>::error;
	using AbstractClassifier<inputType,dataType>::epsilon;
	int numberOfThreads;
	ParallelTrainData<inputType,dataType> * threadsData;
public:
	using AbstractClassifier<T,U>::model;
	void train(TrainData<inputType> & train,AbstractKernel<inputType> & _kernel,
			bool cacheKernel = false) {

		if(cacheKernel) {
			std::cout << "Caching kernel before training is not supported in parallel SMO."
					<< std::endl;
			return;
		}
		model = new TrainedModel<inputType,dataType>
				(train.X,train.Y);
		model->cachedKernel = new Matrix<dataType>(train.X.rows(),train.X.rows());
		std::fill_n(model->cachedKernel->matrixData(),model->cachedKernel->rows()
				*model->cachedKernel->cols(),-1);
		kernel = &_kernel;
		std::cout << "Training model: " << std::endl;
		std::cout << "Number of examples: " <<  model->X.rows() << std::endl;
		std::cout << "Number of features: " <<  model->X.cols() << std::endl;

		errorCache = new Vector<dataType>(train.X.rows());
		dataType fHigh = 0.0;
		int iHigh = 0;
		dataType fLow = 0.0;
		int iLow = 0;
		dataType alphaHighOld = 0.0;
		dataType alphaLowOld = 0.0;
		/**
		 * f(i) = -y(i)
		 */
		for(unsigned int i = 0;i < errorCache->size();++i) {
			(*errorCache)(i) = -1*train.Y(i);
		}
		numberOfThreads = 2;
		int threadDataSize = model->X.rows()/2;
		int * iHighArray = new int[numberOfThreads];
		int * iLowArray = new int[numberOfThreads];
		pthread_t * threadID = new pthread_t[numberOfThreads];
		int * ret = new int[numberOfThreads];
		std::cout << "dlopen" << std::endl;
		threadsData = new ParallelTrainData<inputType,dataType>[numberOfThreads];
		for(int i = 0;i < numberOfThreads;++i) {
			threadsData[i].trainDataSize = threadDataSize;
			threadsData[i].yArray = model->Y.vectorData(threadDataSize*i);
			threadsData[i].alphaArray = model->alphas.vectorData(threadDataSize*i);
			threadsData[i].errorArray = errorCache->vectorData(threadDataSize*i);
			threadsData[i].cost = C;
			threadsData[i].error = error;
		}
		std::cout << threadsData[0].cost << " " << threadsData[0].error << std::endl;
		std::cout << "dlopen" << std::endl;
		void* handle = dlopen("./libcpp.so",RTLD_LAZY);
		assert(handle != nullptr);
	    // reset errors
		dlerror();
		threadFunction function = (threadFunction) dlsym(handle, "findHighLow");
		assert(function != nullptr);
		for(int i = 0;i < numberOfThreads;++i) {
			pthread_create(&(threadID[i]),NULL,function,(void*)(&threadsData[i]));//(void *(*)(void*))function,(void*)&threadsData[i]);
			pthread_join(threadID[i],nullptr);//(void**)&ret[i]);
		}
		iHigh = threadsData[0].iHigh;
		iLow = threadsData[0].iLow;
		for(int i = 1;i < numberOfThreads;++i) {
			if((*errorCache)(iHigh) > (*errorCache)(threadsData[i].iHigh)) {
				iHigh = threadsData[i].iHigh+threadDataSize*i;
			}
			if((*errorCache)(iLow) < (*errorCache)(threadsData[i].iLow)) {
				iLow = threadsData[i].iLow+threadDataSize*i;
			}
		}
		fLow = (*errorCache)(iLow);
		fHigh = (*errorCache)(iHigh);
		//getfHigh(fHigh,iHigh,threadsData,numberOfThreads);
		//getfLow(fLow,iLow,threadsData);
		alphaHighOld = model->alphas(iHigh);
		alphaLowOld = model->alphas(iLow);
		updateAlpha(iHigh,iLow);
		std::cout << "alphaHigh " << alphaHighOld << " -> " << model->alphas(iHigh) << std::endl;
		std::cout << "alphaLow " << alphaLowOld << " -> " << model->alphas(iLow) << std::endl;
		std::cout << "fHigh " << fHigh << " fLow " << fLow << std::endl;
		updateErrorCache(iHigh,iLow,alphaHighOld,alphaLowOld);
		//std::cout << "alphaHigh " << alphaHighOld << " -> " << model->alphas(iHigh) << std::endl;
		//std::cout << "alphaLow " << alphaLowOld << " -> " << model->alphas(iLow) << std::endl;
		//getfLow(fLow,iLow);
		//getfHigh(fHigh,iHigh);
		std::cout << "iHigh " << iHigh <<
			"	fHigh " << fHigh << " iLow "<< iLow << " fLow " << fLow << std::endl;
		//if(0 < -1*error)
		std::cout << "C " << C << " error " << error << std::endl;
		int nr_iter = 0;

		do {// {//(fHigh - fLow) < -1*error) {

			if(nr_iter == 558) {
				std::cout << "558" << std::endl;
			}
			updateErrorCache(iHigh,iLow,alphaHighOld,alphaLowOld);
			//getfHigh(fHigh,iHigh);
			//getfLow(fLow,iLow);

			for(int i = 0;i < numberOfThreads;++i) {
				pthread_create(&(threadID[i]),NULL,function,(void*)(&threadsData[i]));//(void *(*)(void*))function,(void*)&threadsData[i]);
				pthread_join(threadID[i],nullptr);//(void**)&ret[i]);
			}

			iHigh = threadsData[0].iHigh;
			iLow = threadsData[0].iLow;
			for(int i = 1;i < numberOfThreads;++i) {
				if((*errorCache)(iHigh) > (*errorCache)(threadsData[i].iHigh)) {
					iHigh = threadsData[i].iHigh+threadDataSize*i;
				}
				if((*errorCache)(iLow) < (*errorCache)(threadsData[i].iLow)) {
					iLow = threadsData[i].iLow+threadDataSize*i;
				}
			}
			fLow = (*errorCache)(iLow);
			fHigh = (*errorCache)(iHigh);
			//std::cout << "fHigh " << fHigh << " fLow " << fLow << std::endl;
			alphaHighOld = model->alphas(iHigh);
			alphaLowOld = model->alphas(iLow);
			updateAlpha(iHigh,iLow);
			//std::cout << "fHigh " << fHigh << " fLow " << fLow << std::endl;
			//std::cout << "iHigh " << iHigh << " alphaHigh " << alphaHighOld << " -> " << model->alphas(iHigh) << std::endl;
			//std::cout << "iLow " << iLow << " alphaLow " << alphaLowOld << " -> " << model->alphas(iLow) << std::endl;
			//std::cout << "yHigh" << model->Y(iHigh) << " yLow " << model->Y(iLow) << std::endl;
			//std::cout << nr_iter++ << std::endl;
			//getfLow(fLow,iLow);
			//getfHigh(fHigh,iHigh);
			//std::cout << "fHigh " << fHigh << " fLow " << fLow << std::endl;
			//std::cout << fHigh - fLow << std::endl;
		}while(fLow > (fHigh+2*error));
		model->b = (fLow+fHigh)/2;

	    dlclose(handle);
		delete iHighArray;
		delete iLowArray;
		delete threadID;
		delete ret;
	}
	/**
	 * Predicts values on test data.
	 */
	Vector<T> predict(Matrix<T> & X/*Matrix<T> & kernel*/) {

		Vector<T> predicts(X.rows());
		U E;
		for(unsigned int i = 0;i < X.rows();++i) {
			E = 0;
			//for(unsigned int k = 0;k < model->alphas.size;++k)
			//	E += model->alphas(k)*model->Y(k)*	\
			//		kernel->kernelFunction(model->X(k),X(i),model->X.cols);
			for(unsigned int k = 0;k < model->alphas.size();++k)
				E += model->alphas(k)*kernel->kernelFunction(X(i),model->X(k),model->X.cols())
				*model->Y(k);
			std::cout << i << " " << E << std::endl;
			E -= model->b;
			std::cout << i << " " << E << std::endl;
			if(E >= 0) {
				predicts(i) = 1;
			}
			else {
				predicts(i) = 0;
			}
		}
		return predicts;

	}
	/**
	 * Destructor.
	 */
	~ModifiedParallelSMOClassifier() {
		delete model;
		delete errorCache;
		delete threadsData;
	}
private:
	/**
	 * Computes new values of alpha.
	 * @param iHigh index of alphaHigh
	 * @param iLow index of alphaLow
	 */
	void updateAlpha(int iHigh,int iLow) {
		//compute boundaries for alphaLow
		dataType alphaLowLowerBound,alphaLowUpperBound;
		dataType eta = computeKernel(iHigh,iHigh) + computeKernel(iLow,iLow)
				- 2 * computeKernel(iHigh,iLow);
		dataType oldAlphaLow = model->alphas(iLow);
		if(model->Y(iHigh)*model->Y(iLow) < 0) {
			if(model->alphas(iLow)-model->alphas(iHigh) < 0) {
				alphaLowLowerBound = 0;
				alphaLowUpperBound = C + model->alphas(iLow)-model->alphas(iHigh);
			}
			else {
				alphaLowLowerBound = model->alphas(iLow)-model->alphas(iHigh);
				alphaLowUpperBound = C;
			}
		}
		else {
			if(model->alphas(iLow)+model->alphas(iHigh) < C) {
				alphaLowUpperBound = model->alphas(iLow)+model->alphas(iHigh);
				alphaLowLowerBound = 0;
			}
			else {
				alphaLowLowerBound = model->alphas(iLow)+model->alphas(iHigh) -C;
				alphaLowUpperBound = C;
			}
		}
		//compute alphaLow
		if(eta > 0) {
			model->alphas(iLow) += model->Y(iLow)*
					((*errorCache)(iHigh) - (*errorCache)(iLow))/eta;
			if(model->alphas(iLow) < alphaLowLowerBound) {
				model->alphas(iLow) = alphaLowLowerBound;
			}
			else if(model->alphas(iLow) > alphaLowUpperBound) {
				model->alphas(iLow) = alphaLowUpperBound;
			}
		}
		else {
			dataType slope = model->Y(iLow)*
					((*errorCache)(iHigh) - (*errorCache)(iLow));
			dataType delta = slope*(alphaLowUpperBound-alphaLowLowerBound);
			if(delta > 0) {
				if(slope > 0) {
					model->alphas(iLow) = alphaLowUpperBound;
				}
				else {
					model->alphas(iLow) = alphaLowUpperBound;
				}
			}
			else {
				model->alphas(iLow) = oldAlphaLow;
			}
		}
		//compute alphaHigh
		model->alphas(iHigh) += model->Y(iHigh)*model->Y(iLow)
				*(oldAlphaLow - model->alphas(iLow));
	}
	/**
	 * Updates error cache.
	 * @param iHigh index of modified alphaHigh
	 * @param iLow index of modified alphaLow
	 * @param highOld old value of alphaHigh
	 * @param lowOld old value of alphaLow
	 */
	void updateErrorCache(int iHigh,int iLow,dataType highOld,dataType lowOld) {
		Vector<dataType> & error = *errorCache;
		for(unsigned int i = 0;i < error.size();++i) {
			//std::cout << "error(" << i << ") " << error(i) <<
				//	" alph " << model->alphas(i) << " y "

					//<< model->Y(i) << std::endl;
			//std::cout << "KerneliHigh " << computeKernel(iHigh,i)
					//<< " KerneliLow " << computeKernel(iLow,i) << std::endl;
			error(i) += (model->alphas(iHigh) - highOld)*
					model->Y(iHigh)*computeKernel(iHigh,i) +
					(model->alphas(iLow) - lowOld)*
					model->Y(iLow)*computeKernel(iLow,i);
			//std::cout << "error(" << i << ") " << error(i) << std::endl;
		}
	}
	/**
	 * Computes kernel value for training examples a,b.
	 * @param a
	 * @param b
	 * @return kernel value
	 */
	U computeKernel(int a,int b) {
		if(cache != nullptr)
			return (*cache)(a,b);
		if((*model->cachedKernel)(a,b) == -1) {
			(*model->cachedKernel)(a,b) =
					kernel->kernelFunction(model->X(a),model->X(b),model->X.cols());
		}
		return (*model->cachedKernel)(a,b);
	}
};


#endif /* MODIFIEDPARALLELSMOCLASSIFIER_H_ */
