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
	Matrix<dataType> * cachedKernel = nullptr;
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
		cachedKernel = new Matrix<dataType>(train.X.rows(),train.X.rows());
		std::fill_n(cachedKernel->matrixData(),cachedKernel->rows()
				*cachedKernel->cols(),-1);
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
		int threadDataSize = model->X.rows()/numberOfThreads;
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
			threadsData[i].kernel = kernel;
			threadsData[i].X = &model->X;
			threadsData[i].cachedKernel = cachedKernel;
			threadsData[i].threadID = i;
		}
		std::cout << threadsData[0].cost << " " << threadsData[0].error << std::endl;
		std::cout << "dlopen" << std::endl;
		//void* handle = dlopen("./libcpp.so",RTLD_LAZY);
		//assert(handle != nullptr);
	    // reset errors
		//dlerror();
		//threadFunction function = (threadFunction) dlsym(handle, "findHighLow");
		//assert(function != nullptr);
		for(int i = 0;i < numberOfThreads;++i) {
			pthread_create(&(threadID[i]),NULL,(void* (*)(void*))findHighLow,(void*)(&threadsData[i]));//(void *(*)(void*))function,(void*)&threadsData[i]);
			pthread_join(threadID[i],nullptr);//(void**)&ret[i]);
		}
		iHigh = threadsData[0].iHigh;
		iLow = threadsData[0].iLow;
		for(int i = 1;i < numberOfThreads;++i) {
			if((*errorCache)(iHigh) > (*errorCache)(threadsData[i].iHigh+threadDataSize*i)) {
				iHigh = threadsData[i].iHigh+threadDataSize*i;
			}
			if((*errorCache)(iLow) < (*errorCache)(threadsData[i].iLow+threadDataSize*i)) {
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
		//updateErrorCache(iHigh,iLow,alphaHighOld,alphaLowOld);
		//std::cout << "alphaHigh " << alphaHighOld << " -> " << model->alphas(iHigh) << std::endl;
		//std::cout << "alphaLow " << alphaLowOld << " -> " << model->alphas(iLow) << std::endl;
		//getfLow(fLow,iLow);
		//getfHigh(fHigh,iHigh);
		std::cout << "iHigh " << iHigh <<
			"	fHigh " << fHigh << " iLow "<< iLow << " fLow " << fLow << std::endl;
		//if(0 < -1*error)
		std::cout << "C " << C << " error " << error << std::endl;
		int nr_iter = 0;

		do {
			for(int i = 0;i < numberOfThreads;++i) {
				threadsData[i].errorUpdateHigh = model->Y(iHigh)*
						(model->alphas(iHigh) - alphaHighOld);
				threadsData[i].errorUpdateLow = model->Y(iLow)*
						(model->alphas(iLow) - alphaLowOld);
				threadsData[i].iHigh = iHigh;
				threadsData[i].iLow = iLow;
			}
			for(int i = 0;i < numberOfThreads;++i) {
				pthread_create(&(threadID[i]),NULL,(void* (*)(void*))updateErrorCache,(void*)(&threadsData[i]));//(void *(*)(void*))function,(void*)&threadsData[i]);
				pthread_join(threadID[i],nullptr);//(void**)&ret[i]);
			}
			//updateErrorCache(iHigh,iLow,alphaHighOld,alphaLowOld);
			std::cout << "updatefLow " << (*errorCache)(iLow) << " updatefHigh " << (*errorCache)(iHigh) << std::endl;
			for(int i = 0;i < numberOfThreads;++i) {
				pthread_create(&(threadID[i]),NULL,(void* (*)(void*))findHighLow,(void*)(&threadsData[i]));//(void *(*)(void*))function,(void*)&threadsData[i]);
				pthread_join(threadID[i],nullptr);//(void**)&ret[i]);
			}

			iHigh = threadsData[0].iHigh;
			iLow = threadsData[0].iLow;
			for(int i = 1;i < numberOfThreads;++i) {
				if((*errorCache)(iHigh) > (*errorCache)(threadsData[i].iHigh+threadDataSize*i)) {
					iHigh = threadsData[i].iHigh+threadDataSize*i;
				}
				if((*errorCache)(iLow) < (*errorCache)(threadsData[i].iLow+threadDataSize*i)) {
					iLow = threadsData[i].iLow+threadDataSize*i;
				}
			}
			fLow = (*errorCache)(iLow);
			fHigh = (*errorCache)(iHigh);
			alphaHighOld = model->alphas(iHigh);
			alphaLowOld = model->alphas(iLow);
			updateAlpha(iHigh,iLow);
			std::cout << "fHigh " << fHigh << " fLow " << fLow << std::endl;
			std::cout << "iHigh " << iHigh << " alphaHigh " << alphaHighOld << " -> " << model->alphas(iHigh) << std::endl;
			std::cout << "iLow " << iLow << " alphaLow " << alphaLowOld << " -> " << model->alphas(iLow) << std::endl;
			//std::cout << "yHigh" << model->Y(iHigh) << " yLow " << model->Y(iLow) << std::endl;
			//std::cout << nr_iter++ << std::endl;
			//getfLow(fLow,iLow);
			//getfHigh(fHigh,iHigh);
			//std::cout << "fHigh " << fHigh << " fLow " << fLow << std::endl;
			//std::cout << fHigh - fLow << std::endl;
		}while(fLow > (fHigh+2*error));
		model->b = (fLow+fHigh)/2;

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
		delete cachedKernel;
		delete errorCache;
		delete threadsData;
	}
private:
	static void findHighLow(ParallelTrainData<inputType,dataType> * x) {
		inputType fHigh = 1000000000;
		dataType fLow = -1000000000;
		int iLow = 0;
		int iHigh = 0;
		for(unsigned int i = 0;i < x->trainDataSize;++i) {
			if((x->alphaArray[i] > x->error && x->alphaArray[i] < (x->cost-x->error)) ||
					(equalsWithTolerance(x->alphaArray[i],0) && x->yArray[i] > 0) ||
					(equalsWithTolerance(x->alphaArray[i],x->cost) && x->yArray[i] < 0)) {
				if(x->errorArray[i] < fHigh) {
					fHigh = x->errorArray[i];
					iHigh = i;
				}
			}
			if((x->alphaArray[i] > x->error && x->alphaArray[i] < (x->cost-x->error)) ||
					(equalsWithTolerance(x->alphaArray[i],0) && x->yArray[i] < 0) ||
					(equalsWithTolerance(x->alphaArray[i],x->cost) && x->yArray[i] > 0)) {
				if(x->errorArray[i] > fLow) {
					fLow = x->errorArray[i];
					iLow = i;
				}
			}
		}
		x->iHigh = iHigh;
		x->iLow = iLow;
	}
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
	static void updateErrorCache(ParallelTrainData<inputType,dataType> * data) {
		unsigned int size = data->trainDataSize;
		AbstractKernel<inputType> * kernel = data->kernel;
		Matrix<dataType> * X = data->X;
		Matrix<dataType> * cachedKernel = data->cachedKernel;
		unsigned int iHigh = data->iHigh;
		unsigned int iLow = data->iLow;
		unsigned int offset = size*data->threadID;
		dataType diff = 0;
		for(unsigned int i = 0;i < size;++i) {
			diff = 0;
			diff += data->errorUpdateHigh*
					computeKernel(kernel,cachedKernel,*X,iHigh,i+offset);
			diff += data->errorUpdateLow*
					computeKernel(kernel,cachedKernel,*X,iLow,i+offset);
			data->errorArray[i] += diff;
		}
	}
	/**
	 * Updates error cache.
	 * @param iHigh index of modified alphaHigh
	 * @param iLow index of modified alphaLow
	 * @param highOld old value of alphaHigh
	 * @param lowOld old value of alphaLow
	 */
/*
	void updateErrorCache(int iHigh,int iLow,dataType highOld,dataType lowOld) {
		Vector<dataType> & error = *errorCache;
		for(unsigned int i = 0;i < error.size();++i) {
			//std::cout << "error(" << i << ") " << error(i) <<
				//	" alph " << model->alphas(i) << " y "

					//<< model->Y(i) << std::endl;
			//std::cout << "KerneliHigh " << computeKernel(iHigh,i)
					//<< " KerneliLow " << computeKernel(iLow,i) << std::endl;
			error(i) += (model->alphas(iHigh) - highOld)*
					model->Y(iHigh)*computeKernel(kernel,cachedKernel,model->X,iHigh,i) +
					(model->alphas(iLow) - lowOld)*
					model->Y(iLow)*computeKernel(kernel,cachedKernel,model->X,iLow,i);
			//std::cout << "error(" << i << ") " << error(i) << std::endl;
		}
	}*/
	static bool equalsWithTolerance(float first, float second,
			float error = 10e-3) {
		if(std::abs(first-second) <= error) {
			return true;
		}
		else {
			return false;
		}
	}
	dataType computeKernel(unsigned int first,unsigned int second) {
		if((*cachedKernel)(first,second) == -1) {
			(*cachedKernel)(first,second) =
					kernel->kernelFunction(model->X(first),
							model->X(second),model->X.cols());
		}
		return (*cachedKernel)(first,second);
	}
	/**
	 * Computes kernel value for training examples a,b.
	 * @param a
	 * @param b
	 * @return kernel value
	 */
	static dataType computeKernel(AbstractKernel<dataType> * kernel,
			Matrix<dataType> * cachedKernel,Matrix<inputType> & X,
			unsigned int first,unsigned int second) {
		if((*cachedKernel)(first,second) == -1) {
			(*cachedKernel)(first,second) =
					kernel->kernelFunction(X(first),X(second),X.cols());
		}
		return (*cachedKernel)(first,second);
	}
};


#endif /* MODIFIEDPARALLELSMOCLASSIFIER_H_ */
