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
#include <stddef.h>
#include <time.h>
#include "AbstractClassifier.h"
#include "../data/TrainedModel.h"
#include "../data/TrainData.h"
#include "../data/ParallelTrainData.h"
#include "../data/ParallelAsmData.h"
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
	ParallelAsmData * threadsAsmData;
	Matrix<dataType> * cachedKernel = nullptr;
public:
	using AbstractClassifier<T,U>::model;
	void trainAsm(TrainData<inputType> & train,AbstractKernel<inputType> & _kernel) {
		std::cout << "TrainASM" << std::endl;
		model = new TrainedModel<inputType,dataType>
				(train.X,train.Y);
		cachedKernel = new Matrix<dataType>(train.X.rows(),train.X.rows());
		std::fill_n(cachedKernel->matrixData(),cachedKernel->rows()
				*cachedKernel->cols(),-1);
		Matrix<dataType> * cachedKernel2 = new Matrix<dataType>(train.X.rows(),train.X.rows());
		std::fill_n(cachedKernel2->matrixData(),cachedKernel2->rows()
				*cachedKernel2->cols(),-1);
		kernel = &_kernel;
		errorCache = new Vector<dataType>(train.X.rows());
		Vector<dataType> * errorCache2 = new Vector<dataType>(train.X.rows());
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
			(*errorCache2)(i) = -1*train.Y(i);
		}
		void * asmHandle = dlopen("./libasm.so",RTLD_LAZY);
		assert(asmHandle != nullptr);
		dlerror();
		threadFunction findHighLowAsm = (threadFunction) dlsym(asmHandle, "findHighLow");
		assert(findHighLowAsm != nullptr);
		threadFunction updateErrorCacheAsm = (threadFunction) dlsym(asmHandle, "updateErrorCache");
		assert(findHighLowAsm != nullptr);
		threadFunction updateAlphaAsm = (threadFunction) dlsym(asmHandle, "updateAlpha");
		assert(findHighLowAsm != nullptr);
		
		
		numberOfThreads = 2;
		int threadDataSize = model->X.rows()/numberOfThreads;
		int * iHighArray = new int[numberOfThreads];
		int * iLowArray = new int[numberOfThreads];
		pthread_t * threadID = new pthread_t[numberOfThreads];
		int * ret = new int[numberOfThreads];
		threadsData = new ParallelTrainData<inputType,dataType>[numberOfThreads];
		threadsAsmData = new ParallelAsmData[numberOfThreads];
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
		for(int i = 0;i < numberOfThreads;++i) {
			threadsAsmData[i].trainDataSize = threadDataSize;
			threadsAsmData[i].yArray = model->Y.vectorData(threadDataSize*i);
			threadsAsmData[i].alphaArray = model->alphas.vectorData(threadDataSize*i);
			threadsAsmData[i].errorArray = errorCache->vectorData(threadDataSize*i);
			threadsAsmData[i].cost = C;
			threadsAsmData[i].error = error;
			threadsAsmData[i].X = model->X(0);
			threadsAsmData[i].threadID = i;
			threadsAsmData[i].numberOfFeatures = model->X.cols();
			threadsAsmData[i].numberOfTrainExamples = model->X.rows();
			threadsAsmData[i].offset = i*threadDataSize*4;
		}
		//threadsAsmData[0].alphaArray[2] = 0.0;
		//threadsAsmData[1].alphaArray[2] = 0.1;
		std::cout << "ASM " << std::endl;
		//(*findHighLowAsm)(&threadsAsmData[0]);
		for(int i = 0;i < numberOfThreads;++i) {
			pthread_create(&(threadID[i]),NULL,(void* (*)(void*))findHighLowAsm,(void*)(&threadsAsmData[i]));//(void *(*)(void*))function,(void*)&threadsData[i]);
			pthread_join(threadID[i],nullptr);//(void**)&ret[i]);
		}
		/*threadsAsmData[0].trainDataSize = 2000;
		threadsAsmData[1].trainDataSize = 2000;
		threadsAsmData[0].alphaArray[2] = 0.005;
		threadsAsmData[1].alphaArray[2] = 0.098;
		threadsAsmData[0].cost = C;
		threadsAsmData[1].cost = C;
		threadsAsmData[0].error = error;
		threadsAsmData[1].error = error;
		std::cout << "ASM " << std::endl;
		for(int i = 0;i < numberOfThreads;++i) {
			pthread_create(&(threadID[i]),NULL,(void* (*)(void*))findHighLowAsm,(void*)(&threadsAsmData[i]));//(void *(*)(void*))function,(void*)&threadsData[i]);
			pthread_join(threadID[i],nullptr);//(void**)&ret[i]);
		}
		std::cout << threadsAsmData[0].trainDataSize << " " << threadsAsmData[1].trainDataSize << std::endl;
		std::cout << threadsAsmData[0].cost << " " << threadsAsmData[0].error << std::endl;
		std::cout << threadsAsmData[1].cost << " " << threadsAsmData[1].error << std::endl;
		//assert(threadsAsmData[0].trainDataSize == 1 && threadsAsmData[1].trainDataSize == 1);
		return;
		for(int i = 0;i < numberOfThreads;++i) {
			pthread_create(&(threadID[i]),NULL,(void* (*)(void*))findHighLow,(void*)(&threadsData[i]));//(void *(*)(void*))function,(void*)&threadsData[i]);
			pthread_join(threadID[i],nullptr);//(void**)&ret[i]);
		}*/
		iHigh = threadsAsmData[0].iHigh;
		iLow = threadsAsmData[0].iLow;
		for(int i = 1;i < numberOfThreads;++i) {
			if((*errorCache)(iHigh) > (*errorCache)(threadsAsmData[i].iHigh+threadDataSize*i)) {
				iHigh = threadsAsmData[i].iHigh+threadDataSize*i;
			}
			if((*errorCache)(iLow) < (*errorCache)(threadsAsmData[i].iLow+threadDataSize*i)) {
				iLow = threadsAsmData[i].iLow+threadDataSize*i;
			}
		}
		fLow = (*errorCache)(iLow);
		fHigh = (*errorCache)(iHigh);
		//getfHigh(fHigh,iHigh,threadsData,numberOfThreads);
		//getfLow(fLow,iLow,threadsData);
		alphaHighOld = model->alphas(iHigh);
		alphaLowOld = model->alphas(iLow);
		//updateAlpha(iHigh,iLow);
		//threadsAsmData[0].yArray = model->Y.vectorData(0);
		//threadsAsmData[0].alphaArray = model->alphas.vectorData(0);
		//threadsAsmData[0].errorArray = errorCache->vectorData(0);
		for(int i = 0;i < numberOfThreads; ++i) {
			threadsAsmData[i].iHigh = iHigh;
			threadsAsmData[i].iLow = iLow;
			threadsAsmData[i].cachedKernelHigh = (*cachedKernel)(iHigh);
			threadsAsmData[i].cachedKernelLow = (*cachedKernel)(iLow);
		}
		//updateAlpha(iHigh,iLow);
		(*updateAlphaAsm)((void*)&threadsAsmData[0]);
		//threadsAsmData[0].yArray = model->Y.vectorData(threadDataSize*i);
		//threadsAsmData[0].alphaArray = model->alphas.vectorData(threadDataSize*i);
		//threadsAsmData[0].errorArray = errorCache->vectorData(threadDataSize*i);
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
		clock_t startClock = clock();
		do {
			for(int i = 0;i < numberOfThreads;++i) {
				threadsData[i].errorUpdateHigh = model->Y(iHigh)*
						(model->alphas(iHigh) - alphaHighOld);
				threadsData[i].errorUpdateLow = model->Y(iLow)*
						(model->alphas(iLow) - alphaLowOld);
				threadsData[i].iHigh = iHigh;
				threadsData[i].iLow = iLow;
				threadsAsmData[i].errorUpdateHigh = model->Y(iHigh)*
						(model->alphas(iHigh) - alphaHighOld);
				threadsAsmData[i].errorUpdateLow = model->Y(iLow)*
						(model->alphas(iLow) - alphaLowOld);
				//threadsAsmData[i].cachedKernelHigh = (*cachedKernel)(iHigh);
				//threadsAsmData[i].cachedKernelLow = (*cachedKernel)(iLow);
			}
			//for(int i = 0;i < numberOfThreads;++i) {
			//	pthread_create(&(threadID[i]),NULL,(void* (*)(void*))updateErrorCache,(void*)(&threadsData[i]));//(void *(*)(void*))function,(void*)&threadsData[i]);
			//	pthread_join(threadID[i],nullptr);//(void**)&ret[i]);
			//}
			for(int i = 0;i < numberOfThreads;++i) {
				pthread_create(&(threadID[i]),NULL,(void* (*)(void*))updateErrorCacheAsm,(void*)(&threadsAsmData[i]));//(void *(*)(void*))function,(void*)&threadsData[i]);
				pthread_join(threadID[i],nullptr);//(void**)&ret[i]);
			}
			/*
			for(int i = 0;i < errorCache2->size();i++) {
				if(std::abs((*errorCache)(i) -(*errorCache2)(i)) >= 0.01) {
					std::cout << nr_iter << " " << i << " " << (*errorCache)(i) << " " << (*errorCache2)(i) << std::endl;
					std::cout << (*cachedKernel)(iHigh,i) << " " << (*cachedKernel2)(iHigh,i)  << std::endl;
					std::cout << iHigh << std::endl;
					throw std::exception();
				}
			}*/
			//std::cout << "SPEC " << (*errorCache)(0) << " " << (*errorCache)(1) << " " << (*errorCache)(2) << " " << (*errorCache)(3) <<std::endl;
			//std::cout << "SPEC " << (*errorCache)(4) << " " << (*errorCache)(5) << " " << (*errorCache)(6) << " " << (*errorCache)(7) <<std::endl;
			//std::cout << "SPEC " << (*errorCache)(2000) << " " << (*errorCache)(2001) << " " << (*errorCache)(2002) << " "<< (*errorCache)(2003)<< std::endl;
			//updateErrorCache(iHigh,iLow,alphaHighOld,alphaLowOld);
			//std::cout << "updatefLow " << (*errorCache)(iLow) << " updatefHigh " << (*errorCache)(iHigh) << std::endl;
			for(int i = 0;i < numberOfThreads;++i) {
				
				pthread_create(&(threadID[i]),NULL,(void* (*)(void*))findHighLowAsm,(void*)(&threadsAsmData[i]));//pthread_create(&(threadID[i]),NULL,(void* (*)(void*))findHighLow,(void*)(&threadsData[i]));//(void *(*)(void*))function,(void*)&threadsData[i]);
				pthread_join(threadID[i],nullptr);//(void**)&ret[i]);
			}
			
			iHigh = threadsAsmData[0].iHigh;
			iLow = threadsAsmData[0].iLow;
			for(int i = 1;i < numberOfThreads;++i) {
				if((*errorCache)(iHigh) > (*errorCache)(threadsAsmData[i].iHigh+threadDataSize*i)) {
					iHigh = threadsAsmData[i].iHigh+threadDataSize*i;
				}
				if((*errorCache)(iLow) < (*errorCache)(threadsAsmData[i].iLow+threadDataSize*i)) {
					iLow = threadsAsmData[i].iLow+threadDataSize*i;
				}
			}
			for(int i = 0;i < numberOfThreads; ++i) {
				threadsAsmData[i].iHigh = iHigh;
				threadsAsmData[i].iLow = iLow;
				threadsAsmData[i].cachedKernelHigh = (*cachedKernel)(iHigh);
				threadsAsmData[i].cachedKernelLow = (*cachedKernel)(iLow);
			}
			
			/*for(int i = 0;i < numberOfThreads;++i) {
				
				pthread_create(&(threadID[i]),NULL,(void* (*)(void*))findHighLow,(void*)(&threadsData[i]));//pthread_create(&(threadID[i]),NULL,(void* (*)(void*))findHighLow,(void*)(&threadsData[i]));//(void *(*)(void*))function,(void*)&threadsData[i]);
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
			}*/
			fLow = (*errorCache)(iLow);
			fHigh = (*errorCache)(iHigh);
			alphaHighOld = model->alphas(iHigh);
			alphaLowOld = model->alphas(iLow);
			(*updateAlphaAsm)((void*)&threadsAsmData[0]);
			/*
			updateAlpha(iHigh,iLow);
			float newAlphaHigh = model->alphas(iHigh);
			float newAlphaLow = model->alphas(iLow);
			model->alphas(iHigh) = alphaHighOld;
			model->alphas(iLow) = alphaLowOld;
			if(std::abs(newAlphaHigh - model->alphas(iHigh)) >= 0.01 || std::abs(newAlphaLow - model->alphas(iLow)) >= 0.01){
				std::cout << "nr " << nr_iter << " new low " << newAlphaLow << " new low2 " <<  model->alphas(iLow)
					<< " new high " << newAlphaHigh << " new high2 " <<  model->alphas(iHigh) << std::endl;
					std::cout << newAlphaLow << " " << model->alphas(iLow) << " " << std::endl;
					std::cout << std::abs(newAlphaHigh - model->alphas(iHigh)) << " " << std::abs(newAlphaLow != model->alphas(iLow)) << std::endl;
					
					model->alphas(iHigh) = alphaHighOld;
					model->alphas(iLow) = alphaLowOld;
					//throw std::exception();
					(*updateAlphaAsm)((void*)&threadsAsmData[0]);
					
					}*/
			
			
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
		std::cout << ((double)clock()-startClock)/CLOCKS_PER_SEC << std::endl;
		model->b = (fLow+fHigh)/2;

		delete iHighArray;
		delete iLowArray;
		delete threadID;
		delete ret;
	}
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
			//std::cout << i << " " << E << std::endl;
			E -= model->b;
			//std::cout << i << " " << E << std::endl;
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
		delete threadsAsmData;
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
		std::cout << "fHigh: " << (*errorCache)(iHigh) << " fLow: " << (*errorCache)(iLow) << std::endl;
		std::cout << "alphaHigh: " << model->alphas(iHigh) << " alphaLow: " << model->alphas(iLow) << std::endl;
		std::cout << "yHigh: " << model->Y(iHigh) << " yLow: " << model->Y(iLow) << std::endl;
		//compute boundaries for alphaLow
		dataType alphaLowLowerBound,alphaLowUpperBound;
		dataType eta = computeKernel(iHigh,iHigh) + computeKernel(iLow,iLow)
				- 2 * computeKernel(iHigh,iLow);
		dataType oldAlphaLow = model->alphas(iLow);
		std::cout << "eta " << eta <<  std::endl;
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
		std::cout << "low bound: " << alphaLowLowerBound << " upper bound: " << alphaLowUpperBound << std::endl;
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
					model->alphas(iLow) = alphaLowLowerBound;
				}
			}
			else {
				model->alphas(iLow) = oldAlphaLow;
			}
		}
		std::cout << "new alpha low " << model->alphas(iLow) << std::endl;
		//compute alphaHigh
		model->alphas(iHigh) += model->Y(iHigh)*model->Y(iLow)
				*(oldAlphaLow - model->alphas(iLow));
		std::cout << "new alpha high " << model->alphas(iHigh) << std::endl;
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
