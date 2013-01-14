/*
 * AsmComputeErrorTest.cpp
 *
 *  Created on: Jan 07, 2012
 *      Author: mcopik
 */

#include <fstream>
#include <pthread.h>
#include <dlfcn.h>
#include <random>
#include "gtest/gtest.h"
#include "../src/data/ParallelAsmData.h"


typedef void *(*threadFunction)(void*);
class AsmComputeErrorTest : public ::testing::Test {
protected:
	threadFunction computeError;
	ParallelAsmData * threadsAsmData;
	void * asmHandle;
	int numberOfThreads;
	pthread_t * threadID;
	int * ret;
	virtual void SetUp() {
		numberOfThreads = 2;
		asmHandle = dlopen("./../bin/libasm.so",RTLD_LAZY);
		assert(asmHandle != nullptr);
		dlerror();
		computeError  = (threadFunction)dlsym(asmHandle, "updateErrorCache");
		assert(computeError != nullptr);
		threadID = new pthread_t[numberOfThreads];
		ret = new int[numberOfThreads];
		threadsAsmData = new ParallelAsmData[numberOfThreads];
	}
	virtual void TearDown() {
		dlclose(asmHandle);
		delete threadsAsmData;
		delete threadID;
		delete ret;
	}
};

TEST_F(AsmComputeErrorTest,testErrorSimple) {
	float firstArray1[] = {-1,-1,-1,-1};
	float X[] = { 1,1,1,1, 0,1,1,1, 1,0,0,1, 0,0,1,1};
	float cache[] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};		
	for(int i = 0;i < numberOfThreads;++i) {
			threadsAsmData[i].trainDataSize = 4;
			threadsAsmData[i].errorArray = firstArray1;
			threadsAsmData[i].cost = 0.1;
			threadsAsmData[i].error = 0.001;
			threadsAsmData[i].X = X;
			threadsAsmData[i].threadID = i;
			threadsAsmData[i].errorUpdateHigh = 0.5;
			threadsAsmData[i].errorUpdateLow = 0.4;
			threadsAsmData[i].iHigh = 2;
			threadsAsmData[i].iLow = 3;
			threadsAsmData[i].cachedKernelHigh = &cache[4*2];
			threadsAsmData[i].cachedKernelLow = &cache[4*3];
			threadsAsmData[i].numberOfFeatures = 4;
			threadsAsmData[i].numberOfTrainExamples = 4;
	}
	(*computeError)((void*)&threadsAsmData[0]);
	std::cout << std::endl;
}