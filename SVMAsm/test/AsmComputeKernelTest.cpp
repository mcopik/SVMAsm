/*
 * AsmComputeKernelTest.cpp
 *
 *  Created on: Jan 07, 2012
 *      Author: mcopik
 */

#include <fstream>
#include <pthread.h>
#include <dlfcn.h>
#include "gtest/gtest.h"


struct kernelData {
public:
	float * first;
	float * second;
	int size;
	float result;
};


typedef void *(*threadFunction)(void*);
class AsmComputeKernelTest : public ::testing::Test {
protected:
	threadFunction computeKernel;
	kernelData * data;
	void * asmHandle;
	int numberOfThreads;
	pthread_t * threadID;
	int * ret;
	virtual void SetUp() {
		numberOfThreads = 2;
		asmHandle = dlopen("./../bin/libasm.so",RTLD_LAZY);
		assert(asmHandle != nullptr);
		dlerror();
		computeKernel  = (threadFunction)dlsym(asmHandle, "computeLinearKernel");
		assert(computeKernel != nullptr);
		threadID = new pthread_t[numberOfThreads];
		ret = new int[numberOfThreads];
		data = new kernelData[numberOfThreads];
	}
	virtual void TearDown() {
		dlclose(asmHandle);
		delete data;
		delete threadID;
		delete ret;
	}
};

TEST_F(AsmComputeKernelTest,testKernelSimple) {
	float firstArray1[] = {1,1,1,1,1,1,1,1,1,1,1,1};
	float secondArray1[] = {1,1,1,1,1,1,1,1,1,1,1,1};
	for(int i = 0;i < numberOfThreads;++i) {
		data[i].size = 12;
		data[i].first = firstArray1;
		data[i].second = secondArray1;
	}
	for(int i = 0;i < numberOfThreads;++i) {
		pthread_create(&(threadID[i]),NULL,(void* (*)(void*))computeKernel,(void*)(&data[i]));
		pthread_join(threadID[i],nullptr);
	}
	for(int i = 0;i < numberOfThreads;++i) {
		ASSERT_EQ(data[i].result,12);
	}
	float firstArray2[] = {1,2,3,4,5,6,7,8,9,10,11,12};
	float secondArray2[] = {9,8,7,0,0,0,0,0,0,0,0,0};
	for(int i = 0;i < numberOfThreads;++i) {
		data[i].size = 12;
		data[i].first = firstArray2;
		data[i].second = secondArray2;
	}
	for(int i = 0;i < numberOfThreads;++i) {
		pthread_create(&(threadID[i]),NULL,(void* (*)(void*))computeKernel,(void*)(&data[i]));
		pthread_join(threadID[i],nullptr);
	}
	for(int i = 0;i < numberOfThreads;++i) {
		ASSERT_EQ(data[i].result,46);
	}
}