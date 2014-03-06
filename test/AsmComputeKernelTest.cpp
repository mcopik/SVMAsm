/*
 * AsmComputeKernelTest.cpp
 *
 *  Created on: Jan 07, 2012
 *      Author: mcopik
 */

#include <fstream>
#include <pthread.h>
#include <dlfcn.h>
#include <random>
#include "gtest/gtest.h"
#include "../src/data/Matrix.h"

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
	float firstArray3[] = {1,2,3,4,5,6,7,8,9,10,11,12};
	float secondArray3[] = {1,-1,1,-1,1,-1,1,-1,1,-1,1,-1};
	for(int i = 0;i < numberOfThreads;++i) {
		data[i].size = 12;
		data[i].first = firstArray3;
		data[i].second = secondArray3;
	}
	for(int i = 0;i < numberOfThreads;++i) {
		pthread_create(&(threadID[i]),NULL,(void* (*)(void*))computeKernel,(void*)(&data[i]));
		pthread_join(threadID[i],nullptr);
	}
	for(int i = 0;i < numberOfThreads;++i) {
		ASSERT_EQ(data[i].result,-6);
	}
	float firstArray4[] = {1,2,3,4,5,6};
	float secondArray4[] = {1,-1,1,-1,1,-1};
	for(int i = 0;i < numberOfThreads;++i) {
		data[i].size = 6;
		data[i].first = firstArray4;
		data[i].second = secondArray4;
	}
	for(int i = 0;i < numberOfThreads;++i) {
		pthread_create(&(threadID[i]),NULL,(void* (*)(void*))computeKernel,(void*)(&data[i]));
		pthread_join(threadID[i],nullptr);
	}
	for(int i = 0;i < numberOfThreads;++i) {
		ASSERT_EQ(data[i].result,-3);
	}
	float firstArray5[] = {1,2,3,4,5,6,7,8,9,10,11,12,100};
	float secondArray5[] = {1,-1,1,-1,1,-1,1,-1,1,-1,1,-1,1};
	for(int i = 0;i < numberOfThreads;++i) {
		data[i].size = 13;
		data[i].first = firstArray5;
		data[i].second = secondArray5;
	}
	for(int i = 0;i < numberOfThreads;++i) {
		pthread_create(&(threadID[i]),NULL,(void* (*)(void*))computeKernel,(void*)(&data[i]));
		pthread_join(threadID[i],nullptr);
	}
	for(int i = 0;i < numberOfThreads;++i) {
		ASSERT_EQ(data[i].result,94);
	}
	
	float * firstArray6 = nullptr;
	float * secondArray6 = nullptr;
	float result = 0;
	int size;
	std::mt19937	rng;
	uint32_t	seed;
	rng.seed(seed);
	std::uniform_int_distribution<uint32_t> uint_dist(1,1000);	
	for(int j = 0;j < 5;++j) {
		size = 960;//uint_dist(rng);
		result = 0;
		firstArray6 = new float[size];
		secondArray6 = new float[size];
		for(int i = 0;i < size;++i) {
			firstArray6[i] = uint_dist(rng) % 2 == 0 ? 1 : 0;
			secondArray6[i] = uint_dist(rng) %  2== 0 ? 1 : 0;
			result += firstArray6[i]*secondArray6[i];
		}
		for(int i = 0;i < numberOfThreads;++i) {
			data[i].size = size;
			data[i].first = firstArray6;
			data[i].second = secondArray6;
		}
		for(int i = 0;i < numberOfThreads;++i) {
			pthread_create(&(threadID[i]),NULL,(void* (*)(void*))computeKernel,(void*)(&data[i]));
			pthread_join(threadID[i],nullptr);
		}
		for(int i = 0;i < numberOfThreads;++i) {
			EXPECT_FLOAT_EQ(data[i].result,result);
		}
		delete	firstArray6;
		delete	secondArray6;
	}
	
	Matrix<float> X = Matrix<float>::loadMatrix("testDataSpam500/X");
	for(int i = 0;i < numberOfThreads;++i) {
		data[i].size = X.cols();
		data[i].first = X(0);
		data[i].second = X(455);
	}
	for(int i = 0;i < numberOfThreads;++i) {
		pthread_create(&(threadID[i]),NULL,(void* (*)(void*))computeKernel,(void*)(&data[i]));
		pthread_join(threadID[i],nullptr);
	}
	for(int i = 0;i < numberOfThreads;++i) {
		ASSERT_EQ(data[i].result,9);
	}
}