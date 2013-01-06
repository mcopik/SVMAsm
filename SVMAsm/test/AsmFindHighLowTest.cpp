/*
 * AsmTest.cpp
 *
 *  Created on: Dec 31, 2012
 *      Author: mcopik
 */

#include <fstream>
#include <pthread.h>
#include <dlfcn.h>
#include "gtest/gtest.h"
#include "../src/data/ParallelAsmData.h"

typedef void *(*threadFunction)(void*);
class AsmFindHighLowTest : public ::testing::Test {
protected:
	threadFunction findHighLow;
	void * asmHandle;
	ParallelAsmData * threadsAsmData;
	int numberOfThreads;
	pthread_t * threadID;
	int * ret;
	virtual void SetUp() {
		numberOfThreads = 2;
		asmHandle = dlopen("./../bin/libasm.so",RTLD_LAZY);
		assert(asmHandle != nullptr);
		dlerror();
		findHighLow  = (threadFunction)dlsym(asmHandle, "findHighLow");
		assert(findHighLow != nullptr);
		threadsAsmData = new ParallelAsmData[numberOfThreads];
		threadID = new pthread_t[numberOfThreads];
		ret = new int[numberOfThreads];
		for(int i = 0;i < numberOfThreads;++i) {
			threadsAsmData[i].cost = 0.1;
			threadsAsmData[i].error = 0.001;
			threadsAsmData[i].threadID = i;
		}
	}
	virtual void TearDown() {
		//dlclose(asmHandle);
		delete threadsAsmData;
		delete threadID;
		delete ret;
	}
};

TEST_F(AsmFindHighLowTest,testMiddleAlpha) {
	float yArray1[] = {1,1};
	float alphaArray1[] = {0.05,0.05,0.1,0.0};
	float errorArray1[] = {1.13,1.14,0.0,2.5};
	for(int i = 0;i < numberOfThreads;++i) {
		threadsAsmData[i].trainDataSize = 2;
		threadsAsmData[i].yArray = yArray1;
		threadsAsmData[i].alphaArray = alphaArray1;
		threadsAsmData[i].errorArray = errorArray1;
	}
	//(*findHighLow)((void*)&threadsAsmData[0]);
	std::cout << threadsAsmData[0].trainDataSize << std::endl;
	for(int i = 0;i < numberOfThreads;++i) {
		pthread_create(&(threadID[i]),NULL,(void* (*)(void*))findHighLow,(void*)(&threadsAsmData[i]));
		pthread_join(threadID[i],nullptr);
	}
		for(int i = 0;i < numberOfThreads;++i) {
		//ASSERT_EQ(threadsAsmData[0].iHigh,0);
		ASSERT_EQ(threadsAsmData[i].iLow,1);
	}
	float yArray2[] = {1,1,-1,-1};
	float alphaArray2[] = {0.0,0.067,0.098,0.0354};
	float errorArray2[] = {0.095,0.093,0.025,0.013};
	for(int i = 0;i < numberOfThreads;++i) {
		threadsAsmData[i].trainDataSize = 4;
		threadsAsmData[i].yArray = yArray2;
		threadsAsmData[i].alphaArray = alphaArray2;
		threadsAsmData[i].errorArray = errorArray2;
	}
	//(*findHighLow)((void*)&threadsAsmData[0]);
	std::cout << threadsAsmData[0].trainDataSize << std::endl;
	for(int i = 0;i < numberOfThreads;++i) {
		pthread_create(&(threadID[i]),NULL,(void* (*)(void*))findHighLow,(void*)(&threadsAsmData[i]));
		pthread_join(threadID[i],nullptr);
	}
	for(int i = 0;i < numberOfThreads;++i) {
		ASSERT_EQ(threadsAsmData[0].iHigh,3);
		ASSERT_EQ(threadsAsmData[i].iLow,1);
	}
	//lowest and highest error are excluded - alpha too close to 0/cost!
	float yArray3[] = {1,-1,-1,1,1,-1};
	float alphaArray3[] = {0.05,0.001,0.067,0.098,0.0354,0.9999};
	float errorArray3[] = {0.093,0.0,0.085,0.013,0.03,0.1};
	for(int i = 0;i < numberOfThreads;++i) {
		threadsAsmData[i].trainDataSize = 6;
		threadsAsmData[i].yArray = yArray3;
		threadsAsmData[i].alphaArray = alphaArray3;
		threadsAsmData[i].errorArray = errorArray3;
	}
	//(*findHighLow)((void*)&threadsAsmData[0]);
	for(int i = 0;i < numberOfThreads;++i) {
		pthread_create(&(threadID[i]),NULL,(void* (*)(void*))findHighLow,(void*)(&threadsAsmData[i]));
		pthread_join(threadID[i],nullptr);
	}
	for(int i = 0;i < numberOfThreads;++i) {
		ASSERT_EQ(threadsAsmData[0].iHigh,3);
		ASSERT_EQ(threadsAsmData[i].iLow,0);
	}
}


