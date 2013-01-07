/*
 * AsmFindHighLowTest.cpp
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
		dlclose(asmHandle);
		delete threadsAsmData;
		delete threadID;
		delete ret;
	}
};

TEST_F(AsmFindHighLowTest,testMiddleAlpha) {
	float yArray1[] = {1,1,1,1};
	float alphaArray1[] = {0.05,0.05,0.098,0.012};
	float errorArray1[] = {1.1,1.14,1.12,1.13};
	for(int i = 0;i < numberOfThreads;++i) {
		threadsAsmData[i].trainDataSize = 4;
		threadsAsmData[i].yArray = yArray1;
		threadsAsmData[i].alphaArray = alphaArray1;
		threadsAsmData[i].errorArray = errorArray1;
	}
	for(int i = 0;i < numberOfThreads;++i) {
		pthread_create(&(threadID[i]),NULL,(void* (*)(void*))findHighLow,(void*)(&threadsAsmData[i]));
		pthread_join(threadID[i],nullptr);
	}
	for(int i = 0;i < numberOfThreads;++i) {
		ASSERT_EQ(threadsAsmData[i].iHigh,0);
		ASSERT_EQ(threadsAsmData[i].iLow,1);
	}
	float yArray2[] = {1,1,-1,-1};
	float alphaArray2[] = {0.01,0.067,0.098,0.0354};
	float errorArray2[] = {0.095,0.093,0.025,0.013};
	for(int i = 0;i < numberOfThreads;++i) {
		threadsAsmData[i].trainDataSize = 4;
		threadsAsmData[i].yArray = yArray2;
		threadsAsmData[i].alphaArray = alphaArray2;
		threadsAsmData[i].errorArray = errorArray2;
	}
	for(int i = 0;i < numberOfThreads;++i) {
		pthread_create(&(threadID[i]),NULL,(void* (*)(void*))findHighLow,(void*)(&threadsAsmData[i]));
		pthread_join(threadID[i],nullptr);
	}
	for(int i = 0;i < numberOfThreads;++i) {
		ASSERT_EQ(threadsAsmData[i].iHigh,3);
		ASSERT_EQ(threadsAsmData[i].iLow,0);
	}
	float yArray3[] = {1,-1,-1,1,1,-1};
	float alphaArray3[] = {0.05,0.0011,0.067,0.098,0.0354,0.0989999};
	float errorArray3[] = {0.093,0.014,0.085,0.013,0.03,0.094};
	for(int i = 0;i < numberOfThreads;++i) {
		threadsAsmData[i].trainDataSize = 6;
		threadsAsmData[i].yArray = yArray3;
		threadsAsmData[i].alphaArray = alphaArray3;
		threadsAsmData[i].errorArray = errorArray3;
	}
	for(int i = 0;i < numberOfThreads;++i) {
		pthread_create(&(threadID[i]),NULL,(void* (*)(void*))findHighLow,(void*)(&threadsAsmData[i]));
		pthread_join(threadID[i],nullptr);
	}
	for(int i = 0;i < numberOfThreads;++i) {
		ASSERT_EQ(threadsAsmData[i].iHigh,3);
		ASSERT_EQ(threadsAsmData[i].iLow,5);
	}
}


TEST_F(AsmFindHighLowTest,testAlpha0) {
	//check iHigh
	float yArray1[] = {1,1,1,1};
	float alphaArray1[] = {0.0,0.000045,0.0,0.00003};
	float errorArray1[] = {1.13,1.14,0.0,2.5};
	for(int i = 0;i < numberOfThreads;++i) {
		threadsAsmData[i].trainDataSize = 4;
		threadsAsmData[i].yArray = yArray1;
		threadsAsmData[i].alphaArray = alphaArray1;
		threadsAsmData[i].errorArray = errorArray1;
	}
	for(int i = 0;i < numberOfThreads;++i) {
		pthread_create(&(threadID[i]),NULL,(void* (*)(void*))findHighLow,(void*)(&threadsAsmData[i]));
		pthread_join(threadID[i],nullptr);
	}
	for(int i = 0;i < numberOfThreads;++i) {
		ASSERT_EQ(threadsAsmData[i].iHigh,2);
		ASSERT_EQ(threadsAsmData[i].iLow,-1);
	}
	//check iLow
	float yArray2[] = {-1,-1,-1,-1};
	float alphaArray2[] = {5e-6,0.0,0.001,0.000099};
	float errorArray2[] = {0.0003,0.093,0.095,0.025};
	for(int i = 0;i < numberOfThreads;++i) {
		threadsAsmData[i].trainDataSize = 4;
		threadsAsmData[i].yArray = yArray2;
		threadsAsmData[i].alphaArray = alphaArray2;
		threadsAsmData[i].errorArray = errorArray2;
	}
	for(int i = 0;i < numberOfThreads;++i) {
		pthread_create(&(threadID[i]),NULL,(void* (*)(void*))findHighLow,(void*)(&threadsAsmData[i]));
		pthread_join(threadID[i],nullptr);
	}
	for(int i = 0;i < numberOfThreads;++i) {
		ASSERT_EQ(threadsAsmData[i].iHigh,-1);
		ASSERT_EQ(threadsAsmData[i].iLow,2);
	}
	//check iLow and iHigh
	float yArray3[] = {1,-1,1,-1,1,-1};
	float alphaArray3[] = {0.0,0.00999,0.000031312,1e-7,0.00098765,0.00032};
	float errorArray3[] = {2.73,1.15,0.0001,2.5,0.82,0.0};
	for(int i = 0;i < numberOfThreads;++i) {
		threadsAsmData[i].trainDataSize = 6;
		threadsAsmData[i].yArray = yArray3;
		threadsAsmData[i].alphaArray = alphaArray3;
		threadsAsmData[i].errorArray = errorArray3;
	}
	for(int i = 0;i < numberOfThreads;++i) {
		pthread_create(&(threadID[i]),NULL,(void* (*)(void*))findHighLow,(void*)(&threadsAsmData[i]));
		pthread_join(threadID[i],nullptr);
	}
	for(int i = 0;i < numberOfThreads;++i) {
		ASSERT_EQ(threadsAsmData[i].iHigh,2);
		ASSERT_EQ(threadsAsmData[i].iLow,3);
	}
}

TEST_F(AsmFindHighLowTest,testAlphaCost) {
	//check iHigh
	float yArray1[] = {-1,-1,-1,-1};
	float alphaArray1[] = {0.1,0.0991234,0.0999999,0.1};
	float errorArray1[] = {0.0,1.14,0.5,2.5};
	for(int i = 0;i < numberOfThreads;++i) {
		threadsAsmData[i].trainDataSize = 4;
		threadsAsmData[i].yArray = yArray1;
		threadsAsmData[i].alphaArray = alphaArray1;
		threadsAsmData[i].errorArray = errorArray1;
	}
	for(int i = 0;i < numberOfThreads;++i) {
		pthread_create(&(threadID[i]),NULL,(void* (*)(void*))findHighLow,(void*)(&threadsAsmData[i]));
		pthread_join(threadID[i],nullptr);
	}
	for(int i = 0;i < numberOfThreads;++i) {
		//ASSERT_EQ(threadsAsmData[i].iHigh,0);
		//ASSERT_EQ(threadsAsmData[i].iLow,-1);
	}
	//check iLow
	float yArray2[] = {1,1,1,1};
	float alphaArray2[] = {0.1,0.1,0.099,0.0993};
	float errorArray2[] = {1.23,0.093,2.095,0.025};
	for(int i = 0;i < numberOfThreads;++i) {
		threadsAsmData[i].trainDataSize = 4;
		threadsAsmData[i].yArray = yArray2;
		threadsAsmData[i].alphaArray = alphaArray2;
		threadsAsmData[i].errorArray = errorArray2;
	}
	for(int i = 0;i < numberOfThreads;++i) {
		pthread_create(&(threadID[i]),NULL,(void* (*)(void*))findHighLow,(void*)(&threadsAsmData[i]));
		pthread_join(threadID[i],nullptr);
	}
	for(int i = 0;i < numberOfThreads;++i) {
		ASSERT_EQ(threadsAsmData[i].iHigh,-1);
		ASSERT_EQ(threadsAsmData[i].iLow,2);
	}
	//check iLow and iHigh
	float yArray3[] = {1,-1,1,-1,1,-1};
	float alphaArray3[] = {0.1,0.0991,0.1,0.1,0.099999999,0.099331221};
	float errorArray3[] = {1.73,1.15,2.031,2.5,0.82,0.0};
	for(int i = 0;i < numberOfThreads;++i) {
		threadsAsmData[i].trainDataSize = 6;
		threadsAsmData[i].yArray = yArray3;
		threadsAsmData[i].alphaArray = alphaArray3;
		threadsAsmData[i].errorArray = errorArray3;
	}
	for(int i = 0;i < numberOfThreads;++i) {
		pthread_create(&(threadID[i]),NULL,(void* (*)(void*))findHighLow,(void*)(&threadsAsmData[i]));
		pthread_join(threadID[i],nullptr);
	}
	for(int i = 0;i < numberOfThreads;++i) {
		ASSERT_EQ(threadsAsmData[i].iHigh,5);
		ASSERT_EQ(threadsAsmData[i].iLow,2);
	}
}