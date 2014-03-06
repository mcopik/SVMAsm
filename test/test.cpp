/*
 * test.cpp
 *
 *  Created on: Oct 12, 2012
 *      Author: mcopik
 */
#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include <dlfcn.h>
#include <pthread.h>
#include "../kernel/GaussianKernel.h"
#include "../kernel/AbstractKernel.h"
#include "../classifier/SMOClassifier.h"

pthread_t threadID[2];

typedef unsigned int (*function)(char*);
#define ASSERT(left,operator,right) \
{ if(!((left) operator (right))) { 	\
	std::cerr << "ASSERT FAILED: " << \
	#left << #operator << #right << 	\
	" @ " << __FILE__ << " (" << __LINE__	\
	<< "). " << #left << "=" << (left) <<	\
	"; " << #right << "=" << (right) <<		\
	std::endl; } }

template<class T>
Matrix<T> loadMatrix(std::ifstream & file) {
	unsigned int rows,cols;
	file >> rows >> cols;
	Matrix<T> X(rows,cols);
	for(unsigned int i = 0;i < rows;++i) {
		for(unsigned int j = 0;j < cols;++j) {
			file >> X(i,j);
		}
	}
	return std::move(X);
}

template<class T>
Vector<T> loadVector(std::ifstream & file) {
	unsigned int size;
	file >> size;
	Vector<T> X(size);
	for(unsigned int i = 0;i < size;++i) {
		file >> X(i);
	}
	return std::move(X);
}
template<class T>
void testGaussianCached(const char * filenameX, const char * filenameResult,double maxError) {
    GaussianKernel<T> kernel;

	std::ifstream file(filenameX);
	ASSERT(file.is_open(),==,true);
    Matrix<T> X = loadMatrix<T>(file);
    file.close();
    file.open(filenameResult);
    Matrix<T> result = loadMatrix<T>(file);
    file.close();

    double sum = 0.0;
    double temp = 0.0;
    Matrix<T> newResult = kernel.cacheKernel(X,0.1);
    ASSERT(result.rows,==,newResult.rows);
    ASSERT(result.cols,==,newResult.cols);
    for(unsigned int i = 0;i < newResult.rows;++i) {
    	for(unsigned int j = 0;j < newResult.cols;++j) {
    		temp = std::abs(newResult(i,j)-result(i,j));
    		ASSERT(temp,<,maxError);
    		sum += temp;
    	}
    }
    std::cout << "Mean error: " << sum/(newResult.rows*newResult.cols) << std::endl;
}

void testSharedLibrary(const char * path,int type,const char* func) {
	void* handle = dlopen(path, type);
	ASSERT(handle, !=, 0);

    // reset errors
    dlerror();
    function hello = (function) dlsym(handle, func);
    assert(hello != nullptr);
    char c = 0;
    int val = (*hello)(&c);
    ASSERT(val,==,5);
    ASSERT(c,==,'A');
    dlclose(handle);
}

void testSMO() {
    GaussianKernel<float> kernel;
    float a[] = {1, 2, 1};
    float b[] = {0, 4, -1};
    ASSERT(std::abs(kernel.kernelFunction(a,b,3,2) - 0.32465),<,0.00001);
/*
    testGaussianCached<float>("../test/gaussianKernelTest2X",
    		"../test/gaussianKernelTestResult",0.002);
    testGaussianCached<double>("../test/gaussianKernelTestX",
    		"../test/gaussianKernelTestResult",10e-10);
*/

    file.open("../test/data");
    Matrix<float> Z = loadMatrix<float>(file);
    file.close();
    file.open("../test/gaussianKernelTest3Xtest");
    Matrix<float> Xtest = loadMatrix<float>(file);
    file.close();
    file.open("../test/gaussianKernelTest3Ytest");
    Vector<float> Ytest = loadVector<float>(file);
    file.close();
    file.open("../test/gaussianKernelTest3TestKernel");
    Matrix<float> kernelTest = loadMatrix<float>(file);
    file.close();
    std::cout << "X " << Xtest.rows << " " << Xtest.cols;
    std::cout << "Y " << Ytest.size;
    std::cout << "Kernel " << kernelTest.rows << " " << kernelTest.cols;
    TrainData<float> data(X,y);
    SMOClassifier<float,float> classifier;
    classifier.setCachedKernel(Z);
    classifier.train(data,kernel,0.1);
    TrainData<float> test(Xtest,Ytest);
    Vector<float> predicts = classifier.predict(test,kernelTest);
    int counter = 0;
    for(int i = 0;i < test.Y.size;++i) {
    	if(predicts(i) == Ytest(i))
    		counter++;
    }
    std::cout << ((float)counter)/predicts.size << std::endl;
    std::cout << "final b: " << classifier.model->b << std::endl;
}

void * thread(void * arg) {
	pthread_t id = pthread_self();

	for(unsigned long int i = 0;i < (0xFFFFFFFF);++i);
	if(pthread_equal(id,threadID[0])) {
		*((int*)arg) = 1;
	}
	else {
		*((int*)arg) = 2;
	}
	pthread_exit(0);
}

void testThreads() {
	int results[] = {0,0};
	void * ret[2];
	for(int i = 0;i < 2;++i) {
		pthread_create(&(threadID[i]),NULL,&thread,&results[i]);
	}
	pthread_join(threadID[0],(void**)&ret[0]);
	pthread_join(threadID[1],(void**)&ret[1]);
	for(int i = 0;i < 2;++i) {
		std::cout << "Thread " << i << " result " <<
				results[i] <<  std::endl;
	}
}

#undef ASSERT
