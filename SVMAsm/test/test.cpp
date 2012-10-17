/*
 * test.cpp
 *
 *  Created on: Oct 12, 2012
 *      Author: mcopik
 */
#include <cassert>
#include <dlfcn.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include "../kernel/GaussianKernel.h"
#include "../kernel/AbstractKernel.h"
#include "../classifier/SMOClassifier.h"

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

void testGaussianKernel() {
    GaussianKernel<float> kernel;
    float a[] = {1, 2, 1};
    float b[] = {0, 4, -1};
    ASSERT(std::abs(kernel.kernelFunction(a,b,3,2) - 0.32465),<,0.00001);

    testGaussianCached<float>("../test/gaussianKernelTestX",
    		"../test/gaussianKernelTestResult",0.002);
    testGaussianCached<double>("../test/gaussianKernelTestX",
    		"../test/gaussianKernelTestResult",10e-10);
    std::ifstream file("../test/gaussianKernelTest3Y");
    Vector<float> y = loadVector<float>(file);
    file.close();
    file.open("../test/gaussianKernelTest3X");
    Matrix<float> X = loadMatrix<float>(file);
    file.close();
    Matrix<float> d = X.multiplyByTranspose();
    //TrainData<float> data(X,y);
    //SMOClassifier<float,float> classifier;
    //classifier.train(data,kernel,0.1);
    //std::cout << classifier.model->b << std::endl;
    //for(int i = 0;i < classifier.model->alphas.size;++i)
      //  std::cout << classifier.model->alphas(i) << std::endl;
}
#undef ASSERT
