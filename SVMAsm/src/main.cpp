/**
 * @file main.cpp
 * Main function.
 * @author Marcin Copik
 * Silesian University of Technology
 *
 * This code is distributed under the terms of the GNU Library
 * General Public License, either version 3 of the license or, at
 * your option, any later version.
 */
#include <iostream>
#include <dlfcn.h>
#include <fstream>
#include "kernel/LinearKernel.h"
#include "classifier/SMOClassifier.h"
#include "data/Vector.h"
#include "data/Matrix.h"
#include "data/TrainData.h"
//extern void testSharedLibrary(const char*,int,const char*);
//extern void testSMO();
//extern void testThreads();


template<class T>
Vector<T> loadVector(const char * path) {
	std::ifstream file(path);
	unsigned int size;
	file >> size;
	Vector<T> X(size);
	for(unsigned int i = 0;i < size;++i) {
		file >> X(i);
	}
	file.close();
	return std::move(X);
}
int main()
{
    //testSharedLibrary("./libasm.so",RTLD_LAZY,"foo");
    //testSharedLibrary("./libcpp.so",RTLD_LAZY,"foo");
    //testSMO();
    //testThreads();    std::ifstream file("../test/gaussianKernelTest3Y");
	LinearKernel<float> kernel;
    Vector<float> y = loadVector<float>("../test/testDataSpam500/Y");
    Matrix<float> X = loadMatrix<float>("../test/testDataSpam500/X");
    Matrix<float> Xtest = loadMatrix<float>("../test/testDataSpam500/Xtest");
    Vector<float> Ytest = loadVector<float>("../test/testDataSpam500/Ytest");
    TrainData<float> data(X,y);
    SMOClassifier<float,float> classifier;
    classifier.train(data,kernel,true);
    std::cout << "b: " << classifier.model->b << std::endl;
	Vector<float> predicts = classifier.predict(Xtest);
	int counter = 0;
    for(unsigned int i = 0;i < predicts.size;++i)
    	if(predicts(i) == Ytest(i))
    		++counter;
    std::cout << counter << std::endl;
    std::cout << ((float)counter)/predicts.size << std::endl;
    return 0;
}

