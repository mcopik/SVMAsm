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
#include <stdexcept>
#include "kernel/LinearKernel.h"
#include "kernel/GaussianKernel.h"
#include "classifier/SMOClassifier.h"
#include "classifier/ModifiedSMOClassifier.h"
#include "classifier/ModifiedParallelSMOClassifier.h"
#include "data/Vector.h"
#include "data/Matrix.h"
#include "data/TrainData.h"
//extern void testSharedLibrary(const char*,int,const char*);
//extern void testSMO();
//extern void testThreads();


int main()
{
    //testSharedLibrary("./libasm.so",RTLD_LAZY,"foo");
    //testSharedLibrary("./libcpp.so",RTLD_LAZY,"foo");
    //testSMO();
    //testThreads();    std::ifstream file("../test/gaussianKernelTest3Y");
	LinearKernel<float> kernel;
	Vector<float> y = Vector<float>::loadVector("../test/testDataSpamBig/Y");
	Matrix<float> X = Matrix<float>::loadMatrix("../test/testDataSpamBig/X");
	Matrix<float> Xtest = Matrix<float>::loadMatrix("../test/testDataSpamBig/Xtest");
	Vector<float> Ytest = Vector<float>::loadVector("../test/testDataSpamBig/Ytest");
	TrainData<float> data(X,y);
	try {
		ModifiedParallelSMOClassifier<float,float> classifier;
		//kernel.setSigma(0.1);
		classifier.setC(0.1);
		classifier.setError(1e-3);
		classifier.setEpsilon(1e-3);
		classifier.trainAsm(data,kernel);
		std::cout << X.rows() << " " << X.cols() << std::endl;
		std::cout << y.size()  << std::endl;
		std::cout << "b: " << classifier.model->b << std::endl;
		Vector<float> predicts = classifier.predict(Xtest);
		int counter = 0;
		for(unsigned int i = 0;i < predicts.size();++i)
		if(predicts(i) == Ytest(i))
			++counter;
		std::cout << counter << std::endl;
		std::cout << ((float)counter)/predicts.size() << std::endl;
	}
	catch(std::runtime_error & e) {
		std::cout << e.what() << std::endl;
	}
	return 0;
}
