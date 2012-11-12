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
#include "../test/functions.h"
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
	std::ifstream file("trainDataY");
    Vector<float> y = loadVector<float>(file);
    file.close();
    file.open("trainDataX");
    Matrix<float> X = loadMatrix<float>(file);
    file.close();
    TrainData<float> data(X,y);
    SMOClassifier<float,float> classifier;
    classifier.train(data,kernel,0.1,true);
    return 0;
}

