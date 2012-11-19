/*
 * SMOClassifierTest.cpp
 *
 *  Created on: Nov 19, 2012
 *      Author: mcopik
 */
#include <fstream>
#include "gtest/gtest.h"
#include "functions.h"
#include "../src/classifier/SMOClassifier.h"
#include "../src/kernel/GaussianKernel.h"

class SMOClassifierTest : public ::testing::Test {
protected:
	SMOClassifier<float,float> classifier;
	GaussianKernel<float> kernel;
};
template<class T>
void testClassifier(TrainData<T> & data,AbstractClassifier<T,T> & classifier,
		AbstractKernel<T> & kernel,Matrix<T> & Xtest,
		Vector<T> & Ytest,float correctness) {
	classifier.train(data,kernel);
	Vector<T> predicts = classifier.predict(Xtest);
	//int counter = 0;
    //for(unsigned int i = 0;i < predicts.size;++i)
    	//if(predicts(i) == Ytest(i))
    		//++counter;
    //float result = ((float)counter)/predicts.size;
    float result = 0;
	EXPECT_GE(result,correctness);
}
TEST_F(SMOClassifierTest,test51x2) {
	Matrix<float> X = loadMatrix<float>("testData51x2/X");
	Vector<float> Y = loadVector<float>("testData51x2/Y");
	Matrix<float> Xtest = loadMatrix<float>("testData51x2/Xtest");
	Vector<float> Ytest = loadVector<float>("testData51x2/Ytest");
	kernel.setSigma(0.1);
	TrainData<float> data(X,Y);
	testClassifier(data,classifier,kernel,Xtest,Ytest,1.0);
}
