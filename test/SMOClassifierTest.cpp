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
#include "../src/kernel/LinearKernel.h"

class SMOClassifierTest : public ::testing::Test {
protected:
	SMOClassifier<float,float> classifier;
	GaussianKernel<float> gaussianKernel;
	LinearKernel<float> linearKernel;
	virtual void setUp() {
		classifier.setEpsilon(1e-3);
		classifier.setError(1e-3);
		classifier.setC(0.1);
		gaussianKernel.setSigma(0.1);
	}
};
template<class T>
void testClassifier(TrainData<T> & data,AbstractClassifier<T,T> & classifier,
		AbstractKernel<T> & kernel,Matrix<T> & Xtest,
		Vector<T> & Ytest,float correctness) {
	classifier.train(data,kernel);
	Vector<T> predicts = classifier.predict(Xtest);
	int counter = 0;
    for(unsigned int i = 0;i < predicts.size();++i)
    	if(predicts(i) == Ytest(i))
    		++counter;
    float result = ((float)counter)/predicts.size();
	EXPECT_GE(result,correctness);
}

TEST_F(SMOClassifierTest,test51x2) {
	Matrix<float> X = loadMatrix<float>("testData51x2/X");
	Vector<float> Y = loadVector<float>("testData51x2/Y");
	Matrix<float> Xtest = loadMatrix<float>("testData51x2/Xtest");
	Vector<float> Ytest = loadVector<float>("testData51x2/Ytest");
	gaussianKernel.setSigma(0.1);
	TrainData<float> data(X,Y);
	testClassifier(data,classifier,gaussianKernel,Xtest,Ytest,1.0);
}

TEST_F(SMOClassifierTest,testSpam500Gaussian) {

    Vector<float> y = loadVector<float>("../test/testDataSpam500/Y");
    Matrix<float> X = Matrix<float>::loadMatrix("../test/testDataSpam500/X");
    Matrix<float> Xtest = Matrix<float>::loadMatrix("../test/testDataSpam500/Xtest");
    Vector<float> Ytest = loadVector<float>("../test/testDataSpam500/Ytest");
	TrainData<float> data(X,y);
	gaussianKernel.setSigma(0.1);
	testClassifier(data,classifier,gaussianKernel,Xtest,Ytest,0.7);
}

TEST_F(SMOClassifierTest,testSpam500Linear) {

    Vector<float> y = loadVector<float>("../test/testDataSpam500/Y");
    Matrix<float> X = Matrix<float>::loadMatrix("../test/testDataSpam500/X");
    Matrix<float> Xtest = Matrix<float>::loadMatrix("../test/testDataSpam500/Xtest");
    Vector<float> Ytest = loadVector<float>("../test/testDataSpam500/Ytest");
	TrainData<float> data(X,y);
	testClassifier(data,classifier,linearKernel,Xtest,Ytest,0.95);
}
