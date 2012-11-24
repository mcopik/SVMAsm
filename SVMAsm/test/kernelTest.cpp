#include <fstream>
#include "gtest/gtest.h"
#include "kernel/GaussianKernel.h"
#include "data/Matrix.h"
#include "data/Vector.h"
#include "functions.h"

class KernelTest : public ::testing::Test {
protected:
	GaussianKernel<float> kernel;
	GaussianKernel<double> kernelDouble;
};
template<class T>
void testGaussianCached(GaussianKernel<T> kernel,const char * filenameX,
		const char * filenameResult,double maxError) {

    Matrix<T> X = loadMatrix<T>(filenameX);
    Matrix<T> result = loadMatrix<T>(filenameResult);

    double sum = 0.0;
    double temp = 0.0;
    Matrix<T> newResult = kernel.cacheKernel(X);
	EXPECT_EQ(result.rows,newResult.rows);
	EXPECT_EQ(result.cols,newResult.cols);
    for(unsigned int i = 0;i < newResult.rows;++i) {
    	for(unsigned int j = 0;j < newResult.cols;++j) {
    		temp = std::abs(newResult(i,j)-result(i,j));
			EXPECT_LT(temp,maxError);
    		sum += temp;
    	}
    }
}
TEST_F(KernelTest,kernelFunction) {
	kernel.setSigma(2);
	float a[] = {1, 2, 1};
	float b[] = {0, 4, -1};
	EXPECT_LT(std::abs(kernel.kernelFunction(a,b,3) - 0.32465),0.00001)
	<< "Gaussian Kernel differs on testdata";
}

TEST_F(KernelTest,KernelCached51x2) {
	kernel.setSigma(0.1);
    testGaussianCached(kernel,"testData51x2/X",
    		"testData51x2/gaussianKernel",0.002);
    kernelDouble.setSigma(0.1);
    testGaussianCached(kernelDouble,"testData51x2/X",
    		"testData51x2/gaussianKernel",10e-10);
}

TEST_F(KernelTest,KernelCached863x2) {
	kernel.setSigma(0.1);
    testGaussianCached(kernel,"testData863x2/X",
    		"testData863x2/gaussianKernel",0.002);
    kernelDouble.setSigma(0.1);
    testGaussianCached(kernelDouble,"testData863x2/X",
    		"testData863x2/gaussianKernel",10e-10);
}
