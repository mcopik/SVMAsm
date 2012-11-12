/*
 * AllTests.cpp
 *
 *  Created on: Nov 12, 2012
 *      Author: mcopik
 */


#include "gtest/gtest.h"
#include "kernel/GaussianKernel.h"
  class FooTest : public ::testing::Test {
  protected:
	    GaussianKernel<float> kernel;
  };
  TEST_F(FooTest,kernel) {
  float a[] = {1, 2, 1};
  float b[] = {0, 4, -1};
  EXPECT_LT(std::abs(kernel.kernelFunction(a,b,3,2) - 0.32465),0.00001)
  	<< "Gaussian Kernel differs on testdata";
  }

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
