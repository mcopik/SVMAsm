/*
 * AllTests.cpp
 *
 *  Created on: Nov 12, 2012
 *      Author: mcopik
 */


#include "gtest/gtest.h"


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
