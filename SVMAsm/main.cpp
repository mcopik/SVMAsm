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
#include "kernel/GaussianKernel.h"

extern void testSharedLibrary(const char*,int,const char*);
extern void testGaussianKernel();
int main()
{
    testSharedLibrary("./libasm.so",RTLD_LAZY,"foo");
    testSharedLibrary("./libcpp.so",RTLD_LAZY,"foo");
    testGaussianKernel();
    return 0;
}

