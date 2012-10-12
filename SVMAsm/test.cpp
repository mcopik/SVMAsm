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
#include "kernel/GaussianKernel.h"
#include "kernel/AbstractKernel.h"
typedef unsigned int (*function)(char*);
#define ASSERT(left,operator,right) \
{ if(!((left) operator (right))) { 	\
	std::cerr << "ASSERT FAILED: " << \
	#left << #operator << #right << 	\
	" @ " << __FILE__ << " (" << __LINE__	\
	<< "). " << #left << "=" << (left) <<	\
	"; " << #right << "=" << (right) <<		\
	std::endl; } }

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
    GaussianKernel kernel;
    int x[] = {1, 2, 1};
    int y[] = {0, 4, -1};
    ASSERT(std::abs(kernel.kernelFunction(x,y,3,2) - 0.32465),<,0.00001);
}
