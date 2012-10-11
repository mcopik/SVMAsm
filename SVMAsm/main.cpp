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

int main()
{
	void* handle = dlopen("./libasm.so", RTLD_LAZY);
	typedef unsigned int (*func)(char*);

    // reset errors
    dlerror();
    func hello = (func) dlsym(handle, "foo");
    const char *dlsym_error = dlerror();
    if (dlsym_error) {
        std::cout << "Cannot load symbol 'hello': " << dlsym_error <<
            '\n';
        dlclose(handle);
        return 1;
    }
    char c = 0;
    unsigned int a = (*hello)(&c);
    //should print "5 A" to output
    std::cout << "LibAsm: " << a << " " << c << std::endl;
    handle = dlopen("./libcpp.so", RTLD_LAZY);
    hello = (func) dlsym(handle, "foo");
    dlsym_error = dlerror();
    if (dlsym_error) {
        std::cout << "Cannot load symbol 'hello': " << dlsym_error <<
            '\n';
        dlclose(handle);
        return 1;
    }
    std::cout << "LibCpp: " << (*hello)(&c) << " " << c << std::endl;
    dlclose(handle);
    return 0;
}

