/*
 * functions.h
 *
 *  Created on: Nov 12, 2012
 *      Author: mcopik
 */

#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_
#include <fstream>
#include "../src/data/Matrix.h"
#include "../src/data/Vector.h"
#include "gtest/gtest.h"

template<class T>
Matrix<T> loadMatrix(const char * path) {
	std::ifstream file(path);
	assert(file.is_open() == true);
	unsigned int rows,cols;
	file >> rows >> cols;
	Matrix<T> X(rows,cols);
	for(unsigned int i = 0;i < rows;++i) {
		for(unsigned int j = 0;j < cols;++j) {
			file >> X(i,j);
		}
	}
	file.close();
	return std::move(X);
}

template<class T>
Vector<T> loadVector(const char * path) {
	std::ifstream file(path);
	assert(file.is_open() == true);
	unsigned int size;
	file >> size;
	Vector<T> X(size);
	for(unsigned int i = 0;i < size;++i) {
		file >> X(i);
	}
	file.close();
	return std::move(X);
}

#endif /* FUNCTIONS_H_ */
