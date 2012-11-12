/*
 * functions.h
 *
 *  Created on: Nov 12, 2012
 *      Author: mcopik
 */

#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

template<class T>
Matrix<T> loadMatrix(std::ifstream & file) {
	unsigned int rows,cols;
	file >> rows >> cols;
	Matrix<T> X(rows,cols);
	for(unsigned int i = 0;i < rows;++i) {
		for(unsigned int j = 0;j < cols;++j) {
			file >> X(i,j);
		}
	}
	return std::move(X);
}

template<class T>
Vector<T> loadVector(std::ifstream & file) {
	unsigned int size;
	file >> size;
	Vector<T> X(size);
	for(unsigned int i = 0;i < size;++i) {
		file >> X(i);
	}
	return std::move(X);
}

#endif /* FUNCTIONS_H_ */
