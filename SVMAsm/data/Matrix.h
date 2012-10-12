/*
 * Matrix.h
 *
 *  Created on: Oct 12, 2012
 *      Author: mcopik
 */

#ifndef MATRIX_H_
#define MATRIX_H_

struct Matrix {
public:
	/**
	 * data[i][j]
	 * i - row
	 * j - column
	 */
	bool ** data;
	int rows = 0;
	int columns = 0;
};


#endif /* MATRIX_H_ */
