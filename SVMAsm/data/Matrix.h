/*
 * Matrix.h
 *
 *  Created on: Oct 12, 2012
 *      Author: mcopik
 */

#ifndef MATRIX_H_
#define MATRIX_H_
#include <algorithm>
#include <cassert>
#include <iostream>

template<class T>
struct Matrix {
public:
	/**
	 * Number of rows.
	 */
	unsigned int rows = 0;
	/**
	 * Number of columns.
	 */
	unsigned int cols = 0;
	/**
	 * data[i][j]
	 * i - row
	 * j - column
	 */
	T * data;
	/**
	 * Constructor.
	 */
	Matrix(int _rows,int _cols):rows(_rows),cols(_cols),data(new T[_rows*_cols]){

	}
	/**
	 * Copy constructor
	 */
	Matrix(const Matrix & m):rows(m.rows),cols(m.cols),data(new T[m.rows*m.rows]) {
		for(int i = 0;i < rows*cols;++i)
			data[i] = m.data[i];
	}
	/**
	 * Move constructor
	 */
	Matrix(Matrix && m):rows(m.rows),cols(m.cols),data(m.data) {
		m.data = nullptr;
	}
	/**
	 * Assignment move operator.
	 */
	Matrix & operator=(Matrix && m){
		std::swap(m.tab,this->data);
		m.data = nullptr;
		return *this;
	}
	/**
	 * Assignment operator.
	 */
	Matrix & operator=(Matrix & m) {
		delete data;
		rows = m.rows;
		cols = m.cols;
		data = new T[rows*cols];
		for(int i = 0;i < rows*cols;++i)
			data[i] = m.data[i];
		return *this;
	}
	/**
	 * Function operator.
	 * Gives access to elements in matrix(R/W).
	 */
	T & operator() (unsigned row,unsigned col) {
		assert(row < rows && col < cols);
		return data[row*cols+col];
	}
	/**
	 * Function operator.
	 * Gives access to elements in matrix(read only).
	 */
	T operator() (unsigned row,unsigned col) const{
		assert(row < rows && col < cols);
		return data[row*cols+col];
	}
	/**
	 * Destructor.
	 */
	~Matrix() {
		//std::cout << "Delete " << data << std::endl;
		delete[] data;
	}
	/**
	 * Overloaded multiply operator.
	 * Multiply by another matrix.
	 */
	Matrix<T> operator *(Matrix<T> & arg) {
		Matrix<T> retval(rows,arg.cols);
		for(unsigned int i = 0;i < rows;++i) {
			for(unsigned int j = 0;j < arg.cols;++j) {
				retval(i,j) = 0;
				for(unsigned int k = 0;k < cols;++k) {
						retval(i,j) += this->operator()(i,k)*arg(k,j);
				}
			}
		}
		return std::move(retval);
	}
	/**
	 * Multiply matrix by its own transposed matrix.
	 */
	Matrix<T> multiplyByTranspose() {
		Matrix<T> retval(rows,rows);
		for(unsigned int i = 0;i < rows;++i) {
			for(unsigned int j = 0;j < rows;++j) {
				retval(i,j) = 0;
				for(unsigned int k = 0;k < cols;++k) {
						retval(i,j) += this->operator()(i,k)*this->operator()(j,k);
				}
			}
		}
		return std::move(retval);
	}
};


#endif /* MATRIX_H_ */
