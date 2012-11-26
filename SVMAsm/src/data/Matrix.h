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

#define ASSERT(left,operator,right) \
{ if(!((left) operator (right))) { 	\
	std::cerr << "ASSERT FAILED: " << \
	#left << #operator << #right << 	\
	" @ " << __FILE__ << " (" << __LINE__	\
	<< "). " << #left << "=" << (left) <<	\
	"; " << #right << "=" << (right) <<		\
	std::endl; } }

template<class T>
struct Matrix {
	/**
	 * Number of rows.
	 */
	unsigned int m_rows = 0;
	/**
	 * Number of columns.
	 */
	unsigned int m_cols = 0;
	/**
	 * data[i][j]
	 * i - row
	 * j - column
	 */
	T * data = nullptr;
public:
	Matrix():data(nullptr),m_rows(0),m_cols(0){}
	/**
	 * Constructor.
	 */
	Matrix(int _rows,int _cols):m_rows(_rows),m_cols(_cols),data(new T[_rows*_cols]){

	}
	/**
	 * Copy constructor
	 */
	Matrix(const Matrix & m):m_rows(m.rows()),m_cols(m.cols()),data(new T[m.rows*m.cols]) {
		for(int i = 0;i < m_rows*m_cols;++i)
			data[i] = m.data[i];
	}
	/**
	 * Move constructor
	 */
	Matrix(Matrix && m):m_rows(m.rows()),m_cols(m.cols()),data(m.data) {
		m.data = nullptr;
	}
	/**
	 * Assignment move operator.
	 */
	Matrix & operator=(Matrix && m){
		std::swap(m.data,this->data);
		m.data = nullptr;
		return *this;
	}
	/**
	 * Assignment operator.
	 */
	Matrix & operator=(Matrix & m) {
		delete data;
		m_rows = m.rows();
		m_cols = m.cols();
		data = new T[m_rows*m_cols];
		for(int i = 0;i < m_rows*m_cols;++i)
			data[i] = m.data[i];
		return *this;
	}
	/**
	 * Function operator.
	 * Gives access to elements in matrix(R/W).
	 */
	T & operator() (unsigned row,unsigned col) {
		ASSERT(row,<,m_rows);
		ASSERT(col,<,m_cols);
		return data[row*m_cols+col];
	}
	/**
	 * Function operator.
	 * Gives access to elements in matrix(read only).
	 */
	T operator() (unsigned row,unsigned col) const{
		ASSERT(row,<,m_rows);
		ASSERT(col,<,m_cols);
		return data[row*m_cols+col];
	}
	/**
	 * Function operator.
	 * Gives access to row in matrix(read only).
	 */
	T * operator() (unsigned row) const {
		ASSERT(row,<,m_rows);
		return &data[row*m_cols];
	}
	/**
	 * Destructor.
	 */
	~Matrix() {
		//std::cout << "Delete " << data << std::endl;
		delete[] data;
	}
	unsigned int rows() {
		return this->m_rows;
	}
	unsigned int cols() {
		return this->m_cols;
	}
	T * matrixData() {
		return data;
	}
	/**
	 * Overloaded multiply operator.
	 * Multiply by another matrix.
	 */
	Matrix<T> operator *(Matrix<T> & arg) {
		Matrix<T> retval(m_rows,arg.cols());
		for(unsigned int i = 0;i < m_rows;++i) {
			for(unsigned int j = 0;j < arg.cols();++j) {
				retval(i,j) = 0;
				for(unsigned int k = 0;k < m_cols;++k) {
						retval(i,j) += this->operator()(i,k)*arg(k,j);
				}
			}
		}
		return std::move(retval);
	}
	/**
	 * Multiply matrix by its own transposed matrix.
	 */
	Matrix<T> multiplyByTranspose(Matrix<T> & arg) {
		///TODO:
		///optimize - output matrix is symetric!
		Matrix<T> retval(m_rows,m_rows);
		for(unsigned int i = 0;i < m_rows;++i) {
			for(unsigned int j = 0;j < m_rows;++j) {
				retval(i,j) = 0;
				for(unsigned int k = 0;k < m_cols;++k) {
						retval.data[i*m_rows+j] += data[i*m_cols+k]*arg.data[j*m_cols+k];//this->operator()(i,k)*this->operator()(j,k);
				}
			}
		}
		return std::move(retval);
	}
	Matrix<T> multiplyEachByEach(Matrix<T> & arg) {
		Matrix<T> retval(m_rows,m_rows);
		for(unsigned int i = 0;i < m_rows;++i) {
			for(unsigned int j = 0;j < m_rows;++j) {
				retval(i,j) = this->operator()(i,j)*arg(i,j);
			}
		}
		return std::move(retval);
	}
	static Matrix<T> loadMatrix(const char * path) {
		std::ifstream file(path);
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
};

#undef ASSERT
#endif /* MATRIX_H_ */
