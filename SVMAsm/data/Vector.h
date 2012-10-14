/*
 * Vector.h
 *
 *  Created on: Oct 12, 2012
 *      Author: mcopik
 */

#ifndef VECTOR_H_
#define VECTOR_H_
#include <algorithm>
#include <cassert>

template<class T>
struct Vector {
public:
	unsigned int size = 0;
	/**
	 * data[i]
	 * i-th element in
	 */
	T * data;
	/**
	 * Constructor.
	 */
	Vector(int _size):size(_size),data(new T[_size]) {
		//don't want to do anything!
	}
	/**
	 * Copy constructor
	 */
	Vector(const Vector & m):size(m.size),data(new T[m.size]) {
		for(int i = 0;i < size;++i)
			data[i] = m.data[i];
	}
	/**
	 * Move constructor
	 */
	Vector(Vector && m):size(m.size),data(m.data) {
		m.data = nullptr;
	}
	/**
	 * Assignment move operator.
	 */
	Vector & operator=(Vector && m){
		std::swap(m.tab,this->data);
		m.data = nullptr;
		return *this;
	}
	/**
	 * Assignment operator.
	 */
	Vector & operator=(Vector & m) {
		delete data;
		data = new T[size];
		for(int i = 0;i < size;++i)
			data[i] = m.data[i];
		return *this;
	}
	/**
	 * Function operator.
	 * Gives access to elements in vector(R/W).
	 */
	T & operator() (unsigned row) {
		assert(row < size);
		return data[row];
	}
	/**
	 * Function operator.
	 * Gives access to elements in vector(read only).
	 */
	T operator() (unsigned row) const {
		assert(row < size);
		return data[row];
	}
	/**
	 * Destructor.
	 */
	~Vector() {
		delete[] data;
	}
};


#endif /* VECTOR_H_ */
