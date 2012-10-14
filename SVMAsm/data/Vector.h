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

#define ASSERT(left,operator,right) \
{ if(!((left) operator (right))) { 	\
	std::cerr << "ASSERT FAILED: " << \
	#left << #operator << #right << 	\
	" @ " << __FILE__ << " (" << __LINE__	\
	<< "). " << #left << "=" << (left) <<	\
	"; " << #right << "=" << (right) <<		\
	std::endl; } }

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
		//std::cout << "Create v" << std::endl;
	}
	/**
	 * Copy constructor
	 */
	Vector(const Vector & m):size(m.size),data(new T[m.size]) {
		for(int i = 0;i < size;++i)
			data[i] = m.data[i];
		//std::cout << "Create v3" << std::endl;
	}
	/**
	 * Move constructor
	 */
	Vector(Vector && m):size(m.size),data(m.data) {
		m.data = nullptr;
		//std::cout << "move con" << std::endl;
	}
	/**
	 * Assignment move operator.
	 */
	Vector & operator=(Vector && m){
		std::swap(m.tab,this->data);
		m.data = nullptr;
		//std::cout << "move as" << std::endl;
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
		ASSERT(row,<,size);
		return data[row];
	}
	/**
	 * Function operator.
	 * Gives access to elements in vector(read only).
	 */
	T operator() (unsigned row) const {
		ASSERT(row,<,size);
		return data[row];
	}
	Vector<T> multiplyEachByEach(Vector<T> & arg) {
		Vector<T> retval(size);
		for(unsigned int i = 0;i < size;++i) {
			retval(i) = this->operator()(i)*arg(i);
		}
		return std::move(retval);
	}
	/**
	 * Destructor.
	 */
	~Vector() {
		delete[] data;
	}
};

#undef ASSERT
#endif /* VECTOR_H_ */
