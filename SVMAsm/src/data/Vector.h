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
	unsigned int m_size = 0;
	/**
	 * data[i]
	 * i-th element in
	 */
	T * data;
public:
	/**
	 * Constructor.
	 */
	Vector(int _size):m_size(_size),data(new T[_size]()) {
		//don't want to do anything!
		//std::cout << "Create v" << std::endl;
	}
	/**
	 * Copy constructor
	 */
	Vector(const Vector & m):m_size(m.m_size),data(new T[m.size]) {
		for(int i = 0;i < m_size;++i)
			data[i] = m.data[i];
		//std::cout << "Create v3" << std::endl;
	}
	/**
	 * Move constructor
	 */
	Vector(Vector && m):m_size(m.m_size),data(m.data) {
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
		data = new T[m_size];
		for(int i = 0;i < m_size;++i)
			data[i] = m.data[i];
		return *this;
	}
	/**
	 * Function operator.
	 * Gives access to elements in vector(R/W).
	 */
	T & operator() (unsigned row) {
		ASSERT(row,<,m_size);
		if(row >= m_size)
			throw std::exception();
		return data[row];
	}
	/**
	 * Function operator.
	 * Gives access to elements in vector(read only).
	 */
	T operator() (unsigned row) const {
		ASSERT(row,<,m_size);
		if(row >= m_size)
			throw std::exception();
		return data[row];
	}
	unsigned int size() {
		return m_size;
	}
	T * vectorData() {
		return data;
	}
	Vector<T> multiplyEachByEach(Vector<T> & arg) {
		Vector<T> retval(m_size);
		for(unsigned int i = 0;i < m_size;++i) {
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
	static Vector<T> loadVector(const char * path) {
		std::ifstream file(path);
		unsigned int size;
		file >> size;
		Vector<T> X(size);
		for(unsigned int i = 0;i < size;++i) {
			file >> X(i);
		}
		file.close();
		return std::move(X);
	}
};

#undef ASSERT
#endif /* VECTOR_H_ */
