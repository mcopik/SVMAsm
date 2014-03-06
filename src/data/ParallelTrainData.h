#ifndef PARALLEL_TRAIN_DATA_H
#define PARALLEL_TRAIN_DATA_H

#include "../kernel/AbstractKernel.h"

template<class T,class U>
struct ParallelTrainData {
public:
	/**
	 * Size of arrays.
	 */
	unsigned int trainDataSize;
	/**
	 * Thread saves here iHigh value.
	 */
	unsigned int iHigh;
	/**
	 * Thread saves here iLow value.
	 */
	unsigned int iLow;
	/**
	 * Number of this thread.
	 * Determines part of data.
	 */
	unsigned int threadID;
	/**
	 * Cost value.
	 */
	U cost;
	/**
	 * Approximation error value.
	 */
	U error;
	/**
	 * Pointer to array with Y.
	 * Used in finding iHigh,iLow.
	 */
	T * yArray;
	/**
	 * Pointer to array with alpha values.
	 * Used in finding iHigh,iLow.
	 */
	T * alphaArray;
	/**
	 * Pointer to array with error function values.
	 * Used in updating error cache.
	 */
	U * errorArray;
	/**
	 * Pointer to kernel class.
	 */
	AbstractKernel<T> * kernel;
	/**
	 * Training matrix.
	 */
	Matrix<T> * X;
	/**
	 * Matrix with cached kernel.
	 */
	Matrix<T> * cachedKernel;
	/**
	 * Value of model->Y(iHigh)*model->alphas(iHigh) - highOld.
	 */
	U errorUpdateHigh;
	/**
	 * Value of model->Y(iLow)*model->alphas(iLow) - lowOld.
	 */
	U errorUpdateLow;
};



#endif
