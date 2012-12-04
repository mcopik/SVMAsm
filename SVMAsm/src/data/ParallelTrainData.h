#ifndef PARALLEL_TRAIN_DATA_H
#define PARALLEL_TRAIN_DATA_H

#include "../kernel/AbstractKernel.h"

template<class T,class U>
struct ParallelTrainData {
public:
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
	 * Size of arrays.
	 */
	unsigned int trainDataSize;
	/**
	 * Cost value.
	 */
	U cost;
	/**
	 * Approximation error value.
	 */
	U error;
	/**
	 * Thread saves here iHigh value.
	 */
	unsigned int iHigh;
	/**
	 * Thread saves here iLow value.
	 */
	unsigned int iLow;
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
	 * Array with X(iHigh) row.
	 */
	T * xHigh;
	/**
	 * Array with X(iLow) row.
	 */
	T * xLow;
	/**
	 * Number of features in model.
	 */
	int featuresNumber;
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
