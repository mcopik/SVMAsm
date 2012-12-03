#ifndef PARALLEL_TRAIN_DATA_H
#define PARALLEL_TRAIN_DATA_H

#include "../kernel/AbstractKernel.h"

template<class T>
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
	int trainDataSize;
	/**
	 * Cost value.
	 */
	int cost;
	/**
	 * Thread saves here iHigh value.
	 */
	int iHigh;
	/**
	 * Thread saves here iLow value.
	 */
	int iLow;
	/**
	 * Pointer to array with error function values.
	 * Used in updating error cache.
	 */
	T * errorArray;
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
	T errorUpdateHigh;
	/**
	 * Value of model->Y(iLow)*model->alphas(iLow) - lowOld.
	 */
	T erroUpdateLow;
};



#endif
