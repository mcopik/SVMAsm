#ifndef PARALLEL_ASM_DATA_H
#define PARALLEL_ASM_DATA_H

#include "../kernel/AbstractKernel.h"

struct ParallelAsmData {
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
	float cost;
	/**
	 * Approximation error value.
	 */
	float error;
	/**
	 * Pointer to array with Y.
	 * Used in finding iHigh,iLow.
	 */
	float * yArray;
	/**
	 * Pointer to array with alpha values.
	 * Used in finding iHigh,iLow.
	 */
	float * alphaArray;
	/**
	 * Pointer to array with error function values.
	 * Used in updating error cache.
	 */
	float * errorArray;
	/**
	 * Value of model->Y(iHigh)*model->alphas(iHigh) - highOld.
	 */
	float errorUpdateHigh;
	/**
	 * Value of model->Y(iLow)*model->alphas(iLow) - lowOld.
	 */
	float errorUpdateLow;
	/**
	 * Matrix with cached kernel.
	 */
	float * cachedKernelLow;
	float * cachedKernelHigh;
	/**
	 * Training matrix.
	 */
	float * X;
	int	numberOfFeatures;
	int	numberOfTrainExamples;
	int	offset;
};



#endif