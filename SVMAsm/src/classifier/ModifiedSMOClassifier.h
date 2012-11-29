
#ifndef MODIFIED_SMO_H
#define MODIFIED_SMO_H

#include <iostream>
#include <algorithm>
#include "AbstractClassifier.h"
#include "../data/TrainedModel.h"
#include "../data/TrainData.h"
#include "../data/Vector.h"
#include "../data/Matrix.h"
#include <random>
#include <cmath>
#include <fstream>
#include <cstdio>

#define DEBUG_OUTPUT true
#define TIME_OUTPUT false


template<class T,class U>
class ModifiedSMOClassifier: public AbstractClassifier<T,U> {
private:
	typedef T inputType;
	typedef U dataType;
	/**
	 * Pointer to kernel class.
	 */
	AbstractKernel<inputType> * kernel;
	/**
	 * Cached kernel values.
	 */
	Matrix<dataType> * cache = nullptr;
	/**
	 * Cached error function.
	 */
	Vector<dataType> * errorCache = nullptr;
	using AbstractClassifier<inputType,dataType>::C;
	using AbstractClassifier<inputType,dataType>::error;
	using AbstractClassifier<inputType,dataType>::epsilon;
public:
	using AbstractClassifier<T,U>::model;
	void train(TrainData<inputType> & train,AbstractKernel<inputType> & _kernel,
			bool cacheKernel = false) {

		if(cacheKernel) {
			model = new TrainedModel<inputType,dataType>
				(train.X,train.Y,std::move(_kernel.cacheKernel(train.X)));
			kernel = &_kernel;
		}
		else {
			model = new TrainedModel<inputType,dataType>
				(train.X,train.Y);
			kernel = &_kernel;
			if(cache == nullptr) {
				model->cachedKernel = new Matrix<dataType>(train.X.rows(),train.X.rows());
				std::fill_n(model->cachedKernel->matrixData(),model->cachedKernel->rows()
						*model->cachedKernel->cols(),-1);
			}
		}
		std::cout << "Training model: " << std::endl;
		std::cout << "Number of examples: " <<  model->X.rows() << std::endl;
		std::cout << "Number of features: " <<  model->X.cols() << std::endl;

		errorCache = new Vector<dataType>(train.X.rows());
		U fHigh = 0.0;
		int iHigh = 0;
		U fLow = 0.0;
		int iLow = 0;
		dataType alphaHighOld = 0.0;
		dataType alphaLowOld = 0.0;
		/**
		 * f(i) = -y(i)
		 */
		for(unsigned int i = 0;i < errorCache->size();++i) {
			(*errorCache)(i) = -1*train.Y(i);
		}
		getfHigh(fHigh,iHigh);
		getfLow(fLow,iLow);
		alphaHighOld = model->alphas(iHigh);
		alphaLowOld = model->alphas(iLow);
		updateAlpha(iHigh,iLow);
		std::cout << "alphaHigh " << alphaHighOld << " -> " << model->alphas(iHigh) << std::endl;
		std::cout << "alphaLow " << alphaLowOld << " -> " << model->alphas(iLow) << std::endl;
		std::cout << "fHigh " << fHigh << " fLow " << fLow << std::endl;
		//updateErrorCache(iHigh,iLow,alphaHighOld,alphaLowOld);
		//std::cout << "alphaHigh " << alphaHighOld << " -> " << model->alphas(iHigh) << std::endl;
		//std::cout << "alphaLow " << alphaLowOld << " -> " << model->alphas(iLow) << std::endl;
		//getfLow(fLow,iLow);
		//getfHigh(fHigh,iHigh);
		std::cout << "iHigh " << iHigh <<
			"	fHigh " << fHigh << " iLow "<< iLow << " fLow " << fLow << std::endl;
		//if(0 < -1*error)
		std::cout << "C " << C << " error " << error << std::endl;
		int nr_iter = 0;
		do {// {//(fHigh - fLow) < -1*error) {

			if(nr_iter == 558) {
				std::cout << "558" << std::endl;
			}
			updateErrorCache(iHigh,iLow,alphaHighOld,alphaLowOld);
			getfHigh(fHigh,iHigh);
			getfLow(fLow,iLow);
			//std::cout << "fHigh " << fHigh << " fLow " << fLow << std::endl;
			alphaHighOld = model->alphas(iHigh);
			alphaLowOld = model->alphas(iLow);
			updateAlpha(iHigh,iLow);
			//std::cout << "fHigh " << fHigh << " fLow " << fLow << std::endl;
			std::cout << "iHigh " << iHigh << " alphaHigh " << alphaHighOld << " -> " << model->alphas(iHigh) << std::endl;
			std::cout << "iLow " << iLow << " alphaLow " << alphaLowOld << " -> " << model->alphas(iLow) << std::endl;
std::cout << nr_iter++ << std::endl;
			//getfLow(fLow,iLow);
			//getfHigh(fHigh,iHigh);
			//std::cout << "fHigh " << fHigh << " fLow " << fLow << std::endl;
			//std::cout << fHigh - fLow << std::endl;
		}while(fLow > (fHigh+2*error));
		model->b = (fLow+fHigh)/2;
	}
	/**
	 * Predicts values on test data.
	 */
	Vector<T> predict(Matrix<T> & X/*Matrix<T> & kernel*/) {

		Vector<T> predicts(X.rows());
		U E;
		for(unsigned int i = 0;i < X.rows();++i) {
			E = 0;
			//for(unsigned int k = 0;k < model->alphas.size;++k)
			//	E += model->alphas(k)*model->Y(k)*	\
			//		kernel->kernelFunction(model->X(k),X(i),model->X.cols);
			for(unsigned int k = 0;k < model->alphas.size();++k)
				E += model->alphas(k)*kernel->kernelFunction(X(i),model->X(k),model->X.cols())
				*model->Y(k);
			std::cout << i << " " << E << std::endl;
			E -= model->b;
			std::cout << i << " " << E << std::endl;
			if(E >= 0) {
				predicts(i) = 1;
			}
			else {
				predicts(i) = 0;
			}
		}
		return predicts;

	}
	/**
	 * Destructor.
	 */
	~ModifiedSMOClassifier() {
		delete model;
		delete errorCache;
	}
private:
	/**
	 * Finds lowest error function value in I_high set.
	 * I_high set is build from union of three subsets:
	 * all i, where 0 < alpha(i) < C
	 * all i, where alpha(i) = 0 and Y(i) > 0
	 * all i, where alpha(i) = C and Y(i) < 0
	 * @param fHigh contains error value
	 * @param iHigh contains error index
	 */
	void getfHigh(dataType & fHigh,int & iHigh) {
		fHigh = 1000000000;
		iHigh = 0;
		for(unsigned int i = 0;i < model->alphas.size();++i) {
			if((model->alphas(i) > error && model->alphas(i) < (C-error)) ||
					(equalsWithTolerance(model->alphas(i),0) && model->Y(i) > 0) ||
					(equalsWithTolerance(model->alphas(i),C) && model->Y(i) < 0)) {
				if((*errorCache)(i) < fHigh) {
					fHigh = (*errorCache)(i);
					iHigh = i;
				}
			}
		}
	}
	/**
	 * Finds highest error function value in I_low set.
	 * I_low set is build from union of three subsets:
	 * all i, where 0 < alpha(i) < C
	 * all i, where alpha(i) = 0 and Y(i) < 0
	 * all i, where alpha(i) = C and Y(i) > 0
	 * @param fLow contains error value
	 * @param iLow contains error index
	 */
	void getfLow(U & fLow,int & iLow) {
		fLow = -1000000000;
		iLow = 0;
		for(unsigned int i = 0;i < model->alphas.size();++i) {
			if((model->alphas(i) > error && model->alphas(i) < (C-error)) ||
					(equalsWithTolerance(model->alphas(i),0) && model->Y(i) < 0) ||
					(equalsWithTolerance(model->alphas(i),C) && model->Y(i) > 0)) {
				if((*errorCache)(i) > fLow) {
					fLow = (*errorCache)(i);
					iLow = i;
				}
			}
		}
	}
	/**
	 * Computes new values of alpha.
	 * @param iHigh index of alphaHigh
	 * @param iLow index of alphaLow
	 */
	void updateAlpha(int iHigh,int iLow) {
		//compute boundaries for alphaLow
		dataType alphaLowLowerBound,alphaLowUpperBound;
		dataType eta = computeKernel(iHigh,iHigh) + computeKernel(iLow,iLow)
				- 2 * computeKernel(iHigh,iLow);
		dataType oldAlphaLow = model->alphas(iLow);
		if(model->Y(iHigh)*model->Y(iLow) < 0) {
			if(model->alphas(iLow)-model->alphas(iHigh) < 0) {
				alphaLowLowerBound = 0;
				alphaLowUpperBound = C + model->alphas(iLow)-model->alphas(iHigh);
			}
			else {
				alphaLowLowerBound = model->alphas(iLow)-model->alphas(iHigh);
				alphaLowUpperBound = C;
			}
		}
		else {
			if(model->alphas(iLow)+model->alphas(iHigh) < C) {
				alphaLowUpperBound = model->alphas(iLow)+model->alphas(iHigh);
				alphaLowLowerBound = 0;
			}
			else {
				alphaLowLowerBound = model->alphas(iLow)+model->alphas(iHigh) -C;
				alphaLowUpperBound = C;
			}
		}
		//compute alphaLow
		if(eta > 0) {
			model->alphas(iLow) += model->Y(iLow)*
					((*errorCache)(iHigh) - (*errorCache)(iLow))/eta;
			if(model->alphas(iLow) < alphaLowLowerBound) {
				model->alphas(iLow) = alphaLowLowerBound;
			}
			else if(model->alphas(iLow) > alphaLowUpperBound) {
				model->alphas(iLow) = alphaLowUpperBound;
			}
		}
		else {
			dataType slope = model->Y(iLow)*
					((*errorCache)(iHigh) - (*errorCache)(iLow));
			dataType delta = slope*(alphaLowUpperBound-alphaLowLowerBound);
			if(delta > 0) {
				if(slope > 0) {
					model->alphas(iLow) = alphaLowUpperBound;
				}
				else {
					model->alphas(iLow) = alphaLowUpperBound;
				}
			}
			else {
				model->alphas(iLow) = oldAlphaLow;
			}
		}
		//compute alphaHigh
		model->alphas(iHigh) += model->Y(iHigh)*model->Y(iLow)
				*(oldAlphaLow - model->alphas(iLow));
	}
	/**
	 * Updates error cache.
	 * @param iHigh index of modified alphaHigh
	 * @param iLow index of modified alphaLow
	 * @param highOld old value of alphaHigh
	 * @param lowOld old value of alphaLow
	 */
	void updateErrorCache(int iHigh,int iLow,dataType highOld,dataType lowOld) {
		Vector<dataType> & error = *errorCache;
		for(unsigned int i = 0;i < error.size();++i) {
			//std::cout << "error(" << i << ") " << error(i) <<
				//	" alph " << model->alphas(i) << " y "

					//<< model->Y(i) << std::endl;
			//std::cout << "KerneliHigh " << computeKernel(iHigh,i)
					//<< " KerneliLow " << computeKernel(iLow,i) << std::endl;
			error(i) += (model->alphas(iHigh) - highOld)*
					model->Y(iHigh)*computeKernel(iHigh,i) +
					(model->alphas(iLow) - lowOld)*
					model->Y(iLow)*computeKernel(iLow,i);
			//std::cout << "error(" << i << ") " << error(i) << std::endl;
		}
	}
	/**
	 * Checks if value first lies in interval <-error+second,second+error>.
	 * @param first
	 * @param second
	 * @param error
	 * @return true if |first-second| <= error
	 */
	bool equalsWithTolerance(dataType first, dataType second,
			dataType error = 10e-3) {
		if(std::abs(first-second) <= error) {
			return true;
		}
		else {
			return false;
		}
	}
	/**
	 * Computes kernel value for training examples a,b.
	 * @param a
	 * @param b
	 * @return kernel value
	 */
	U computeKernel(int a,int b) {
		if(cache != nullptr)
			return (*cache)(a,b);
		if((*model->cachedKernel)(a,b) == -1) {
			(*model->cachedKernel)(a,b) =
					kernel->kernelFunction(model->X(a),model->X(b),model->X.cols());
		}
		return (*model->cachedKernel)(a,b);
	}
};

#endif
