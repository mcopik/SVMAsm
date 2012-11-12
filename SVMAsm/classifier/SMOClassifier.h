/*
 * SMOClassifier.h
 *
 *  Created on: Oct 12, 2012
 *      Author: mcopik
 */

#ifndef SMOCLASSIFIER_H_
#define SMOCLASSIFIER_H_
#include <iostream>
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

template<class T>
Matrix<T> loadMatrix2(std::ifstream & file) {
	unsigned int rows,cols;
	file >> rows >> cols;
	Matrix<T> X(rows,cols);
	for(unsigned int i = 0;i < rows;++i) {
		for(unsigned int j = 0;j < cols;++j) {
			file >> X(i,j);
		}
	}
	return std::move(X);
}

template<class T,class U>
class SMOClassifier: public AbstractClassifier<T,U> {
public:
	using AbstractClassifier<T,U>::model;
	AbstractKernel<T> * kernel;
	int counterKernel = 0;
	double sigma = 0;
	Matrix<U> * cache = nullptr;
	/**
	 * Trains model on given data.
	 */
	void train(TrainData<T> & train,AbstractKernel<T> & _kernel,double _sigma,bool cacheKernel = false) {

		if(cacheKernel) {
			model = new TrainedModel<T,U>(train.X,train.Y,std::move(_kernel.cacheKernel(train.X,_sigma)));
		}
		else {
			model = new TrainedModel<T,U>(train.X,train.Y);
			sigma = _sigma;
			kernel = &_kernel;
			if(cache == nullptr) {
				model->cachedKernel.rows = train.X.rows;
				model->cachedKernel.cols = train.X.rows;
				model->cachedKernel.data = new U[train.X.rows*train.X.rows]();
			}
		}
		std::cout << "Training model: " << std::endl;
		std::cout << "Number of examples: " <<  model->X.rows << std::endl;
		std::cout << "Number of features: " <<  model->X.cols << std::endl;

		/**
		 * Error cache.
		 */
		Vector<U> E(model->X.rows);
		int maxPasses = 5;
		int pass = 0;
		U C = 0.1;
		int numChangedAlphas = 0;
		U error = 1e-3;
		/**
		 * C++11 random
		 */
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<float> randDistr;

		unsigned int j = 0;
		U alpha_i_old = 0;
		U alpha_j_old = 0;
		U L,H;
		U eta;
		U b1,b2;
		U kernelVal = 0;
		int c = 0;
		int t = 0;
		while(pass < maxPasses) {
			numChangedAlphas = 0;
			//for each training example
			if(TIME_OUTPUT)
				c = clock();
			for(unsigned int i = 0;i < model->X.rows;++i) {

				if(TIME_OUTPUT)
					t=clock();
				if(model->alphas(i) == 0 || model->alphas(i) == C) {
					E(i) = 0;
					for(unsigned int k = 0;k < model->alphas.size;++k)
						if(model->alphas(k) != 0)
							E(i) += model->alphas(k)*model->Y(k)*computeKernel(k,i);
					E(i) = E(i) - model->Y(i)+ model->b;
				}
				if(TIME_OUTPUT)
					std::cout << "Compute E(i): " << ((double)(clock()-t))/CLOCKS_PER_SEC << std::endl;

				U temp = E(i)*model->Y(i);
				if( (temp < (error*-1) && model->alphas(i) < C) ||
						(temp > error && model->alphas(i) > 0) ) {
					if(TIME_OUTPUT)
						t=clock();
					kernelVal = computeKernel(i,j);
					if(TIME_OUTPUT)
						std::cout << "Compute kernelVal: " << ((double)(clock()-t))/CLOCKS_PER_SEC << std::endl;

					/**
					 * Choose j.
					 */
					j = floor(model->X.rows *(((double)rand()) / (RAND_MAX)));//randDistr(gen));
		            while (j == i)
			            j = floor(model->X.rows *(((double)rand()) / (RAND_MAX)));//randDistr(gen));

					if(TIME_OUTPUT)
						t=clock();
					/**
					 * Compute error on j if necessary.
					 */
					if(model->alphas(j) == 0 || model->alphas(j) == C) {
						E(j) = model->b-model->Y(j);
						for(unsigned int k = 0;k < model->alphas.size;++k)
							if(model->alphas(k) > 0)
								E(j) += model->alphas(k)*model->Y(k)*computeKernel(k,j);
					}
		            if(TIME_OUTPUT)
		            	std::cout << "Compute E(j): " << ((double)(clock()-t))/CLOCKS_PER_SEC << std::endl;

		            alpha_i_old = model->alphas(i);
		            alpha_j_old = model->alphas(j);

		            /**
		             * Equation (14)
		             */
		            if(model->Y(i) == model->Y(j)) {
		                L = max(0, model->alphas(j) + model->alphas(i) - C);
		                H = min(C, model->alphas(j) + model->alphas(i));
		            }
		            /**
		             * Equation (13)
		             */
		            else {
		                L = max(0, model->alphas(j) - model->alphas(i));
		                H = min(C, C + model->alphas(j) - model->alphas(i));
		            }
		            if (L == H){
		                continue;
		            }
		            /**
		             * Equation (15)
		             */
					if(TIME_OUTPUT)
						t=clock();
		            eta = 2 * kernelVal - computeKernel(i,i)
				            		- computeKernel(j,j);
		            if (eta >= 0)
						continue;
		            if(TIME_OUTPUT)
		            	std::cout << "Compute eta: " << ((double)(clock()-t))/CLOCKS_PER_SEC << std::endl;
		            /**
		             * Equation(16)
		             */
		            model->alphas(j) = model->alphas(j) - (model->Y(j) * (E(i) - E(j))) / eta;
		            /**
		             * Equation(17)
		             */
		            model->alphas(j) = min (H, model->alphas(j));
		            model->alphas(j) = max (L, model->alphas(j));

		            if (std::abs(model->alphas(j) - alpha_j_old) < error) {
		                model->alphas(j) = alpha_j_old;
		                continue;
		            }
		            /**
		             * Equation(18)
		             */
					model->alphas(i) = model->alphas(i) + model->Y(i)*model->Y(j) \
						*(alpha_j_old - model->alphas(j));
					/**
					 * Equation(20)
					 */
					b1 = model->b - E(i) \
							- model->Y(i) * (model->alphas(i) - alpha_i_old) *  kernelVal \
							- model->Y(j) * (model->alphas(j) - alpha_j_old) *  kernelVal;
					/**
					 * Equation(21)
					 */
					b2 = model->b - E(j) \
							- model->Y(i) * (model->alphas(i) - alpha_i_old) *  kernelVal- \
							- model->Y(j) * (model->alphas(j) - alpha_j_old) *  computeKernel(j,j);
					U old_bias = model->b;
					if (0 < model->alphas(i) && model->alphas(i) < C)
					   model->b = b1;
					else if (0 < model->alphas(j) && model->alphas(j) < C)
						model->b = b2;
					else
					   model->b = (b1+b2)/2;
					/**
					 * Update error cache.
					 */
					for (unsigned int k = 0; k < E.size; ++k) {
						if (0 !=  model->alphas(k) && C != model->alphas(k)){
							E(k) += model->Y(i) * (model->alphas(i) - alpha_i_old) * computeKernel( i,k)
							+ model->Y(j) * (model->alphas(j) - alpha_j_old) * computeKernel(j,k)+ (model->b - old_bias);
						}
					}
					++numChangedAlphas;
					if(TIME_OUTPUT)
						std::cout << "Loop end: " << ((double)(clock()-c))/CLOCKS_PER_SEC << std::endl;
				}
			}
			pass = numChangedAlphas == 0 ? pass + 1 : 0;
			if(TIME_OUTPUT)
				std::cout << "loop time: " << ((double)(clock()-c))/CLOCKS_PER_SEC << std::endl;
			if(DEBUG_OUTPUT) {
				std::cout << "numChanged: " << numChangedAlphas << std::endl;
				std::cout << "PASS: " << pass << std::endl;
				std::cout << "KERNEL COMPUTE: " << counterKernel << std::endl;
			}
		}
	}
	/**
	 * Predicts values on test data.
	 */
	Vector<T> predict(TrainData<T> & test,Matrix<T> & kernel) {
		Vector<T> predicts(test.X.rows);
		U E;
		/*
		for(int i = 0;i < test.X.rows;++i) {
			E = 0;
			for(unsigned int k = 0;k < model->alphas.size;++k)
				if(model->alphas(k) != 0)
					E += model->alphas(k)*model->Y(k)*kernel(i,k);
			//E = E - model->Y(i)+ model->b;
			if(E >= 0) {
				predicts(i) = 1;
			}
			else {
				predicts(i) = 0;
			}
		}*/

		for(int i = 0;i < test.X.rows;++i) {
			for(unsigned int k = 0;k < model->alphas.size;++k)
				E += model->alphas(k)*test.X(i,k) + model->b;
			if(E >= 0) {
				predicts(i) = 1;
			}
			else {
				predicts(i) = 0;
			}
		}
		return predicts;
	}
	void setCachedKernel(Matrix<U> & _cache) {
		cache = &_cache;
	}
	/**
	 * Destructor.
	 */
	~SMOClassifier() {
		delete model;
	}
private:
	U computeE(int row) {
		U result = model->b - model->Y(row);
		Vector<U> kernel(model->cachedKernel.rows);
		for(unsigned int i = 0; i < kernel.size;++i) {
			kernel(i) = computeKernel(i,row);
		}
		Vector<U> vec = (model->alphas.multiplyEachByEach(model->Y)).multiplyEachByEach(kernel);
		for(unsigned int i = 0;i < kernel.size;++i) {
			result += vec(i);
		}
		return result;
	}
	/*
	U computeE(int row,AbstractKernel<U> & kernel,double sigma) {
		U result = model->b - model->Y(row);
		Vector<U> kernel_vec(model->X.rows);
		U kernelVal = kernel.kernelFunction(1,0,sigma);
		Vector<U> temp(model->X.cols);
		U temp2;
		for(unsigned int i = 0; i < model->X.rows;++i) {
			for(unsigned int j = 0;j < model->X.cols;++j)
				temp(j) = model->X(row,j) - model->X(i,j);
			temp2 = 0;
			for(unsigned int j = 0;j < model->X.cols;++j)
				temp2 += pow(temp(j),2.0);
			kernel_vec(i) = pow(kernelVal,temp2);
		}
		Vector<U> vec = (model->alphas.multiplyEachByEach(model->Y)).multiplyEachByEach(kernel_vec);
		for(unsigned int i = 0;i < kernel_vec.size;++i) {
			result += vec(i);
		}
		return result;
	}*/
	U max(U arg1,U arg2) {
		if(arg1 >= arg2) {
			return arg1;
		}
		else {
			return arg2;
		}
	}
	U min(U arg1,U arg2) {
		return arg1 >= arg2 ? arg2 : arg1;
	}
	U computeKernel(int a,int b) {
		if(cache != nullptr)
			return (*cache)(a,b);

		//only for gaussian!
		if(a == b)
			return 1;
		if(model->cachedKernel(a,b) == 0) {
			U kernelVal = kernel->kernelFunction(1,0,sigma);
			U temp2;
			temp2 = 0;
			for(unsigned int j = 0;j < model->X.cols;++j)
				temp2 += pow(model->X(a,j) - model->X(b,j),2.0);
			model->cachedKernel(a,b) = pow(kernelVal,temp2);
			counterKernel++;
		}
		return model->cachedKernel(a,b);
	}
};



#endif /* SMOCLASSIFIER_H_ */
