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
	/**
	 * Trains model on given data.
	 */
	void train(TrainData<T> & train,AbstractKernel<T> & kernel,double sigma) {
		model = new TrainedModel<T,U>(train.X,train.Y);
				//std::move(kernel.cacheKernel(train.X,sigma)));
		Vector<U> E(model->X.rows);
		std::cout << "Training model: " << std::endl;
		std::cout << "Number of examples: " <<  model->X.rows << std::endl;
		std::cout << "Number of features: " <<  model->X.cols << std::endl;



		int maxPasses = 5;
		int pass = 0;
		int C = 1;
		int numChangedAlphas = 0;
		U error = 1e-3;
		U temp;
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<float> randDistr;
		unsigned int j = 0;
		U alpha_i_old = 0;
		U alpha_j_old = 0;
		U L,H;
		U eta;
		U b1,b2;
		int counter =0;
		U kernelVal = 0;
		while(pass < maxPasses) {
			numChangedAlphas = 0;
			//for each training example
			for(unsigned int i = 0;i < model->X.rows;++i) {
				E(i) = computeE(i,kernel,sigma);

				temp = E(i)*model->Y(i);
				if( (temp < (error*-1) && model->alphas(i) < C) ||
						(temp > error && model->alphas(i) > 0) ) {
					kernelVal = kernel.kernelFunction((U)i,(U)j,sigma);
					j = floor(model->X.rows *(((double)rand()) / (RAND_MAX)));
		            while (j == i)
			            j = floor(model->X.rows *(((double)rand()) / (RAND_MAX)));
		            E(j) = computeE(j,kernel,sigma);
		            alpha_i_old = model->alphas(i);
		            alpha_j_old = model->alphas(j);

		            if(model->Y(i) == model->Y(j)) {
		                L = max(0, model->alphas(j) + model->alphas(i) - C);
		                H = min(C, model->alphas(j) + model->alphas(i));
		            }
		            else {
		                L = max(0, model->alphas(j) - model->alphas(i));
		                H = min(C, C + model->alphas(j) - model->alphas(i));
		            }
		            if (L == H){
		                continue;
		            }
/*
		            eta = 2 * model->cachedKernel(i,j) - model->cachedKernel(i,i)
		            		- model->cachedKernel(j,j);*/
		            eta = 2 * kernelVal - kernel.kernelFunction((U)i,(U)i,sigma);
				            		- kernel.kernelFunction((U)j,(U)j,sigma);;
					if (eta >= 0){
						continue;
					}

		            model->alphas(j) = model->alphas(j) - (model->Y(j) * (E(i) - E(j))) / eta;

		            model->alphas(j) = min (H, model->alphas(j));
		            model->alphas(j) = max (L, model->alphas(j));

		            if (std::abs(model->alphas(j) - alpha_j_old) < error) {
		                model->alphas(j) = alpha_j_old;
		                continue;
		            }

					model->alphas(i) = model->alphas(i) + model->Y(i)*model->Y(j) \
						*(alpha_j_old - model->alphas(j));
					/*
					b1 = model->b - E(i) \
						- model->Y(i) * (model->alphas(i) - alpha_i_old) *  model->cachedKernel(i,j) \
						- model->Y(j) * (model->alphas(j) - alpha_j_old) *  model->cachedKernel(i,j);
					b2 = model->b - E(j) \
						- model->Y(i) * (model->alphas(i) - alpha_i_old) *  model->cachedKernel(i,j) \
						- model->Y(j) * (model->alphas(j) - alpha_j_old) *  model->cachedKernel(j,j);
					*/
					b1 = model->b - E(i) \
							- model->Y(i) * (model->alphas(i) - alpha_i_old) *  kernelVal \
							- model->Y(j) * (model->alphas(j) - alpha_j_old) *  kernelVal;
					b2 = model->b - E(j) \
							- model->Y(i) * (model->alphas(i) - alpha_i_old) *  kernelVal- \
							- model->Y(j) * (model->alphas(j) - alpha_j_old) *  kernel.kernelFunction(j,j,sigma);
					if (0 < model->alphas(i) && model->alphas(i) < C)
					   model->b = b1;
					else if (0 < model->alphas(j) && model->alphas(j) < C)
						model->b = b2;
					else
					   model->b = (b1+b2)/2;

					++numChangedAlphas;
				}
			}
			printf("numChanged: %d\n",numChangedAlphas);
			pass = numChangedAlphas == 0 ? pass + 1 : 0;
			printf("PASS: %d\n",pass);
		}
		/**
		idx = alphas > 0;
		model->X= X(idx,:);
		model->y= Y(idx);
		*/

		//model->w = ((alphas.*Y)'*X)';
	}
	/**
	 * Predicts values on test data.
	 */
	void predict(TrainData<T> & test) {

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
			kernel(i) = model->cachedKernel(i,row);
		}
		Vector<U> vec = (model->alphas.multiplyEachByEach(model->Y)).multiplyEachByEach(kernel);
		for(unsigned int i = 0;i < kernel.size;++i) {
			result += vec(i);
		}
		return result;
	}
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
	}
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
};



#endif /* SMOCLASSIFIER_H_ */
