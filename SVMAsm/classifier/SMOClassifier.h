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

template<class T,class U>
class SMOClassifier: public AbstractClassifier<T,U> {
public:
	using AbstractClassifier<T,U>::model;
	/**
	 * Trains model on given data.
	 */
	void train(TrainData<T> & train,AbstractKernel<T> & kernel,double sigma) {
		model = new TrainedModel<T,U>(train.X,train.Y,
				std::move(kernel.cacheKernel(train.X,sigma)));
		Vector<U> E(model->X.rows);
		std::cout << "Training model: " << std::endl;
		std::cout << "Number of examples: " <<  model->X.rows << std::endl;
		std::cout << "Number of features: " <<  model->X.cols << std::endl;
		int maxPasses = 5;
		int pass = 0;
		int C = 1;
		int numChangedAlphas = 0;
		double error = 10e-3;
		U temp;
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<float> randDistr;
		unsigned int j = 0;
		int alpha_i_old = 0;
		int alpha_j_old = 0;
		U L,H;
		U eta;
		U b1,b2;
		while(pass < maxPasses) {
			numChangedAlphas = 0;
			//for each training example
			for(unsigned int i = 0;i < model->X.rows;++i) {
				E(i) = computeE(i);
				std::cout << "E(" << i << ") " << E(i) << std::endl;
				temp = E(i)*model->Y(i);
				if( (temp < -error && model->alphas(i) < C) ||
						(temp > error && model->alphas(i) > 0) ) {
		            j = floor(model->X.rows *(((double)rand()) / (RAND_MAX)));//randDistr(gen));
		            while (j == i)
			            j = floor(model->X.rows *(((double)rand()) / (RAND_MAX)));// randDistr(gen));
		            std::cout << i << " " << j << std::endl;
		            E(j) = computeE(j);
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
		            std::cout << "L " << L << "H " << H << std::endl;
		            if (L == H)
		                continue;

		            eta = 2 * model->cachedKernel(i,j) - model->cachedKernel(i,i)
		            		- model->cachedKernel(j,j);
					if (eta >= 0)
						continue;

		            model->alphas(j) = model->alphas(j) - (model->Y(j) * (E(i) - E(j))) / eta;

		            model->alphas(j) = min (H, model->alphas(j));
		            model->alphas(j) = max (L, model->alphas(j));

		            if (std::abs(model->alphas(j) - alpha_j_old) < error) {
		                model->alphas(j) = alpha_j_old;
		                continue;
		            }

					model->alphas(i) = model->alphas(i) + model->Y(i)*model->Y(j) \
						*(alpha_j_old - model->alphas(j));

					b1 = model->b - E(i) \
						- model->Y(i) * (model->alphas(i) - alpha_i_old) *  model->cachedKernel(i,j) \
						- model->Y(j) * (model->alphas(j) - alpha_j_old) *  model->cachedKernel(i,j);
					b2 = model->b - E(j) \
						- model->Y(i) * (model->alphas(i) - alpha_i_old) *  model->cachedKernel(i,j) \
						- model->Y(j) * (model->alphas(j) - alpha_j_old) *  model->cachedKernel(j,j);

					if (0 < model->alphas(i) && model->alphas(i) < C)
					   model->b = b1;
					else if (0 < model->alphas(j) && model->alphas(j) < C)
					   model->b = b2;
					else
					   model->b = (b1+b2)/2;

					++numChangedAlphas;
					std::cout << "b: " << model->b << std::endl;
				}
			}

			pass = numChangedAlphas == 0 ? pass + 1 : 0;
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
