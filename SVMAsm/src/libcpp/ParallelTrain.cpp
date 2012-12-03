/*
 * ParallelTrain.cpp
 *
 *  Created on: Dec 3, 2012
 *      Author: mcopik
 */

#include "../data/ParallelTrainData.h"



bool equalsWithTolerance(float first, float second,
		float error = 10e-3) {
	if(std::abs(first-second) <= error) {
		return true;
	}
	else {
		return false;
	}
}

extern "C" void findHighLow(ParallelTrainData<float,float> * x) {
	float fHigh = 1000000000;
	float fLow = -1000000000;
	int iLow = 0;
	int iHigh = 0;
	for(unsigned int i = 0;i < x->trainDataSize;++i) {
		if((x->alphaArray[i] > x->error && x->alphaArray[i] < (x->cost-x->error)) ||
				(equalsWithTolerance(x->alphaArray[i],0) && x->yArray[i] > 0) ||
				(equalsWithTolerance(x->alphaArray[i],x->cost) && x->yArray[i] < 0)) {
			if(x->errorArray[i] < fHigh) {
				//std::cout << "High in " << i << " " << x->errorArray[i] << " " << x->alphaArray[i] << " " << x->yArray[i] << std::endl;
				fHigh = x->errorArray[i];
				iHigh = i;
			}
		}
		if((x->alphaArray[i] > x->error && x->alphaArray[i] < (x->cost-x->error)) ||
				(equalsWithTolerance(x->alphaArray[i],0) && x->yArray[i] < 0) ||
				(equalsWithTolerance(x->alphaArray[i],x->cost) && x->yArray[i] > 0)) {
			if(x->errorArray[i] > fLow) {
				std::cout << "Low in " << i << " " << x->errorArray[i] << " " << x->alphaArray[i] << " " << x->yArray[i] << std::endl;
				fLow = x->errorArray[i];
				iLow = i;
			}
		}
	}
	x->iHigh = iHigh;
	x->iLow = iLow;
}
//template void findHighLow<float,float>(ParallelTrainData<float,float> * x);

