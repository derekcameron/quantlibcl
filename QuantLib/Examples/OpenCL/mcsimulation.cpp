/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*!
 Copyright (C) 2010 William Gross

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#include "mcsimulation.hpp"
#include <math.h>

void initializeMersenneTwister(const mt_params_stripped* mtParams, mt_state* mtState) {
    //Initialize current state
    mtState->mt[0] = mtParams->seed;
    for (mtState->iState = 1; mtState->iState < MT_NN; mtState->iState++)
        mtState->mt[mtState->iState] = (1812433253U * (mtState->mt[mtState->iState - 1] ^ (mtState->mt[mtState->iState - 1] >> 30)) + mtState->iState) & MT_WMASK;

	mtState->iState = 0;
    mtState->mti1 = mtState->mt[0];
}

inline float vanillaCallPayoff(const float S, const float X) {
    float payoff = S - X;
    return (payoff > 0) ? payoff : 0;
}

inline float vanillaPutPayoff(const float S, const float X) {
	float payoff = X - S;
	return (payoff > 0) ? payoff : 0;
}

inline void MersenneTwister_Bulk(float* randOutput,
			      const size_t count,
			      const mt_params_stripped* mtParams,
			      mt_state* mtState)
{
	//The following variables are not part of the Mersenne Twister's state
	uint32_t x;
	uint32_t mti, mtiM;

    for (uint32_t iOut = 0; iOut < count; iOut++) {
        mtState->iState1 = mtState->iState + 1;
        mtState->iStateM = mtState->iState + MT_MM;
        if(mtState->iState1 >= MT_NN) mtState->iState1 -= MT_NN;
        if(mtState->iStateM >= MT_NN) mtState->iStateM -= MT_NN;
        mti = mtState->mti1;
        mtState->mti1 = mtState->mt[mtState->iState1];
        mtiM = mtState->mt[mtState->iStateM];

	    // MT recurrence
        x = (mti & MT_UMASK) | (mtState->mti1 & MT_LMASK);
	    x = mtiM ^ (x >> 1) ^ ((x & 1) ? mtParams->matrix_a : 0);

        mtState->mt[mtState->iState] = x;
        mtState->iState = mtState->iState1;

        //Tempering transformation
        x ^= (x >> MT_SHIFT0);
        x ^= (x << MT_SHIFTB) & mtParams->mask_b;;
        x ^= (x << MT_SHIFTC) & mtParams->mask_c;
        x ^= (x >> MT_SHIFT1);

        //Convert to (0, 1] float and write to global memory
        randOutput[iOut] = ((float)x + 1.0f) / 4294967296.0f;
    }
}

// An implementation of Peter J. Acklam's inverse cumulative normal distribution function
// Generates one random variable
inline float AcklamInvCND(float input) {
    const float a1 = -39.6968302866538f;
    const float a2 = 220.946098424521f;
    const float a3 = -275.928510446969f;
    const float a4 = 138.357751867269f;
    const float a5 = -30.6647980661472f;
    const float a6 = 2.50662827745924f;
    const float b1 = -54.4760987982241f;
    const float b2 = 161.585836858041f;
    const float b3 = -155.698979859887f;
    const float b4 = 66.8013118877197f;
    const float b5 = -13.2806815528857f;
    const float c1 = -7.78489400243029E-03f;
    const float c2 = -0.322396458041136f;
    const float c3 = -2.40075827716184f;
    const float c4 = -2.54973253934373f;
    const float c5 = 4.37466414146497f;
    const float c6 = 2.93816398269878f;
    const float d1 = 7.78469570904146E-03f;
    const float d2 = 0.32246712907004f;
    const float d3 = 2.445134137143f;
    const float d4 = 3.75440866190742f;
    const float p_low = 0.02425f;
    const float p_high = 1.0f - p_low;
    float z, R;

	if((input) <= 0 || (input) >= 1.0)
		(input) = (float)0x7FFFFFFF;

	if((input) < p_low){
		z = (float)sqrt(-2.0 * log((input)));
		z = (((((c1 * z + c2) * z + c3) * z + c4) * z + c5) * z + c6) /
			((((d1 * z + d2) * z + d3) * z + d4) * z + 1.0f);
	}else{
		if((input) > p_high){
			z = (float)sqrt(-2.0 * log(1.0 - (input)));
			z = -(((((c1 * z + c2) * z + c3) * z + c4) * z + c5) * z + c6) /
				 ((((d1 * z + d2) * z + d3) * z + d4) * z + 1.0f);
		} else {
			z = (input) - 0.5f;
			R = z * z;
			z = (((((a1 * R + a2) * R + a3) * R + a4) * R + a5) * R + a6) * z /
				(((((b1 * R + b2) * R + b3) * R + b4) * R + b5) * R + 1.0f);
		}
	}

	return z;
}

// An implementation of Peter J. Acklam's inverse cumulative normal distribution function
// Generates bufferLength random variables
inline void AcklamInvCND_Bulk(float* randOutput, const size_t bufferLength) {
    const float a1 = -39.6968302866538f;
    const float a2 = 220.946098424521f;
    const float a3 = -275.928510446969f;
    const float a4 = 138.357751867269f;
    const float a5 = -30.6647980661472f;
    const float a6 = 2.50662827745924f;
    const float b1 = -54.4760987982241f;
    const float b2 = 161.585836858041f;
    const float b3 = -155.698979859887f;
    const float b4 = 66.8013118877197f;
    const float b5 = -13.2806815528857f;
    const float c1 = -7.78489400243029E-03f;
    const float c2 = -0.322396458041136f;
    const float c3 = -2.40075827716184f;
    const float c4 = -2.54973253934373f;
    const float c5 = 4.37466414146497f;
    const float c6 = 2.93816398269878f;
    const float d1 = 7.78469570904146E-03f;
    const float d2 = 0.32246712907004f;
    const float d3 = 2.445134137143f;
    const float d4 = 3.75440866190742f;
    const float p_low = 0.02425f;
    const float p_high = 1.0f - p_low;
    float z, R;
	float* target;

	for(uint32_t i = 0; i < bufferLength; i++) {
		target = &(randOutput[i]);

		if((*target) <= 0 || (*target) >= 1.0)
			(*target) = (float)0x7FFFFFFF;

		if((*target) < p_low){
			z = (float)sqrt(-2.0 * log((*target)));
			z = (((((c1 * z + c2) * z + c3) * z + c4) * z + c5) * z + c6) /
				((((d1 * z + d2) * z + d3) * z + d4) * z + 1.0f);
		}else{
			if((*target) > p_high){
				z = (float)sqrt(-2.0 * log(1.0 - (*target)));
				z = -(((((c1 * z + c2) * z + c3) * z + c4) * z + c5) * z + c6) /
					 ((((d1 * z + d2) * z + d3) * z + d4) * z + 1.0f);
			} else {
				z = (*target) - 0.5f;
				R = z * z;
				z = (((((a1 * R + a2) * R + a3) * R + a4) * R + a5) * R + a6) * z /
					(((((b1 * R + b2) * R + b3) * R + b4) * R + b5) * R + 1.0f);
			}
		}

		(*target) = z;
    }
}

inline float generateLognormalPath(OpenCL_Option* option, const uint32_t timeSteps, const mt_params_stripped* mtParams, mt_state* mtState) {
	const float dt = option->T / (float)timeSteps;

	float path = option->S;

	for(uint32_t i = 0; i < timeSteps; i++) {
		float rv;
		MersenneTwister_Bulk(&rv, 1, mtParams, mtState);
		rv = AcklamInvCND(rv);
		path *= exp(
					(option->R - (0.5f * option->V * option->V)) * dt +
					(option->V * sqrt(dt) * rv)
					);
	}

	return path;
}

void valueOptions(OpenCL_Option* d_Options, const uint32_t numberOfOptions, const uint32_t pathsPerOption, const uint32_t timeStepsPerOption, mt_params_stripped* d_MT) {
	for(size_t globalID = 0; globalID < numberOfOptions; globalID++) {
		OpenCL_Option* option = &d_Options[globalID];						//Pointer to this thread's allocated buffer
		const mt_params_stripped* mtParams = &d_MT[globalID];		//Pointer to this thread's mtParams struct
		mt_state mtState;													//Stores the Mersenne Twister state between calls to the MersenneTwister function

		initializeMersenneTwister(mtParams, &mtState);

		OptionAccumulator accum;											//Stores the sum of all generated paths
		float discountFactor = exp(-option->R * option->T);					//Discount factor

		accum.expectedCallValueSum = 0.0f;
		accum.expectedPutValueSum = 0.0f;

		//Iterate
		for(uint32_t j = 0; j < pathsPerOption; j++) {
			float futureUnderlyingValue = generateLognormalPath(option, timeStepsPerOption, mtParams, &mtState);
			accum.expectedCallValueSum += vanillaCallPayoff(futureUnderlyingValue, option->X);
			accum.expectedPutValueSum += vanillaPutPayoff(futureUnderlyingValue, option->X);
		}

		option->callValue = accum.expectedCallValueSum / (float)pathsPerOption * discountFactor;
		option->putValue = accum.expectedPutValueSum / (float)pathsPerOption * discountFactor;
	}
}
