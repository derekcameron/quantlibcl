/*
 * Copyright 1993-2009 NVIDIA Corporation.  All rights reserved.
 *
 * NVIDIA Corporation and its licensors retain all intellectual property and 
 * proprietary rights in and to this software and related documentation. 
 * Any use, reproduction, disclosure, or distribution of this software 
 * and related documentation without an express license agreement from
 * NVIDIA Corporation is strictly prohibited.
 *
 * Please refer to the applicable NVIDIA end user license agreement (EULA) 
 * associated with this source code for terms and conditions that govern 
 * your use of this NVIDIA software.
 * 
 */
 
 /*
 # Architecture details
 # --------------------
 # Each kernel calculates the call and put values for a set of vanilla European options
 #
 # At kernel execution, the following parameters are passed to each kernel
 # # uint32_t numberOfOptions:  Represents the number of options this thread will compute
 # # uint32_t timeSteps:  The number of time steps to simulate (i.e. the dimensionality of each generated path)
 # # struct* optionsParams:  Pointer to an array of option parameter structs, each containing:
 # # ## float S:  Underlying spot price at time t=0
 # # ## float X:  Option strike price
 # # ## float r:  Continuously compounded risk-free rate
 # # ## float T:  Time at option maturity
 # # ## V:  Volatility of the underlying asset spot price
 # # results:  Pointer to a buffer for storing the results.  The buffer will store one array of the following type for each option:
 # # ## float callValue:  The present value of a call option
 # # ##
 # # ## float putValue:  The present value of a put option
 # # ## 
 # # 
 # Prior to kernel execution, the host program allocates the following OpenCL device buffers:
 # 
 # Summary of algorithm:
 */
 
#define          MT_MM 9
#define          MT_NN 19
#define       MT_WMASK 0xFFFFFFFFU
#define       MT_UMASK 0xFFFFFFFEU
#define       MT_LMASK 0x1U
#define      MT_SHIFT0 12
#define      MT_SHIFTB 7
#define      MT_SHIFTC 15
#define      MT_SHIFT1 18
#define PI 3.14159265358979f

typedef struct{
  uint matrix_a;
  uint mask_b;
  uint mask_c;
  uint seed;
} mt_struct_stripped;

typedef struct {
    int iState, iState1, iStateM;
    uint mti1;
    uint mt[MT_NN];
} mt_state;

typedef struct {
	float expectedCallValueSum;
	float expectedPutValueSum;
	float callValueErrorEstimate;
	float putValueErrorEstimate;
} OptionAccumulator;

typedef struct {
	float S;			//Underlying spot price at t=0
	float X;			//Exercise price
	float V;			//Volatility
	float R;			//Risk-free rate
	float T;			//Time to maturity
	float callValue;	//Value of a call option
	float putValue;		//Value of a put option
} Option;

void initializeMersenneTwister(__global const mt_struct_stripped* mtParams, mt_state* mtState) {
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
			      __global const mt_struct_stripped* mtParams,
			      mt_state* mtState)
{
	//The following variables are not part of the Mersenne Twister's state
	uint x;
	uint mti, mtiM;
	
    for (uint iOut = 0; iOut < count; iOut++) {
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
    const float a1 = -39.6968302866538;
    const float a2 = 220.946098424521;
    const float a3 = -275.928510446969;
    const float a4 = 138.357751867269;
    const float a5 = -30.6647980661472;
    const float a6 = 2.50662827745924;
    const float b1 = -54.4760987982241;
    const float b2 = 161.585836858041;
    const float b3 = -155.698979859887;
    const float b4 = 66.8013118877197;
    const float b5 = -13.2806815528857;
    const float c1 = -7.78489400243029E-03;
    const float c2 = -0.322396458041136;
    const float c3 = -2.40075827716184;
    const float c4 = -2.54973253934373;
    const float c5 = 4.37466414146497;
    const float c6 = 2.93816398269878;
    const float d1 = 7.78469570904146E-03;
    const float d2 = 0.32246712907004;
    const float d3 = 2.445134137143;
    const float d4 = 3.75440866190742;
    const float p_low = 0.02425;
    const float p_high = 1.0 - p_low;
    float z, R;
		
	if((input) <= 0 || (input) >= 1.0)
		(input) = (float)0x7FFFFFFF;

	if((input) < p_low){
		z = sqrt(-2.0 * log((input)));
		z = (((((c1 * z + c2) * z + c3) * z + c4) * z + c5) * z + c6) /
			((((d1 * z + d2) * z + d3) * z + d4) * z + 1.0);
	}else{
		if((input) > p_high){
			z = sqrt(-2.0 * log(1.0 - (input)));
			z = -(((((c1 * z + c2) * z + c3) * z + c4) * z + c5) * z + c6) /
				 ((((d1 * z + d2) * z + d3) * z + d4) * z + 1.0);
		} else {
			z = (input) - 0.5;
			R = z * z;
			z = (((((a1 * R + a2) * R + a3) * R + a4) * R + a5) * R + a6) * z /
				(((((b1 * R + b2) * R + b3) * R + b4) * R + b5) * R + 1.0);
		}
	}

	return z;
}

// An implementation of Peter J. Acklam's inverse cumulative normal distribution function
// Generates bufferLength random variables
inline void AcklamInvCND_Bulk(float* randOutput, const size_t bufferLength) {
    const float a1 = -39.6968302866538;
    const float a2 = 220.946098424521;
    const float a3 = -275.928510446969;
    const float a4 = 138.357751867269;
    const float a5 = -30.6647980661472;
    const float a6 = 2.50662827745924;
    const float b1 = -54.4760987982241;
    const float b2 = 161.585836858041;
    const float b3 = -155.698979859887;
    const float b4 = 66.8013118877197;
    const float b5 = -13.2806815528857;
    const float c1 = -7.78489400243029E-03;
    const float c2 = -0.322396458041136;
    const float c3 = -2.40075827716184;
    const float c4 = -2.54973253934373;
    const float c5 = 4.37466414146497;
    const float c6 = 2.93816398269878;
    const float d1 = 7.78469570904146E-03;
    const float d2 = 0.32246712907004;
    const float d3 = 2.445134137143;
    const float d4 = 3.75440866190742;
    const float p_low = 0.02425;
    const float p_high = 1.0 - p_low;
    float z, R;
	float* target;

	for(uint i = 0; i < bufferLength; i++) {
		target = &(randOutput[i]);
		
		if((*target) <= 0 || (*target) >= 1.0)
			(*target) = (float)0x7FFFFFFF;

		if((*target) < p_low){
			z = sqrt(-2.0 * log((*target)));
			z = (((((c1 * z + c2) * z + c3) * z + c4) * z + c5) * z + c6) /
				((((d1 * z + d2) * z + d3) * z + d4) * z + 1.0);
		}else{
			if((*target) > p_high){
				z = sqrt(-2.0 * log(1.0 - (*target)));
				z = -(((((c1 * z + c2) * z + c3) * z + c4) * z + c5) * z + c6) /
					 ((((d1 * z + d2) * z + d3) * z + d4) * z + 1.0);
			} else {
				z = (*target) - 0.5;
				R = z * z;
				z = (((((a1 * R + a2) * R + a3) * R + a4) * R + a5) * R + a6) * z /
					(((((b1 * R + b2) * R + b3) * R + b4) * R + b5) * R + 1.0);
			}
		}

		(*target) = z;
    }
}

inline float generateLognormalPath(__global Option* option, const uint timeSteps, __global const mt_struct_stripped* mtParams, mt_state* mtState) {
	const float dt = option->T / (float)timeSteps;
	
	float path = option->S;
	
	for(uint i = 0; i < timeSteps; i++) {
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

__kernel void valueOptions(__global Option* d_Options, const uint numberOfOptions, const uint pathsPerOption, const uint timeStepsPerOption, __global mt_struct_stripped* d_MT) {
	const size_t globalID = get_global_id(0);							//This thread's id
    __global Option* option = &d_Options[globalID];						//Pointer to this thread's allocated buffer
    __global const mt_struct_stripped* mtParams = &d_MT[globalID];		//Pointer to this thread's mtParams struct
    mt_state mtState;													//Stores the Mersenne Twister state between calls to the MersenneTwister function
	
    initializeMersenneTwister(mtParams, &mtState);
    
    OptionAccumulator accum;											//Stores the sum of all generated paths
	float discountFactor = exp(-option->R * option->T);					//Discount factor
	
	accum.expectedCallValueSum = 0.0f;
	accum.expectedPutValueSum = 0.0f;
	
    //Iterate
	for(uint j = 0; j < pathsPerOption; j++) {
		float futureUnderlyingValue = generateLognormalPath(option, timeStepsPerOption, mtParams, &mtState);
		accum.expectedCallValueSum += vanillaCallPayoff(futureUnderlyingValue, option->X);
		accum.expectedPutValueSum += vanillaPutPayoff(futureUnderlyingValue, option->X);
	}
	
	option->callValue = accum.expectedCallValueSum / (float)pathsPerOption * discountFactor;
	option->putValue = accum.expectedPutValueSum / (float)pathsPerOption * discountFactor;
}
/*
// Pseudocode for 1-factor MonteCarlo simulation code, assuming no control-variation path pricer
// adapted from Quantlib : mcsimulation.hpp

__kernel void MC_Calculate(const float requiredTolerance) {
	float order, errorEstimate;
//	unsigned __int64 sambleNumber, nextBatch, minSamples = 1023, maxSamples = 2147483647;

//	sampleNumber = getCurrentSampleNumber();
//	if(sampleNumber <minSamples) {
//		generateSamples(minSamples - sampleNumber);
//		sampleNumber = getCurrentSampleNumber();
//	}
//	errorEstimate = generateErrorEstimate();
//	while(errorEstimate > requiredTolerance) {
//		order = (errorEstimate * errorEstimate) / (requiredTolerance * requiredTolerance);
//		nextBatch = max((double)sampleNumber * order * 0.8 - (double)sampleNumber, (double)minSamples);
//		nextBatch = min(nextBatch,maxSamples-sampleNumber);
//		sampleNumber += nextBatch;
//		generateSamples(nextBatch);
//		errorEstimate = generateErrorEstimate();
//	}

//	return averageOfSamples();
}

*/
