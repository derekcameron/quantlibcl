/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
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

/*! \file mcsimulation.hpp
    \brief Supporting header file for mcsimulation.cl
*/

#ifndef quantlib_opencl_example_mcsimulation_hpp
#define quantlib_opencl_example_mcsimulation_hpp

#include "mt_params.h"

#define      MT_SHIFT0 12
#define      MT_SHIFTB 7
#define      MT_SHIFTC 15
#define      MT_SHIFT1 18
#define PI 3.14159265358979f

typedef struct {
	float S;			//Underlying spot price at t=0
	float X;			//Exercise price
	float V;			//Volatility
	float R;			//Risk-free rate
	float T;			//Time to maturity
	float callValue;	//Value of a call option
	float putValue;		//Value of a put option
} OpenCL_Option;

typedef struct {
    int iState, iState1, iStateM;
    uint32_t mti1;
    uint32_t mt[MT_NN];
} mt_state;

typedef struct {
	float expectedCallValueSum;
	float expectedPutValueSum;
	float callValueErrorEstimate;
	float putValueErrorEstimate;
} OptionAccumulator;

void initializeMersenneTwister(const mt_params_stripped* mtParams, mt_state* mtState);
inline float vanillaCallPayoff(const float S, const float X);
inline float vanillaPutPayoff(const float S, const float X);
inline void MersenneTwister_Bulk(float* randOutput,
			      const size_t count,
			      const mt_params_stripped* mtParams,
			      mt_state* mtState);
inline float AcklamInvCND(float input);
inline void AcklamInvCND_Bulk(float* randOutput, const size_t bufferLength);
inline float generateLognormalPath(OpenCL_Option* option, const uint32_t timeSteps, const mt_params_stripped* mtParams, mt_state* mtState);
void valueOptions(OpenCL_Option* d_Options, const uint32_t numberOfOptions, const uint32_t pathsPerOption, const uint32_t timeStepsPerOption, mt_params_stripped* d_MT);
#endif
