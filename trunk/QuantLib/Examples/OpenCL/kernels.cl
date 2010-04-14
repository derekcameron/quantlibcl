/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

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

// Testing kernels
__kernel void oclTest1(__global unsigned long* out) {
	unsigned long tid = get_global_id(0);
	out[tid] = tid;
}

// Pseudocode for a samples generator
// adapted from MonteCarloModel<MC,RNG,S>::addSamples(Size samples)
// in QuantLib : montecarlomodel.hpp

void generateSamples(ulong samples, float* ptr) {
	for(ulong i = 0; i < samples; i++) {
		
		// Generate one path
		float path = pathGenerator();
		
		// Calculate a value based on that path
		float price = calculateAssetValue(path);
		
		
	}
}

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