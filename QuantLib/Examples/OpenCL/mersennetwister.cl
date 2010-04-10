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

typedef struct {
	//Crucial data (required for MT generation)
	unsigned int aaa;
	unsigned int maskB, maskC;
	unsigned int seed;

	//Non-crucial data (can be discarded after initialization)
	int mm,nn,rr,ww;
	unsigned int *state;
	unsigned int wmask,umask,lmask;
	int shift0, shift1, shiftB, shiftC;
} mt_params;

// An implementation of Matsumoto and Nishimura's Mersenne pseudorandom number generator
// with period 2^19937-1

#define MT_RNG_COUNT 4096
#define MT_NN 624
#define MT_MM 397
#define MT_WMASK 0xFFFFFFFFUL
#define MT_UMASK 0x80000000UL
#define MT_LMASK 0x7FFFFFFFUL
#define TEMPERING_SHIFT_U(x)  (x >> 11)
#define TEMPERING_SHIFT_S(x)  (x << 7)
#define TEMPERING_SHIFT_T(x)  (x << 15)
#define TEMPERING_SHIFT_L(x)  (x >> 18)

/*!	OpenCL kernel for generating pseudorandom numbers using Mersenne Twister 19937
	\param dest Pointer to a memory address where the generated random sequence will be stored
	\param params Pointer to an array containing the MT parameters for each thread
	\param randN The number of random variables this thread will calculate before returning
*/

__kernel void MersenneTwister(__global float* dest, __global mt_params* params, int randN)
{
    int globalID = get_global_id(0);

    int iState, iState1, iStateM, iOut;
    unsigned int mti, mti1, mtiM, x;
    unsigned int mt[MT_NN];
    unsigned int matrix_a, mask_b, mask_c; 

    //Load bit-vector Mersenne Twister parameters
    matrix_a = params[globalID].aaa;
    mask_b   = params[globalID].maskB;
    mask_c   = params[globalID].maskC;
        
    //Initialize current state
    mt[0] = params[globalID].seed;
    for (iState = 1; iState < MT_NN; iState++)
        mt[iState] = (1812433253UL * (mt[iState - 1] ^ (mt[iState - 1] >> 30)) + iState) & MT_WMASK;

    iState = 0;
    mti1 = mt[0];

    for (iOut = 0; iOut < randN; iOut++) {
        iState1 = iState + 1;
        iStateM = iState + MT_MM;
        if(iState1 >= MT_NN) iState1 -= MT_NN;
        if(iStateM >= MT_NN) iStateM -= MT_NN;
        mti  = mti1;
        mti1 = mt[iState1];
        mtiM = mt[iStateM];

	    // MT recurrence
        x = (mti & MT_UMASK) | (mti1 & MT_LMASK);
	    x = mtiM ^ (x >> 1) ^ ((x & 1) ? matrix_a : 0);

        mt[iState] = x;
        iState = iState1;

        //Tempering transformation
		x ^= TEMPERING_SHIFT_U(x);
		x ^= TEMPERING_SHIFT_S(x) & mask_b;
		x ^= TEMPERING_SHIFT_T(x) & mask_c;
		x ^= TEMPERING_SHIFT_L(x);
		
        //Convert to (0, 1] float and write to global memory
        dest[globalID + iOut * MT_RNG_COUNT] = ((float)x + 1.0f) / 4294967296.0f;
    }
}
