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
#define TEMPERING_SHIFT_U 11
#define TEMPERING_SHIFT_S 7
#define TEMPERING_SHIFT_T 15
#define TEMPERING_SHIFT_L 18

/*!	OpenCL kernel for generating pseudorandom numbers using Mersenne Twister 19937
	\param dest Pointer to a memory address where the generated random sequence will be stored
	\param params Pointer to an array containing the MT parameters for each thread
	\param randN The number of random variables this thread will calculate before returning
*/