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

// An implementation of Peter J. Acklam's inverse cumulative normal distribution function
float acklamInvCND(float p){
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

    if(p <= 0 || p >= 1.0)
        return (float)0x7FFFFFFF;

    if(p < p_low){
        z = sqrt(-2.0 * log(p));
        z = (((((c1 * z + c2) * z + c3) * z + c4) * z + c5) * z + c6) /
            ((((d1 * z + d2) * z + d3) * z + d4) * z + 1.0);
    }else{
        if(p > p_high){
            z = sqrt(-2.0 * log(1.0 - p));
            z = -(((((c1 * z + c2) * z + c3) * z + c4) * z + c5) * z + c6) /
                 ((((d1 * z + d2) * z + d3) * z + d4) * z + 1.0);
        }else{
            z = p - 0.5;
            R = z * z;
            z = (((((a1 * R + a2) * R + a3) * R + a4) * R + a5) * R + a6) * z /
                (((((b1 * R + b2) * R + b3) * R + b4) * R + b5) * R + 1.0);
        }
    }

    return z;
}
