/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003 Ferdinando Ametrano
 Copyright (C) 2007, 2008 StatPro Italia srl
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

/*! \file mceuropeanengineOCL.hpp
    \brief OpenCL Monte Carlo European option engine
*/

#ifndef quantlib_montecarlo_european_engine_opencl_hpp
#define quantlib_montecarlo_european_engine_opencl_hpp

#include <ql/pricingengines/vanilla/mcvanillaengine.hpp>
#include <ql/processes/blackscholesprocess.hpp>
#include <ql/termstructures/volatility/equityfx/blackconstantvol.hpp>
#include <ql/termstructures/volatility/equityfx/blackvariancecurve.hpp>
#include <ql/utilities/opencl.hpp>

namespace QuantLib {

    //! Monte Carlo European engine factory
    template <class RNG = PseudoRandom, class S = Statistics>
    class MakeMCEuropeanEngineOCL {
      public:
        MakeMCEuropeanEngineOCL(
                    const boost::shared_ptr<GeneralizedBlackScholesProcess>&);
        // named parameters
        MakeMCEuropeanEngineOCL& withSteps(Size steps);
        MakeMCEuropeanEngineOCL& withStepsPerYear(Size steps);
        MakeMCEuropeanEngineOCL& withBrownianBridge(bool b = true);
        MakeMCEuropeanEngineOCL& withSamples(Size samples);
        MakeMCEuropeanEngineOCL& withAbsoluteTolerance(Real tolerance);
        MakeMCEuropeanEngineOCL& withMaxSamples(Size samples);
        MakeMCEuropeanEngineOCL& withSeed(BigNatural seed);
        MakeMCEuropeanEngineOCL& withAntitheticVariate(bool b = true);
		MakeMCEuropeanEngineOCL& withOclDevice(boost::shared_ptr<OclDevice> oclDevice);
        // conversion to pricing engine
        operator boost::shared_ptr<PricingEngine>() const;
      private:
        boost::shared_ptr<GeneralizedBlackScholesProcess> process_;
		boost::shared_ptr<OclDevice> oclDevice_;
        bool antithetic_;
        Size steps_, stepsPerYear_, samples_, maxSamples_;
        Real tolerance_;
        bool brownianBridge_;
        BigNatural seed_;
    };

    // inline definitions

    template <class RNG, class S>
    inline MakeMCEuropeanEngineOCL<RNG,S>::MakeMCEuropeanEngineOCL(
             const boost::shared_ptr<GeneralizedBlackScholesProcess>& process)
    : process_(process), antithetic_(false),
      steps_(Null<Size>()), stepsPerYear_(Null<Size>()),
      samples_(Null<Size>()), maxSamples_(Null<Size>()),
      tolerance_(Null<Real>()), brownianBridge_(false), seed_(0) {}

    template <class RNG, class S>
    inline MakeMCEuropeanEngineOCL<RNG,S>&
    MakeMCEuropeanEngineOCL<RNG,S>::withSteps(Size steps) {
        steps_ = steps;
        return *this;
    }

    template <class RNG, class S>
    inline MakeMCEuropeanEngineOCL<RNG,S>&
    MakeMCEuropeanEngineOCL<RNG,S>::withStepsPerYear(Size steps) {
        stepsPerYear_ = steps;
        return *this;
    }

    template <class RNG, class S>
    inline MakeMCEuropeanEngineOCL<RNG,S>&
    MakeMCEuropeanEngineOCL<RNG,S>::withSamples(Size samples) {
        QL_REQUIRE(tolerance_ == Null<Real>(),
                   "tolerance already set");
        samples_ = samples;
        return *this;
    }

    template <class RNG, class S>
    inline MakeMCEuropeanEngineOCL<RNG,S>&
    MakeMCEuropeanEngineOCL<RNG,S>::withAbsoluteTolerance(Real tolerance) {
        QL_REQUIRE(samples_ == Null<Size>(),
                   "number of samples already set");
        QL_REQUIRE(RNG::allowsErrorEstimate,
                   "chosen random generator policy "
                   "does not allow an error estimate");
        tolerance_ = tolerance;
        return *this;
    }

    template <class RNG, class S>
    inline MakeMCEuropeanEngineOCL<RNG,S>&
    MakeMCEuropeanEngineOCL<RNG,S>::withMaxSamples(Size samples) {
        maxSamples_ = samples;
        return *this;
    }

    template <class RNG, class S>
    inline MakeMCEuropeanEngineOCL<RNG,S>&
    MakeMCEuropeanEngineOCL<RNG,S>::withSeed(BigNatural seed) {
        seed_ = seed;
        return *this;
    }

    template <class RNG, class S>
    inline MakeMCEuropeanEngineOCL<RNG,S>&
    MakeMCEuropeanEngineOCL<RNG,S>::withBrownianBridge(bool brownianBridge) {
        brownianBridge_ = brownianBridge;
        return *this;
    }

    template <class RNG, class S>
    inline MakeMCEuropeanEngineOCL<RNG,S>&
    MakeMCEuropeanEngineOCL<RNG,S>::withAntitheticVariate(bool b) {
        antithetic_ = b;
        return *this;
    }

	template <class RNG, class S>
	inline MakeMCEuropeanEngineOCL<RNG,S>&
	MakeMCEuropeanEngineOCL<RNG,S>::withOclDevice(boost::shared_ptr<OclDevice> oclDevice) {
		oclDevice_ = oclDevice;
		return *this;
	}

    template <class RNG, class S>
    inline
    MakeMCEuropeanEngineOCL<RNG,S>::operator boost::shared_ptr<PricingEngine>()
                                                                      const {
        QL_REQUIRE(steps_ != Null<Size>() || stepsPerYear_ != Null<Size>(),
                   "number of steps not given");
        QL_REQUIRE(steps_ == Null<Size>() || stepsPerYear_ == Null<Size>(),
                   "number of steps overspecified");

        return boost::shared_ptr<PricingEngine>(new
            MCEuropeanEngine<RNG,S>(process_,
                                    steps_,
                                    stepsPerYear_,
                                    brownianBridge_,
                                    antithetic_,
                                    samples_, tolerance_,
                                    maxSamples_,
                                    seed_));
    }

}


#endif
