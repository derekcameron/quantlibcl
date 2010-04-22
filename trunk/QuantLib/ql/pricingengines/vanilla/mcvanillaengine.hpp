/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2003 Ferdinando Ametrano
 Copyright (C) 2003, 2004, 2005, 2007 StatPro Italia srl
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

/*! \file mcvanillaengine.hpp
    \brief Monte Carlo vanilla option engine
*/

#ifndef quantlib_mcvanilla_engine_hpp
#define quantlib_mcvanilla_engine_hpp

#include <ql/pricingengines/mcsimulation.hpp>
#include <ql/instruments/vanillaoption.hpp>

namespace QuantLib {

    //! Pricing engine for vanilla options using Monte Carlo simulation
    /*! \ingroup vanillaengines */
    template <template <class> class MC, class RNG,
              class S = Statistics, class Inst = VanillaOption, template <template <class> class,class,class> class SIM = McSimulation>
    class MCVanillaEngine : public Inst::engine,
                            public SIM<MC,RNG,S> {
      public:
        void calculate() const {
            SIM<MC,RNG,S>::calculate(requiredTolerance_,
                                              requiredSamples_,
                                              maxSamples_);
            this->results_.value = this->mcModel_->sampleAccumulator().mean();
            if (RNG::allowsErrorEstimate)
            this->results_.errorEstimate =
                this->mcModel_->sampleAccumulator().errorEstimate();
        }
      protected:
        typedef typename SIM<MC,RNG,S>::path_generator_type
            path_generator_type;
        typedef typename SIM<MC,RNG,S>::path_pricer_type
            path_pricer_type;
        typedef typename SIM<MC,RNG,S>::stats_type
            stats_type;
        typedef typename SIM<MC,RNG,S>::result_type
            result_type;
        // constructor
        MCVanillaEngine(const boost::shared_ptr<StochasticProcess>&,
                        Size timeSteps,
                        Size timeStepsPerYear,
                        bool brownianBridge,
                        bool antitheticVariate,
                        bool controlVariate,
                        Size requiredSamples,
                        Real requiredTolerance,
                        Size maxSamples,
                        BigNatural seed);
        // McSimulation implementation
        TimeGrid timeGrid() const;
        boost::shared_ptr<path_generator_type> pathGenerator() const {

            Size dimensions = process_->factors();
            TimeGrid grid = this->timeGrid();
            typename RNG::rsg_type generator =
                RNG::make_sequence_generator(dimensions*(grid.size()-1),seed_);
            return boost::shared_ptr<path_generator_type>(
                   new path_generator_type(process_, grid,
                                           generator, brownianBridge_));
        }
        result_type controlVariateValue() const;
        // data members
        boost::shared_ptr<StochasticProcess> process_;
        Size timeSteps_, timeStepsPerYear_;
        Size requiredSamples_, maxSamples_;
        Real requiredTolerance_;
        bool brownianBridge_;
        BigNatural seed_;
    };


    // template definitions

    template <template <class> class MC, class RNG, class S, class Inst, template <template <class> class,class,class> class SIM>
    inline MCVanillaEngine<MC,RNG,S,Inst,SIM>::MCVanillaEngine(
                          const boost::shared_ptr<StochasticProcess>& process,
                          Size timeSteps,
                          Size timeStepsPerYear,
                          bool brownianBridge,
                          bool antitheticVariate,
                          bool controlVariate,
                          Size requiredSamples,
                          Real requiredTolerance,
                          Size maxSamples,
                          BigNatural seed)
    : SIM<MC,RNG,S>(antitheticVariate, controlVariate),
      process_(process), timeSteps_(timeSteps),
      timeStepsPerYear_(timeStepsPerYear),
      requiredSamples_(requiredSamples), maxSamples_(maxSamples),
      requiredTolerance_(requiredTolerance),
      brownianBridge_(brownianBridge), seed_(seed) {
        QL_REQUIRE(timeSteps != Null<Size>() ||
                   timeStepsPerYear != Null<Size>(),
                   "no time steps provided");
        QL_REQUIRE(timeSteps == Null<Size>() ||
                   timeStepsPerYear == Null<Size>(),
                   "both time steps and time steps per year were provided");
        QL_REQUIRE(timeSteps != 0,
                   "timeSteps must be positive, " << timeSteps <<
                   " not allowed");
        QL_REQUIRE(timeStepsPerYear != 0,
                   "timeStepsPerYear must be positive, " << timeStepsPerYear <<
                   " not allowed");
        this->registerWith(process_);
    }

    template <template <class> class MC, class RNG, class S, class Inst, template <template <class> class,class,class> class SIM>
    inline typename MCVanillaEngine<MC,RNG,S,Inst,SIM>::result_type
    MCVanillaEngine<MC,RNG,S,Inst,SIM>::controlVariateValue() const {

        boost::shared_ptr<PricingEngine> controlPE =
            this->controlPricingEngine();
        QL_REQUIRE(controlPE,
                   "engine does not provide "
                   "control variation pricing engine");

        typename Inst::arguments* controlArguments =
                dynamic_cast<typename Inst::arguments*>(
                                                   controlPE->getArguments());

        QL_REQUIRE(controlArguments, "engine is using inconsistent arguments");

        *controlArguments = this->arguments_;
        controlPE->calculate();

        const typename Inst::results* controlResults =
                dynamic_cast<const typename Inst::results*>(
                                                     controlPE->getResults());
        QL_REQUIRE(controlResults,
                   "engine returns an inconsistent result type");

        return result_type(controlResults->value);
    }


    template <template <class> class MC, class RNG, class S, class Inst, template <template <class> class,class,class> class SIM>
    inline TimeGrid MCVanillaEngine<MC,RNG,S,Inst,SIM>::timeGrid() const {
        Date lastExerciseDate = this->arguments_.exercise->lastDate();
        Time t = process_->time(lastExerciseDate);
        if (this->timeSteps_ != Null<Size>()) {
            return TimeGrid(t, this->timeSteps_);
        } else if (this->timeStepsPerYear_ != Null<Size>()) {
            Size steps = static_cast<Size>(this->timeStepsPerYear_*t);
            return TimeGrid(t, std::max<Size>(steps, 1));
        } else {
            QL_FAIL("time steps not specified");
        }
    }

}


#endif