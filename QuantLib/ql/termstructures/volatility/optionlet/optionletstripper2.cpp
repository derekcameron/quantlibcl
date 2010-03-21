/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2007 Giorgio Facchinetti

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

#include <ql/termstructures/volatility/optionlet/optionletstripper2.hpp>
#include <ql/termstructures/volatility/optionlet/optionletstripper1.hpp>
#include <ql/termstructures/volatility/optionlet/strippedoptionletadapter.hpp>
#include <ql/termstructures/volatility/optionlet/spreadedoptionletvol.hpp>
#include <ql/termstructures/volatility/capfloor/capfloortermvolcurve.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/math/solvers1d/brent.hpp>
#include <ql/instruments/makecapfloor.hpp>
#include <ql/pricingengines/capfloor/blackcapfloorengine.hpp>
#include <ql/indexes/iborindex.hpp>


namespace QuantLib {

    OptionletStripper2::OptionletStripper2(
            const boost::shared_ptr<OptionletStripper1>& optionletStripper1,
            const Handle<CapFloorTermVolCurve>& atmCapFloorTermVolCurve)
    : OptionletStripper(optionletStripper1->termVolSurface(),
                        optionletStripper1->iborIndex()),
      stripper1_(optionletStripper1),
      atmCapFloorTermVolCurve_(atmCapFloorTermVolCurve),
      dc_(stripper1_->termVolSurface()->dayCounter()),
      nOptionExpiries_(atmCapFloorTermVolCurve->optionTenors().size()),
      atmCapFloorStrikes_(nOptionExpiries_),
      atmCapFloorPrices_(nOptionExpiries_),
      spreadsVolImplied_(nOptionExpiries_),
      caps_(nOptionExpiries_),
      maxEvaluations_(10000),
      accuracy_(1.e-6) {
        registerWith(stripper1_);
        registerWith(atmCapFloorTermVolCurve_);

        QL_REQUIRE(dc_ == atmCapFloorTermVolCurve->dayCounter(),
                   "different day counters provided");
     }

    void OptionletStripper2::performCalculations() const {

        //// optionletStripper data
        optionletDates_ = stripper1_->optionletFixingDates();
        optionletPaymentDates_ = stripper1_->optionletPaymentDates();
        optionletAccrualPeriods_ = stripper1_->optionletAccrualPeriods();
        optionletTimes_ = stripper1_->optionletFixingTimes();
        atmOptionletRate_ = stripper1_->atmOptionletRates();
        for (Size i=0; i<optionletTimes_.size(); ++i) {
            optionletStrikes_[i] = stripper1_->optionletStrikes(i);
            optionletVolatilities_[i] = stripper1_->optionletVolatilities(i);
        }

        // atmCapFloorTermVolCurve data
        const std::vector<Period>& optionExpiriesTenors =
                                    atmCapFloorTermVolCurve_->optionTenors();
        const std::vector<Time>& optionExpiriesTimes =
                                    atmCapFloorTermVolCurve_->optionTimes();

        for (Size j=0; j<nOptionExpiries_; ++j) {
            caps_[j] = MakeCapFloor(CapFloor::Cap,
                                    optionExpiriesTenors[j],
                                    iborIndex_,
                                    Null<Real>(),
                                    0*Days);
            atmCapFloorStrikes_[j] =
                caps_[j]->atmRate(**iborIndex_->forwardingTermStructure());
            Volatility atmOptionVol =
                atmCapFloorTermVolCurve_->volatility(optionExpiriesTimes[j],
                                                     atmCapFloorStrikes_[j]);
            boost::shared_ptr<BlackCapFloorEngine> engine(new
                    BlackCapFloorEngine(iborIndex_->forwardingTermStructure(),
                                        atmOptionVol, dc_));
            caps_[j]->setPricingEngine(engine);
            atmCapFloorPrices_[j] = caps_[j]->NPV();
        }

        spreadsVolImplied_ = spreadsVolImplied();

        StrippedOptionletAdapter adapter(stripper1_);

        Volatility unadjustedVol, adjustedVol;
        for (Size j=0; j<nOptionExpiries_; ++j) {
            for (Size i=0; i<optionletVolatilities_.size(); ++i) {
                if (i<=caps_[j]->floatingLeg().size()) {
                    unadjustedVol = adapter.volatility(optionletTimes_[i],
                                                       atmCapFloorStrikes_[j]);
                    adjustedVol = unadjustedVol + spreadsVolImplied_[j];

                    // insert adjusted volatility
                    std::vector<Rate>::const_iterator previous =
                        std::lower_bound(optionletStrikes_[i].begin(),
                                         optionletStrikes_[i].end(),
                                         atmCapFloorStrikes_[j]);
                    Size insertIndex = previous - optionletStrikes_[i].begin();

                    optionletStrikes_[i].insert(
                                optionletStrikes_[i].begin() + insertIndex,
                                atmCapFloorStrikes_[j]);
                    optionletVolatilities_[i].insert(
                                optionletVolatilities_[i].begin() + insertIndex,
                                adjustedVol);
                }
            }
        }
    }

    std::vector<Volatility> OptionletStripper2::spreadsVolImplied() const {

        Brent solver;
        std::vector<Volatility> result(nOptionExpiries_);
        Volatility guess = 0.0001, minSpread = -0.1, maxSpread = 0.1;
        for (Size j=0; j<nOptionExpiries_; ++j) {
            ObjectiveFunction f(stripper1_, caps_[j], atmCapFloorPrices_[j]);
            solver.setMaxEvaluations(maxEvaluations_);
            Volatility root = solver.solve(f, accuracy_, guess,
                                           minSpread, maxSpread);
            result[j] = root;
        }
        return result;
    }

    std::vector<Volatility> OptionletStripper2::spreadsVol() const {
        calculate();
        return spreadsVolImplied_;
    }

    std::vector<Rate> OptionletStripper2::atmCapFloorStrikes() const{
        calculate();
        return atmCapFloorStrikes_;
    }

    std::vector<Real> OptionletStripper2::atmCapFloorPrices() const {
        calculate();
        return atmCapFloorPrices_;
    }

//===========================================================================//
//                 OptionletStripper2::ObjectiveFunction                     //
//===========================================================================//

    OptionletStripper2::ObjectiveFunction::ObjectiveFunction(
            const boost::shared_ptr<OptionletStripper1>& optionletStripper1,
            const boost::shared_ptr<CapFloor>& cap,
            Real targetValue):
       stripper1_(optionletStripper1),
       cap_(cap),
       targetValue_(targetValue) {}

    Real OptionletStripper2::ObjectiveFunction::operator()(Volatility s) const
    {
        boost::shared_ptr<OptionletVolatilityStructure> adapter(new
            StrippedOptionletAdapter(stripper1_));

        boost::shared_ptr<SimpleQuote> spreadQuote(new SimpleQuote(s));

        boost::shared_ptr<OptionletVolatilityStructure> spreadedAdapter(new
            SpreadedOptionletVolatility(Handle<OptionletVolatilityStructure>(adapter),
                                 Handle<Quote>(spreadQuote)));

        boost::shared_ptr<BlackCapFloorEngine> engine(new
            BlackCapFloorEngine(
                stripper1_->iborIndex()->forwardingTermStructure(),
                Handle<OptionletVolatilityStructure>(spreadedAdapter)));

        cap_->setPricingEngine(engine);
        return cap_->NPV()-targetValue_;
    }
}