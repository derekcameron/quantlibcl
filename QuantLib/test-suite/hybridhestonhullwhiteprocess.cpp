/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2007, 2008, 2009 Klaus Spanderen

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

#include "utilities.hpp"
#include "hybridhestonhullwhiteprocess.hpp"

#include <ql/time/schedule.hpp>
#include <ql/time/calendars/target.hpp>
#include <ql/quotes/simplequote.hpp>
#include <ql/instruments/europeanoption.hpp>
#include <ql/instruments/impliedvolatility.hpp>
#include <ql/processes/blackscholesprocess.hpp>
#include <ql/processes/hybridhestonhullwhiteprocess.hpp>
#include <ql/math/randomnumbers/rngtraits.hpp>
#include <ql/math/optimization/simplex.hpp>
#include <ql/math/optimization/levenbergmarquardt.hpp>
#include <ql/math/statistics/generalstatistics.hpp>
#include <ql/math/statistics/sequencestatistics.hpp>
#include <ql/math/statistics/incrementalstatistics.hpp>
#include <ql/math/matrixutilities/svd.hpp>
#include <ql/time/daycounters/actual360.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>
#include <ql/methods/montecarlo/multipathgenerator.hpp>
#include <ql/termstructures/yield/zerocurve.hpp>
#include <ql/termstructures/volatility/equityfx/blackconstantvol.hpp>
#include <ql/models/equity/hestonmodelhelper.hpp>
#include <ql/models/shortrate/onefactormodels/hullwhite.hpp>
#include <ql/pricingengines/vanilla/analytichestonengine.hpp>
#include <ql/pricingengines/vanilla/analyticeuropeanengine.hpp>
#include <ql/pricingengines/vanilla/mchestonhullwhiteengine.hpp>
#include <ql/pricingengines/vanilla/analyticbsmhullwhiteengine.hpp>
#include <ql/pricingengines/vanilla/analytichestonhullwhiteengine.hpp>
#include <ql/experimental/finitedifferences/fdhestonhullwhitevanillaengine.hpp>

#include <boost/bind.hpp>

using namespace QuantLib;
using namespace boost::unit_test_framework;

void HybridHestonHullWhiteProcessTest::testBsmHullWhiteEngine() {
    BOOST_MESSAGE("Testing European option pricing for a BSM process"
                  " with one-factor Hull-White model...");

    SavedSettings backup;

    DayCounter dc = Actual365Fixed();

    const Date today = Date::todaysDate();
    const Date maturity = today + Period(20, Years);

    Settings::instance().evaluationDate() = today;

    const Handle<Quote> spot(
                         boost::shared_ptr<Quote>(new SimpleQuote(100.0)));
    boost::shared_ptr<SimpleQuote> qRate(new SimpleQuote(0.04));
    const Handle<YieldTermStructure> qTS(flatRate(today, qRate, dc));
    boost::shared_ptr<SimpleQuote> rRate(new SimpleQuote(0.0525));
    const Handle<YieldTermStructure> rTS(flatRate(today, rRate, dc));
    boost::shared_ptr<SimpleQuote> vol(new SimpleQuote(0.25));
    const Handle<BlackVolTermStructure> volTS(flatVol(today, vol, dc));

    // FLOATING_POINT_EXCEPTION
    boost::shared_ptr<HullWhite> hullWhiteModel(
        new HullWhite(Handle<YieldTermStructure>(rTS), 0.00883, 0.00526));

    boost::shared_ptr<BlackScholesMertonProcess> stochProcess(
                      new BlackScholesMertonProcess(spot, qTS, rTS, volTS));

    boost::shared_ptr<Exercise> exercise(new EuropeanExercise(maturity));

    Real fwd = spot->value()*qTS->discount(maturity)/rTS->discount(maturity);
    boost::shared_ptr<StrikedTypePayoff> payoff(new
                                      PlainVanillaPayoff(Option::Call, fwd));

    EuropeanOption option(payoff, exercise);

    const Real tol = 1e-8;
    const Real corr[] = {-0.75, -0.25, 0.0, 0.25, 0.75 };
    const Volatility expectedVol[] = { 0.217064577, 0.243995801,
                                       0.256402830, 0.268236596, 0.290461343 };

    for (Size i=0; i < LENGTH(corr); ++i) {
        boost::shared_ptr<PricingEngine> bsmhwEngine(
                         new AnalyticBSMHullWhiteEngine(corr[i], stochProcess,
                                                        hullWhiteModel));

        option.setPricingEngine(bsmhwEngine);
        const Real npv = option.NPV();

        const Handle<BlackVolTermStructure> compVolTS(
                                        flatVol(today, expectedVol[i], dc));

        boost::shared_ptr<BlackScholesMertonProcess> bsProcess(
                    new BlackScholesMertonProcess(spot, qTS, rTS, compVolTS));
        boost::shared_ptr<PricingEngine> bsEngine(
                                       new AnalyticEuropeanEngine(bsProcess));

        EuropeanOption comp(payoff, exercise);
        comp.setPricingEngine(bsEngine);

        Volatility impliedVol =
            comp.impliedVolatility(npv, bsProcess, 1e-10, 100);

        if (std::fabs(impliedVol - expectedVol[i]) > tol) {
            BOOST_FAIL("Failed to reproduce implied volatility"
                       << "\n    calculated: " << impliedVol
                       << "\n    expected  : " << expectedVol[i]);
        }
        if (std::fabs((comp.NPV() - npv)/npv) > tol) {
            BOOST_FAIL("Failed to reproduce NPV"
                       << "\n    calculated: " << npv
                       << "\n    expected  : " << comp.NPV());
        }
        if (std::fabs(comp.delta() - option.delta()) > tol) {
            BOOST_FAIL("Failed to reproduce NPV"
                       << "\n    calculated: " << npv
                       << "\n    expected  : " << comp.NPV());
        }
        if (std::fabs((comp.gamma() - option.gamma())/npv) > tol) {
            BOOST_FAIL("Failed to reproduce NPV"
                       << "\n    calculated: " << npv
                       << "\n    expected  : " << comp.NPV());
        }
        if (std::fabs((comp.theta() - option.theta())/npv) > tol) {
            BOOST_FAIL("Failed to reproduce NPV"
                       << "\n    calculated: " << npv
                       << "\n    expected  : " << comp.NPV());
        }
        if (std::fabs((comp.vega() - option.vega())/npv) > tol) {
            BOOST_FAIL("Failed to reproduce NPV"
                       << "\n    calculated: " << npv
                       << "\n    expected  : " << comp.NPV());
        }
   }
}

void HybridHestonHullWhiteProcessTest::testCompareBsmHWandHestonHW() {
    BOOST_MESSAGE("Comparing European option pricing for a BSM process"
                  " with one-factor Hull-White model...");

    SavedSettings backup;

    DayCounter dc = Actual365Fixed();

    const Date today = Date::todaysDate();

    Settings::instance().evaluationDate() = today;

    const Handle<Quote> spot(
                         boost::shared_ptr<Quote>(new SimpleQuote(100.0)));
    std::vector<Date> dates;
    std::vector<Time> times;
    std::vector<Rate> rates, divRates;

    for (Size i=0; i <= 40; ++i) {
        dates.push_back(today+Period(i, Years));
        // FLOATING_POINT_EXCEPTION
        rates.push_back(0.01 + 0.0002*std::exp(std::sin(i/4.0)));
        divRates.push_back(0.02 + 0.0001*std::exp(std::sin(i/5.0)));
        times.push_back(dc.yearFraction(today, dates.back()));
    }

    const Handle<Quote> s0(boost::shared_ptr<Quote>(new SimpleQuote(100)));
    const Handle<YieldTermStructure> rTS(
       boost::shared_ptr<YieldTermStructure>(new ZeroCurve(dates, rates, dc)));
    const Handle<YieldTermStructure> qTS(
       boost::shared_ptr<YieldTermStructure>(
                                          new ZeroCurve(dates, divRates, dc)));

    boost::shared_ptr<SimpleQuote> vol(new SimpleQuote(0.25));
    const Handle<BlackVolTermStructure> volTS(flatVol(today, vol, dc));

    boost::shared_ptr<BlackScholesMertonProcess> bsmProcess(
                      new BlackScholesMertonProcess(spot, qTS, rTS, volTS));

    boost::shared_ptr<HestonProcess> hestonProcess(
                   new HestonProcess(rTS, qTS, spot,
                                     vol->value()*vol->value(), 1.0,
                                     vol->value()*vol->value(), 1e-4, 0.0));

    boost::shared_ptr<HestonModel> hestonModel(new HestonModel(hestonProcess));

    boost::shared_ptr<HullWhite> hullWhiteModel(
        new HullWhite(Handle<YieldTermStructure>(rTS), 0.01, 0.01));

    boost::shared_ptr<PricingEngine> bsmhwEngine(
             new AnalyticBSMHullWhiteEngine(0.0, bsmProcess, hullWhiteModel));

    boost::shared_ptr<PricingEngine> hestonHwEngine(
          new AnalyticHestonHullWhiteEngine(hestonModel, hullWhiteModel, 192));


    const Real tol = 1e-6;
    const Real strike[] = { 0.25, 0.5, 0.75, 0.8, 0.9,
                            1.0, 1.1, 1.2, 1.5, 2.0, 4.0 };
    const Size maturity[] = { 1, 2, 3, 5, 10, 15, 20, 25, 30 };
    const Option::Type types[] = { Option::Put, Option::Call };

    for (Size i=0; i < LENGTH(types); ++i) {
        for (Size j=0; j < LENGTH(strike); ++j) {
            for (Size l=0; l < LENGTH(maturity); ++l) {
                const Date maturityDate = today + Period(maturity[l], Years);

                boost::shared_ptr<Exercise> exercise(
                                         new EuropeanExercise(maturityDate));

                Real fwd = strike[j] * spot->value()
                    *qTS->discount(maturityDate)/rTS->discount(maturityDate);

                boost::shared_ptr<StrikedTypePayoff> payoff(new
                                          PlainVanillaPayoff(types[i], fwd));

                EuropeanOption option(payoff, exercise);

                option.setPricingEngine(bsmhwEngine);
                const Real calculated = option.NPV();

                option.setPricingEngine(hestonHwEngine);
                const Real expected = option.NPV();

                if (std::fabs(calculated-expected) > calculated*tol &&
                    std::fabs(calculated-expected) > tol) {
                    BOOST_ERROR("Failed to reproduce npvs"
                                << "\n    calculated: " << calculated
                                << "\n    expected  : " << expected
                                << "\n    strike    : " << strike[j]
                                << "\n    maturity  : " << maturity[l]
                                << "\n    type      : "
                                << ((types[i] == Option::Put)? "Put" : "Call")
                                );
                }
            }
        }
   }
}

void HybridHestonHullWhiteProcessTest::testZeroBondPricing() {
    BOOST_MESSAGE("Testing Monte-Carlo zero bond pricing...");

    SavedSettings backup;

    DayCounter dc = Actual360();
    const Date today = Date::todaysDate();

    Settings::instance().evaluationDate() = today;

    // construct a strange yield curve to check drifts and discounting
    // of the joint stochastic process

    std::vector<Date> dates;
    std::vector<Time> times;
    std::vector<Rate> rates;

    dates.push_back(today);
    rates.push_back(0.02);
    times.push_back(0.0);
    for (Size i=120; i < 240; ++i) {
        dates.push_back(today+Period(i, Months));
        rates.push_back(0.02 + 0.0002*std::exp(std::sin(i/8.0)));
        times.push_back(dc.yearFraction(today, dates.back()));
    }

    const Date maturity = dates.back() + Period(10, Years);
    dates.push_back(maturity);
    rates.push_back(0.04);
    times.push_back(dc.yearFraction(today, dates.back()));

    const Handle<Quote> s0(boost::shared_ptr<Quote>(new SimpleQuote(100)));

    const Handle<YieldTermStructure> ts(
       boost::shared_ptr<YieldTermStructure>(new ZeroCurve(dates, rates, dc)));
    const Handle<YieldTermStructure> ds(flatRate(today, 0.0, dc));

    const boost::shared_ptr<HestonProcess> hestonProcess(
            new HestonProcess(ts, ds, s0, 0.02, 1.0, 0.2, 0.5, -0.8));
    const boost::shared_ptr<HullWhiteForwardProcess> hwProcess(
                   new HullWhiteForwardProcess(ts, 0.05, 0.05));
    hwProcess->setForwardMeasureTime(dc.yearFraction(today, maturity));
    const boost::shared_ptr<HullWhite> hwModel(new HullWhite(ts, 0.05, 0.05));

    const boost::shared_ptr<HybridHestonHullWhiteProcess> jointProcess(
        new HybridHestonHullWhiteProcess(hestonProcess, hwProcess, -0.4));

    TimeGrid grid(times.begin(), times.end()-1);

    typedef PseudoRandom::rsg_type rsg_type;
    typedef MultiPathGenerator<rsg_type>::sample_type sample_type;

    BigNatural seed = 1234;
    rsg_type rsg = PseudoRandom::make_sequence_generator(
                              jointProcess->factors()*(grid.size()-1), seed);

    MultiPathGenerator<rsg_type> generator(jointProcess, grid, rsg, false);

    const Size m = 90;
    std::vector<GeneralStatistics> zeroStat(m);
    std::vector<GeneralStatistics> optionStat(m);

    const Size nrTrails = 5000;
    const Size optionTenor = 24;
    const DiscountFactor strike = 0.5;

    std::vector<DiscountFactor> tmpZero(m);
    std::vector<DiscountFactor> tmpOption(m);

    for (Size i=0; i < nrTrails; ++i) {
        const bool antithetic = (i%2)==0 ? false : true;
        sample_type path = (!antithetic) ? generator.next()
                                         : generator.antithetic();

        for (Size j=1; j < m; ++j) {
            const Time t = grid[j];            // zero end and option maturity
            const Time T = grid[j+optionTenor];// maturity of zero bond
                                               // of option

            Array states(3);
            Array optionStates(3);
            for (Size k=0; k < jointProcess->size(); ++k) {
                states[k]       = path.value[k][j];
                optionStates[k] = path.value[k][j+optionTenor];
            }

            const DiscountFactor zeroBond
                = 1.0/jointProcess->numeraire(t, states);
            const DiscountFactor zeroOption = zeroBond
                * std::max(0.0, hwModel->discountBond(t, T, states[2])-strike);

            if (antithetic) {
                zeroStat[j].add(0.5*(tmpZero[j] + zeroBond));
                optionStat[j].add(0.5*(tmpOption[j] + zeroOption));
            }
            else {
                tmpZero[j] = zeroBond;
                tmpOption[j] = zeroOption;
            }
        }
    }

    for (Size j=1; j < m; ++j) {
        const Time t = grid[j];
        Real calculated = zeroStat[j].mean();
        Real error = zeroStat[j].errorEstimate();
        Real expected = ts->discount(t);

        if (std::fabs(calculated - expected) > 5*error) {
            BOOST_ERROR("Failed to reproduce expected zero bond prices"
                        << "\n   t:          " << t
                        << "\n   calculated: " << calculated
                        << "\n   error:      " << error
                        << "\n   expected:   " << expected);
        }

        const Time T = grid[j+optionTenor];

        calculated = optionStat[j].mean();
        error = optionStat[j].errorEstimate();
        expected = hwModel->discountBondOption(Option::Call, strike, t, T);

        if (std::fabs(calculated - expected) > 5*error) {
            BOOST_ERROR("Failed to reproduce expected zero bond option prices"
                        << "\n   t:          " << t
                        << "\n   T:          " << T
                        << "\n   calculated: " << calculated
                        << "\n   error:      " << error
                        << "\n   expected:   " << expected);
        }
    }
}

void HybridHestonHullWhiteProcessTest::testMcVanillaPricing() {
    BOOST_MESSAGE("Testing Monte-Carlo vanilla option pricing...");

    SavedSettings backup;

    DayCounter dc = Actual360();
    const Date today = Date::todaysDate();

    Settings::instance().evaluationDate() = today;

    // construct a strange yield curve to check drifts and discounting
    // of the joint stochastic process

    std::vector<Date> dates;
    std::vector<Time> times;
    std::vector<Rate> rates, divRates;

    for (Size i=0; i <= 40; ++i) {
        dates.push_back(today+Period(i, Years));
        // FLOATING_POINT_EXCEPTION
        rates.push_back(0.03 + 0.0003*std::exp(std::sin(i/4.0)));
        divRates.push_back(0.02 + 0.0001*std::exp(std::sin(i/5.0)));
        times.push_back(dc.yearFraction(today, dates.back()));
    }

    const Date maturity = today + Period(20, Years);

    const Handle<Quote> s0(boost::shared_ptr<Quote>(new SimpleQuote(100)));
    const Handle<YieldTermStructure> rTS(
       boost::shared_ptr<YieldTermStructure>(new ZeroCurve(dates, rates, dc)));
    const Handle<YieldTermStructure> qTS(
       boost::shared_ptr<YieldTermStructure>(
                                          new ZeroCurve(dates, divRates, dc)));
    boost::shared_ptr<SimpleQuote> vol(new SimpleQuote(0.25));
    const Handle<BlackVolTermStructure> volTS(flatVol(today, vol, dc));

    const boost::shared_ptr<BlackScholesMertonProcess> bsmProcess(
              new BlackScholesMertonProcess(s0, qTS, rTS, volTS));
    const boost::shared_ptr<HestonProcess> hestonProcess(
              new HestonProcess(rTS, qTS, s0, 0.0625, 0.5, 0.0625, 1e-5, 0.3));
    const boost::shared_ptr<HullWhiteForwardProcess> hwProcess(
              new HullWhiteForwardProcess(rTS, 0.01, 0.01));
    hwProcess->setForwardMeasureTime(dc.yearFraction(today, maturity));

    const Real tol = 0.05;
    const Real corr[] = {-0.9, -0.5, 0.0, 0.5, 0.9 };
    const Real strike[] = { 100 };

    for (Size i=0; i < LENGTH(corr); ++i) {
        for (Size j=0; j < LENGTH(strike); ++j) {
            boost::shared_ptr<HybridHestonHullWhiteProcess> jointProcess(
                new HybridHestonHullWhiteProcess(hestonProcess,
                                                 hwProcess, corr[i]));

            boost::shared_ptr<StrikedTypePayoff> payoff(
                               new PlainVanillaPayoff(Option::Put, strike[j]));
            boost::shared_ptr<Exercise> exercise(
                               new EuropeanExercise(maturity));

            VanillaOption optionHestonHW(payoff, exercise);
            boost::shared_ptr<PricingEngine> engine =
                MakeMCHestonHullWhiteEngine<PseudoRandom>(jointProcess)
                .withSteps(1)
                .withAntitheticVariate()
                .withControlVariate()
                .withAbsoluteTolerance(tol)
                .withSeed(42);

            optionHestonHW.setPricingEngine(engine);

            const boost::shared_ptr<HullWhite> hwModel(
                        new HullWhite(Handle<YieldTermStructure>(rTS),
                                      hwProcess->a(), hwProcess->sigma()));

            VanillaOption optionBsmHW(payoff, exercise);
            optionBsmHW.setPricingEngine(
                boost::shared_ptr<PricingEngine>(
                           new AnalyticBSMHullWhiteEngine(corr[i], bsmProcess,
                                                          hwModel)));

            const Real calculated = optionHestonHW.NPV();
            const Real error      = optionHestonHW.errorEstimate();
            const Real expected   = optionBsmHW.NPV();

            if (  (corr[i] != 0.0 && std::fabs(calculated - expected) > 3*error)
                ||(corr[i] == 0.0 &&  std::fabs(calculated - expected)> 1e-4)) {
                BOOST_ERROR("Failed to reproduce BSM-HW vanilla prices"
                        << "\n   corr:       " << corr[i]
                        << "\n   strike:     " << strike[j]
                        << "\n   calculated: " << calculated
                        << "\n   error:      " << error
                        << "\n   expected:   " << expected);
            }
        }
    }
}


void HybridHestonHullWhiteProcessTest::testMcPureHestonPricing() {
    BOOST_MESSAGE("Testing Monte-Carlo Heston option pricing...");

    SavedSettings backup;

    DayCounter dc = Actual360();
    const Date today = Date::todaysDate();

    Settings::instance().evaluationDate() = today;

    // construct a strange yield curve to check drifts and discounting
    // of the joint stochastic process

    std::vector<Date> dates;
    std::vector<Time> times;
    std::vector<Rate> rates, divRates;

    for (Size i=0; i <= 100; ++i) {
        dates.push_back(today+Period(i, Months));
        // FLOATING_POINT_EXCEPTION
        rates.push_back(0.02 + 0.0002*std::exp(std::sin(i/10.0)));
        divRates.push_back(0.02 + 0.0001*std::exp(std::sin(i/20.0)));
        times.push_back(dc.yearFraction(today, dates.back()));
    }

    const Date maturity = today + Period(2, Years);

    const Handle<Quote> s0(boost::shared_ptr<Quote>(new SimpleQuote(100)));
    const Handle<YieldTermStructure> rTS(
       boost::shared_ptr<YieldTermStructure>(new ZeroCurve(dates, rates, dc)));
    const Handle<YieldTermStructure> qTS(
       boost::shared_ptr<YieldTermStructure>(
                                          new ZeroCurve(dates, divRates, dc)));

    const boost::shared_ptr<HestonProcess> hestonProcess(
              new HestonProcess(rTS, qTS, s0, 0.08, 1.5, 0.0625, 0.5, -0.8));
    const boost::shared_ptr<HullWhiteForwardProcess> hwProcess(
              new HullWhiteForwardProcess(rTS, 0.1, 1e-8));
    hwProcess->setForwardMeasureTime(dc.yearFraction(
                                        today, maturity+Period(1, Years)));

    const Real tol = 0.001;
    const Real corr[] = { -0.45, 0.45, 0.25 };
    const Real strike[] = { 100, 75, 50, 150 };

    for (Size i=0; i < LENGTH(corr); ++i) {
        for (Size j=0; j < LENGTH(strike); ++j) {
            boost::shared_ptr<HybridHestonHullWhiteProcess> jointProcess(
                new HybridHestonHullWhiteProcess(
                        hestonProcess, hwProcess,
                        corr[i], HybridHestonHullWhiteProcess::Euler));

            boost::shared_ptr<StrikedTypePayoff> payoff(
                               new PlainVanillaPayoff(Option::Put, strike[j]));
            boost::shared_ptr<Exercise> exercise(
                               new EuropeanExercise(maturity));

            VanillaOption optionHestonHW(payoff, exercise);
            VanillaOption optionPureHeston(payoff, exercise);
            optionPureHeston.setPricingEngine(
                boost::shared_ptr<PricingEngine>(
                    new AnalyticHestonEngine(
                          boost::shared_ptr<HestonModel>(
                                           new HestonModel(hestonProcess)))));

            Real expected   = optionPureHeston.NPV();

            optionHestonHW.setPricingEngine(
                MakeMCHestonHullWhiteEngine<PseudoRandom>(jointProcess)
                .withSteps(2)
                .withAntitheticVariate()
                .withControlVariate()
                .withAbsoluteTolerance(tol)
                .withSeed(42));

            Real calculated = optionHestonHW.NPV();
            Real error      = optionHestonHW.errorEstimate();

            if (   std::fabs(calculated - expected) > 3*error
                && std::fabs(calculated - expected) > tol) {
                BOOST_ERROR("Failed to reproduce pure heston vanilla prices"
                        << "\n   corr:       " << corr[i]
                        << "\n   strike:     " << strike[j]
                        << "\n   calculated: " << calculated
                        << "\n   error:      " << error
                        << "\n   expected:   " << expected);
            }
        }
    }
}


void HybridHestonHullWhiteProcessTest::testAnalyticHestonHullWhitePricing() {
    BOOST_MESSAGE("Testing analytic Heston Hull-White option pricing...");

    SavedSettings backup;

    DayCounter dc = Actual360();
    const Date today = Date::todaysDate();

    Settings::instance().evaluationDate() = today;

    // construct a strange yield curve to check drifts and discounting
    // of the joint stochastic process

    std::vector<Date> dates;
    std::vector<Time> times;
    std::vector<Rate> rates, divRates;

    for (Size i=0; i <= 40; ++i) {
        dates.push_back(today+Period(i, Years));
        // FLOATING_POINT_EXCEPTION
        rates.push_back(0.03 + 0.0001*std::exp(std::sin(i/4.0)));
        divRates.push_back(0.02 + 0.0002*std::exp(std::sin(i/3.0)));
        times.push_back(dc.yearFraction(today, dates.back()));
    }

    const Date maturity = today + Period(5, Years);
    const Handle<Quote> s0(boost::shared_ptr<Quote>(new SimpleQuote(100)));
    const Handle<YieldTermStructure> rTS(
       boost::shared_ptr<YieldTermStructure>(new ZeroCurve(dates, rates, dc)));
    const Handle<YieldTermStructure> qTS(
       boost::shared_ptr<YieldTermStructure>(
                                          new ZeroCurve(dates, divRates, dc)));

    const boost::shared_ptr<HestonProcess> hestonProcess(
            new HestonProcess(rTS, qTS, s0, 0.08, 1.5, 0.0625, 0.5, -0.8));
    const boost::shared_ptr<HestonModel> hestonModel(
                                            new HestonModel(hestonProcess));

    const boost::shared_ptr<HullWhiteForwardProcess> hwFwdProcess(
              new HullWhiteForwardProcess(rTS, 0.01, 0.01));
    hwFwdProcess->setForwardMeasureTime(dc.yearFraction(today, maturity));
    const boost::shared_ptr<HullWhite> hullWhiteModel(new HullWhite(
                               rTS, hwFwdProcess->a(), hwFwdProcess->sigma()));

    const Real tol = 0.001;
    const Real strike[] = { 80, 120 };
    const Option::Type types[] = { Option::Put, Option::Call };

    for (Size i=0; i < LENGTH(types); ++i) {
        for (Size j=0; j < LENGTH(strike); ++j) {
            boost::shared_ptr<HybridHestonHullWhiteProcess> jointProcess(
                new HybridHestonHullWhiteProcess(
                        hestonProcess, hwFwdProcess, 0.0,
                        HybridHestonHullWhiteProcess::Euler));

            boost::shared_ptr<StrikedTypePayoff> payoff(
                               new PlainVanillaPayoff(types[i], strike[j]));
            boost::shared_ptr<Exercise> exercise(
                               new EuropeanExercise(maturity));

            VanillaOption optionHestonHW(payoff, exercise);
            optionHestonHW.setPricingEngine(
                    MakeMCHestonHullWhiteEngine<PseudoRandom>(jointProcess)
                    .withSteps(1)
                    .withAntitheticVariate()
                    .withControlVariate()
                    .withAbsoluteTolerance(tol)
                    .withSeed(42));

            VanillaOption optionPureHeston(payoff, exercise);
            optionPureHeston.setPricingEngine(
                boost::shared_ptr<PricingEngine>(
                    new AnalyticHestonHullWhiteEngine(hestonModel,
                                                      hullWhiteModel, 192)));

            Real calculated = optionHestonHW.NPV();
            Real error      = optionHestonHW.errorEstimate();
            Real expected   = optionPureHeston.NPV();

            if (   std::fabs(calculated - expected) > 3*error
                && std::fabs(calculated - expected) > tol) {
                BOOST_ERROR("Failed to reproduce hw heston vanilla prices"
                        << "\n   strike:     " << strike[j]
                        << "\n   calculated: " << calculated
                        << "\n   error:      " << error
                        << "\n   expected:   " << expected);
            }
        }
    }
}

void HybridHestonHullWhiteProcessTest::testCallableEquityPricing() {
    BOOST_MESSAGE("Testing the pricing of a callable equity product...");

    SavedSettings backup;

    /*
       For the definition of the example product see
       Alexander Giese, On the Pricing of Auto-Callable Equity
       Structures in the Presence of Stochastic Volatility and
       Stochastic Interest Rates .
       http://workshop.mathfinance.de/2006/papers/giese/slides.pdf
    */

    const Size maturity = 7;
    DayCounter dc = Actual365Fixed();
    const Date today = Date::todaysDate();

    Settings::instance().evaluationDate() = today;

    Handle<Quote> spot(boost::shared_ptr<Quote>(new SimpleQuote(100.0)));
    boost::shared_ptr<SimpleQuote> qRate(new SimpleQuote(0.04));
    Handle<YieldTermStructure> qTS(flatRate(today, qRate, dc));
    boost::shared_ptr<SimpleQuote> rRate(new SimpleQuote(0.04));
    Handle<YieldTermStructure> rTS(flatRate(today, rRate, dc));

    const boost::shared_ptr<HestonProcess> hestonProcess(
            new HestonProcess(rTS, qTS, spot, 0.0625, 1.0,
                              0.24*0.24, 1e-4, 0.0));
    // FLOATING_POINT_EXCEPTION
    const boost::shared_ptr<HullWhiteForwardProcess> hwProcess(
            new HullWhiteForwardProcess(rTS, 0.00883, 0.00526));
    hwProcess->setForwardMeasureTime(
                      dc.yearFraction(today, today+Period(maturity+1, Years)));

    const boost::shared_ptr<HybridHestonHullWhiteProcess> jointProcess(
        new HybridHestonHullWhiteProcess(hestonProcess, hwProcess, -0.4));

    Schedule schedule(today, today + Period(maturity, Years),
                      Period(1, Years), TARGET(),
                      Following, Following,
                      DateGeneration::Forward, false);

    std::vector<Time> times(maturity+1);
    std::transform(schedule.begin(), schedule.end(), times.begin(),
                   boost::bind(&Actual365Fixed::yearFraction,
                               dc, today, _1, Date(), Date()));

    for (Size i=0; i<=maturity; ++i)
        times[i] = static_cast<Time>(i);

    TimeGrid grid(times.begin(), times.end());

    std::vector<Real> redemption(maturity);
    for (Size i=0; i < maturity; ++i) {
        redemption[i] = 1.07 + 0.03*i;
    }

    typedef PseudoRandom::rsg_type rsg_type;
    typedef MultiPathGenerator<rsg_type>::sample_type sample_type;

    BigNatural seed = 42;
    rsg_type rsg = PseudoRandom::make_sequence_generator(
                              jointProcess->factors()*(grid.size()-1), seed);

    MultiPathGenerator<rsg_type> generator(jointProcess, grid, rsg, false);
    GeneralStatistics stat;

    Real antitheticPayoff=0;
    const Size nrTrails = 40000;
    for (Size i=0; i < nrTrails; ++i) {
        const bool antithetic = (i%2)==0 ? false : true;

        sample_type path = antithetic ? generator.antithetic()
                                      : generator.next();

        Real payoff=0;
        for (Size j=1; j <= maturity; ++j) {
            if (path.value[0][j] > spot->value()) {
                Array states(3);
                for (Size k=0; k < 3; ++k) {
                    states[k] = path.value[k][j];
                }
                payoff = redemption[j-1]
                    / jointProcess->numeraire(grid[j], states);
                break;
            }
            else if (j == maturity) {
                Array states(3);
                for (Size k=0; k < 3; ++k) {
                    states[k] = path.value[k][j];
                }
                payoff = 1.0 / jointProcess->numeraire(grid[j], states);
            }
        }

        if (antithetic){
            stat.add(0.5*(antitheticPayoff + payoff));
        }
        else {
            antitheticPayoff = payoff;
        }
    }

    const Real expected = 0.938;
    const Real calculated = stat.mean();
    const Real error = stat.errorEstimate();

    if (std::fabs(expected - calculated) > 3*error) {
        BOOST_ERROR("Failed to reproduce auto-callable equity structure price"
                    << "\n   calculated: " << calculated
                    << "\n   error:      " << error
                    << "\n   expected:   " << expected);
    }
}

void HybridHestonHullWhiteProcessTest::testDiscretizationError() {
    BOOST_MESSAGE("Testing the discretization error of the "
                  "Heston Hull-White process...");

    SavedSettings backup;

    DayCounter dc = Actual360();
    const Date today = Date::todaysDate();

    Settings::instance().evaluationDate() = today;

    // construct a strange yield curve to check drifts and discounting
    // of the joint stochastic process

    std::vector<Date> dates;
    std::vector<Time> times;
    std::vector<Rate> rates, divRates;

    for (Size i=0; i <= 31; ++i) {
        dates.push_back(today+Period(i, Years));
        // FLOATING_POINT_EXCEPTION
        rates.push_back(0.04 + 0.0001*std::exp(std::sin(double(i))));
        divRates.push_back(0.04 + 0.0001*std::exp(std::sin(double(i))));
        times.push_back(dc.yearFraction(today, dates.back()));
    }

    const Date maturity = today + Period(10, Years);
    const Volatility v = 0.25;

    const Handle<Quote> s0(boost::shared_ptr<Quote>(new SimpleQuote(100)));
    const boost::shared_ptr<SimpleQuote> vol(new SimpleQuote(v));
    const Handle<BlackVolTermStructure> volTS(flatVol(today, vol, dc));
    const Handle<YieldTermStructure> rTS(
       boost::shared_ptr<YieldTermStructure>(new ZeroCurve(dates, rates, dc)));
    const Handle<YieldTermStructure> qTS(
       boost::shared_ptr<YieldTermStructure>(
                                          new ZeroCurve(dates, divRates, dc)));

    const boost::shared_ptr<BlackScholesMertonProcess> bsmProcess(
                          new BlackScholesMertonProcess(s0, qTS, rTS, volTS));

    const boost::shared_ptr<HestonProcess> hestonProcess(
           new HestonProcess(rTS, qTS, s0, v*v, 1, v*v, 1e-6, -0.4));

    const boost::shared_ptr<HullWhiteForwardProcess> hwProcess(
              new HullWhiteForwardProcess(rTS, 0.01, 0.01));
    hwProcess->setForwardMeasureTime(20.1472222222222222);

    const Real tol = 0.05;
    const Real corr[] = {-0.85, 0.5 };
    const Real strike[] = { 50, 100, 125 };

    for (Size i=0; i < LENGTH(corr); ++i) {
        for (Size j=0; j < LENGTH(strike); ++j) {
            boost::shared_ptr<StrikedTypePayoff> payoff(
                               new PlainVanillaPayoff(Option::Put, strike[j]));
            boost::shared_ptr<Exercise> exercise(
                               new EuropeanExercise(maturity));

            VanillaOption optionBsmHW(payoff, exercise);
            const boost::shared_ptr<HullWhite> hwModel(new HullWhite(
                               rTS, hwProcess->a(), hwProcess->sigma()));
            optionBsmHW.setPricingEngine(
                boost::shared_ptr<PricingEngine>(
                         new AnalyticBSMHullWhiteEngine(corr[i], bsmProcess,
                                                        hwModel)));

            Real expected = optionBsmHW.NPV();

            VanillaOption optionHestonHW(payoff, exercise);
            boost::shared_ptr<HybridHestonHullWhiteProcess> jointProcess(
                new HybridHestonHullWhiteProcess(hestonProcess,
                                                 hwProcess, corr[i]));
            optionHestonHW.setPricingEngine(
                    MakeMCHestonHullWhiteEngine<PseudoRandom>(jointProcess)
                    .withSteps(1)
                    .withAntitheticVariate()
                    .withAbsoluteTolerance(tol)
                    .withSeed(42));

            Real calculated = optionHestonHW.NPV();
            Real error      = optionHestonHW.errorEstimate();

            if ((   std::fabs(calculated - expected) > 3*error
                 && std::fabs(calculated - expected) > 1e-5)) {
                BOOST_ERROR("Failed to reproduce discretization error"
                        << "\n   corr:       " << corr[i]
                        << "\n   strike:     " << strike[j]
                        << "\n   calculated: " << calculated
                        << "\n   error:      " << error
                        << "\n   expected:   " << expected);
            }
        }
    }
}

void HybridHestonHullWhiteProcessTest::testFdmHestonHullWhiteEngine() {
    BOOST_MESSAGE("Testing the FDM Heston Hull-White engine...");

    SavedSettings backup;

    const Date today = Date(28, March, 2004);
    Settings::instance().evaluationDate() = today;
    const Date exerciseDate = Date(28, March, 2012);
    DayCounter dc = Actual365Fixed();
    
    Handle<Quote> s0(boost::shared_ptr<Quote>(new SimpleQuote(100.0)));

    const Handle<YieldTermStructure> rTS(flatRate(0.05, dc));
    const Handle<YieldTermStructure> qTS(flatRate(0.02, dc));

    const Volatility vol = 0.30;
    const Handle<BlackVolTermStructure> volTS(flatVol(vol, dc));

    const Real v0 = vol*vol;
    boost::shared_ptr<HestonProcess> hestonProcess(
        new HestonProcess(rTS, qTS, s0, v0, 1.0, v0, 0.000001, 0.0));

    boost::shared_ptr<BlackScholesMertonProcess> stochProcess(
                      new BlackScholesMertonProcess(s0, qTS, rTS, volTS));

    boost::shared_ptr<HullWhiteProcess> hwProcess(
                              new HullWhiteProcess(rTS, 0.00883, 0.01));
    boost::shared_ptr<HullWhite> hwModel(
                    new HullWhite(rTS, hwProcess->a(), hwProcess->sigma()));
    
    boost::shared_ptr<Exercise> exercise(new EuropeanExercise(exerciseDate));
    const Real corr[] = {-0.85, 0.5 };
    const Real strike[] = { 75, 120, 160 };

    for (Size i=0; i < LENGTH(corr); ++i) {
        for (Size j=0; j < LENGTH(strike); ++j) {
            boost::shared_ptr<StrikedTypePayoff> payoff(
                              new PlainVanillaPayoff(Option::Call, strike[j]));
            VanillaOption option(payoff, exercise);
            
            option.setPricingEngine(boost::shared_ptr<PricingEngine>(
                new FdHestonHullWhiteVanillaEngine(
                     boost::shared_ptr<HestonModel>(
                             new HestonModel(hestonProcess)),
                                       hwProcess, corr[i], 50, 200, 10, 15)));    
            const Real calculated = option.NPV();
            
            option.setPricingEngine(boost::shared_ptr<PricingEngine>(
                new AnalyticBSMHullWhiteEngine(corr[i],
                                               stochProcess, hwModel)));
            const Real expected = option.NPV();

            const Real tol = 0.01;
            if (std::fabs(calculated - expected) > tol) {
                 BOOST_ERROR("Failed to reproduce analytic values"
                         << "\n   corr:       " << corr[i]
                         << "\n   strike:     " << strike[j]
                         << "\n   calculated: " << calculated
                         << "\n   expected:   " << expected);
            }
        }
    }
}

namespace {
    class HestonHullWhiteCorrelationConstraint : public Constraint {
      private:
        class Impl : public Constraint::Impl {
          public:
            Impl(Real equityShortRateCorr)
            : equityShortRateCorr_(equityShortRateCorr) {}

            bool test(const Array& params) const {
                const Real rho = params[3];

                return (  square<Real>()(rho)
                        + square<Real>()(equityShortRateCorr_) <= 1.0);
            }
          private:
            const Real equityShortRateCorr_;
        };
      public:
        HestonHullWhiteCorrelationConstraint(Real equityShortRateCorr)
        : Constraint(boost::shared_ptr<Constraint::Impl>(
             new HestonHullWhiteCorrelationConstraint::Impl(
                                                     equityShortRateCorr))) {}
    };
}

void HybridHestonHullWhiteProcessTest::testHestonHullWhiteCalibration() {
    BOOST_MESSAGE("Testing the Heston Hull-White calibration...");

    // Calibration of a hybrid Heston-Hull-White model using 
    // the finite difference HestonHullWhite pricing engine
    //
    // Imput surface is based on a Heston-Hull-White model with
    // Hull-White: a = 0.00883, \sigma = 0.00631
    // Heston    : \nu = 0.12, \kappa = 2.0, 
    //             \theta = 0.09, \sigma = 0.5, \rho=-0.75
    
    SavedSettings backup;

    const DayCounter dc = Actual365Fixed();
    const Calendar calendar = TARGET();
    const Date today = Date(28, March, 2004);
    Settings::instance().evaluationDate() = today;

    const Handle<YieldTermStructure> rTS(flatRate(0.05, dc));

    // assuming, that the Hull-White process is already calibrated 
    // on a given set of pure interest rate calibration instruments.
    boost::shared_ptr<HullWhiteProcess> hwProcess(
                              new HullWhiteProcess(rTS, 0.00883, 0.00631));
    boost::shared_ptr<HullWhite> hullWhiteModel(
                    new HullWhite(rTS, hwProcess->a(), hwProcess->sigma()));
                    
    const Handle<YieldTermStructure> qTS(flatRate(0.02, dc));
    Handle<Quote> s0(boost::shared_ptr<Quote>(new SimpleQuote(100.0)));

    // starting point the the pure Heston calibration
    const Real start_v0    = 0.2*0.2;
    const Real start_theta = start_v0;
    const Real start_kappa = 0.5;
    const Real start_sigma = 0.25;
    const Real start_rho   = -0.5;
    
    boost::shared_ptr<HestonProcess> hestonProcess(
        new HestonProcess(rTS, qTS, s0, start_v0, start_kappa, 
                          start_theta, start_sigma, start_rho));
    boost::shared_ptr<HestonModel> analyticHestonModel
                                            (new HestonModel(hestonProcess));
    boost::shared_ptr<PricingEngine> analyticHestonEngine(
                         new AnalyticHestonEngine(analyticHestonModel, 164));
    boost::shared_ptr<HestonModel> fdmHestonModel
                                            (new HestonModel(hestonProcess));
  
    const Real equityShortRateCorr = -0.5;
    
    const Real strikes[]    = { 50, 75, 90, 100, 110, 125, 150, 200 };
    const Time maturities[] = { 1/12., 3/12., 0.5, 1.0, 2.0, 3.0, 5.0, 7.5, 10};
    
    const Volatility vol[] = {
        0.482627,0.407617,0.366682,0.340110,0.314266,0.280241,0.252471,0.325552,
        0.464811,0.393336,0.354664,0.329758,0.305668,0.273563,0.244024,0.244886,
        0.441864,0.375618,0.340464,0.318249,0.297127,0.268839,0.237972,0.225553,
        0.407506,0.351125,0.322571,0.305173,0.289034,0.267361,0.239315,0.213761,
        0.366761,0.326166,0.306764,0.295279,0.284765,0.270592,0.250702,0.222928,
        0.345671,0.314748,0.300259,0.291744,0.283971,0.273475,0.258503,0.235683,
        0.324512,0.303631,0.293981,0.288338,0.283193,0.276248,0.266271,0.250506,
        0.311278,0.296340,0.289481,0.285482,0.281840,0.276924,0.269856,0.258609,
        0.303219,0.291534,0.286187,0.283073,0.280239,0.276414,0.270926,0.262173
    };
    
    std::vector<boost::shared_ptr<CalibrationHelper> > options;
    
    for (Size i=0; i < LENGTH(maturities); ++i) {
        const Period maturity((int)(maturities[i]*12.0+0.5), Months);
        boost::shared_ptr<Exercise> exercise(
                                        new EuropeanExercise(today + maturity));

        for (Size j=0; j < LENGTH(strikes); ++j) {
            boost::shared_ptr<StrikedTypePayoff> payoff(
                             new PlainVanillaPayoff(Option::Call, strikes[j]));
            RelinkableHandle<Quote> v(boost::shared_ptr<Quote>(
                                   new SimpleQuote(vol[i*LENGTH(strikes)+j])));
            options.push_back(boost::shared_ptr<CalibrationHelper>(
                new HestonModelHelper(maturity, calendar,s0->value(), 
                                      strikes[j], v, rTS, qTS,
                                      CalibrationHelper::PriceError)));
            const Real marketValue = options.back()->marketValue();
            
            // Improve the quality of the starting point 
            // for the full Heston-Hull-White calibration 
            boost::shared_ptr<SimpleQuote> volQuote(new SimpleQuote);
            boost::shared_ptr<GeneralizedBlackScholesProcess> bsProcess =
                detail::ImpliedVolatilityHelper::clone(
                    boost::shared_ptr<GeneralizedBlackScholesProcess>(
                        new GeneralizedBlackScholesProcess(
                            s0, qTS, rTS, Handle<BlackVolTermStructure>(
                                                    flatVol(v->value(), dc)))),
                        volQuote);
            
            VanillaOption dummyOption(payoff, exercise);

            boost::shared_ptr<PricingEngine> bshwEngine(
                new AnalyticBSMHullWhiteEngine(equityShortRateCorr, 
                                               bsProcess, hullWhiteModel));
            
            Volatility vt = detail::ImpliedVolatilityHelper::calculate(
                dummyOption, *bshwEngine, *volQuote, 
                marketValue, 1e-8, 100, 0.0001, 10);
            
            v.linkTo(boost::shared_ptr<Quote>(new SimpleQuote(vt)));
            
            options.back()->setPricingEngine(
                boost::shared_ptr<PricingEngine>(analyticHestonEngine));
        }
    }    
        
    HestonHullWhiteCorrelationConstraint corrConstraint(equityShortRateCorr);
    LevenbergMarquardt om(1e-6, 1e-8, 1e-8);
    analyticHestonModel->calibrate(options, om, 
                                   EndCriteria(400, 40, 1.0e-8, 1.0e-4, 1.0e-8),
                                   corrConstraint);

    options.clear();
    fdmHestonModel->setParams(analyticHestonModel->params());

    for (Size i=0; i < LENGTH(maturities); ++i) {
        const Size tGrid = static_cast<Size>(std::max(10.0, maturities[i]*10.0));
        boost::shared_ptr<FdHestonHullWhiteVanillaEngine> engine(
            new FdHestonHullWhiteVanillaEngine(fdmHestonModel, hwProcess, 
                                               equityShortRateCorr, 
                                               tGrid, 61, 13, 9, true));
        
        engine->enableMultipleStrikesCaching(
                     std::vector<Real>(strikes, strikes + LENGTH(strikes)));
        
        const Period maturity((int)(maturities[i]*12.0+0.5), Months);
        
        for (Size j=0; j < LENGTH(strikes); ++j) {
            // multiple strikes engine works best if the first option
            // per maturity has the average strike (because the first option
            // is priced first during the calibration and the first pricing
            // is used to calculate the prices for all strikes
            const Size js = (j + (LENGTH(strikes)-1)/2) % LENGTH(strikes);
 
            boost::shared_ptr<StrikedTypePayoff> payoff(
                             new PlainVanillaPayoff(Option::Call, strikes[js]));
            Handle<Quote> v(boost::shared_ptr<Quote>(
                                   new SimpleQuote(vol[i*LENGTH(strikes)+js])));
            options.push_back(boost::shared_ptr<CalibrationHelper>(
                new HestonModelHelper(maturity, calendar, s0->value(), 
                                      strikes[js], v, rTS, qTS,
                                      CalibrationHelper::PriceError)));
            
            options.back()->setPricingEngine(engine);
        }
    }    

    LevenbergMarquardt vm(1e-6, 1e-2, 1e-2);
    fdmHestonModel->calibrate(options, vm, 
                              EndCriteria(400, 40, 1.0e-8, 1.0e-4, 1.0e-8),
                              corrConstraint);
    
    const Real relTol = 0.01;
    const Real expected_v0    =  0.12;
    const Real expected_kappa =  2.0;
    const Real expected_theta =  0.09;
    const Real expected_sigma =  0.5;
    const Real expected_rho   = -0.75;
    
    if (std::fabs(fdmHestonModel->v0() - expected_v0)/expected_v0 > relTol) {
         BOOST_ERROR("Failed to reproduce Heston-Hull-White model"
                 << "\n   v0 calculated: " << fdmHestonModel->v0()
                 << "\n   v0 expected  : " << expected_v0
                 << "\n   relatove tol : " << relTol);
    }
    if (std::fabs(fdmHestonModel->theta() - expected_theta)/expected_theta 
                                                                    > relTol) {
         BOOST_ERROR("Failed to reproduce Heston-Hull-White model"
                 << "\n   theta calculated: " << fdmHestonModel->theta()
                 << "\n   theta expected  : " << expected_theta
                 << "\n   relatove tol    : " << relTol);
    }
    if (std::fabs(fdmHestonModel->kappa() - expected_kappa)/expected_kappa 
                                                                    > relTol) {
        BOOST_ERROR("Failed to reproduce Heston-Hull-White model"
                << "\n   kappa calculated: " << fdmHestonModel->kappa()
                << "\n   kappa expected  : " << expected_kappa
                << "\n   relatove tol    : " << relTol);
    }
    if (std::fabs(fdmHestonModel->sigma() - expected_sigma)/expected_sigma 
                                                                    > relTol) {
       BOOST_ERROR("Failed to reproduce Heston-Hull-White model"
               << "\n   sigma calculated: " << fdmHestonModel->sigma()
               << "\n   sigma expected  : " << expected_sigma
               << "\n   relatove tol    : " << relTol);
    }
    if (std::fabs(fdmHestonModel->rho() - expected_rho)/expected_rho > relTol) {
         BOOST_ERROR("Failed to reproduce Heston-Hull-White model"
                 << "\n   rho calculated: " << fdmHestonModel->rho()
                 << "\n   rho expected  : " << expected_rho
                 << "\n   relatove tol  : " << relTol);
    }
}
    
    
test_suite* HybridHestonHullWhiteProcessTest::suite() {
    test_suite* suite = BOOST_TEST_SUITE("Hybrid Heston-HullWhite tests");

    // FLOATING_POINT_EXCEPTION
    suite->add(QUANTLIB_TEST_CASE(
        &HybridHestonHullWhiteProcessTest::testBsmHullWhiteEngine));
    // FLOATING_POINT_EXCEPTION
    suite->add(QUANTLIB_TEST_CASE(
        &HybridHestonHullWhiteProcessTest::testCompareBsmHWandHestonHW));
    suite->add(QUANTLIB_TEST_CASE(
        &HybridHestonHullWhiteProcessTest::testZeroBondPricing));
    // FLOATING_POINT_EXCEPTION
    suite->add(QUANTLIB_TEST_CASE(
        &HybridHestonHullWhiteProcessTest::testMcVanillaPricing));
    // FLOATING_POINT_EXCEPTION
    suite->add(QUANTLIB_TEST_CASE(
        &HybridHestonHullWhiteProcessTest::testMcPureHestonPricing));
    // FLOATING_POINT_EXCEPTION
    suite->add(QUANTLIB_TEST_CASE(
      &HybridHestonHullWhiteProcessTest::testAnalyticHestonHullWhitePricing));
    suite->add(QUANTLIB_TEST_CASE(
        &HybridHestonHullWhiteProcessTest::testCallableEquityPricing));
    suite->add(QUANTLIB_TEST_CASE(
        &HybridHestonHullWhiteProcessTest::testDiscretizationError));
    suite->add(QUANTLIB_TEST_CASE(
        &HybridHestonHullWhiteProcessTest::testFdmHestonHullWhiteEngine));
    suite->add(QUANTLIB_TEST_CASE(
        &HybridHestonHullWhiteProcessTest::testHestonHullWhiteCalibration));
    
    return suite;
}
