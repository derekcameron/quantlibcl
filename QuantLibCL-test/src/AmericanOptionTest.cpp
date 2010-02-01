/*
 * AmericanOptionTest.cpp
 *
 *  Created on: Jan 31, 2010
 *      Author: William Gross
 */

#include "AmericanOptionTest.hpp"
#include "utilities.hpp"

#include <ql/pricingengines/vanilla/all.hpp>
#include <ql/time/calendars/target.hpp>

#include <boost/shared_ptr.hpp>

#include <iostream>

using namespace QuantLib;

AmericanOptionTest::AmericanOptionTest() {
	// TODO Auto-generated constructor stub

}

AmericanOptionTest::~AmericanOptionTest() {
	// TODO Auto-generated destructor stub
}

void AmericanOptionTest::testLongstaffSchwartzEngine()
{
    //Define our option type (let's make it a call option)
	Option::Type type(Option::Call);

	//Define our Monte Carlo seed value (of course it would be 42)
    Size mcSeed = 42;

	//Define the value of the underlying
	Real underlying = 360;

	//Define the option strike price
    Real strike = 40;

	//Define the option's settlement date
	Date settlementDate(17, May, 1998);

	//Define the option's maturity
    Date maturity(17, May, 1999);

	//Define the dividend yield
	Spread dividendYield = 0.00;

	//Define our risk free rate
    Rate riskFreeRate = 0.06;

	//Define our calendar
    Calendar calendar = TARGET();

    //Define our day counter
	DayCounter dayCounter = Actual365Fixed();

    //Define our expected volatility
	Volatility volatility = 0.20;

	//Define the characteristics of our option's exercise
    boost::shared_ptr<Exercise> americanExercise(
                                     new AmericanExercise(settlementDate,
                                                          maturity));

	//Define a handle to the underlying
    Handle<Quote> underlyingH(
	            boost::shared_ptr<Quote>(new SimpleQuote(underlying)));

    //Define a handle to our flat term structure
    Handle<YieldTermStructure> flatTermStructure(
        boost::shared_ptr<YieldTermStructure>(
            new FlatForward(settlementDate, riskFreeRate, dayCounter)));


    Handle<YieldTermStructure> flatDividendTS(
    		boost::shared_ptr<YieldTermStructure>(
    				new FlatForward(settlementDate, dividendYield, dayCounter)));

    Handle<BlackVolTermStructure> flatVolTS(
    		boost::shared_ptr<BlackVolTermStructure>(
    				new BlackConstantVol(settlementDate, calendar, volatility, dayCounter)));

    boost::shared_ptr<StrikedTypePayoff> payoff(
                                    new PlainVanillaPayoff(type, strike));

    boost::shared_ptr<BlackScholesMertonProcess> bsmProcess(
			new BlackScholesMertonProcess(underlyingH, flatDividendTS,
					flatTermStructure, flatVolTS));

    //Define our option
    VanillaOption americanOption(payoff, americanExercise);

	std::string method = "MC (Longstaff Schwartz)";
	boost::shared_ptr<PricingEngine> mcengine;
	mcengine = MakeMCAmericanEngine<PseudoRandom>(bsmProcess)
			.withSteps(100)
			.withAntitheticVariate()
			.withCalibrationSamples(4096)
			.withAbsoluteTolerance(0.02)
			.withSeed(mcSeed);

    americanOption.setPricingEngine(mcengine);

	BOOST_TEST_MESSAGE(
			"Option NPV (QuantLib method) = "<< americanOption.NPV()
			);
}

boost::unit_test_framework::test_suite* AmericanOptionTest::suite() {
    boost::unit_test_framework::test_suite* suite = BOOST_TEST_SUITE("American option tests");

    suite->add( BOOST_TEST_CASE(&AmericanOptionTest::testLongstaffSchwartzEngine) );
    return suite;
}
