/*
 * AmericanOptionTest.hpp
 *
 *  Created on: Jan 31, 2010
 *      Author: inzite
 */

#ifndef AMERICANOPTIONTEST_HPP_
#define AMERICANOPTIONTEST_HPP_

#include <boost/test/unit_test.hpp>

class AmericanOptionTest {
public:
	AmericanOptionTest();
	virtual ~AmericanOptionTest();
	static void testLongstaffSchwartzEngine();
	static boost::unit_test_framework::test_suite* suite();
};

#endif /* AMERICANOPTIONTEST_HPP_ */
