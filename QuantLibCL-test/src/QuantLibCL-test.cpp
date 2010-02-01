/*
 ============================================================================
 Name        : QuantLibCL-test.cpp
 Author      : William Gross
 Version     :
 Copyright   : 
 Description : QuantLibCL test suite
 ============================================================================
 */
#include <boost/test/included/unit_test.hpp>

#include "AmericanOptionTest.hpp"

using namespace boost::unit_test;

test_suite* init_unit_test_suite( int argc, char* argv[] )
{
//    test_suite* ts1 = &AmericanOptionTest::suite();

//    ts1->add( BOOST_TEST_CASE( &test_case1 ) );

	framework::master_test_suite().add( AmericanOptionTest::suite() );

    /*    ts1->add( BOOST_TEST_CASE( &test_case2 ) );

    test_suite* ts2 = BOOST_TEST_SUITE( "test_suite2" );
    ts2->add( BOOST_TEST_CASE( &test_case3 ) );
    ts2->add( BOOST_TEST_CASE( &test_case4 ) );

    framework::master_test_suite().add( ts2 );*/

    return 0;
}
