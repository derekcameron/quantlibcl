bin_PROGRAMS = test_suite
test_suite_SOURCES = \
	AmericanOptionTest.cpp \
	QuantLibCL-test.cpp
	
test_suite_LDFLAGS = -static -L/home/inzite/workspace/QuantLibCL/ql/.libs -lQuantLib

AM_CPPFLAGS = ${BOOST_UNIT_TEST_DEFINE}
AM_LDFLAGS = -l${BOOST_UNIT_TEST_LIB}