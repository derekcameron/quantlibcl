dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(QuantLibCL-test, 1.0)
AC_CANONICAL_SYSTEM

AC_CONFIG_SRCDIR([src/QuantLibCL-test.cpp])

AC_CONFIG_HEADERS([src/config.hpp])

AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_MACRO_DIR([m4])

# Help the compiler find external software

AC_ARG_WITH([boost-include],
            AC_HELP_STRING([--with-boost-include=INCLUDE_PATH],
                           [Supply the location of Boost header files]),
            [ql_boost_include_path="`cd ${withval} 2>/dev/null && pwd`"],
            [ql_boost_include_path=""])
if test [ -n "$ql_boost_include_path" ] ; then
   AC_SUBST([BOOST_INCLUDE],["-I${ql_boost_include_path}"])
   AC_SUBST([CPPFLAGS],["${CPPFLAGS} -I${ql_boost_include_path}"])
fi
AC_ARG_WITH([boost-lib],
            AC_HELP_STRING([--with-boost-lib=LIB_PATH],
                           [Supply the location of Boost libraries]),
            [ql_boost_lib_path="`cd ${withval} 2>/dev/null && pwd`"],
            [ql_boost_lib_path=""])
if test [ -n "$ql_boost_lib_path" ] ; then
   AC_SUBST([BOOST_LIB],["-L${ql_boost_lib_path}"])
   AC_SUBST([LDFLAGS],["${LDFLAGS} -L${ql_boost_lib_path}"])
fi
AC_ARG_WITH([quantlib-include],
            AC_HELP_STRING([--with-quantlib-include=INCLUDE_PATH],
                           [Supply the location of QuantLib header files]),
            [quantlib_include_path="`cd ${withval} 2>/dev/null && pwd`"],
            [quantlib_include_path=""])
if test [ -n "$quantlib_include_path" ] ; then
   AC_SUBST([QUANTLIB_INCLUDE],["-I${quantlib_include_path}"])
   AC_SUBST([CPPFLAGS],["${CPPFLAGS} -I${quantlib_include_path}"])
fi
AC_ARG_WITH([quantlib-lib],
            AC_HELP_STRING([--with-quantlib-lib=LIB_PATH],
                           [Supply the location of QuantLib libraries]),
            [quantlib_lib_path="`cd ${withval} 2>/dev/null && pwd`"],
            [quantlib_lib_path=""])
if test [ -n "$quantlib_lib_path" ] ; then
   AC_SUBST([QUANTLIB_LIB],["-L${quantlib_lib_path}"])
   AC_SUBST([LDFLAGS],["${LDFLAGS} -L${quantlib_lib_path}"])
fi

AC_PROG_CC
AC_PROG_CXX
AC_PROG_CPP
AC_PROG_LIBTOOL

AM_CONDITIONAL(ENABLE_STATIC, test "${enable_static}" == "yes")

AC_LANG([C++])

# Check for Boost components
QLTEST_CHECK_BOOST
AM_CONDITIONAL(BOOST_UNIT_TEST_FOUND, test "x${BOOST_UNIT_TEST_LIB}" != "x")

# Check for QuantLib
QLTEST_CHECK_QUANTLIB

AC_CONFIG_FILES(Makefile src/Makefile)
AC_OUTPUT

