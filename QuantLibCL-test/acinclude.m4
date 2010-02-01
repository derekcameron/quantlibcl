
# QLTEST_CHECK_BOOST_DEVEL
# --------------------
# Check whether the Boost headers are available
AC_DEFUN([QLTEST_CHECK_BOOST_DEVEL],
[AC_MSG_CHECKING([for Boost development files])
 AC_TRY_COMPILE(
    [@%:@include <boost/version.hpp>
     @%:@include <boost/shared_ptr.hpp>
     @%:@include <boost/assert.hpp>
     @%:@include <boost/current_function.hpp>],
    [],
    [AC_MSG_RESULT([yes])],
    [AC_MSG_RESULT([no])
     AC_MSG_ERROR([Boost development files not found])
    ])
])

# QLTEST_CHECK_BOOST_LIB
# --------------------
# Check whether the Boost library is available
#AC_DEFUN([QLTEST_CHECK_BOOST_DEVEL],
#[AC_MSG_CHECKING([for Boost library file])
# AC_CHECK_LIB(boost_lib, avl_create, [], [
#        echo "Error! You need to have the boost libraries."
#        exit -1
#        ])
#])

# QLTEST_CHECK_BOOST_VERSION
# ----------------------
# Check whether the Boost installation is up to date
AC_DEFUN([QLTEST_CHECK_BOOST_VERSION],
[AC_MSG_CHECKING([Boost version])
 AC_REQUIRE([QLTEST_CHECK_BOOST_DEVEL])
 AC_TRY_COMPILE(
    [@%:@include <boost/version.hpp>],
    [@%:@if BOOST_VERSION < 103100
     @%:@error too old
     @%:@endif],
    [AC_MSG_RESULT([yes])],
    [AC_MSG_RESULT([no])
     AC_MSG_ERROR([outdated Boost installation])
    ])
])

# QLTEST_CHECK_BOOST_UBLAS
# --------------------
# Check whether the Boost headers are available
AC_DEFUN([QLTEST_CHECK_BOOST_UBLAS],
[AC_MSG_CHECKING([for Boost::uBLAS support])
 AC_TRY_COMPILE(
    [@%:@include <boost/numeric/ublas/vector_proxy.hpp>
     @%:@include <boost/numeric/ublas/triangular.hpp>
     @%:@include <boost/numeric/ublas/lu.hpp>],
    [],
    [AC_MSG_RESULT([yes])],
    [AC_MSG_RESULT([no])
     AC_MSG_WARN([Some functionality will be disabled.])
     AC_DEFINE([QLTEST_NO_UBLAS_SUPPORT],[],
               [Define this if your compiler does not support Boost::uBLAS.])
    ])
])

# QLTEST_CHECK_BOOST_UNIT_TEST
# ------------------------
# Check whether the Boost unit-test framework is available
AC_DEFUN([QLTEST_CHECK_BOOST_UNIT_TEST],
[AC_MSG_CHECKING([for Boost unit-test framework])
 AC_REQUIRE([AC_PROG_CC])
 ql_original_LIBS=$LIBS
 ql_original_CXXFLAGS=$CXXFLAGS
 CC_BASENAME=`basename $CC`
 CC_VERSION=`echo "__GNUC__ __GNUC_MINOR__" | $CC -E -x c - | tail -n 1 | $SED -e "s/ //"`
 for boost_lib in boost_unit_test_framework-$CC_BASENAME$CC_VERSION \
                  boost_unit_test_framework-$CC_BASENAME \
                  boost_unit_test_framework \
                  boost_unit_test_framework-mt-$CC_BASENAME$CC_VERSION \
                  boost_unit_test_framework-$CC_BASENAME$CC_VERSION-mt \
                  boost_unit_test_framework-x$CC_BASENAME$CC_VERSION-mt \
                  boost_unit_test_framework-mt-$CC_BASENAME \
                  boost_unit_test_framework-$CC_BASENAME-mt \
                  boost_unit_test_framework-mt ; do
     LIBS="$ql_original_LIBS -l$boost_lib"
     # 1.33.1 or 1.34 static
     CXXFLAGS="$ql_original_CXXFLAGS"
     boost_unit_found=no
     AC_LINK_IFELSE(
         [@%:@include <boost/test/unit_test.hpp>
          using namespace boost::unit_test_framework;
          test_suite*
          init_unit_test_suite(int argc, char** argv)
          {
              return (test_suite*) 0;
          }
         ],
         [boost_unit_found=$boost_lib
          boost_defines=""
          break],
         [])
     # 1.34 shared
     CXXFLAGS="$ql_original_CXXFLAGS -DBOOST_TEST_MAIN -DBOOST_TEST_DYN_LINK"
     boost_unit_found=no
     AC_LINK_IFELSE(
         [@%:@include <boost/test/unit_test.hpp>
          using namespace boost::unit_test_framework;
          test_suite*
          init_unit_test_suite(int argc, char** argv)
          {
              return (test_suite*) 0;
          }
         ],
         [boost_unit_found=$boost_lib
          boost_defines="-DBOOST_TEST_DYN_LINK"
          break],
         [])
 done
 LIBS="$ql_original_LIBS"
 CXXFLAGS="$ql_original_CXXFLAGS"
 if test "$boost_unit_found" = no ; then
     AC_MSG_RESULT([no])
     AC_SUBST([BOOST_UNIT_TEST_LIB],[""])
     AC_SUBST([BOOST_UNIT_TEST_MAIN_CXXFLAGS],[""])
     AC_MSG_WARN([Boost unit-test framework not found.])
     AC_MSG_WARN([The test suite will be disabled.])
 else
     AC_MSG_RESULT([yes])
     AC_SUBST([BOOST_UNIT_TEST_LIB],[$boost_lib])
     AC_SUBST([BOOST_UNIT_TEST_MAIN_CXXFLAGS],[$boost_defines])
 fi
])

# QLTEST_CHECK_BOOST_TEST_STREAM
# --------------------------
# Check whether Boost unit-test stream accepts std::fixed
AC_DEFUN([QLTEST_CHECK_BOOST_TEST_STREAM],
[AC_MSG_CHECKING([whether Boost unit-test streams work])
 AC_REQUIRE([AC_PROG_CC])
 AC_TRY_COMPILE(
    [@%:@include <boost/test/unit_test.hpp>
     @%:@include <iomanip>],
    [BOOST_ERROR("foo " << std::fixed << 42.0);],
    [AC_MSG_RESULT([yes])
     AC_SUBST(BOOST_UNIT_TEST_DEFINE,[-DQLTEST_WORKING_BOOST_STREAMS])],
    [AC_MSG_RESULT([no])
     AC_SUBST(BOOST_UNIT_TEST_DEFINE,[""])
    ])
])

# QLTEST_CHECK_BOOST
# ------------------------
# Boost-related tests
AC_DEFUN([QLTEST_CHECK_BOOST],
[AC_REQUIRE([QLTEST_CHECK_BOOST_DEVEL])
 AC_REQUIRE([QLTEST_CHECK_BOOST_VERSION])
 AC_REQUIRE([QLTEST_CHECK_BOOST_UBLAS])
 AC_REQUIRE([QLTEST_CHECK_BOOST_UNIT_TEST])
 AC_REQUIRE([QLTEST_CHECK_BOOST_TEST_STREAM])
])

# QLTEST_CHECK_QUANTLIB_DEVEL
# --------------------
# Check whether the QuantLib headers are available
AC_DEFUN([QLTEST_CHECK_QUANTLIB_DEVEL],
[AC_MSG_CHECKING([for QuantLib header files])
 AC_TRY_COMPILE(
    [@%:@include <ql/quantlib.hpp>],
    [],
    [AC_MSG_RESULT([yes])],
    [AC_MSG_RESULT([no])
     AC_MSG_ERROR([QuantLib header files not found])
    ])
])

# QLTEST_CHECK_QUANTLIB_LIB
# --------------------
# Check whether the QuantLib library is available
#AC_DEFUN([QLTEST_CHECK_QUANTLIB_DEVEL],
#[AC_MSG_CHECKING([for QuantLib header files])
# AC_TRY_COMPILE(
#    [@%:@include <ql/quantlib.hpp>],
#    [],
#    [AC_MSG_RESULT([yes])],
#    [AC_MSG_RESULT([no])
#     AC_MSG_ERROR([QuantLib header files not found])
#    ])
#])


# QLTEST_CHECK_QUANTLIB
# ------------------------
# QuantLib-related tests
AC_DEFUN([QLTEST_CHECK_QUANTLIB],
[AC_REQUIRE([QLTEST_CHECK_QUANTLIB_DEVEL])
# AC_REQUIRE([QLTEST_CHECK_QUANTLIB_LIB])
])