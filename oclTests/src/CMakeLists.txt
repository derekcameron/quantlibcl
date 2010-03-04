SET(OCL_UTILS_INCLUDE_DIR "../../oclUtils/inc/")
SET(OCL_UTILS_LINK_DIR "../../oclUtils/lib/")

INCLUDE_DIRECTORIES(${OCL_UTILS_INCLUDE_DIR})
LINK_DIRECTORIES(${OCL_UTILS_LINK_DIR})

SET(EXTRA_LIBS ${EXTRA_LIBS} "oclUtils" "OpenCL")

ADD_EXECUTABLE (test1 test1.cpp)
TARGET_LINK_LIBRARIES (test1 ${EXTRA_LIBS})

SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

# make sure we have Boost
find_package (Boost)
INCLUDE_DIRECTORIES ( ${Boost_INCLUDE_DIRS} )
#this line shouldn't be necessary, because we don't (yet) use any linkable parts of Boost
LINK_DIRECTORIES ( ${Boost_LIBRARY_DIRS} )