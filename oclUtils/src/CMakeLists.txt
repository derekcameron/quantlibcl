add_library(oclUtils oclUtils.cpp)

INCLUDE_DIRECTORIES( "../inc" )

SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

# make sure we have Boost
FIND_PACKAGE(Boost)

INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIRS} )
#this line shouldn't be necessary, because we don't (yet) use any linkable parts of Boost
LINK_DIRECTORIES( ${Boost_LIBRARY_DIRS} )