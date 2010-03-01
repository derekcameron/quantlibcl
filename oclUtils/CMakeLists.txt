cmake_minimum_required (VERSION 2.6)
project (OpenCL_Utilities) 

# add our subdirectories
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(inc)

SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)