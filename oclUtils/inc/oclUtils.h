/*
 * oclUtils.h
 *
 *  Created on: Feb 28, 2010
 *      Author: inzite
 */

#ifndef OCLUTILS_H_
#define OCLUTILS_H_

#include <boost/exception/all.hpp>
#include <boost/shared_ptr.hpp>
#include <CL/cl.h>
//#include <CL/cl_ext.h>

typedef boost::error_info<struct tag_ocl_exception,std::string> oclExceptionInfo;
struct oclException : virtual boost::exception, virtual std::exception { };

cl_int oclGetPlatformID(cl_platform_id* clSelectedPlatformID);

inline void oclCheckError(cl_int iSample, cl_int iReference)
{
    if (iReference != iSample)
    	//throw oclException();
    	throw oclException() << oclExceptionInfo("Error in oclCheckError");
}

#endif /* OCLUTILS_H_ */
