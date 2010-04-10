/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2005, 2006 StatPro Italia srl
 Copyright (C) 2010 William Gross

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

/*! \file opencl.hpp
    \brief utilities for working with the underlying OpenCL libraries and hardware
*/

#ifndef quantlib_opencl_hpp
#define quantlib_opencl_hpp

//Only include these functions if we're actually using OpenCL
#ifndef QUANTLIB_DISABLE_OPENCL

#define __CL_ENABLE_EXCEPTIONS // Enable OpenCL exceptions
#define __NO_STD_VECTOR // Use cl::vector and cl::string and 
#define __NO_STD_STRING  // not the STL versions

#if defined(__APPLE__) || defined(__MACOSX)
#include <OpenCL/cl.hpp>
#else
#include <CL/cl.hpp>
#endif


namespace QuantLib {

	class OclDevice {
	public:
		// constructor
		OclDevice(
			// any required components should be passed in here
			const cl::Context context,
			const cl::vector<cl::Device> devices,
			const cl::Device device,
			const cl::CommandQueue commandQueue
		);
		//! Load and compile OpenCL sources on this device
		unsigned int loadSources(const cl::Program::Sources &sources);
		unsigned int loadKernel(const size_t programIndex, const char* kernelName);
		unsigned int allocateBuffer(void* ptr, const size_t bufferSize);
		unsigned int launchKernel(const unsigned int kernelHandle, const unsigned int numberOfThreads, const unsigned int localWorkSize = 128);
		void wait(const unsigned int eventHandle);
		void readBuffer(const unsigned int bufferHandle, void* dest);
		cl::Buffer buffer(unsigned int i);
		template <class T> void setKernelArg(unsigned int kernelHandle, cl_uint argIndex, T arg);
	private:
		cl::Context context_;
		cl::vector<cl::Device> devices_;
		cl::Device device_;
		cl::CommandQueue commandQueue_;
		cl::vector<cl::Program> programs_;
		cl::vector<cl::Buffer> buffers_;
		cl::vector<size_t> bufferSizes_;
		cl::vector<cl::Kernel> kernels_;
		cl::vector<cl::Event> events_;
	};

	// constructor
	inline OclDevice::OclDevice(
		const cl::Context context,
		const cl::vector<cl::Device> devices,
		const cl::Device device,
		const cl::CommandQueue commandQueue
		) {
		context_ = context;
		devices_ = devices;
		device_ = device;
		commandQueue_ = commandQueue;
	}

	inline unsigned int OclDevice::loadSources(const cl::Program::Sources &sources) {
		programs_.push_back(cl::Program(context_, sources));
		programs_.back().build(devices_,"");
		return programs_.size() - 1;
	}

	inline unsigned int OclDevice::loadKernel(const size_t programIndex, const char* kernelName) {
		// Load the kernel
		kernels_.push_back(cl::Kernel(programs_[programIndex], kernelName));
		return kernels_.size() - 1;
	}

	inline unsigned int OclDevice::allocateBuffer(void* ptr, const size_t bufferSize) {
		buffers_.push_back(cl::Buffer(context_, CL_MEM_USE_HOST_PTR, bufferSize, ptr));
		bufferSizes_.push_back(bufferSize);
		return buffers_.size() - 1;
	}

	inline unsigned int OclDevice::launchKernel(const unsigned int kernelHandle, const unsigned int numberOfThreads, const unsigned int localWorkSize) {
		events_.push_back(cl::Event());
		commandQueue_.enqueueNDRangeKernel(
				kernels_[kernelHandle],
				cl::NullRange,
				cl::NDRange(numberOfThreads),
				cl::NDRange(localWorkSize),
				NULL,
				&events_.back()
			);

		return events_.size() - 1;
	}

	inline void OclDevice::wait(const unsigned int eventHandle) {
		events_[eventHandle].wait();
	}

	inline void OclDevice::readBuffer(const unsigned int bufferHandle, void* dest) {
		commandQueue_.enqueueReadBuffer(
				buffers_[bufferHandle],
				CL_TRUE,
				0,
				bufferSizes_[bufferHandle],
				dest
				);
	}

	inline cl::Buffer OclDevice::buffer(unsigned int i) {
		return buffers_[i];
	}

	template <class T> inline void OclDevice::setKernelArg(unsigned int kernelHandle, cl_uint argIndex, T arg) {
		kernels_[kernelHandle].setArg(argIndex, arg);
	}

	class MakeOclDevice {
	public:
        MakeOclDevice();
        // named parameters
        MakeOclDevice& withDeviceType(cl_int deviceType);
        //MakeOclDevice& withOtherParameter(Parameter someOtherParameter);
        // conversion to OclDevice
        operator boost::shared_ptr<OclDevice>() const;
	private:
		cl_int deviceType_;
    };

	// inline definitions

	inline MakeOclDevice::MakeOclDevice() {}

	inline MakeOclDevice::operator boost::shared_ptr<OclDevice>() const {
		// check for any missing requirements here
        /*QL_REQUIRE(steps_ == Null<Size>() || stepsPerYear_ == Null<Size>(),
                   "number of steps overspecified");*/
		
		//Find the available platforms and select one
		cl::vector<cl::Platform> allPlatforms;
		cl::Platform::get(&allPlatforms);
		cl::Platform targetPlatform = allPlatforms[0];

		//Get a context matching deviceType_
		cl_context_properties cprops[3] = { CL_CONTEXT_PLATFORM, (cl_context_properties)targetPlatform(), 0 };
		cl::Context context = cl::Context(deviceType_, cprops);

		//And a list of all devices matching that context
		cl::vector<cl::Device> devices = context.getInfo<CL_CONTEXT_DEVICES>();

		//Lastly, create a command queue for device 0
		cl::CommandQueue commandQueue = cl::CommandQueue(context,devices[0]);

		//And if all the creation goes according to plan, return a shared pointer to the created object
		return boost::shared_ptr<OclDevice>(
			new OclDevice(
				//Insert any required components here
				context,
				devices,
				devices[0],
				commandQueue
			)
		);
    }

	inline MakeOclDevice& MakeOclDevice::withDeviceType(cl_int deviceType) {
		deviceType_ = deviceType;
		return *this;
	}
}

//End of 'ifdef QUANTLIB_DISABLE_OPENCL'
#endif

#endif
