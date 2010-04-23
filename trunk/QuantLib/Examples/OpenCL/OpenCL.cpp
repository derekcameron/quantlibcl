/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*!
 Copyright (C) 2005, 2006, 2007, 2009 StatPro Italia srl
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

// the only header you need to use QuantLib
#include <ql/quantlib.hpp>

#ifdef BOOST_MSVC
/* Uncomment the following lines to unmask floating-point
   exceptions. Warning: unpredictable results can arise...

   See http://www.wilmott.com/messageview.cfm?catid=10&threadid=9481
   Is there anyone with a definitive word about this?
*/
// #include <float.h>
// namespace { unsigned int u = _controlfp(_EM_INEXACT, _MCW_EM); }
#endif

#include "dcmt.hpp"
#include "mersennetwister.hpp"
#include "NVIDIA_mersennetwister.hpp"
#include "mcsimulation.hpp"
#include <boost/timer.hpp>
#include <boost/cstdint.hpp>
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace QuantLib;
using std::ifstream;

#if defined(QL_ENABLE_SESSIONS)
namespace QuantLib {

    Integer sessionId() { return 0; }

}
#endif

#define OCL_THREAD_TOTAL 16384

// Test 1 - Launch parallel threads that write consecutive integers into an array
void test1(boost::shared_ptr<OclDevice> ocldevice1) {
	std::cout << "Running OpenCL test 1 with " << OCL_THREAD_TOTAL << " threads..." << std::endl;
	
	// Load and compile sources
	std::ifstream file1("kernels.cl");
	std::string kernels(std::istreambuf_iterator<char>(file1), (std::istreambuf_iterator<char>()));
	file1.close();
	cl::Program::Sources sources1(1, std::make_pair(kernels.c_str(), kernels.length()+1));
	size_t sourcesIndex = ocldevice1->loadSources(sources1);

   	boost::uint64_t outH[OCL_THREAD_TOTAL];
	for(int i = 0; i < OCL_THREAD_TOTAL; i++)
		outH[i] = 0;

	//Allocate buffers and launch threads on the device
	unsigned int bufferHandle = ocldevice1->allocateBuffer(&outH, sizeof(outH));
	unsigned int kernelHandle = ocldevice1->loadKernel(sourcesIndex,"oclTest1");
	ocldevice1->setKernelArg(kernelHandle, 0, *(ocldevice1->buffer(bufferHandle)));
	unsigned int eventHandle = ocldevice1->launchKernel(kernelHandle, OCL_THREAD_TOTAL);

	// Wait for OpenCL execution to complete
	ocldevice1->wait(eventHandle);

	// Read the result
	ocldevice1->readBuffer(bufferHandle, &outH);

	bool passed = 1;
	std::cout << "Test complete:  ";
	for(int i = 0; i < OCL_THREAD_TOTAL; i++)
	{
		if (i != outH[i])
		{
			passed = 0;
			std::cout << "failed" << std::endl << std::endl;
			break;
		}
	}
	if(passed)
		std::cout << "passed" << std::endl << std::endl;
	//End of test 1
}

// Test 2 - Generate (0,1) uniformly distributed random variables using the Thomas-Luk algorithm
void test2(boost::shared_ptr<OclDevice> ocldevice1) {
	std::ifstream file2("thomaslukrng.cl");
	std::string thomaslukrng(std::istreambuf_iterator<char>(file2), (std::istreambuf_iterator<char>()));
	file2.close();
	cl::Program::Sources sources2(1, std::make_pair(thomaslukrng.c_str(), thomaslukrng.length()+1));
	size_t sources2Index = ocldevice1->loadSources(sources2);
}

int test3FindParameters(mt_params*** mtpp, int desiredCount) {
	int count;
	*mtpp = get_mt_parameters_st(32,521,0,desiredCount-1,42, &count);
	return count;
}

// Test 3 - Generate (0,1) uniformly distributed random variables using the Mersenne Twister algorithm
void test3(boost::shared_ptr<OclDevice> ocldevice1) {
	mt_params **mtpp;
	int paramsCount = test3FindParameters(&mtpp, 16);

	for (int i = 0; i < paramsCount; i++) {
		std::cout << "Item " << i << std::endl;
		std::cout << "aaa = " << mtpp[i]->aaa << std::endl;
		std::cout << "maskB = " << mtpp[i]->maskB << std::endl;
		std::cout << "maskC = " << mtpp[i]->maskC << std::endl;
		std::cout << "seed = " << mtpp[i]->seed << std::endl << std::endl;
	}

	std::ifstream file3("mersennetwister.cl");
	std::string string3(std::istreambuf_iterator<char>(file3), (std::istreambuf_iterator<char>()));
	file3.close();
	cl::Program::Sources sources3(1, std::make_pair(string3.c_str(), string3.length()+1));
	size_t sources3Index = ocldevice1->loadSources(sources3);
}

void loadMersenneTwisterParams(const char *fname, 
	       const unsigned int seed, 
	       mt_params_stripped *mtpp,
	       const size_t size)
{
    FILE* fd = 0;
    #ifdef _WIN32
        // open the file for binary read
        errno_t err;
        if ((err = fopen_s(&fd, fname, "rb")) != 0)
    #else
        // open the file for binary read
        if ((fd = fopen(fname, "rb")) == 0)
    #endif
        {
            if(fd)
            {
                fclose (fd);
            }
			throw std::runtime_error("Error opening mt_params raw data file in test 4");
        }
  
    for (unsigned int i = 0; i < size; i++)
        fread(&mtpp[i], sizeof(mt_params_stripped), 1, fd);
    fclose(fd);

    for(unsigned int i = 0; i < size; i++)
        mtpp[i].seed = seed;
}

// Test 4 - Generate (0,1) uniformly distributed random variables using the GPL'd NVIDIA Mersenne Twister algorithm
void test4(boost::shared_ptr<OclDevice> ocldevice1) {
	
	//const parameters
	const int seed = 777;

	// Allocate space for the result
	boost::shared_array<float> h_RandGPU(new float[NVIDIA_nRand]);
	const size_t size_h_RandGPU = sizeof(float[NVIDIA_nRand]);
	// Allocate space for dynamic creation parameters
	boost::shared_array<mt_params_stripped> h_mtParams(new mt_params_stripped[NVIDIA_MT_RNG_COUNT]);
	const size_t size_h_mtParams = sizeof(mt_params_stripped[NVIDIA_MT_RNG_COUNT]);

	loadMersenneTwisterParams("data/MersenneTwister.dat", seed, h_mtParams.get(), NVIDIA_MT_RNG_COUNT);

	std::ifstream file4("NVIDIA_mersennetwister.cl");
	std::string string4(std::istreambuf_iterator<char>(file4), (std::istreambuf_iterator<char>()));
	file4.close();
	cl::Program::Sources sources4(1, std::make_pair(string4.c_str(), string4.length()+1));
	unsigned int programHandle = ocldevice1->loadSources(sources4);

	//Allocate device buffers
	unsigned int d_RandGPU = ocldevice1->allocateBuffer(h_RandGPU.get(), size_h_RandGPU);
	unsigned int d_mtParams = ocldevice1->allocateBuffer(h_mtParams.get(), size_h_mtParams);
	
	//Load the kernel, set arguments, and launch it
	unsigned int kernelHandle = ocldevice1->loadKernel(programHandle,"MersenneTwister");
	ocldevice1->setKernelArg(kernelHandle, 0, *(ocldevice1->buffer(d_RandGPU)));
	ocldevice1->setKernelArg(kernelHandle, 1, *(ocldevice1->buffer(d_mtParams)));
	ocldevice1->setKernelArg(kernelHandle, 2, NVIDIA_RVs_PER_THREAD);
	unsigned int eventHandle = ocldevice1->launchKernel(kernelHandle, NVIDIA_MT_RNG_COUNT);

	// Wait for OpenCL execution to complete
	ocldevice1->wait(eventHandle);

	// Copy the result from the device buffer (d_RandGPU) to the host buffer (h_RandGPU)
	ocldevice1->readBuffer(d_RandGPU, h_RandGPU.get());

	//Calculate the average and print it to the screen
	float sum = 0.0;
	for(uint32_t i = 0; i < NVIDIA_nRand; i++)
		sum += h_RandGPU[i];

	float avg = sum / NVIDIA_nRand;

	std::cout << "Test 4 average = " << avg << std::endl;

	//All allocated OpenCL objects should be released automatically here
	//But for some reason one of our buffer objects is causing an exception at program termination
	//WTF?!?!?
}

// Test 5 - Generate normally distributed random variables using the GPL'd NVIDIA Mersenne Twister algorithm
// and a Box Muller transform
void test5(boost::shared_ptr<OclDevice> ocldevice) {
	
	//const parameters
	const int seed = 777;

	// Allocate space for the result
	boost::shared_array<float> h_RandGPU(new float[NVIDIA_nRand]);
	const size_t size_h_RandGPU = sizeof(float[NVIDIA_nRand]);
	// Allocate space for dynamic creation parameters
	boost::shared_array<mt_params_stripped> h_mtParams(new mt_params_stripped[NVIDIA_MT_RNG_COUNT]);
	const size_t size_h_mtParams = sizeof(mt_params_stripped[NVIDIA_MT_RNG_COUNT]);

	loadMersenneTwisterParams("data/MersenneTwister.dat", seed, h_mtParams.get(), NVIDIA_MT_RNG_COUNT);

	std::ifstream kernel_file("NVIDIA_mersennetwister.cl");
	std::string kernel_string(std::istreambuf_iterator<char>(kernel_file), (std::istreambuf_iterator<char>()));
	kernel_file.close();
	cl::Program::Sources sources(1, std::make_pair(kernel_string.c_str(), kernel_string.length()+1));
	unsigned int programHandle = ocldevice->loadSources(sources);

	//Allocate device buffers
	unsigned int d_RandGPU = ocldevice->allocateBuffer(h_RandGPU.get(), size_h_RandGPU);
	unsigned int d_mtParams = ocldevice->allocateBuffer(h_mtParams.get(), size_h_mtParams);
	
	//Load the Mersenne Twister kernel, set arguments, and launch it
	unsigned int mersenneTwisterKernelHandle = ocldevice->loadKernel(programHandle,"MersenneTwister");
	ocldevice->setKernelArg(mersenneTwisterKernelHandle, 0, *(ocldevice->buffer(d_RandGPU)));
	ocldevice->setKernelArg(mersenneTwisterKernelHandle, 1, *(ocldevice->buffer(d_mtParams)));
	ocldevice->setKernelArg(mersenneTwisterKernelHandle, 2, NVIDIA_RVs_PER_THREAD);
	unsigned int mersenneTwisterKernelEventHandle = ocldevice->launchKernel(mersenneTwisterKernelHandle, NVIDIA_MT_RNG_COUNT);

	//Load the Box Muller kernel, set arguments, and launch it
	unsigned int boxMullerKernelHandle = ocldevice->loadKernel(programHandle,"BoxMuller");
	ocldevice->setKernelArg(boxMullerKernelHandle, 0, *(ocldevice->buffer(d_RandGPU)));
	ocldevice->setKernelArg(boxMullerKernelHandle, 1, NVIDIA_RVs_PER_THREAD);
	unsigned int boxMullerKernelEventHandle = ocldevice->launchKernel(boxMullerKernelHandle, NVIDIA_MT_RNG_COUNT);

	// Wait for OpenCL execution to complete
	ocldevice->wait(boxMullerKernelEventHandle );

	// Copy the result from the device buffer (d_RandGPU) to the host buffer (h_RandGPU)
	ocldevice->readBuffer(d_RandGPU, h_RandGPU.get());

	//Calculate the average and print it to the screen
	float sum = 0.0;
	for(uint32_t i = 0; i < NVIDIA_nRand; i++)
		sum += h_RandGPU[i];

	float avg = sum / NVIDIA_nRand;

	std::cout << "Test 5 average = " << avg << std::endl;

	//All allocated OpenCL objects should be released automatically here
	//But for some reason one of our buffer objects is causing an exception at program termination
	//WTF?!?!?
}

// Test 6 - Compute the values of 32 options
void test6(boost::shared_ptr<OclDevice> ocldevice) {
	
	//const parameters
	const int seed = 42;
	const uint32_t numberOfOptions = 32;
	const uint32_t numberOfThreads = numberOfOptions;
	const uint32_t numberOfPaths = 1000000;
	const uint32_t timeStepsPerPath = 1000;

	// Allocate space for the result
	boost::shared_array<OpenCL_Option> h_Options(new OpenCL_Option[numberOfOptions]);
	const size_t size_h_Options = sizeof(OpenCL_Option[numberOfOptions]);
	// Allocate space for dynamic creation parameters
	boost::shared_array<mt_params_stripped> h_mtParams(new mt_params_stripped[numberOfThreads]);
	const size_t size_h_mtParams = sizeof(mt_params_stripped);

	//Give our option some values
	h_Options[0].X = 20.0f;
	h_Options[0].S = 50.0f;
	h_Options[0].V = 0.3f;
	h_Options[0].R = 0.03f;
	h_Options[0].T = 20.0f;

	//Give our option some values
	h_Options[1].X = 20.0f;
	h_Options[1].S = 50.0f;
	h_Options[1].V = 0.3f;
	h_Options[1].R = 0.03f;
	h_Options[1].T = 20.0f;

	loadMersenneTwisterParams("data/MersenneTwister.dat", seed, h_mtParams.get(), numberOfThreads);

	std::ifstream kernel_file("mcsimulation.cl");
	std::string kernel_string(std::istreambuf_iterator<char>(kernel_file), (std::istreambuf_iterator<char>()));
	kernel_file.close();
	cl::Program::Sources sources(1, std::make_pair(kernel_string.c_str(), kernel_string.length()+1));
	unsigned int programHandle = ocldevice->loadSources(sources);

	//Allocate device buffers
	unsigned int d_Options = ocldevice->allocateBuffer(h_Options.get(), size_h_Options);
	unsigned int d_mtParams = ocldevice->allocateBuffer(h_mtParams.get(), size_h_mtParams);
	
	//Load the kernel, set arguments, and launch it
	unsigned int kernelHandle = ocldevice->loadKernel(programHandle,"valueOptions");
	ocldevice->setKernelArg(kernelHandle, 0, *(ocldevice->buffer(d_Options)));
	ocldevice->setKernelArg(kernelHandle, 1, numberOfOptions);
	ocldevice->setKernelArg(kernelHandle, 2, numberOfPaths);	//Number of paths to generate
	ocldevice->setKernelArg(kernelHandle, 3, timeStepsPerPath);	//Timesteps per path
	ocldevice->setKernelArg(kernelHandle, 4, *(ocldevice->buffer(d_mtParams)));
	unsigned int kernelEventHandle = ocldevice->launchKernel(kernelHandle, numberOfThreads, 32);

	// Wait for OpenCL execution to complete
	ocldevice->wait(kernelEventHandle);

	// Copy the result from the device buffer (d_RandGPU) to the host buffer (h_RandGPU)
	ocldevice->readBuffer(d_Options, h_Options.get());

	std::cout << "Option 1 call value = " << h_Options[0].callValue << std::endl;
	std::cout << "Option 1 put value = " << h_Options[0].putValue << std::endl;

	std::cout << "Option 2 call value = " << h_Options[1].callValue << std::endl;
	std::cout << "Option 2 put value = " << h_Options[1].putValue << std::endl;
}

// Test 7 - Calculate the values of 65536 randomly generated options using OpenCL
// Then calculate the option values without OpenCL.  Compare the computation times.
void test7(boost::shared_ptr<OclDevice> ocldevice) {

	//const parameters
	const int seed = 42;
	const uint32_t numberOfOptions = 65536;
	const uint32_t numberOfThreads = 4096;
	const uint32_t numberOfPaths = 10000;
	const uint32_t timeStepsPerPath = 1;

	// Number of days to maturity for each option
	boost::shared_array<int> daysToMaturity(new int[numberOfOptions]);

	//Create a timer
	boost::timer t0;
	double clTime, scalarTime, quantlibTime;

	// Allocate space for the result
	boost::shared_array<OpenCL_Option> h_Options(new OpenCL_Option[numberOfOptions]);
	const size_t size_h_Options = sizeof(OpenCL_Option[numberOfOptions]);
	// Allocate space for dynamic creation parameters
	boost::shared_array<mt_params_stripped> h_mtParams(new mt_params_stripped[numberOfThreads]);
	const size_t size_h_mtParams = sizeof(mt_params_stripped);

	//Give our option some values
	for(uint32_t i = 0; i < numberOfOptions; i++) {
		daysToMaturity[i] = 1 + (int)(60.0 *(rand() / (RAND_MAX + 1.0)));	//An integer between 1 and 60
		h_Options[i].X = (float)(99.0f * (rand() / (RAND_MAX + 1.0)) + 1.0f);	//A number between 1.0 and 100.0
		h_Options[i].S = (float)(99.0f * (rand() / (RAND_MAX + 1.0)) + 1.0f);	//A number between 1.0 and 100.0
		h_Options[i].V = (float)(0.4f * (rand() / (RAND_MAX + 1.0)) + 0.1f);	//A number between 0.0 and 0.5
		h_Options[i].R = (float)(0.09f * (rand() / (RAND_MAX + 1.0)) + 0.01f);	//A number between 0.01 and 0.1
		h_Options[i].T = daysToMaturity[i] / 365.0f;	//A number of days between 1 and 60
	}

	loadMersenneTwisterParams("data/MersenneTwister.dat", seed, h_mtParams.get(), numberOfThreads);

	std::ifstream kernel_file("mcsimulation.cl");
	std::string kernel_string(std::istreambuf_iterator<char>(kernel_file), (std::istreambuf_iterator<char>()));
	kernel_file.close();
	cl::Program::Sources sources(1, std::make_pair(kernel_string.c_str(), kernel_string.length()+1));
	unsigned int programHandle = ocldevice->loadSources(sources);

	//Allocate device buffers
	unsigned int d_Options = ocldevice->allocateBuffer(h_Options.get(), size_h_Options);
	unsigned int d_mtParams = ocldevice->allocateBuffer(h_mtParams.get(), size_h_mtParams);

	//Load the kernel, set arguments, and launch it
	unsigned int kernelHandle = ocldevice->loadKernel(programHandle,"valueOptions");
	ocldevice->setKernelArg(kernelHandle, 0, *(ocldevice->buffer(d_Options)));
	ocldevice->setKernelArg(kernelHandle, 1, numberOfOptions);
	ocldevice->setKernelArg(kernelHandle, 2, numberOfPaths);	//Number of paths to generate
	ocldevice->setKernelArg(kernelHandle, 3, timeStepsPerPath);	//Timesteps per path
	ocldevice->setKernelArg(kernelHandle, 4, *(ocldevice->buffer(d_mtParams)));

	std::cout << "Launching OpenCL simulation..." << std::endl;

	t0.restart();	//restart our timer

	unsigned int kernelEventHandle = ocldevice->launchKernel(kernelHandle, numberOfThreads, 32);

	// Wait for OpenCL execution to complete
	ocldevice->wait(kernelEventHandle);

	clTime = t0.elapsed();

	std::cout << "Completed OpenCL simulation in " << clTime << " seconds" << std::endl;

	// Copy the result from the device buffer (d_RandGPU) to the host buffer (h_RandGPU)
	ocldevice->readBuffer(d_Options, h_Options.get());

    //Run the scalar simulation
    std::cout << "Launching scalar CPU-based simulation..." << std::endl;

    t0.restart();	//restart our timer
    valueOptions(h_Options.get(), numberOfOptions, numberOfPaths, timeStepsPerPath, h_mtParams.get());

    scalarTime = t0.elapsed();
    std::cout << "Completed scalar CPU-based simulation in " << scalarTime << " seconds" << std::endl;

	// BEGIN OF QuantLib TEST!!!
	boost::shared_array< boost::shared_ptr<PricingEngine> > mcEngine(new boost::shared_ptr<PricingEngine>[numberOfOptions]);

	Calendar calendar = TARGET();
    Date todaysDate(31, March, 1998);
    Date settlementDate(31, March, 1998);

    boost::shared_array<Real> npv(new Real[numberOfOptions]);	//Somewhere to store the results
    boost::shared_array< boost::shared_ptr<VanillaOption> > optionsArray(new boost::shared_ptr<VanillaOption>[numberOfOptions]);

    for(uint32_t i = 0; i < numberOfOptions; i++) {
    	 // set up dates
    	        Calendar calendar = TARGET();
    	        Settings::instance().evaluationDate() = todaysDate;

    	        boost::shared_ptr<Exercise> europeanExercise;

    	    	if(daysToMaturity[i] < 31) {
    	    		Date maturity(daysToMaturity[i], April, 1998);
    				boost::shared_ptr<Exercise> tmp(new EuropeanExercise(maturity));
    				europeanExercise = tmp;
    	    	}
    	    	else {
    				Date maturity(daysToMaturity[i]-30, May, 1998);
    				boost::shared_ptr<Exercise> tmp(new EuropeanExercise(maturity));
    				europeanExercise = tmp;
    	    	}

    	        Option::Type type(Option::Put);
				Real underlying = h_Options[i].S;
				Real strike = h_Options[i].X;
				Spread dividendYield = 0.00;
				Rate riskFreeRate = h_Options[i].R;
				Volatility volatility = h_Options[i].V;
				DayCounter dayCounter = Actual365Fixed();

    	        Handle<Quote> underlyingH(
    	            boost::shared_ptr<Quote>(new SimpleQuote(underlying)));

    	        // bootstrap the yield/dividend/vol curves
    	        Handle<YieldTermStructure> flatTermStructure(
    	            boost::shared_ptr<YieldTermStructure>(
    	                new FlatForward(settlementDate, riskFreeRate, dayCounter)));
    	        Handle<YieldTermStructure> flatDividendTS(
    	            boost::shared_ptr<YieldTermStructure>(
    	                new FlatForward(settlementDate, dividendYield, dayCounter)));
    	        Handle<BlackVolTermStructure> flatVolTS(
    	            boost::shared_ptr<BlackVolTermStructure>(
    	                new BlackConstantVol(settlementDate, calendar, volatility,
    	                                     dayCounter)));
    	        boost::shared_ptr<StrikedTypePayoff> payoff(
    	                                        new PlainVanillaPayoff(type, strike));
    	        boost::shared_ptr<BlackScholesMertonProcess> bsmProcess(
    	                 new BlackScholesMertonProcess(underlyingH, flatDividendTS,
    	                                               flatTermStructure, flatVolTS));

    	        // options
    	        boost::shared_ptr<VanillaOption> option(new VanillaOption(payoff, europeanExercise));

    	        // Monte Carlo Method: MC (crude)
    	        boost::shared_ptr<PricingEngine> mcengine1;
    	        mcengine1 = MakeMCEuropeanEngine<PseudoRandom>(bsmProcess)
    	            .withSteps(timeStepsPerPath)
    	            .withSamples(numberOfPaths)
    	            .withSeed(seed);
    	        option->setPricingEngine(mcengine1);
    	        // Real errorEstimate = europeanOption.errorEstimate();
    	        optionsArray[i] = option;
    	            }

    //Run the QuantLib simulation
    std::cout << "Launching QuantLib simulation..." << std::endl;

    t0.restart();	//restart our timer
    for(uint32_t j = 0; j < numberOfOptions; j++) {
    	//Do calculations here
        npv[j] = optionsArray[j]->NPV();
    }

    quantlibTime = t0.elapsed();
    std::cout << "Completed QuantLib simulation in " << quantlibTime << " seconds" << std::endl;

    //Print first 50 results to screen
	for(uint32_t j = 0; j < (numberOfOptions > 50 ? 50 : numberOfOptions); j++) {
		std::cout << "Option " << j << ":  OpenCL put value = " << h_Options[j].putValue << ", non-OpenCL put value = " << npv[j] << std::endl;
	}
}

int main(int, char* []) {

    try {
        boost::timer timer;
        std::cout << std::endl;

		boost::shared_ptr<OclDevice> ocldevice1;
		ocldevice1 = MakeOclDevice()
			.withDeviceType(CL_DEVICE_TYPE_GPU);

		// Run the tests
		//test1(ocldevice1);
		//test2(ocldevice1);
		//test3(ocldevice1);
		//test4(ocldevice1);
		//test5(ocldevice1);
		//test6(ocldevice1);
		test7(ocldevice1);

        // set up dates
        Calendar calendar = TARGET();
        Date todaysDate(15, May, 1998);
        Date settlementDate(17, May, 1998);
        Settings::instance().evaluationDate() = todaysDate;

        // our options
        Option::Type type(Option::Put);
        Real underlying = 36;
        Real strike = 40;
        Spread dividendYield = 0.00;
        Rate riskFreeRate = 0.06;
        Volatility volatility = 0.20;
        Date maturity(17, May, 1999);
        DayCounter dayCounter = Actual365Fixed();

        std::cout << "Option type = "  << type << std::endl;
        std::cout << "Maturity = "        << maturity << std::endl;
        std::cout << "Underlying price = "        << underlying << std::endl;
        std::cout << "Strike = "                  << strike << std::endl;
        std::cout << "Risk-free interest rate = " << io::rate(riskFreeRate)
                  << std::endl;
        std::cout << "Dividend yield = " << io::rate(dividendYield)
                  << std::endl;
        std::cout << "Volatility = " << io::volatility(volatility)
                  << std::endl;
        std::cout << std::endl;
        std::string method;
        std::cout << std::endl ;

        // write column headings
        Size widths[] = { 35, 14, 14, 14 };
        std::cout << std::setw(widths[0]) << std::left << "Method"
                  << std::setw(widths[1]) << std::left << "European"
                  << std::setw(widths[2]) << std::left << "Bermudan"
                  << std::setw(widths[3]) << std::left << "American"
                  << std::endl;

        boost::shared_ptr<Exercise> europeanExercise(
                                         new EuropeanExercise(maturity));

        Handle<Quote> underlyingH(
            boost::shared_ptr<Quote>(new SimpleQuote(underlying)));

        // bootstrap the yield/dividend/vol curves
        Handle<YieldTermStructure> flatTermStructure(
            boost::shared_ptr<YieldTermStructure>(
                new FlatForward(settlementDate, riskFreeRate, dayCounter)));
        Handle<YieldTermStructure> flatDividendTS(
            boost::shared_ptr<YieldTermStructure>(
                new FlatForward(settlementDate, dividendYield, dayCounter)));
        Handle<BlackVolTermStructure> flatVolTS(
            boost::shared_ptr<BlackVolTermStructure>(
                new BlackConstantVol(settlementDate, calendar, volatility,
                                     dayCounter)));
        boost::shared_ptr<StrikedTypePayoff> payoff(
                                        new PlainVanillaPayoff(type, strike));
        boost::shared_ptr<BlackScholesMertonProcess> bsmProcess(
                 new BlackScholesMertonProcess(underlyingH, flatDividendTS,
                                               flatTermStructure, flatVolTS));

        // options
        VanillaOption europeanOption(payoff, europeanExercise);

		// Analytic formulas:

        // Monte Carlo Method: MC (crude)
        Size timeSteps = 10;
        method = "MC (crude)";
        Size mcSeed = 42;
        boost::shared_ptr<PricingEngine> mcengine1;
        mcengine1 = MakeMCEuropeanEngine<PseudoRandom>(bsmProcess)
            .withSteps(timeSteps)
            .withAbsoluteTolerance(0.1)
            .withSeed(mcSeed);
        europeanOption.setPricingEngine(mcengine1);
        // Real errorEstimate = europeanOption.errorEstimate();
        Real npv = europeanOption.NPV();
        std::cout << std::setw(widths[0]) << std::left << method
                  << std::fixed
                  << std::setw(widths[1]) << std::left << npv
                  << std::setw(widths[2]) << std::left << "N/A"
                  << std::setw(widths[3]) << std::left << "N/A"
                  << std::endl;

		// OpenCL Monte Carlo Method: MC (crude)
        timeSteps = 10;
        method = "OpenCL MC (crude)";
        Size mcSeed2 = 42;
		boost::shared_ptr<PricingEngine> mcengine2;
		mcengine2 = MakeMCEuropeanEngine<PseudoRandom,Statistics,McSimulationCl>(bsmProcess)
            .withSteps(timeSteps)
            .withAbsoluteTolerance(0.1)
            .withSeed(mcSeed2);
        europeanOption.setPricingEngine(mcengine2);
        // Real errorEstimate = europeanOption.errorEstimate();
        std::cout << std::setw(widths[0]) << std::left << method
                  << std::fixed
                  << std::setw(widths[1]) << std::left << europeanOption.NPV()
                  << std::setw(widths[2]) << std::left << "N/A"
                  << std::setw(widths[3]) << std::left << "N/A"
                  << std::endl;

        // End test
        Real seconds = timer.elapsed();
        Integer hours = int(seconds/3600);
        seconds -= hours * 3600;
        Integer minutes = int(seconds/60);
        seconds -= minutes * 60;
        std::cout << " \nRun completed in ";
        if (hours > 0)
            std::cout << hours << " h ";
        if (hours > 0 || minutes > 0)
            std::cout << minutes << " m ";
        std::cout << std::fixed << std::setprecision(0)
                  << seconds << " s\n" << std::endl;

        return 0;

    } 
	catch (cl::Error& e) {
        std::cerr << e.what() << "(" << e.err() << ")" << std::endl;
		return 1;
    } catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "unknown error" << std::endl;
        return 1;
    }
}
