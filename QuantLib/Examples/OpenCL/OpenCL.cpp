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

#define OCL_THREAD_TOTAL 1024

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
	unsigned int kernelHandle = ocldevice1->loadKernel(sourcesIndex,"oclTest1", 1, bufferHandle);
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

void test4LoadParameters(const char *fname, 
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
			throw std::exception("Error opening mt_params raw data file in test 4");
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

	test4LoadParameters("data/MersenneTwister.dat", seed, h_mtParams.get(), NVIDIA_MT_RNG_COUNT);

	for (int i = 0; i < NVIDIA_MT_RNG_COUNT; i++) {
		std::cout << "Item " << i << std::endl;
		std::cout << "aaa = " << h_mtParams[i].matrix_a << std::endl;
		std::cout << "maskB = " << h_mtParams[i].mask_b << std::endl;
		std::cout << "maskC = " << h_mtParams[i].mask_c << std::endl;
		std::cout << "seed = " << h_mtParams[i].seed << std::endl << std::endl;
	}

	std::ifstream file4("NVIDIA_mersennetwister.cl");
	std::string string4(std::istreambuf_iterator<char>(file4), (std::istreambuf_iterator<char>()));
	file4.close();
	cl::Program::Sources sources4(1, std::make_pair(string4.c_str(), string4.length()+1));
	size_t sources4Index = ocldevice1->loadSources(sources4);

	//Allocate device buffers
	unsigned int d_RandGPU = ocldevice1->allocateBuffer(h_RandGPU.get(), size_h_RandGPU);
	unsigned int d_mtParams = ocldevice1->allocateBuffer(h_mtParams.get(), size_h_mtParams);
	
	//Load and launch the kernel
	//Our error is in the next line...the loadKernel command assumes the only arguments to kernels are buffers
	//But in our case, we need to modify our code to support other argument types (like ints)
	unsigned int kernelHandle = ocldevice1->loadKernel(sources4Index,"MersenneTwister", 1, ocldevice1->buffer(d_RandGPU), ocldevice1->buffer(d_mtParams), NVIDIA_RVs_PER_THREAD);
	unsigned int eventHandle = ocldevice1->launchKernel(kernelHandle, NVIDIA_MT_RNG_COUNT);
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
		test3(ocldevice1);
		test4(ocldevice1);

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

        std::vector<Date> exerciseDates;
        for (Integer i=1; i<=4; i++)
            exerciseDates.push_back(settlementDate + 3*i*Months);

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
        Size timeSteps = 1;
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
        timeSteps = 1;
        method = "OpenCL MC (crude)";
        Size mcSeed2 = 42;
		boost::shared_ptr<PricingEngine> mcengine2;
		mcengine2 = MakeMCEuropeanEngineOCL<PseudoRandom>(bsmProcess)
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