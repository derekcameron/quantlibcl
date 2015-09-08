To use the MCEuropeanEngine, one simply makes a call to the [MakeMCEuropeanEngine() method](http://quantlib.org/reference/class_quant_lib_1_1_make_m_c_european_engine.html).  The MakeMCEuropeanEngine() method takes one argument, a boost shared pointer to an instance of the [GeneralizedBlackScholesProcess class](http://quantlib.org/reference/class_quant_lib_1_1_generalized_black_scholes_process.html) (or one of its inherited classes).

## Usage ##
```
#include <ql/quantlib.hpp>
#include <iostream>
#include <iomanip>

// Set the number of time steps
timeSteps = 1;

// Create a seed for our random number BSM process
Size mcSeed = 42;

// Create a boost shared pointer to our MC engine
boost::shared_ptr<PricingEngine> mcengine1;

mcengine1 = MakeMCEuropeanEngine<PseudoRandom>(bsmProcess)
            .withSteps(timeSteps)
            .withAbsoluteTolerance(0.02)
            .withSeed(mcSeed);

europeanOption.setPricingEngine(mcengine1);

// Real errorEstimate = europeanOption.errorEstimate();
std::cout << "NPV = " << europeanOption.NPV() << endl;
```

All the magic is done in the call to europeanOption.NPV().  Execution follows the following path:

Instrument::NPV() -> Instrument::calculate() -> LazyObject::calculate() -> Instrument::performCalculations() -> MCVanillaEngine::calculate() -> McSimulation<MC,RNG,S>::value()

The McSimulation<MC,RNG,S>::value() method is where the actual Monte Carlo simulation is performed, and all the dozens of method calls up to this point were simply preparing us for this step.  The pseudocode for this method is:

```
ARGUMENTS:  tolerance, maxSamples, minSamples

1. Determine a number of samples to compute (sampleNumber)
2. If sampleNumber is less than minSamples, generate samples until sampleNumber >= minSamples
3. Estimate the maximum error in our sample set
4. Loop until maximum error <= tolerance
     a.  samplesNeeded = max(currentNumberOfSamples*currentSampleError^2/tolerance^2*0.8-currentNumberOfSamples, minSamples)
     b.  If samplesNeeded > maxSamples, samplesNeeded = maxSamples - currentNumberOfSamples
     c.  Generate samplesNeeded new samples
     d.  Recalculate the error of the sample set
5.  value = mean of the sample set
6.  Return value
```

## Characteristics of a MonteCarloModel ##
  * Components include path generator, path pricer, and sample accumulator