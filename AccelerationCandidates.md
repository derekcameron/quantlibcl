# Overview #

By their very nature, many of the quantitative problems facing the financial sector today are suited to Monte Carlo simulation, which is generally parallelizable and can be accelerated on modern GPU hardware.  Among the pricing engines implemented in QuantLib, the following appear to be good candidates for GPU acceleration:

**Asian option engines**:
  * [MCDiscreteArithmeticAPEngine](http://quantlib.org/reference/class_quant_lib_1_1_m_c_discrete_arithmetic_a_p_engine.html) - Monte Carlo pricing engine for discrete arithmetic average price Asian

**Vanilla option engines**:
  * [MCAmericanEngine](MCAmericanEngine.md) - A Monte Carlo American option pricing engine
  * [MCVanillaEngine](http://quantlib.org/reference/class_quant_lib_1_1_m_c_vanilla_engine.html) - A pricing engine for vanilla options using Monte Carlo simulation
  * [MCEuropeanEngine](http://quantlib.org/reference/class_quant_lib_1_1_m_c_european_engine.html) - A European option pricing engine using Monte Carlo simulation