# Introduction #

This project explores methods to offload certain quantitative finance computations to a GPU.

# Project details #

The project's main goals are as follows:
  * Identify algorithms currently in use in QuantLib which can be effectively offloaded to a modern GPU.
  * Port existing algorithms to a GPU platform.
  * Evaluate performance of original and GPU implementations.
  * Evaluate stability and accuracy of original and GPU implementations.
  * Evaluate the suitability of single-precision vs. double-precision computation.
  * Compare compute performance on various GPU platforms

While actual performance gains have not yet been determined, the goal of this project is to realize a minimum 100x increase in computation speed while preserving the stability and accuracy of the original algorithms.

# Technical specifications #
  * Target architecture:  Nvidia Fermi architecture (Compute Capability 3.0)
    * Secondary architecture:  Nvidia Compute Capability 1.1, 2.0
    * Secondary architecture:  AMD ATI Stream
  * Target language:  OpenCL
    * Secondary language:  CUDA