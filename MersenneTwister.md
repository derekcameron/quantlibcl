# Design considerations #
From NVIDIA's CUDA Programming Guide:
> _"A warp executes one common instruction at a time, so full efficiency is realized when all 32 threads of a warp agree on their execution path. If threads of a warp diverge via a data-dependent conditional branch, the warp serially executes each branch path taken, disabling threads that are not on that path, and when all paths complete, the threads converge back to the same execution path."_

# Design assumptions #
  * 32 threads per warp
  * 

# Details #