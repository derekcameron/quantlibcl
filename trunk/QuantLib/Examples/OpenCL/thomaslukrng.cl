/*
Copyright (c) 2009, Imperial College London
Copyright (c) 2010, William Gross
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/



//////////////////////////////////////////////////////////////////////////////////////
// Private constants
// By definition, in OpenCL "__constant" variables are global

__constant unsigned int WarpBuffered_Q[3*32]={
  2,26,11,6,12,8,13,31,16,4,25,3,21,23,1,14,17,7,9,20,27,19,30,22,18,5,24,29,15,0,28,10,
  8,20,31,21,3,25,26,29,17,30,7,0,11,2,22,12,19,18,27,5,10,1,24,16,4,15,23,13,14,28,9,6,
  31,30,28,24,10,9,27,13,23,22,3,2,1,21,16,20,25,29,5,11,15,17,7,8,26,14,4,18,6,12,19,0
};
__constant unsigned int WarpBuffered_Z0=19;
__constant unsigned int WarpBuffered_Z1[32]={
  9,16,12,16,11,16,13,11,13,16,16,13,9,16,14,9,15,9,15,12,14,15,11,9,15,8,9,16,8,9,12,13};
__constant unsigned int WarpBuffered_SHMEM_WORDS=32;
__constant unsigned int WarpBuffered_GMEM_WORDS=32;


////////////////////////////////////////////////////////////////////////////////////////
// Public functions
/*
__device__ void WarpBuffered_LoadState(const unsigned *seed, unsigned *regs, unsigned *shmem)
{
  unsigned offset=threadIdx.x % 32;  unsigned base=threadIdx.x-offset;
  // setup constants
  regs[0]=WarpBuffered_Z1[offset];
  regs[1]=base + WarpBuffered_Q[0][offset];
  regs[2]=base + WarpBuffered_Q[1][offset];
  regs[3]=base + WarpBuffered_Q[2][offset];
  // Setup state
  unsigned stateOff=blockDim.x * blockIdx.x * 2 + threadIdx.x * 2;
  shmem[threadIdx.x]=seed[stateOff];
  regs[4]=seed[stateOff+1];
}*/
/*
__device__ void WarpBuffered_SaveState(const unsigned *regs, const unsigned *shmem, unsigned *seed)
{
  unsigned stateOff=blockDim.x * blockIdx.x * 2 + threadIdx.x * 2;
  seed[stateOff] = shmem[threadIdx.x];
  seed[stateOff+1]=regs[4];
}*/

/*__device__ unsigned WarpBuffered_Generate(unsigned int *regs, unsigned int *shmem)
{
  unsigned int t0=shmem[regs[1]], t1=shmem[regs[2]];
  unsigned int res=(t0<<WarpBuffered_Z0) ^ (t1>>regs[0]) ^ regs[4];
  regs[4] = shmem[regs[3]];
  
  shmem[threadIdx.x]=res;
  return t0+t1;
};
*/