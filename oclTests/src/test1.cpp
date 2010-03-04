//============================================================================
// Name        : test1.cpp
// Author      : William Gross
// Version     :
// Copyright   : Public Domain
// Description : Hello World in C, Ansi-style
//============================================================================

#include <CL/cl.h>
#include <stdlib.h>
#include <oclUtils.h>

float randFloat(float low, float high){
    float t = (float)rand() / (float)RAND_MAX;
    return (1.0f - t) * low + t * high;
}

int main(int argc, const char **argv)
{
    cl_platform_id cpPlatform;       //OpenCL platform
    cl_device_id cdDevice;           //OpenCL devices
    cl_context      cxGPUContext;    //OpenCL context
    cl_command_queue cqCommandQueue; //OpenCL command que
    cl_mem                           //OpenCL memory buffer objects
        d_Call,
        d_Put,
        d_S,
        d_X,
        d_T;

    cl_int ciErrNum;

    float
        *h_CallCPU,
        *h_PutCPU,
        *h_CallGPU,
        *h_PutGPU,
        *h_S,
        *h_X,
        *h_T;

    const unsigned int   optionCount = 4000000;
    const float                    R = 0.02f;
    const float                    V = 0.30f;

    //Allocating and initializing host memory
    h_CallCPU = (float *)malloc(optionCount * sizeof(float));
    h_PutCPU  = (float *)malloc(optionCount * sizeof(float));
    h_CallGPU = (float *)malloc(optionCount * sizeof(float));
    h_PutGPU  = (float *)malloc(optionCount * sizeof(float));
    h_S       = (float *)malloc(optionCount * sizeof(float));
    h_X       = (float *)malloc(optionCount * sizeof(float));
    h_T       = (float *)malloc(optionCount * sizeof(float));

    for(unsigned int i = 0; i < optionCount; i++)
    {
    	h_CallCPU[i] = -1.0f;
    	h_PutCPU[i]  = -1.0f;
    	h_S[i]       = randFloat(5.0f, 30.0f);
    	h_X[i]       = randFloat(1.0f, 100.0f);
    	h_T[i]       = randFloat(0.25f, 10.0f);
    }

    // Get the NVIDIA platform
    ciErrNum = oclGetPlatformID(&cpPlatform);
    oclCheckError(ciErrNum, CL_SUCCESS);
}
