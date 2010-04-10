#ifndef NVIDIA_MERSENNETWISTER_HPP
#define NVIDIA_MERSENNETWISTER_H

// Total # of RNG threads to spawn
#define NVIDIA_MT_RNG_COUNT 4096

// # of rvs that each thread will generate
#define NVIDIA_RVs_PER_THREAD 5860

// Total # of rvs that will be generated
#define NVIDIA_nRand (NVIDIA_MT_RNG_COUNT * NVIDIA_RVs_PER_THREAD)

typedef struct {
  unsigned int matrix_a;
  unsigned int mask_b;
  unsigned int mask_c;
  unsigned int seed;
} mt_params_stripped;

#endif