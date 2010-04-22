#ifndef mt_params_h
#define mt_params_h

#include <boost/cstdint.hpp>

//Period parameters
#define MT_NN 19
#define MT_MM 9
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define MT_WMASK 0xFFFFFFFFU
#define MT_UMASK 0xFFFFFFFEU /* most significant w-r bits */
#define MT_LMASK 0x1U /* least significant r bits */
#define TEMPERING_MASK_B 0x9d2c5680UL
#define TEMPERING_MASK_C 0xefc60000UL
#define TEMPERING_SHIFT_U(x)  (x >> 12)
#define TEMPERING_SHIFT_S(x)  (x << 7)
#define TEMPERING_SHIFT_T(x)  (x << 15)
#define TEMPERING_SHIFT_L(x)  (x >> 18)
#define SSS 7
#define TTT 15
/* #define S00 11 */
#define S00 12
#define S01 18

typedef struct {
	//Crucial data (required for MT generation)
	uint32_t aaa;
	uint32_t maskB, maskC;
	uint32_t seed;

	//Non-crucial data (can be discarded after initialization)
	int mm,nn,rr,ww;
	uint32_t *state;
	uint32_t wmask,umask,lmask;
	int shift0, shift1, shiftB, shiftC;
} mt_params;

typedef struct{
  uint32_t matrix_a;
  uint32_t mask_b;
  uint32_t mask_c;
  uint32_t seed;
} mt_params_stripped;

#endif