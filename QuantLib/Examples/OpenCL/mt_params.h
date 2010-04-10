#ifndef mt_params_h
#define mt_params_h

//Period parameters
#define MT_NN 624
#define MT_MM 397
#define MATRIX_A UINT32_C(0x9908b0df)   /* constant vector a */
#define MT_UMASK UINT32_C(0x80000000) /* most significant w-r bits */
#define MT_LMASK UINT32_C(0x7fffffff) /* least significant r bits */
#define TEMPERING_MASK_B UINT32_C(0x9d2c5680)
#define TEMPERING_MASK_C UINT32_C(0xefc60000)
#define TEMPERING_SHIFT_U(x)  (x >> 11)
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

#endif