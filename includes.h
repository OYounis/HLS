#define AP_INT_MAX_W 4096

#include <ap_int.h>
#include <hls_math.h>

//N can be 128, 256, 512, 1024
#define SIZE_N 128
#define MONT_IMP 1
#define CONST_N

#ifdef CONST_N
	#define R_PPRIME
	#define R_SQUARED 12
#endif


	#if (SIZE_N == 128)
    	#define R	3.4028237e+38//2^n
    	#define R_prime	2.938735877055739 //

	#elif (SIZE_N == 256)
    	#define R	1.1579209e+77//2^n
    	#define R_prime	8.636168555094478//

	#elif (SIZE_N == 512)
    	#define R	1.340781e+154//2^n
    	#define R_prime	7.458340731200255//

	#else //(SIZE_N == 1024)
    	#define R	1.797693134 E+308//2^n
    	#define R_prime	5.562684646268309//
	#endif
