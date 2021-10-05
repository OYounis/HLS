#define AP_INT_MAX_W 4096

#include <ap_int.h>
#include <hls_math.h>

//N can be 128, 256, 512, 1024
#define SIZE_N 128
#define MONT_IMP 2
//#define CONST_N

#ifdef CONST_N
	#define R_CONST
	#define R_CONST2 12
#else
	#if (SIZE_N == 128)
    	#define R_CONST	 	3.4028237e+38	//2^n
    	#define R_CONST2 	1.1579209e+77	//2^2n
        
	#elif (SIZE_N == 256)
    	#define R_CONST		1.1579209e+77	//2^n
    	#define R_CONST2	1.340781e+154	//2^2n

	#elif (SIZE_N == 512)
    	#define R_CONST		1.340781e+154	//2^n
    	#define R_CONST2	1.797693134e+308//2^2n

	#else //(SIZE_N == 1024)
    	#define R_CONST		1.797693134e+308//2^n
    	#define R_CONST 	3.231700607e+616//2^2n
	#endif
#endif
