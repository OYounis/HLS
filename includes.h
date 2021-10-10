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


//#pragma HLS unroll

ap_uint<SIZE_N> barrett_reduction(
	ap_uint<SIZE_N+1> value,
	ap_uint<SIZE_N> modulus
	){

	ap_uint<SIZE_N> temp_1;
	ap_uint<SIZE_N+1> temp_2;
	ap_uint<2*SIZE_N+1> temp_3;

	ap_uint<SIZE_N> result;

	temp_1 = ap_uint<SIZE_N+1>(R_CONST) >> (SIZE_N - 1);
	temp_2 = ap_uint<2*SIZE_N+1>(R_CONST2) / modulus;

	temp_3 = temp_1 * temp_2;

	result = temp_3 >> (SIZE_N - 1);
	return result;
}
//TODO: REvise these
ap_uint<SIZE_N*2> barrett_reduction2(
	ap_uint<SIZE_N*2+1> value,
	ap_uint<SIZE_N*2> modulus
	){

	ap_uint<SIZE_N*2> temp_1;
	ap_uint<(SIZE_N*2)+1> temp_2;
	ap_uint<4*SIZE_N+1> temp_3;

	ap_uint<SIZE_N*2> result;

	temp_1 = ap_uint<SIZE_N*2+1>(R_CONST) >> (SIZE_N*2 - 1);
	temp_2 = ap_uint<2*SIZE_N*2+1>(R_CONST2) / modulus;

	temp_3 = temp_1 * temp_2;

	result = temp_3 >> (SIZE_N - 1);
	return result;
}

#if(MONT_IMP == 1)
    ap_uint<SIZE_N> rad2_mont_mult(
    	/*Source: A Scalable Architecture for Montgomery Multiplication*/
        ap_uint<SIZE_N> multiplicand,
        ap_uint<SIZE_N> multiplier,

        ap_uint<SIZE_N> modulus
    ){
        ap_uint<SIZE_N> result;
        ap_uint<SIZE_N> temp_1;
        ap_uint<SIZE_N> temp_2;

        result = 0;

        for(int i = 0; i < SIZE_N; i++){

            temp_1 = multiplicand[i] ? multiplier : ap_uint<SIZE_N>(0);
            temp_2 = result + temp_1;

            if(temp_2[0] == 0){ //even
                result = temp_2 >> 1;   //divide by 2
            } else{
                result = (temp_2 + modulus) >> 1; //divide by 2
            }
        }

        return result;
    }

#elif(MONT_IMP == 2)
    ap_uint<SIZE_N> rad2_mont_mult(
        ap_uint<SIZE_N> multiplicand,
        ap_uint<SIZE_N> multiplier,

        ap_uint<SIZE_N> modulus
    ){
        ap_uint<SIZE_N> result;
        ap_uint<1> temp_1;

        ap_uint<SIZE_N> temp_2;
        ap_uint<SIZE_N> temp_3;

        result = 0;

        for(int i = 0; i < SIZE_N; i++){
            temp_1 = result[0] ^ (multiplicand[i] & multiplier[0]);
            temp_2 = temp_1? modulus : ap_uint<SIZE_N>(0);
            temp_3 = multiplicand[i] ? multiplier : ap_uint<SIZE_N>(0);

            result = (result + temp_3 + temp_2) >> 1;
        }

        return result;
    }
#endif

ap_uint<SIZE_N> mont_mod_mult(

#ifndef CONST_N
	ap_uint<SIZE_N> R_SQUARED,
#endif

	ap_uint<SIZE_N> multiplicand,
	ap_uint<SIZE_N> multiplier,
	ap_uint<SIZE_N> modulus
){
	ap_uint<SIZE_N> temp_1;
	ap_uint<SIZE_N> result;

	temp_1 = rad2_mont_mult(multiplicand, multiplier, modulus);
	result = rad2_mont_mult(temp_1, R_SQUARED, modulus);

    return result;
}

/*********************************************************************************************************/
/*********************************************************************************************************/

//(a*b mod n)
ap_uint<SIZE_N> mod_mult(
	ap_uint<SIZE_N> a,
	ap_uint<SIZE_N> b,
	ap_uint<SIZE_N> modulus
){
	ap_uint<SIZE_N*2> temp;
	ap_uint<SIZE_N*2> quotient;
	ap_uint<SIZE_N*2> remainder;
	ap_uint<SIZE_N> result;

	temp     = (a * b);
	quotient = temp / modulus;
	remainder = temp - quotient * modulus;

	return result = remainder(SIZE_N-1,0);
}

/*********************************************************************************************************/
/*********************************************************************************************************/

/*	SQUARE-AND-MULTPLY EXPONENTIATION	*/

ap_uint<SIZE_N> sqrmul_exp(
#ifndef CONST_N
	ap_uint<SIZE_N> R_SQUARED,
#endif
	ap_uint<SIZE_N> base,
	ap_uint<SIZE_N> exponent,

	ap_uint<SIZE_N> modulus
){
	ap_uint<SIZE_N> result;
	ap_uint<SIZE_N> temp;

	if(exponent == 0) return ap_uint<SIZE_N>(1);

	temp = base;

	if(exponent[0] == 1)result = base;
	else result = 1;

	for(int i = 0; i < SIZE_N; i++){
		temp = mont_mod_mult(R_SQUARED, temp, temp, modulus);
		exponent = exponent >> 1;
		if(exponent & 1) result = mont_mod_mult(R_SQUARED, result, temp, modulus);
	}

	return result;
}

/*	MONTGOMREY EXPONENTIATION	*/
ap_uint<SIZE_N> mont_exp(
	ap_uint<SIZE_N> base_prime,
	ap_uint<SIZE_N> exponent,
	ap_uint<SIZE_N> i_prime,
	ap_uint<SIZE_N> modulus
){

	ap_uint<SIZE_N> z;
	ap_uint<SIZE_N> result;

	z = base_prime;
	result = i_prime;

	for(int i = 0; i < SIZE_N; i++){
		if(exponent[i] == 1)
			result = rad2_mont_mult(result, z, modulus);

		z = rad2_mont_mult(z, z, modulus);
	}

	result = rad2_mont_mult(result, z, modulus);
	result = rad2_mont_mult(result, 1, modulus);

	return result;
}

ap_uint<SIZE_N/2> gcd(ap_uint<SIZE_N/2> n1, ap_uint<SIZE_N/2> n2)
{
	ap_uint<SIZE_N/2> hcf;
	// swapping variables n1 and n2 if n2 is greater than n1.
	if ( n2 > n1) {
		ap_uint<SIZE_N/2> temp = n2;
		n2 = n1;
		n1 = temp;
	}

	for (int i = 1; i <=  n2; ++i){
		if (n1 % i == 0 && n2 % i ==0){
			hcf = i;
		}
	}

	return hcf;
}

ap_uint<SIZE_N> modInverse(ap_uint<SIZE_N> a, ap_uint<SIZE_N> m)
{
	ap_uint<SIZE_N> m0 = m;
	ap_uint<SIZE_N> y = 0, x = 1;

	if (m == 1)
		return 0;

	while (a > 1) {
		ap_uint<SIZE_N> q = a / m;
		ap_uint<SIZE_N> t = m;
		m = a % m, a = t;
		t = y;
		y = x - q * y;
		x = t;
	}

	if (x < 0)
		x += m0;

	return x;
}

//LFSR
ap_uint<SIZE_N> pseudo_random(ap_uint<SIZE_N> seed, ap_uint<1> load) {
  static ap_uint<SIZE_N> lfsr;

  if (load ==1 )
    lfsr = seed;

  bool b_32 = lfsr.get_bit(32-32);
  bool b_22 = lfsr.get_bit(32-22);
  bool b_2 = lfsr.get_bit(32-2);
  bool b_1 = lfsr.get_bit(32-1);
  bool new_bit = b_32 ^ b_22 ^ b_2 ^ b_1;

  lfsr = lfsr >> 1;
  lfsr.set_bit(31, new_bit);

  return lfsr;
}
