#include "includes.h"

//#pragma HLS unroll

ap_uint<SIZE_N> barrett_reduction(
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

ap_uint<SIZE_N> square_r(
	ap_uint<SIZE_N> modulus
	){

	ap_uint<SIZE_N> r_reduced;

	r_reduced = barrett_reduction(modulus);

	return barrett_reduction(r_reduced * r_reduced);
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

ap_uint<SIZE_N> paillier_encrypt(
	ap_uint<SIZE_N> i_prime,
	ap_uint<SIZE_N> square_r,

	ap_uint<SIZE_N> key_g,
	ap_uint<SIZE_N> message,

	ap_uint<SIZE_N> r,
	ap_uint<SIZE_N> n,

	ap_uint<SIZE_N> n_squared
){

	ap_uint<SIZE_N> temp_1;
	ap_uint<SIZE_N> key_g_prime;

	ap_uint<SIZE_N> temp_2;
	ap_uint<SIZE_N> r_prime;

	ap_uint<SIZE_N> result;

	key_g_prime = mont_mod_mult(square_r, key_g, R_CONST, n_squared);
	temp_1 = mont_exp(key_g_prime, message, i_prime, n_squared);

	r_prime = mont_mod_mult(square_r, r, R_CONST, n_squared);
	temp_2 = mont_exp(r_prime, n, i_prime, n_squared);

	result = mont_mod_mult(square_r, temp_1, temp_2, n_squared);
	return result;
}
