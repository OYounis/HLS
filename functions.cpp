#include "includes.h"
//#pragma HLS unroll

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
	ap_uint<SIZE_N> R_PRIME,
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
		temp = mont_mod_mult(temp, temp, modulus);
		exponent = exponent >> 1;
		if(exponent & 1) result = mont_mod_mult(result, temp, modulus);
	}

	return result;
}

/*	MONTGOMREY EXPONENTIATION	*/
ap_uint<SIZE_N> mont_exp(
	ap_uint<SIZE_N> base_prime,
	ap_uint<SIZE_N> exponent,
	ap_uint<SIZE_N> i_prime,
	ap_uint<SIZE_N> r_prime,
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


/*********************************************************************************************************/
/*********************************************************************************************************/


ap_uint<SIZE_N> mexp(
	ap_uint<SIZE_N> base,
	ap_uint<SIZE_N> exponent,

	ap_uint<SIZE_N> modulus
){
	//return sqrmul_exp(base, exponent, modulus);
    //return rad2_mont_mult(base, exponent, modulus);
    return mont_mod_mult(modulus,base, exponent);
	//return mont_exp(11, exponent,9, 18, modulus);
	//return mont_exp(base, exponent, modulus);
	//return mul(base, exponent, modulus);
}
