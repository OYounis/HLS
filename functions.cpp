#include "includes.h"

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

ap_uint<SIZE_N> paillier_decrypt(
	ap_uint<SIZE_N> i_prime,
	ap_uint<SIZE_N> square_r,

	ap_uint<SIZE_N> cipher,

	ap_uint<SIZE_N> key_lambda,
	ap_uint<SIZE_N> key_mu,

	ap_uint<SIZE_N> n,
	ap_uint<SIZE_N> n_squared
){
	ap_uint<SIZE_N> cipher_prime;
	ap_uint<SIZE_N> temp_1;

	ap_uint<SIZE_N> l_result;

	///c^lambda mod n^2
	cipher_prime = mont_mod_mult(square_r, cipher, R_CONST, n_squared);
	temp_1 = mont_exp(cipher_prime, key_lambda, i_prime, n_squared);

	//L = (x-1) / n
	l_result = (temp_1 - 1) / n;

	return mont_mod_mult(square_r, l_result, key_mu, n);
}

void  paillier_key_gen(
		ap_uint<SIZE_N> i_prime,
		ap_uint<SIZE_N> square_r,

		ap_uint<SIZE_N/2> p,
		ap_uint<SIZE_N/2> q,

		ap_uint<SIZE_N> *n,
		ap_uint<SIZE_N> *g,

		ap_uint<SIZE_N> *lambda,
		ap_uint<SIZE_N> *mu
){

	ap_uint<SIZE_N> n_temp;
	ap_uint<SIZE_N> g_temp;

	ap_uint<SIZE_N> lambda_temp;
	ap_uint<SIZE_N> mu_temp;

	n_temp = p*q;

	ap_uint<SIZE_N> n_squared;
	n_squared = n_temp * n_temp;

	lambda_temp = ((p-1)/gcd(p-1, q-1))*(q-1);


	for(int i = 0; i < (n_squared) - 1; i++){
		if(gcd(n_squared, i) == 1) g_temp = i;	//TODO: g should be randomly chosen from a set that meets the conditions
	}


	ap_uint<SIZE_N> g_prime;
	g_prime = mont_mod_mult(square_r, g_temp, R_CONST, n_squared);

	ap_uint<SIZE_N> temp_1;
	temp_1 = mont_exp(g_prime, lambda_temp, i_prime, n_squared);

	ap_uint<SIZE_N> temp_2;
	temp_2 = (temp_1 - 1) / n_temp;
	mu_temp = modInverse(temp_2, n_temp);

	*n      = n_temp;
	*g      = g_temp;
	*lambda = lambda_temp;
	*mu     = mu_temp;
}

void paillier(
	//enc = 1, dec = 2, noop = 0
	ap_uint<2> operation,

	ap_uint<SIZE_N*2> data_in,

	ap_uint<SIZE_N/2> p,
	ap_uint<SIZE_N/2> q,

	ap_uint<SIZE_N*2> *data_out
){

	ap_uint<SIZE_N> n;
	ap_uint<SIZE_N*2> n_squared;
	ap_uint<SIZE_N> lambda;
	ap_uint<SIZE_N> index;
	ap_uint<SIZE_N> g;
	ap_uint<SIZE_N> g_prime;
	ap_uint<SIZE_N> temp_1;
	ap_uint<SIZE_N> temp_2;
	ap_uint<SIZE_N> mu;

	ap_uint<SIZE_N> R_REDUCED;
	ap_uint<SIZE_N> R_SQUARED;

	ap_uint<SIZE_N*2> R2_REDUCED;
	ap_uint<SIZE_N*2> R2_SQUARED;

	/********************************************/

	ap_uint<SIZE_N*2> temp_3;

	ap_uint<SIZE_N*2> temp_4;
	ap_uint<SIZE_N> r_prime;

	ap_uint<SIZE_N> result;

	ap_uint<SIZE_N> cipher_prime;

	ap_uint<SIZE_N> l_result;

	if(operation == 3){	//if enc or dec
		//key Generation
		n = p * q;
		R_REDUCED = barrett_reduction(R_CONST ,n); //i_prime
		R_SQUARED = barrett_reduction(R_CONST ,n); //

		n_squared = n * n;
		R2_REDUCED = barrett_reduction2(R_CONST2, n_squared);
		R2_SQUARED = barrett_reduction2(R_CONST2, n_squared);

		lambda = ((p-1)/gcd(p-1, q-1))*(q-1);

		index =  pseudo_random(n, 1);

		for(int i = 0; i < n; i++){
			if(gcd(n_squared, i) == 1){
				if(i >= index) g = i;
			}
		}

		g_prime = mont_mod_mult(R2_SQUARED, g, R_CONST, n_squared); /*||||||||*/

		temp_1 = mont_exp(g_prime, lambda, R2_REDUCED, n_squared);

		temp_2 = (temp_1 - 1) / n;
		mu = modInverse(temp_2, n);

	}else if(operation == 1){
		temp_3 = mont_exp(g_prime, data_in, R2_REDUCED, n_squared);

		r_prime = pseudo_random(0, 0);
		temp_4 = mont_exp(r_prime, n, R2_REDUCED, n_squared);

		*data_out = mont_mod_mult(R2_SQUARED, temp_3, temp_2, n_squared);

	} else if(operation == 2){
		///c^lambda mod n^2
		cipher_prime = mont_mod_mult(R2_SQUARED, data_in, R_CONST, n_squared);
		temp_3 = mont_exp(cipher_prime, lambda, R2_REDUCED, n_squared);

		//L = (x-1) / n
		l_result = (temp_1 - 1) / n;

		*data_out = mont_mod_mult(R_SQUARED, l_result, mu, n);
	}
}
