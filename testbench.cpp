#include <ap_int.h>
#include <stdio.h>
#include <iostream>

#define n_size 128

ap_uint<n_size> mexp(ap_uint<n_size> base, ap_uint<n_size> exponent, ap_uint<n_size> modulus);

int  main(){
	ap_uint<n_size> base = 14;
	ap_uint<n_size> exponent = 25
	ap_uint<n_size> modulus = 23;
	ap_uint<n_size> result;


	result = mexp(base, exponent, modulus);

	std::cout << result << std::endl;

	return 0;
}
