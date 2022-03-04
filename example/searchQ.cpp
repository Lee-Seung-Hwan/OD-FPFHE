#include "fpfhecontext.h"
//#include "time.h" 
#include <chrono>
#include <iostream>
//#include "lwe.h"
#include <cmath>

using namespace lbcrypto;
using namespace std;

int main() {	
	std::cout << "Q check starts!" << std::endl;
	
	NativeInteger test_n = 12289;
	bool test_res = MillerRabinPrimalityTest(test_n);
	std::cout << "Miler Test result is "<<test_res << std::endl;
	
	uint32_t start_bitsQ1 = 20;
	uint32_t search_gapQ1 = 20;
	uint32_t end_bitsQ1 = start_bitsQ1 + search_gapQ1;
	
	//uint32_t prime_gap = 10;
	//uint32_t prime_gap = 9;
	uint32_t prime_gap = 17;

	NativeInteger P = 1;
	NativeInteger Q = 1;
	uint64_t iter_num = 1;
	for (uint64_t i = 0; i < end_bitsQ1 - start_bitsQ1 - 1; i++) {
		iter_num = iter_num * 2;
	}
	for (uint32_t i = 0; i < start_bitsQ1 +prime_gap; i++){
		P = P * 2;
	}
	for (uint32_t i = 0; i < (start_bitsQ1); i++){
		Q =Q * 2;
	}
	std::cout<<"P is " << P << std::endl;
	bool result;
	uint32_t P_bits;
	uint32_t Q_bits;
	for (uint64_t i = 0; i < iter_num; i ++) { 
		result = MillerRabinPrimalityTest(P*(2*i-1)+1) && MillerRabinPrimalityTest(Q*(2*i-1)+1);

		if (result) {
			P_bits = (uint32_t)std::ceil(log((P*(2*i-1)+1).ConvertToDouble()) / log(static_cast<double>(2)));
			Q_bits = (uint32_t)std::ceil(log((Q*(2*i-1)+1).ConvertToDouble()) / log(static_cast<double>(2)));

			std::cout << "We Find it !!!" << std::endl
				<< "multiplicative odd is " << (2*i-1) << std::endl
				<< "P bits is "<< P_bits <<", and value is " << P*(2*i-1)+1 << std::endl
				<< "Q bits is "<< Q_bits <<", and value is " << Q*(2*i-1)+1 << std::endl;
		}
	}


	// Overflow Test
	uint64_t Q1_nat = 1451355348664321;
	NativeInteger Q1 = NativeInteger(Q1_nat);
	BigInteger Q1_big = BigInteger(1451355348664321);

	uint32_t scalar = 1;
	scalar <<= 12;
	//scalar *= 8;
	//scalar *=
	Q1.MulEq(scalar);
	Q1_nat *= scalar;
	Q1_big.MulEq(scalar);
	std::cout << "Mul eq result is " << Q1
		<< ", native mul res is " << Q1_nat << std::endl
		<< ", Big mul res is " << Q1_big << std::endl;
	/*
	NativeInteger Zero = NativeInteger(0);
	Zero.SubEq(5);
	Zero.AddEq(5);
	std::cout << "minus Zero is " << Zero << std::endl;
	*/

	return 0;
}
