#include "fpfhecontext.h"
#include <cmath>
#include <time.h>
#include <limits>
#include <float.h>
#include <omp.h>
//#include <double.h>
typedef std::numeric_limits<double> dbl;

using namespace lbcrypto;
using namespace std;


int main() {
	std::cout << "Amortized test" << std::endl;
	std::cout <<omp_get_max_threads() << std::endl;
	std::cout << "Nested condition: " << omp_get_nested() << std::endl;
	omp_set_nested(0);
	cout.precision(dbl::max_digits10);
	// MAYBE GINX
	auto cc = FPFHEContext();
	bool PARALLELs = true;
	omp_set_num_threads(20);
	

	//cc.GenerateFHEContext(FPFHE, PARALLELs);
	//cc.GenerateFHEContext(STD128D, PARALLELs);
	cc.GenerateFHEContext(STD160D, PARALLELs);
	//cc.GenerateFHEContext(STD192D, PARALLELs);
	//cc.GenerateFHEContext(STD128S, PARALLELs);
    //cc.GenerateFHEContext(STD160S, PARALLELs);
	//cc.GenerateFHEContext(STD192S, PARALLELs);

    auto Rsk = cc.RKeyGen();
	cc.BTKeyGen(Rsk);
	//auto OFUF = cc.MakeOverFlowDetectCT();
	uint32_t nums = 1000;	
    

	double msg2 = 2.0;
	//vector<FPRLWECiphertext> ind1(nums);
	//vector<FPRLWECiphertext> ind2(nums);
	//vector<FPRLWECiphertext> outad(nums);
	//vector<FPRLWECiphertext> outpd(nums);
	//vector<std::shared_ptr<FPOFUF>> OFUF_arr2(nums);

	std::shared_ptr<vector<FPRLWECiphertext>> ind1 = 
	 std::make_shared<vector<FPRLWECiphertext>> (nums);
	std::shared_ptr<vector<FPRLWECiphertext>> ind2 = 
	 std::make_shared<vector<FPRLWECiphertext>> (nums);
	std::shared_ptr<vector<FPRLWECiphertext>> outad = 
	 std::make_shared<vector<FPRLWECiphertext>> (nums);
	std::shared_ptr<vector<FPRLWECiphertext>> outpd = 
	 std::make_shared<vector<FPRLWECiphertext>> (nums);
	std::shared_ptr<vector<std::shared_ptr<FPOFUF>>> OFUF_arr2 = 
	 std::make_shared<vector<std::shared_ptr<FPOFUF>>> (nums);
	



	#pragma omp parallel for collapse(1)
	for (uint32_t i = 0; i < nums; i++) {
		(*ind1)[i] = cc.EncryptDouble(Rsk, msg2);		
		(*ind2)[i] = cc.EncryptDouble(Rsk, msg2);		
		(*OFUF_arr2)[i] = cc.MakeOverFlowDetectCT();
	}
	std::cout << "Double Encryption is end" << std::endl ;  
	for (uint32_t i = 0; i < nums; i++) {
		if ( cc.DecryptDouble(Rsk, (*ind1)[i]) != msg2) {
			std::cout << i << "'st double1 is not " << msg2 <<", but" << 
			cc.DecryptDouble(Rsk, (*ind1)[i]) << std::endl;
		}
		if ( cc.DecryptDouble(Rsk, (*ind2)[i]) != msg2) {
			std::cout << i << "'st double2 is not " << msg2 <<", but" << 
			cc.DecryptDouble(Rsk, (*ind2)[i]) << std::endl;
		}
	}
	std::cout << "Pre check is ends" << std::endl;


	
	// Add Double
	auto add_startd = std::chrono::system_clock::now();
	#pragma omp parallel for collapse(1)
	for (uint32_t i = 0; i < nums; i++) {
		(*outad)[i] = cc.Add((*ind1)[i], (*ind2)[i], (*OFUF_arr2)[i]);			
	} 
	std::chrono::duration<double> add_endd = std::chrono::system_clock::now() - add_startd;
	std::cout << "Double Add time is " <<add_endd.count() << std::endl ;  

	
	// Prouduct Double
	auto product_startd = std::chrono::system_clock::now();
	#pragma omp parallel for collapse(1)
	for (uint32_t i = 0; i < nums; i++) {
		(*outpd)[i] = cc.Product((*ind1)[i], (*ind2)[i], (*OFUF_arr2)[i]);	
	} 
	std::chrono::duration<double> product_endd = std::chrono::system_clock::now() - product_startd;
	std::cout << "Double product time is " <<product_endd.count() << std::endl ;  

	// Correctness check
	double ansd = 4.0;
	for (uint32_t i = 0; i < nums; i++) {
		if ( cc.DecryptDouble(Rsk, (*outad)[i]) != ansd) {
			std::cout << i << "'st double adding answer is not " << ansd <<", but" << 
			cc.DecryptDouble(Rsk, (*outad)[i]) << std::endl;
		}
		if ( cc.DecryptDouble(Rsk, (*outpd)[i]) != ansd) {
			std::cout << i << "'st double mul answer is not " << ansd <<", but" << 
			cc.DecryptDouble(Rsk, (*outpd)[i]) << std::endl;
		}
	}
	std::cout << "Decrypt check is ends" << std::endl;


/////////////////////////////////////////////////////////////////////

    //vector<FPRLWECiphertext> inf1(nums);
	//vector<FPRLWECiphertext> inf2(nums);
	//vector<FPRLWECiphertext> outaf(nums);
	//vector<FPRLWECiphertext> outpf(nums);
	//vector<std::shared_ptr<FPOFUF>> OFUF_arr(nums);
	std::shared_ptr<vector<FPRLWECiphertext>> inf1 = 
	 std::make_shared<vector<FPRLWECiphertext>> (nums);
	std::shared_ptr<vector<FPRLWECiphertext>> inf2 = 
	 std::make_shared<vector<FPRLWECiphertext>> (nums);
	std::shared_ptr<vector<FPRLWECiphertext>> outaf = 
	 std::make_shared<vector<FPRLWECiphertext>> (nums);
	std::shared_ptr<vector<FPRLWECiphertext>> outpf = 
	 std::make_shared<vector<FPRLWECiphertext>> (nums);
	std::shared_ptr<vector<std::shared_ptr<FPOFUF>>> OFUF_arr = 
	 std::make_shared<vector<std::shared_ptr<FPOFUF>>> (nums);
	
	
	//std::make_shared<vector<FPRLWECiphertext>> (inf1(nums));
	//std::make_shared<vector<FPRLWECiphertext>> (inf2(nums));
	//std::make_shared<vector<FPRLWECiphertext>> (outaf(nums));
	//std::make_shared<vector<FPRLWECiphertext>> (outpf(nums));
	//std::make_shared<vector<std::shared_ptr<FPOFUF>>> (OFUF_arr(nums));
	
	
	float msg = 2.0;	
	#pragma omp parallel for collapse(1)
	for (uint32_t i = 0; i < nums; i++) {
		(*inf1)[i] = cc.EncryptFloat(Rsk, msg);		
		(*inf2)[i] = cc.EncryptFloat(Rsk, msg);		
		(*OFUF_arr)[i] = cc.MakeOverFlowDetectCT();
	}
	std::cout << "Float Encryption is end" << std::endl ;  
	for (uint32_t i = 0; i < nums; i++) {
		if ( cc.DecryptFloat(Rsk, (*inf1)[i]) != msg) {
			std::cout << i << "'st float1 is not " << msg <<", but" << 
			cc.DecryptFloat(Rsk, (*inf1)[i]) << std::endl;
		}
		if ( cc.DecryptFloat(Rsk, (*inf2)[i]) != msg) {
			std::cout << i << "'st float2 is not " << msg <<", but" << 
			cc.DecryptFloat(Rsk, (*inf2)[i]) << std::endl;
		}
	}

	std::cout << "Pre check is ends" << std::endl;


	// Add Float
	auto add_startf = std::chrono::system_clock::now();
	#pragma omp parallel for collapse(1)
	for (uint32_t i = 0; i < nums; i++) {
		(*outaf)[i] = cc.Add((*inf1)[i], (*inf2)[i], (*OFUF_arr)[i]);		
	} 
	std::chrono::duration<double> add_endf = std::chrono::system_clock::now() - add_startf;
	std::cout << "Float Add time is " <<add_endf.count() << std::endl ;  

	
	// Prouduct Float
	auto product_startf = std::chrono::system_clock::now();
	#pragma omp parallel for collapse(1)
	for (uint32_t i = 0; i < nums; i++) {
		//std::cout << " index is " << i << std::endl;
		(*outpf)[i] = cc.Product((*inf1)[i], (*inf2)[i], (*OFUF_arr)[i]);		
	} 
	std::chrono::duration<double> product_endf = std::chrono::system_clock::now() - product_startf;
	std::cout << "Float Product time is " <<product_endf.count() << std::endl ;  

	// Correctness check
	float ans = 4.0;
	for (uint32_t i = 0; i < nums; i++) {
		if ( cc.DecryptFloat(Rsk, (*outaf)[i]) != ans) {
			std::cout << i << "'st float adding answer is not " << ans <<", but" << 
			cc.DecryptFloat(Rsk, (*outaf)[i]) << std::endl;
		}
		if ( cc.DecryptFloat(Rsk, (*outpf)[i]) != ans) {
			std::cout << i << "'st float mul answer is not " << ans <<", but" << 
			cc.DecryptFloat(Rsk, (*outpf)[i]) << std::endl;
		}
	}
	std::cout << "Decrypt check is ends" << std::endl;

	return 0;
}
