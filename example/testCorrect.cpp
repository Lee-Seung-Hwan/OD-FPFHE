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
	cout.precision(dbl::max_digits10);
	// MAYBE GINX
	auto cc = FPFHEContext();
	bool PARALLELs = true;
	//bool PARALLELs = false;
	//omp_set_num_threads(10);
	
  if (PARALLELs == false ){
		omp_set_num_threads(1);
	}
	cc.GenerateFHEContext(TOY, PARALLELs);
	//cc.GenerateFHEContext(STD128D, PARALLELs);
	//cc.GenerateFHEContext(STD160D, PARALLELs);
	//cc.GenerateFHEContext(STD192D, PARALLELs);
	//cc.GenerateFHEContext(STD128S, PARALLELs);
	//cc.GenerateFHEContext(STD160S, PARALLELs);
	//cc.GenerateFHEContext(STD192S, PARALLELs);



	// Make Secret Key
	//auto sk = cc.KeyGen();
	
	auto Rsk = cc.RKeyGen();
	cc.BTKeyGen(Rsk);
	
	auto OFUF = cc.MakeOverFlowDetectCT();
	//std::shared_ptr<FPOFUF> OFUF;
	
	
	std::cout << "Half================================================" << std::endl;
	
	float test_val1h = -0.0917635142362524246263723462452467245236235235325325;
	float test_val2h = 0.06246724724637586345346;
	
	float test_val3h = 24.235268723623733298623623723532446;
	float test_val4h = -5.324663335297272335623523446;
	
	auto ct1h = cc.EncryptHalf(Rsk, test_val1h);
	auto ct2h = cc.EncryptHalf(Rsk, test_val2h);
	auto ct3h = cc.EncryptHalf(Rsk, test_val3h);
	auto ct4h = cc.EncryptHalf(Rsk, test_val4h);
	
	float res1h = cc.DecryptHalf(Rsk, ct1h);
	float res2h = cc.DecryptHalf(Rsk, ct2h);
	float res3h = cc.DecryptHalf(Rsk, ct3h);
	float res4h = cc.DecryptHalf(Rsk, ct4h);
	std::cout << "Orign msg1 is             " << test_val1h << std::endl;
	std::cout << "Decrypt Half1 value is  " << res1h << std::endl;
	std::cout << "Orign msg2 is             " << test_val2h << std::endl;
	std::cout << "Decrypt Half2 value is  " << res2h << std::endl;
	std::cout << "Orign msg3 is             " << test_val3h << std::endl;
	std::cout << "Decrypt Half1 value is  " << res3h << std::endl;
	std::cout << "Orign msg4 is             " << test_val4h << std::endl;
	std::cout << "Decrypt Half2 value is  " << res4h << std::endl;
	
	std::cout << "================================================" << std::endl;
	
	
	
	std::cout << "Float================================================" << std::endl;
	

	float test_val1f = -0.0000000000000027915;
	float test_val2f =  0.0000000000083867;
	float test_val3f =  182634000000000.0;
	float test_val4f =      -6278270000.0;
	
	auto ct1f = cc.EncryptFloat(Rsk, test_val1f);
	auto ct2f = cc.EncryptFloat(Rsk, test_val2f);
	auto ct3f = cc.EncryptFloat(Rsk, test_val3f);
	auto ct4f = cc.EncryptFloat(Rsk, test_val4f);
	
	float res1f = cc.DecryptFloat(Rsk, ct1f);
	float res2f = cc.DecryptFloat(Rsk, ct2f);
	float res3f = cc.DecryptFloat(Rsk, ct3f);
	float res4f = cc.DecryptFloat(Rsk, ct4f);
	std::cout << "Orign msg1 is             " << test_val1f << std::endl;
	std::cout << "Decrypt float1  value is  " << res1f << std::endl;
	std::cout << "Orign msg2 is             " << test_val2f << std::endl;
	std::cout << "Decrypt float2  value is  " << res2f << std::endl;
	std::cout << "Orign msg3 is             " << test_val3f << std::endl;
	std::cout << "Decrypt float1  value is  " << res3f << std::endl;
	std::cout << "Orign msg4 is             " << test_val4f << std::endl;
	std::cout << "Decrypt float2  value is  " << res4f << std::endl;
	
	std::cout << "================================================" << std::endl;
	
	auto add_startf = std::chrono::system_clock::now();
	auto ct5f = cc.Add(ct1f, ct2f, OFUF);
	std::chrono::duration<double> add_endf = std::chrono::system_clock::now() - add_startf;
	std::cout << "Add time is " <<add_endf.count() << std::endl ;
	float ans5f = test_val1f + test_val2f;
	std::cout << "Orign Adding answer is	" << ans5f << std::endl;
	float res5f = cc.DecryptFloat(Rsk, ct5f);
	std::cout << "Decrypt Addin result is	" << res5f << std::endl;
	cc.CheckOFUF(Rsk,OFUF);
	


	auto sub_startf = std::chrono::system_clock::now();
	auto ct6f = cc.Sub(ct3f, ct4f, OFUF);
	std::chrono::duration<double> sub_endf = std::chrono::system_clock::now() - sub_startf;
	std::cout << "Sub time is " << sub_endf.count() << std::endl ;  
	float ans6f = test_val3f-test_val4f;
	std::cout << "Orign Sub answer is	" << ans6f << std::endl;
	float res6f = cc.DecryptFloat(Rsk, ct6f);
	std::cout << "Decrypt Sub is		" << res6f << std::endl;
	cc.CheckOFUF(Rsk,OFUF);
	


	auto product_startf = std::chrono::system_clock::now();
	auto ct7f = cc.Product(ct5f, ct6f, OFUF);
	std::chrono::duration<double> product_endf = std::chrono::system_clock::now() - product_startf;
	std::cout << "Product time is " <<product_endf.count() << std::endl ;  
	float ans7f = ans5f * ans6f;
	std::cout << "Orign Prod answer is 	" << ans7f << std::endl;
	float res7f = cc.DecryptFloat(Rsk, ct7f);
	std::cout << "Decrypt Prod Result is  " << res7f << std::endl;
	cc.CheckOFUF(Rsk,OFUF);
	
	auto product_startf2 = std::chrono::system_clock::now();
	auto ct8f = cc.Product(ct7f, ct7f, OFUF);
	std::chrono::duration<double> product_endf2 = std::chrono::system_clock::now() - product_startf2;
	std::cout << "Product time^2 is " <<product_endf2.count() << std::endl ;  
	float ans8f = ans7f * ans7f;
	std::cout << "Orign Pr^2 answer is 	" << ans8f << std::endl;
	float res8f = cc.DecryptFloat(Rsk, ct8f);
	std::cout << "Decrypt Pr^2 Result is  " << res8f << std::endl;
	cc.CheckOFUF(Rsk,OFUF);
	






	
	std::cout << "================================================" << std::endl;
	
								
	double test_val1 = -0.00000000000000000000000000000009176351423625429;
	double test_val2 =  0.0000000000000000000000062467247246375865;
	double test_val3 = 24523526872362373000000.0;
	double test_val4 =     -543246633352972740.0;
	
	auto ct1 = cc.EncryptDouble(Rsk, test_val1);
	auto ct2 = cc.EncryptDouble(Rsk, test_val2);
	auto ct3 = cc.EncryptDouble(Rsk, test_val3);
	auto ct4 = cc.EncryptDouble(Rsk, test_val4);
	
	double res1 = cc.DecryptDouble(Rsk, ct1);
	double res2 = cc.DecryptDouble(Rsk, ct2);
	double res3 = cc.DecryptDouble(Rsk, ct3);
	double res4 = cc.DecryptDouble(Rsk, ct4);
	std::cout << "Orign msg1 is             " << test_val1 << std::endl;
	std::cout << "Decrypt double1 value is  " << res1 << std::endl;
	std::cout << "Orign msg2 is             " << test_val2 << std::endl;
	std::cout << "Decrypt double2 value is  " << res2 << std::endl;
	std::cout << "Orign msg3 is             " << test_val3 << std::endl;
	std::cout << "Decrypt double1 value is  " << res3 << std::endl;
	std::cout << "Orign msg4 is             " << test_val4 << std::endl;
	std::cout << "Decrypt double2 value is  " << res4 << std::endl;
	
	std::cout << "================================================" << std::endl;
	
	auto add_start = std::chrono::system_clock::now();
	auto ct5 = cc.Add(ct1, ct2, OFUF);
	std::chrono::duration<double> add_end = std::chrono::system_clock::now() - add_start;
	std::cout << "Add time is " <<add_end.count() << std::endl ;
	double ans5 = test_val1 + test_val2;
	std::cout << "Orign Adding answer is	" << ans5 << std::endl;
	double res5 = cc.DecryptDouble(Rsk, ct5);
	std::cout << "Decrypt Addin result is	" << res5 << std::endl;
	cc.CheckOFUF(Rsk,OFUF);
	

	auto sub_start = std::chrono::system_clock::now();
	auto ct6 = cc.Sub(ct3, ct4, OFUF);
	std::chrono::duration<double> sub_end = std::chrono::system_clock::now() - sub_start;
	std::cout << "Sub time is " << sub_end.count() << std::endl ;  
	double ans6 = test_val3-test_val4;
	std::cout << "Orign Sub answer is	" << ans6 << std::endl;
	double res6 = cc.DecryptDouble(Rsk, ct6);
	std::cout << "Decrypt Sub is		" << res6 << std::endl;
	cc.CheckOFUF(Rsk,OFUF);
	


	auto product_start = std::chrono::system_clock::now();
	auto ct7 = cc.Product(ct5, ct6, OFUF);
	std::chrono::duration<double> product_end = std::chrono::system_clock::now() - product_start;
	std::cout << "Product time is " <<product_end.count() << std::endl ;  
	double ans7 = ans5 * ans6;
	std::cout << "Orign Prod answer is 	" << ans7 << std::endl;
	double res7 = cc.DecryptDouble(Rsk, ct7);
	std::cout << "Decrypt Prod Result is  " << res7 << std::endl;
	cc.CheckOFUF(Rsk,OFUF);
	

	auto product_start2 = std::chrono::system_clock::now();
	auto ct8 = cc.Product(ct7, ct7, OFUF);
	std::chrono::duration<double> product_end2 = std::chrono::system_clock::now() - product_start2;
	std::cout << "Product time^2 is " <<product_end2.count() << std::endl ;  
	double ans8 = ans7 * ans7;
	std::cout << "Orign Pr^2 answer is 	" << ans8 << std::endl;
	double res8 = cc.DecryptDouble(Rsk, ct8);
	std::cout << "Decrypt Pr^2 Result is  " << res8 << std::endl;
	cc.CheckOFUF(Rsk,OFUF);
	
	std::cout << "================================================" << std::endl;
	std::cout << "====OF UF CHECK Code!!!!!!!================================================" << std::endl;
	std::cout << "Float================================================" << std::endl;
	
	
	

	float max_f = FLT_MAX /4;
	//float max_mf = - max_f;
	float min_f = FLT_MIN;
	//max_f /= 4;
	//min_f *= 1;
	float one_f =   1000.0;
	float one_overf = 1.0/1000;
	
	/*
	std::cout << "Max Float  is             " << max_f << std::endl;
	std::cout << "Maxm Float  is             " << max_mf << std::endl;
	std::cout << "Min Float  is             " << min_f << std::endl;
	*/
	
	auto ct_maxf = cc.EncryptFloat(Rsk, max_f);
	auto ct_minf = cc.EncryptFloat(Rsk, min_f);
	auto ct_onef = cc.EncryptFloat(Rsk, one_f);
	auto ct_oneoverf = cc.EncryptFloat(Rsk, one_overf);


	float res_maxf = cc.DecryptFloat(Rsk, ct_maxf);
	float res_minf = cc.DecryptFloat(Rsk, ct_minf);
	float res_maxf_prod = res_maxf * 1000.0;
	float res_minf_prod = res_minf / 1000;


	//float res_onef = cc.DecryptFloat(Rsk, ct_onef);
	//float res_oneoverf = cc.DecryptFloat(Rsk, ct_oneoverf);
	std::cout << "Orign max float is        " << max_f << std::endl;
	std::cout << "Decrypt max float is      " << res_maxf << std::endl;
	std::cout << "Product max * 1000 is      " << res_maxf_prod << std::endl;
	/*
	std::cout << "Orign one_f is             " << one_f << std::endl;
	std::cout << "Decrypt one_f value is  " << res_onef << std::endl;
	std::cout << "Orign one_overf is             " << one_overf << std::endl;
	std::cout << "Decrypt one_overf value is  " << res_oneoverf << std::endl;
	*/
	auto ct5_OF1 = cc.Product(ct_maxf, ct_onef, OFUF);
	float max_dec = 	cc.DecryptFloat(Rsk, ct5_OF1);
	
	std::cout << " Max * 1000 Decrypt: is    " << 
		max_dec << std::endl;
	cc.CheckOFUF(Rsk,OFUF);
	
	std::cout << "Orign min float is        " << min_f << std::endl;
	std::cout << "Decrypt min float is      "	<< res_minf << std::endl;
	std::cout << "Product min * 1/1000 is    " << res_minf_prod << std::endl;
	

	auto ct5_UF2 = cc.Product(ct_minf, ct_oneoverf, OFUF);
	float min_dec = 	cc.DecryptFloat(Rsk, ct5_UF2);
	 
	std::cout << "Min * (1/1000) Decrypt: is    " << 
	min_dec << std::endl;
	cc.CheckOFUF(Rsk,OFUF);
	

	
	std::cout << "DOUBLe================================================" << std::endl;
	

	double max_d = DBL_MAX /4;
	//double max_md = - max_d;
	double min_d = DBL_MIN;
	//max_f /= 4;
	//min_f *= 1;
	double one_d =   1000.0;
	double one_overd = 1.0/1000;
	
	/*
	std::cout << "Max Double is             " << max_d << std::endl;
	std::cout << "Maxm Double  is             " << max_md << std::endl;
	std::cout << "Min Double  is             " << min_d << std::endl;
	*/
	
	auto ct_maxd = cc.EncryptDouble(Rsk, max_d);
	auto ct_mind = cc.EncryptDouble(Rsk, min_d);
	auto ct_oned = cc.EncryptDouble(Rsk, one_d);
	auto ct_oneoverd = cc.EncryptDouble(Rsk, one_overd);


	double res_maxd = cc.DecryptDouble(Rsk, ct_maxd);
  double res_mind = cc.DecryptDouble(Rsk, ct_mind);
	double res_maxd_prod = res_maxd * 1000;
	double res_mind_prod = res_mind / 1000;


	//double res_oned = cc.DecryptDouble(Rsk, ct_oned);
	//double res_oneoverd = cc.DecryptDouble(Rsk, ct_oneoverd);
	std::cout << "Orign max double is       " << max_d << std::endl;
	std::cout << "Decrypt max double is     " << res_maxd << std::endl;
	std::cout << "Product max * 100 is      " << res_maxd_prod << std::endl;
	
	/*
	std::cout << "Orign one_d is             " << one_d << std::endl;
	std::cout << "Decrypt one_d value is  " << res_oned << std::endl;
	std::cout << "Orign one_overd is             " << one_overd << std::endl;
	std::cout << "Decrypt one_overd value is  " << res_oneoverd << std::endl;
	*/
	auto ct5_OF1d = cc.Product(ct_maxd, ct_oned, OFUF);
	double max_resd = cc.DecryptDouble(Rsk, ct5_OF1d);
	
	std::cout << "Decrypt Max * 1000: OF is  " << 
		max_resd << std::endl;
		cc.CheckOFUF(Rsk,OFUF);
	

	std::cout << "Orign min double is       " << min_d << std::endl;
	std::cout << "Decrypt min double is     " << res_mind << std::endl;
	std::cout << "Product min * 1/1000 is   " << res_mind_prod << std::endl;
	
	auto ct5_UF2d = cc.Product(ct_mind, ct_oneoverd, OFUF);
	double min_resd = cc.DecryptDouble(Rsk, ct5_UF2d);
	
	std::cout << "Decrypt Min * (1/1000), is  " << 
		min_resd << std::endl;
	cc.CheckOFUF(Rsk,OFUF);
	


	//double res4 = cc.DecryptDouble(Rsk, ct4);
	//std::cout << "Relu 2 value is  " << res4 << std::endl;

	//uint64_t two = 1234567;
	//auto ct2 = cc.EncryptUInt64(Rsk, two);
	
	//auto ct_mul1 = cc.ProductUOrign(ct1, ct1);
	//auto ct_mul_only = cc.ProductOnly(ct1, ct1);
	
	//auto ct_mul2 = cc.ProductUOrign(ct1, ct2);
	//auto ct_mul3 = cc.ProductUOrign(ct_mul1, ct_mul2);


	//uint64_t res_msg1 = cc.DecryptUInt64(Rsk, ct_mul1);
	//uint64_t res_msg2 = cc.DecryptUInt64(Rsk, ct_mul2);
	//uint64_t res_msg3 = cc.DecryptUInt64(Rsk, ct_mul3);

	//std::cout<< "CT mul decrypt 1 is " <<std::endl << res_msg1 << std::endl;
	//std::cout<< "CT mul decrypt 1 is " <<std::endl << res_msg2 << std::endl;
	//std::cout<< "CT mul decrypt 1 is " <<std::endl << res_msg3 << std::endl;
	
	
	//std::cout << "================== Key Gen is End!!!!!==================" << std::endl;

	// Orign Test
	// 64 bit Encryption
	
	
	//uint64_t one = 1;
	//auto ct1 = cc.EncryptUInt64(Rsk, one);
	//uint64_t two = 1234567;
	//auto ct2 = cc.EncryptUInt64(Rsk, two);
	//auto ct_mul = cc.ProductOnly(ct1, ct2);
	//uint64_t res_msg1 = cc.DecryptUInt64(Rsk, ct2);
	//std::cout<< "CT mul decrypt 1 is " <<std::endl << res_msg1 << std::endl;
	
	


	/*
	uint64_t one = 1;
	auto ct1 = cc.EncryptUInt64(Rsk, one);
	uint64_t two = 1234567;
	auto ct2 = cc.EncryptUInt64(Rsk, two);
	
	auto ct_mul1 = cc.ProductUOrign(ct1, ct1);
	//auto ct_mul_only = cc.ProductOnly(ct1, ct1);
	
	auto ct_mul2 = cc.ProductUOrign(ct1, ct2);
	auto ct_mul3 = cc.ProductUOrign(ct_mul1, ct_mul2);


	uint64_t res_msg1 = cc.DecryptUInt64(Rsk, ct_mul1);
	uint64_t res_msg2 = cc.DecryptUInt64(Rsk, ct_mul2);
	uint64_t res_msg3 = cc.DecryptUInt64(Rsk, ct_mul3);

	std::cout<< "CT mul decrypt 1 is " <<std::endl << res_msg1 << std::endl;
	std::cout<< "CT mul decrypt 1 is " <<std::endl << res_msg2 << std::endl;
	std::cout<< "CT mul decrypt 1 is " <<std::endl << res_msg3 << std::endl;
	*/
	
	//uint64_t res_msg0 = cc.DecryptUInt64(Rsk, ct1);
	//std::cout<< "CT mul decrypt res is " <<std::endl << res_msg0 << std::endl;
	//uint64_t res_msg11 = cc.DecryptUInt64(Rsk, ct_mul_only, true);
	//std::cout<< "CT mul only is " <<std::endl << res_msg11 << std::endl;

	//auto ct_mul2 = cc.ProductUOrign(ct_mul1, ct_mul1);
	//uint64_t res_msg2 = cc.DecryptUInt64(Rsk, ct_mul1);
	//std::cout<< "CT mul decrypt 1 is " <<std::endl << res_msg2 << std::endl;
	

	//std::cout<< "CT mul decrypt 1234567 is " <<std::endl << res_msg2 << std::endl;
	//std::cout<< "CT mul decrypt 1234567 is " <<std::endl << res_msg3 << std::endl;
	

	// Boot Orign
	/*
	std::cout << "================================================" << std::endl;
	uint64_t one = 12345678;
	auto ct1 = cc.EncryptUInt64(Rsk,one);
		
	uint64_t two = 1;
	auto ct2 = cc.EncryptUInt64(Rsk,two);
	
	auto ct_boot = cc.ProductUOrign(ct1,ct2);
	//auto ct_boot2 = cc.BootstrapOrign(ct_boot);
	uint64_t res_msg = cc.DecryptUInt64(Rsk, ct_boot);
	std::cout<< "CT 1 decrypt is " <<std::endl << res_msg << std::endl;
	*/
	
	//auto cc_mul = cc.Product(ct1, ct2);
	//uint64_t res_msg = cc.DecryptUInt64(Rsk,cc_mul);
	//std::cout<< "CT 1 decrypt is " <<std::endl << res_msg << std::endl;
	
	//int64_t res_msg_64_add = cc.DecryptInt64(Rsk,ct_64_1);
	//int64_t res_msg_64_add2 = cc.DecryptInt64(Rsk,ct_64_1);
	//std::cout << "Decrypt Result is " << res_msg_64_add << std::endl;
	//std::cout << "Decrypt Result is " << res_msg_64_add2 << std::endl;
	



	/*
	auto ct_64_1_bt = cc.BootstrapOrign(cc.Product(ct_64_1,ct_64_1));
	std::cout<< "CT 1 decrypt is " <<std::endl << cc.DecryptUInt64(Rsk, ct_64_1_bt)<< std::endl;
	
	uint64_t msg_64 = 12345678;
	auto ct_64 = cc.EncryptUInt64(Rsk, msg_64);	
	auto ct_64_2_bt = cc.BootstrapOrign(cc.Product(ct_64,ct_64_1_bt));
	std::cout<< "CT 12345678 decrypt is " <<std::endl << cc.DecryptUInt64(Rsk, ct_64_2_bt)<< std::endl;
	for (uint32_t i = 0 ; i < 10; i++) {
		ct_64_2_bt = cc.BootstrapOrign(cc.Product(ct_64_1_bt,ct_64_2_bt));
	std::cout << i << "'st bootstrap res :" 
		<< "CT 12345678 decrypt is " <<std::endl << cc.DecryptUInt64(Rsk, ct_64_2_bt)<< std::endl;
	}
	
	*/
	// Bootstrap TEST Arithmetic 
	
	
	/*
	
	std::cout << "================================================" << std::endl;
	

	
	int64_t msg64_1 = - ((1 << 13) + 1622);
	auto ct64_1 = cc.EncryptInt64(Rsk, msg64_1);	
	int64_t msg64_2 = (1 << 12) - 994;
	auto ct64_2 = cc.EncryptInt64(Rsk, msg64_2);	
	

	std::cout 
		<< "Msg 1 is                " << (msg64_1) << std::endl
		<< "Msg 2 is                " << (msg64_2) << std::endl
		<< "Add correct answer is   " << (msg64_1  +  msg64_2) << std::endl;
	auto ct_add = cc.Add(ct64_1,ct64_2, Rsk);	
	int64_t res_msg_64_add = cc.DecryptInt64(Rsk,ct_add);
	std::cout 
		<< "an adding ciphertext is " << res_msg_64_add << std::endl;

	std::cout << "================================================" << std::endl;
	
	std::cout 
		<< "Msg 1 is                " << (msg64_1) << std::endl
		<< "Msg 2 is                " << (msg64_2) << std::endl
		<< "Sub correct answer is   " << (msg64_1  -  msg64_2) << std::endl;
	auto ct_sub = cc.Sub(ct64_1,ct64_2, Rsk);	
	int64_t res_msg_64_sub = cc.DecryptInt64(Rsk,ct_sub);
	std::cout 
		<< "an subs ciphertext is	" << res_msg_64_sub << std::endl;
	
	std::cout << "================================================" << std::endl;



	// First Product
	int64_t mul_res = (msg64_1 * msg64_1) - (msg64_2 * msg64_2);
	std::cout 
		<< "Msg 1 is                " << (msg64_1 + msg64_2) << std::endl
		<< "Msg 2 is                " << (msg64_1 - msg64_2) << std::endl
		<< "Prod correct answer is  " << mul_res << std::endl;
	


	auto ct_mul = cc.Product(ct_add,ct_sub);	
	
	// TMP
	// int64_t msg64_11 = - 6712;
	// auto ct64_11 = cc.EncryptInt64(Rsk, msg64_11);	
	// int64_t msg64_22 = -12916;
	// auto ct64_22 = cc.EncryptInt64(Rsk, msg64_22);
	// auto ct_mul = cc.Product(ct64_11,ct64_22);	
	
	int64_t res_msg_64 = cc.DecryptInt64(Rsk,ct_mul);
	std::cout
		<< "Product of ciphertex is	" << res_msg_64 << std::endl;
		
	std::cout << "================================================" << std::endl;

	
	// TMP2
	// Second Product
	std::cout 
		<< "Msg is		        " << mul_res << std::endl
		<< "Pow correct answer is   " << mul_res * mul_res << std::endl;
	
	
	
	// int64_t msg64_111 = 86692192;
	// auto ct64_111 = cc.EncryptInt64(Rsk, msg64_111);	
	
	//auto ct_mul2 = cc.Product(ct64_111,ct64_111);	
	
	auto ct_mul2 = cc.Product(ct_mul,ct_mul);	
	int64_t res_msg_64_2 = cc.DecryptInt64(Rsk,ct_mul2);
	std::cout 
		<< "Pow ciphertext is	" <<  res_msg_64_2 << std::endl;
	
	std::cout << "================================================" << std::endl;

	*/

	//Parallel
	/*
	uint32_t N_FP = (1 <<12);
  vector<NativeInteger> moduliQ(2);
  vector<NativeInteger> rootsQ(2);
	vector<shared_ptr<ILNativeParams>> tmp_params(2);
	
	//moduliQ[0] = 44965627529527297;
  //moduliQ[1] = 10977936408577;	
	//moduliQ[0] = 1952732650930177;
  //moduliQ[1] = 476741369857;
	moduliQ[0] = 1811939329;
  moduliQ[1] = 1769473;
  rootsQ[0] = RootOfUnity<NativeInteger>(N_FP*2,moduliQ[0]);
  rootsQ[1] = RootOfUnity<NativeInteger>(N_FP*2,moduliQ[1]);
  ChineseRemainderTransformFTT<NativeVector>::PreCompute(rootsQ, 2*N_FP, moduliQ);
	
	tmp_params[0] = std::make_shared<ILNativeParams>( 2 * N_FP, moduliQ[0], rootsQ[0]);
  tmp_params[1] = std::make_shared<ILNativeParams>( 2 * N_FP, moduliQ[1], rootsQ[1]);
   
	const std::shared_ptr<ILDCRTParams<BigInteger>> DCRTParamset = std::make_shared     <ILDCRTParams<BigInteger>>(N_FP*2,moduliQ,rootsQ);
	DCRTPoly::DugType dug;
	DCRTPoly::TugType tug;
	//NativePoly::DugType dug_nat;
	DiscreteUniformGeneratorImpl<NativeVector> dug_nat;
	dug_nat.SetModulus(moduliQ[0]);

	const std::shared_ptr<ILNativeParams> NativeParamset = std::make_shared<ILNativeParams>(N_FP*2,moduliQ[0],rootsQ[0]);
	
	uint32_t nums = 50000;
	vector<DCRTPoly>		DC_vec(nums);
	vector<DCRTPoly>		DC_vec_OZ(nums);
	vector<NativePoly>	NT_vec(nums);
	vector<NativePoly>	NT_vec_OZ(nums);
	//#pragma omp parallel for
	{for (uint32_t i = 0 ; i <nums; i++) {
		DC_vec[i] = DCRTPoly(dug, DCRTParamset, Format::COEFFICIENT);
		DC_vec_OZ[i] = DC_vec[i];
		NT_vec[i] = NativePoly(dug_nat, NativeParamset, Format::COEFFICIENT);
		NT_vec_OZ[i] = NT_vec[i];	
	}}
	std::cout << "Gen is end" << std::endl;
	#pragma omp parallel for
	{for (uint32_t i = 0; i < nums; i++) {
		DC_vec[i].SetFormat(Format::EVALUATION);
		//DC_vec[i].SetFormat(Format::COEFFICIENT);
	
	}}
	#pragma omp parallel for
	{for (uint32_t i = 0; i < nums; i++) {
		NT_vec[i].SetFormat(Format::EVALUATION);
		//NT_vec[i].SetFormat(Format::COEFFICIENT);
	}}
	
	uint32_t correct_DC = 0;
	uint32_t correct_NT = 0;

	for (uint32_t i = 0; i < nums; i++) {
		DC_vec[i].SetFormat(Format::COEFFICIENT);
		NT_vec[i].SetFormat(Format::COEFFICIENT);


		correct_DC += ( DC_vec_OZ[i] == DC_vec[i]) ? 1 : 0;
		correct_NT += ( NT_vec_OZ[i] == NT_vec[i]) ? 1 : 0;
	}

	std::cout << " DC count is " << correct_DC << ", NT count is " << correct_NT << std::endl;
	*/

	/*
	DCRTPoly tmp = DCRTPoly(dug, DCRTParamset, Format::COEFFICIENT);
	DCRTPoly tmp2 =  DCRTPoly(DCRTParamset, Format::COEFFICIENT, true);
 	
	clock_t startd1, endd1, startd2, endd2;
	startd1 = clock();
	
	DCRTPoly rounding = cc.DivideRounding(&tmp, Format::COEFFICIENT);
	endd1 = clock();
	startd2 = clock();
	for (uint32_t i = 0; i < N_FP; i++) {
		tmp2.GetElementW(0)[i]= tmp[i].DivideAndRound(BigInteger(moduliQ[1])).ConvertToInt();
	}
	endd2 = clock();
	

	
	std::cout << rounding << std::endl;
	std::cout << tmp2 << std::endl;
	std::cout << "Nat Time is " << std::endl;
	printf("%f\n ", (double) (endd1 - startd1));
	std::cout << ", Big Time is " << std::endl;
	printf("%f\n", (double) (endd2 - startd2));
	*/


	
	// GD Test
	/*
	uint32_t base = 5;
	uint32_t rm = 3;
	
	DCRTPoly tmp = DCRTPoly(dug, DCRTParamset, Format::COEFFICIENT);
	tmp.SetFormat(Format::EVALUATION);
	
	DCRTPoly SK = DCRTPoly(tug, DCRTParamset, Format::EVALUATION);
	vector<DCRTPoly> SK_vec = SK.PowersOfBase(base);
	DCRTPoly result = tmp * SK;
	
	tmp.SetFormat(Format::COEFFICIENT);
	// Origns
	
	clock_t startg1, endg1, startg2, endg2;

	startg1 = clock();
	vector<DCRTPoly> Uniform_vec1 = tmp.BaseDecompose(base, true);
	endg1 = clock();

	startg2 = clock();
	vector<DCRTPoly> Uniform_vec2 = cc.SGDLevel2(tmp, base, 0, Format::EVALUATION);
	endg2 = clock();
	
	std::cout << "Big Time is " << std::endl;
	printf("%f\n ", (double) (endg1 - startg1));
	std::cout << ", Nat Time is " << std::endl;
	printf("%f\n", (double) (endg2 - startg2));


	// Checking
	DCRTPoly res1 = SK_vec[0] * Uniform_vec1[0];
	for (uint32_t i = 1; i < Uniform_vec1.size(); i++) {
		res1 += SK_vec[i] * Uniform_vec1[i];
	}
	DCRTPoly res2 = SK_vec[0] * Uniform_vec2[0];
	for (uint32_t i = 1; i < Uniform_vec2.size(); i++) {
		res2 += SK_vec[i] * Uniform_vec2[i];
	}

	bool equal_res1 = (result == res1);
	bool equal_res2 = (result == res2);

	std::cout << "res1 is " << (equal_res1) << ", res2 is "<< equal_res2 << std::endl; 
	
	vector<DCRTPoly> Uniform_vec3 = cc.SGDLevel2(tmp, base, rm, Format::EVALUATION);
	bool equal_res3 = true;
	for (uint32_t i = 0; i < Uniform_vec3.size(); i++) {
		equal_res3 &= (Uniform_vec3[i] == Uniform_vec2[i+rm]);
	}
	std::cout << "Remove res are " <<  equal_res3 << std::endl;
	*/

	/*
	std::cout << rounding << std::endl;
	std::cout << tmp2 << std::endl;
	std::cout << "Nat Time is " << std::endl;
	printf("%f\n ", (double) (endd1 - startd1));
	std::cout << ", Big Time is " << std::endl;
	printf("%f\n", (double) (endd2 - startd2));
	*/

	//NativePoly & aa = tmp.GetElementW(0);
	//tmp.GetElementW(0)[1] = 1;
	//std::cout << tmp << std::endl;
	//std::cout << tmp[1] << std::endl;
	


	//std::shared_ptr<FPDCRTParams> Params = std::make_shared<FPDCRTParams> (FPDCRTParams(tmp_params,moduliQ,N_FP)); 
	
	// Power of Basis & GD Test 
	/*
	uint32_t rm_num = 4;
	uint32_t basis = 5;
	
	Format test_fm = Format::EVALUATION;
	//Format test_fm = Format::COEFFICIENT;


	FPDCRTPoly SK(DISTType::Tug ,  Params, test_fm); 
	std::cout << "SK is " << std::endl << SK << std::endl;
	FPDCRTPoly Uniform(DISTType::Dug ,  Params, test_fm); 
	std::cout << "Uniform is " << std::endl << Uniform << std::endl;

	// Make Locals...
	
	vector<unique_ptr<FPDCRTPoly>> SK_vec1 = SK.PowerOfBase(basis, 0, test_fm);
	std::cout << "Power Is Over1" <<  std::endl;
	vector<unique_ptr<FPDCRTPoly>> SK_vec2 = SK.PowerOfBase(basis, rm_num, test_fm);
	std::cout << "Power Is Over2" <<  std::endl;

	// PLZ...
	uint64_t Tot_nums = (Params->GetQBit() / basis);
	if (Params->GetQBit() % basis > 1) {Tot_nums++;}
	// Takes vals
	vector<unique_ptr<FPDCRTPoly>> Uniform_GD1(Tot_nums);
	for (uint32_t i = 0; i < Tot_nums; i++) {
		Uniform_GD1[i] = make_unique<FPDCRTPoly>(Params, Format::COEFFICIENT, true);
	}

	Uniform.SignedGadgetDecompositionLevel2(basis, 0, &Uniform_GD1, test_fm);
	vector<unique_ptr<FPDCRTPoly>> Uniform_GD2 = Uniform.SignedGadgetDecompositionLevel2(basis, rm_num, test_fm);
	vector<unique_ptr<FPDCRTPoly>> Uniform_GD3 = Uniform.SignedGadgetDecompositionLevel2(basis, 0, test_fm);
	
	std::cout << "SGD Is Over" <<  std::endl;


	SK.SetFormat(Format::EVALUATION);
	Uniform.SetFormat(Format::EVALUATION);
	FPDCRTPoly muls = SK * Uniform;
	
	SK.SetFormat(test_fm);
	Uniform.SetFormat(test_fm);

	// Product Preparation
	SK_vec1[0]->SetFormat(Format::EVALUATION);
	Uniform_GD1[0]->SetFormat(Format::EVALUATION);
	Uniform_GD3[0]->SetFormat(Format::EVALUATION);
	

	//FPDCRTPoly mul_chunk; 
	//((*SK_vec1[0])*(*Uniform_GD1[0]));
	FPDCRTPoly mul_chunk = ((*SK_vec1[0])*(*Uniform_GD1[0]));
	SK_vec1[0]->SetFormat(test_fm);
	Uniform_GD1[0]->SetFormat(test_fm);

	for (uint32_t i = 1; i < SK_vec1.size(); i++) {
		SK_vec1[i]->SetFormat(Format::EVALUATION);
		Uniform_GD1[i]->SetFormat(Format::EVALUATION);
		Uniform_GD3[i]->SetFormat(Format::EVALUATION);
		
		mul_chunk += ((*SK_vec1[i])*(*Uniform_GD1[i]));
	
		SK_vec1[i]->SetFormat(test_fm);
		Uniform_GD1[i]->SetFormat(test_fm);
	}

	
	bool isEqual = (muls == mul_chunk);
	if (isEqual) {
		std::cout <<"SGD is Success!!!!!!!" << std::endl;
	} else {
		std::cout <<"SGD is Fail!!!!!!!!!!" << std::endl;
	}

	std::cout << "Numbers are " << SK_vec1.size() << " and " << Uniform_GD1.size() << std::endl;

	bool GD_bool = true;
	bool PW_bool = true;
	for (uint32_t i = 0; i < Uniform_GD2.size(); i++) {
		GD_bool = GD_bool && (Uniform_GD1[i+rm_num] == Uniform_GD2[i]);
		PW_bool = PW_bool && (SK_vec1[i+rm_num] == SK_vec2[i]);
	}
	
	std::cout << "Again SK is " << std::endl << SK << std::endl;
	std::cout << "Again Uniform is " << std::endl << Uniform << std::endl;

	*/



	/* GD Test */

	//uint32_t rm_num = 4;
	/*
	uint32_t basis = 5;
	
	Format test_fm = Format::COEFFICIENT;
	
	FPDCRTPoly Uniform(DISTType::Dug ,  Params, test_fm); 
	std::cout << "Uniform is " << std::endl << Uniform << std::endl;

	vector<FPDCRTPoly> Uniform_GD1 = Uniform.SignedGadgetDecompositionLevel2(basis, 0, test_fm);
	
	
	FPDCRTPoly Add_vals; 
	Add_vals = Uniform_GD1[Uniform_GD1.size() - 1];

	for (int32_t i = Uniform_GD1.size() - 2;  i >= 0; i--) {
		Add_vals *= (1 << basis);
		Add_vals += Uniform_GD1[i];
	}
	

	std::cout << "Orign is " << std::endl;
	std::cout << Uniform << std::endl;
	std::cout << " Reconstruction is " << std::endl;
	std::cout << Add_vals << std::endl;
	
	

	bool isEqual = (Uniform == Add_vals);

	if (isEqual) {
		std::cout <<"SGD is Success!!!!!!!" << std::endl;
	} else {
		std::cout <<"SGD is Fail!!!!!!!!!!" << std::endl;
	}
	*/

	/* Scalar Multiplication Test */
	/*
	FPDCRTPoly Tmp_Poly1(DISTType::Tug ,  Params, Format::EVALUATION); 
	std::cout << Tmp_Poly1 << std::endl;
	
	Tmp_Poly1 *= 3;
	std::cout << Tmp_Poly1 << std::endl;
	*/
	/* Scalar Multiplication Test End */


	//std::cout <<"Q depth is " << Params->GetQDepth() << std::endl;
 
	//std::shared_ptr<FPDCRTPoly> Tmp_Poly1 = std::make_shared<FPDCRTPoly> (FPDCRTPoly(DISTType::Dgg,1.0 ,  Params, Format::COEFFICIENT)); 
		//std::shared_ptr<FPDCRTPoly> Tmp_Poly2 = std::make_shared<FPDCRTPoly> (FPDCRTPoly(DISTType::Dgg,1.0 ,  Params, Format::COEFFICIENT)); 
	

 
	//std::shared_ptr<FPDCRTPoly> Tmp_Poly1 = std::make_shared<FPDCRTPoly> (FPDCRTPoly(DISTType::Dgg,1.0 ,  Params, Format::COEFFICIENT)); 
	//std::shared_ptr<FPDCRTPoly> Tmp_Poly2 = std::make_shared<FPDCRTPoly> (FPDCRTPoly(DISTType::Dgg,1.0 ,  Params, Format::COEFFICIENT)); 
	
  /* Rounding Test */
	/*
	FPDCRTPoly Tmp_Poly1(DISTType::Dug ,  Params, Format::COEFFICIENT); 
	uint32_t Nums = 1 << 6;
	clock_t start1, end1, start2, end2;
	start1 = clock();
	for (uint32_t i = 0; i < Nums; i++) {
		FPDCRTPoly &tmps1 = Tmp_Poly1.MakeDivideRoundingBig(1);
		tmps1 += tmps1;
	}
	end1 = clock();
	start2 = clock();
	for (uint32_t i = 0; i < Nums; i++) {
		FPDCRTPoly &tmps2 = Tmp_Poly1.MakeDivideRounding(1);
		tmps2 += tmps2;
	}
	end2 = clock();
	std::cout << "Big Time is " << std::endl;
	printf("%f\n ", (double) (end1 - start1));
	std::cout << ", Nat Time is " << std::endl;
	printf("%f\n", (double) (end2 - start2));
	*/
	/* Rounding Test End*/


/*
	std::cout << "Big Dividings is " << std::endl  
		<< tmps1 << std::endl
		<< "Nat Dividings is " << std::endl  
		<< tmps2 << std::endl;
	*/
	//Tmp_Poly3 = Tmp_Poly2 + Tmp_Poly1;
	//Tmp_Poly3.SetFormat(Format::EVALUATION);
	//std::cout<< "After Ops, Tmp_Poly3 is " << std::endl << Tmp_Poly3 << std::endl;
	//std::cout<< "After Ops, Tmp_Poly1 is " << std::endl << Tmp_Poly1 << std::endl;
	//std::cout<< "After Ops, Tmp_Poly2 is " << std::endl << Tmp_Poly2 << std::endl;
	
	//Tmp_Poly1 += Tmp_Poly2;
	//FPDCRTPoly Tmp_Poly3 = Tmp_Poly1 + Tmp_Poly2;
	//std::cout << Tmp_Poly->GetElementAtIndex(0) << std::endl;
	//std::cout << Tmp_Poly->GetElementAtIndex(1) << std::endl;
	





	// = is ....







	// C++ 코딩 테스트
	
	/*
	std::cout << "before fix it" << std::endl;
	std::cout << Tmp_Poly->GetElementAtIndex(0) << std::endl;

	// 새로운 객체 생성
	NativePoly tmps = (Tmp_Poly->GetElementAtIndex(0));
	
	// Pointer를 빼오고 원래 값은 소멸
	//NativePoly tmps = std::move((*Tmp_Poly)[0]);
	std::cout << "transfered tmps is " << std::endl;
	std::cout << tmps << std::endl;
	
	
	tmps[0] = 0;
	std::cout << "Before TmpPoly Address addrss is " << std::endl;
	std::cout << &(Tmp_Poly->GetElementAtIndex(0)) << std::endl;
	std::cout << "Tems Address addrss is " << std::endl;
	std::cout << &(tmps) << std::endl;
	
	//직접 수정 & 값 복사하는 것임
	//(*Tmp_Poly)[0] = tmps;
	// 소유권을 넘겨주고 싶을 땐 다음처럼 쓰자.
	//(*Tmp_Poly)[0] = std::move(tmps);
	// 주소를 바로 넘겨줘보자.
	//Tmp_Poly->SetElementAddressAtIndex(0, tmps);


	//NativeVector tmpss = tmps.GetValues();
	//NativeInteger tmp = NativeInteger(0);
	


	std::cout << "after fix it tmps" << std::endl;
	std::cout << tmps << std::endl;
	std::cout << "remaining val is " << std::endl;
	std::cout << Tmp_Poly->GetElementAtIndex(0) << std::endl;
	//std::cout << "tmps addrss is " << std::endl;
	//std::cout << &tmps << std::endl;
	std::cout << "TmpPoly Address addrss is " << std::endl;
	std::cout << &(Tmp_Poly->GetElementAtIndex(0)) << std::endl;
	std::cout << &(Tmp_Poly->GetElementAtIndex(0)) << std::endl;
	std::cout << &(Tmp_Poly->GetElementAtIndex(0)) << std::endl;
	

	*/
	//개발 체크 OK
	



	//FPDCRTPoly *Tmp_Poly(Params);
//const std::shared_ptr<ILDCRTParams<BigInteger>> DCRTParamset = std::make_shared<ILDCRTParams<BigInteger>>(N_FP*2,moduliQ,rootsQ);
   // Maybe?
    //ChineseRemainderTransformFTT<NativeVector>::PreCompute(rootsQ, N_FP, moduliQ);



  // Product
	/*
	uint64_t msg_64 = 2;
	auto ct_64 = cc.EncryptUInt64(Rsk, msg_64);	
	uint64_t msg_642 = 1;
	auto ct_642 = cc.EncryptUInt64(Rsk, msg_642);	
	ct_64 = cc.Product(ct_64,ct_642);	
	auto ct_mul = cc.BootstrapOrign(ct_64);
	//auto ct_mul_BT2 = cc.BootstrapMul(ct_mul2);
	uint64_t res_msg_64_2 = cc.DecryptUInt64(Rsk,ct_mul);
	std::cout 
		<< "Product 1 & BT Orign  Decyption Result of 64 bit Encryption is \n"
		<< res_msg_64_2 << std::endl;
	*/

	//auto reduct_ct = cc.DivideReduction(ct_mul);
	/*
	uint64_t product_decrypt = cc.DecryptUInt64(Rsk,ct_mul);
	std::cout 
		<< "Product  Decyption Result of 64 bit Encryption is \n"
		<< product_decrypt << std::endl;
	*/
	/*
	uint64_t res_bt_before = cc.CheckDecryptBeforeBS(Rsk, ct_mul, 2);
	std::cout 
		<< " BTtt before modulus down  Decyption Result of 64 bit Encryption is \n"
		<< res_bt_before << std::endl;
	*/

//uint64_t res_msg_64 = cc.BSDecryptUInt64_TEST(ct_mul_BT);
	//std::cout 
	//	<< "Product 1 & BT  Decyption Result of 64 bit Encryption is \n"
	//	<< res_msg_64 << std::endl;



	
	/*
	NativePoly tmp1 = Rsk->GetElement().GetElementAtIndex(0);
	tmp1.SetFormat(Format::COEFFICIENT);

	std::cout << "SK is " << std::endl
		<< tmp1 << std::endl;
	*/
	//	cout.precision(64);
	//	std::cout 
	//		<< "(Product 1) * Orign Result of 64 bit Encryption is \n"
	//		<< res_msg_64 << std::endl;



//	auto ct_mul2 = cc.BootstrapOrign(cc.DivideReduction(cc.Product(Bootstrapped_ct,Bootstrapped_ct)));
//	uint64_t res_msg2_64 = cc.DecryptUInt64(Rsk, ct_mul2);
	
//	std::cout 
//		<< "Decyption Result of 64 bit Encryption 2 is \n"
//		<< res_msg2_64 << std::endl;
	



	/* Rounding Test Clear */
	/*
	uint64_t mm = 3;
	DCRTPoly msg_polys = cc.EncodingUInt64(mm);
	
	msg_polys += cc.MakeDG() * 100000;
	uint64_t dec_mm = cc.DecodingUInt64(msga_polys);
	
	msg_polys.SetFormat(Format::EVALUATION);

	DCRTPoly msg_muls = (msg_polys * msg_polys);
	msg_polys.SetFormat(Format::COEFFICIENT);
	msg_muls.SetFormat(Format::COEFFICIENT);
	

	DCRTPoly re_msg = cc.DivideRounding(msg_muls);
	
	uint64_t re_msg_int = cc.DecodingUInt64(re_msg);
	std::cout
		
		<< "Orign MSG" << std::endl
		<< msg_polys << std::endl
		<< "Mul MSG" << std::endl
		<< msg_muls << std::endl
		<< "After Rounding Mul" << std::endl
		<< re_msg << std::endl


		<< "Encoding msg is " << std::endl
		<< mm << std::endl
		<< "msg should be " << std::endl
		<< dec_mm << std::endl
		<< "and result is " << std::endl
		<< re_msg_int << std::endl;

		*/
	

	// After ....
	/*
	//const std::shared_ptr<FPRLWECiphertextImpl> ct_mul = cc.Product(ct_64,ct_64);
	

	auto ct_mul = cc.Product(ct_64,ct_64);
	
	//auto reduct_ct = cc.DivideReduction(ct_mul);
	//uint64_t res_msg_64 = cc.DecryptUInt64(Rsk,ct_64);
	uint64_t res_msg_64 = cc.DecryptUInt64(Rsk,ct_mul, true);
	//uint64_t res_msg_64 = cc.DecryptUInt64(Rsk,reduct_ct);
	// Sample Extract
	cout.precision(64);
	std::cout 
		<< "Result of 64 bit Encryption \n"
		<< res_msg_64 << std::endl;
	*/
	

	/* one-bits
	DCRTPoly tmp_a = ct_64->GetA();
	vector<DCRTPoly> tmp_a_decomp = tmp_a.BaseDecompose(1,false);

	std::cout
		<< "BaseDecom" << std::endl
		<< tmp_a_decomp << std::endl;
	std::cout << "\nEXIT!\n";
	*/



	
	/*
	std::cout << "Q check starts!" << std::endl;
	
	NativeInteger test_n = 12289;
	bool test_res = MillerRabinPrimalityTest(test_n);
	std::cout << "Miler Test result is "<<test_res << std::endl;
	
	uint32_t start_bits = 46;
	uint32_t search_gap = 10;
	uint32_t end_bits = start_bits + search_gap;
	uint32_t prime_gap = 12;

	NativeInteger P = 1;
	NativeInteger Q = 1;
	uint64_t iter_num = 1;
	for (uint64_t i = 0; i < end_bits - start_bits - 1; i++) {
		iter_num = iter_num * 2;
	}
	for (uint32_t i = 0; i < start_bits; i++){
		P = P * 2;
	}
	for (uint32_t i = 0; i < (start_bits - prime_gap); i++){
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
	
*/
	return 0;
}
