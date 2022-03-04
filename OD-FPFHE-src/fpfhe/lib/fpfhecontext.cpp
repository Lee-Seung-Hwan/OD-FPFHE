// @file binfhecontext.cpp - Implementation file for Boolean Circuit FHE context
// class
// @author TPOC: contact@palisade-crypto.org
// @copyright Copyright (c) 2019, Duality Technologies Inc.
// All rights reserved.
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution. THIS SOFTWARE IS
// PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
// EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "fpfhecontext.h"
#include "fprlwe.h"
#include <cmath>
//#include <math.h>
using namespace std;     // SH

namespace lbcrypto {
	bool DEBUGSH = false;

	/* Dont use it */
	void FPFHEContext::GenerateFHEContext(uint32_t n, uint32_t N,
                                          const uint32_t Packing_N,
																					const NativeInteger &q,
                                          const NativeInteger &Q, double std, double r_std,
																					double gsw_std,
                                          uint32_t baseKS, uint32_t scaling, uint32_t baseC,  uint32_t baseG,
                                          uint32_t baseR, uint32_t GD_base, FPFHEMETHOD method) {


		/*
		vector<NativeInteger> moduliQ(2);
		vector<NativeInteger> rootsQ(2);
		moduliQ[0] = 753227037677715457;
		moduliQ[1] = 718333280257;
		
		rootsQ[0] = RootOfUnity<NativeInteger>(2048,moduliQ[0]);
		rootsQ[1] = RootOfUnity<NativeInteger>(2048,moduliQ[1]);
		const std::shared_ptr<ILDCRTParams<BigInteger>> DCRTParamset = std::make_shared<ILDCRTParams<BigInteger>>(2048,moduliQ,rootsQ);
		//auto lweparams = std::make_shared<FPLWECryptoParams>(n, N, Packing_N, q, Q, std, baseKS, scaling);
  
		auto rlweparams = std::make_shared<FPRLWECryptoParams>(moduliQ, DCRTParamset, Packing_N, 2, 2,  r_std);
		
		//m_params = std::make_shared<FPRingGSWCryptoParams>(lweparams, rlweparams, baseC, baseG, baseR, 4, gsw_std, method);
	*/
	}

	void FPFHEContext::GenerateFHEContext(FPFHEPARAMSET set,
		bool parallel_flag) {

		shared_ptr<FPLWECryptoParams> lweparams;
		shared_ptr<FPRLWECryptoParams> rlweparams;
		
		uint32_t N_FP;
		uint32_t K_FP;
		uint32_t N_GSW;
		uint32_t K_GSW;
		uint32_t n_PKS;
		vector<NativeInteger> moduliQ(2);
		vector<NativeInteger> rootsQ(2);
		vector<NativeInteger> moduliQ_GSW(2);
		vector<NativeInteger> rootsQ_GSW(2);
		
		std::shared_ptr<ILDCRTParams<BigInteger>> DCRTParamset;
		std::shared_ptr<ILDCRTParams<BigInteger>> DCRTGSWParamset;
		
		double FP_std;
		double KS_std;
	  double BK_std;
	  double EV_std;
		double PKS_std;
	



		//Encoding
		uint32_t Encode_base_bit; // bit 4
		uint32_t C_base_bit;
		uint32_t carry_bit;
	

		//GSW Params Gadget
		uint32_t GD_base_bit; // 128
		uint32_t l_rm_GD;

		//KS Key Switching
		uint32_t KS_base_bit; // 128
		uint32_t l_rm_KS;

		//EV Evaluation bits
		uint32_t EV_base_bit; // 128
		uint32_t l_rm_EV;
	
		//PKS
		uint32_t PKS_base_bit; // Must 1 !!!
		uint32_t l_rm_PKS;
		uint32_t q_bit;
		uint32_t msg_bit;
		
		// sparse h
		uint32_t h_FP;
		uint32_t h_GSW;
		uint32_t h_PKS;

		// scalle_bit
		uint32_t exp_bit;
		
		// OMP Setting
		//omp_set_max_active_levels(2);
		bool PARALLEL_FLAG = parallel_flag;
		if (PARALLEL_FLAG) {
			std::cout << "Number of usable threads are " << omp_get_max_threads() << std::endl;
			std::cout << "And current threads are " << omp_get_num_threads() << std::endl;
		}

		switch (set) {
			case TOY:
				//uint32_t n_LWE = 4096;
				#ifdef WITH_INTEL_HEXL
				N_FP = (1 << 6);
				#else
				N_FP = (1 << 8);
				#endif
				K_FP = 1;
				N_GSW = (1 << 8);
				K_GSW = 1;
				n_PKS = 8;
				
				// TOY Modulus
				moduliQ[0] = 236223201281;
				moduliQ[1] = 230686721;
			
				moduliQ_GSW[0] = moduliQ[0];
				moduliQ_GSW[1] = moduliQ[1];
				
				// Calculation
				// FP Paramset 
				rootsQ[0] = RootOfUnity<NativeInteger>(N_FP*2,moduliQ[0]);
				rootsQ[1] = RootOfUnity<NativeInteger>(N_FP*2,moduliQ[1]);
				DCRTParamset = std::make_shared<ILDCRTParams<BigInteger>>(N_FP*2,moduliQ,rootsQ);
				ChineseRemainderTransformFTT<NativeVector>::PreCompute(rootsQ, 2*N_FP, moduliQ);
				
				rootsQ_GSW[0] = RootOfUnity<NativeInteger>(N_GSW*2,moduliQ[0]);
				rootsQ_GSW[1] = RootOfUnity<NativeInteger>(N_GSW*2,moduliQ[1]);
				DCRTGSWParamset = std::make_shared<ILDCRTParams<BigInteger>>(N_GSW*2,moduliQ_GSW,rootsQ_GSW);
				ChineseRemainderTransformFTT<NativeVector>::PreCompute(rootsQ_GSW, 2*N_GSW, moduliQ_GSW);
				
				FP_std		= 0.398942280;
				KS_std		= 0.398942280;
				BK_std		= 0.398942280;
				EV_std		= 0.398942280;
				PKS_std	  = 0.398942280;
				//Encoding
				Encode_base_bit = 2; // bit 4
				C_base_bit = 4;
				carry_bit = 2;
				//GSW Params Gadget
				GD_base_bit = 9; // 128
				l_rm_GD = 1;

				//KS Key Switching
				KS_base_bit = 9; // 128
				l_rm_KS = 1;

				//EV Evaluation bits
				EV_base_bit = 6; // 128
				l_rm_EV = 1;
	
				//PKS
				PKS_base_bit = 1; // Must 1 !!!
				l_rm_PKS = 0;
				q_bit = 18;
				msg_bit = C_base_bit + carry_bit;

				// scalle_bit
				exp_bit = q_bit;




				// sparse h
				h_FP = 10000000;
				h_GSW = 10000000;
				h_PKS = n_PKS;
				lweparams = std::make_shared<FPLWECryptoParams>(n_PKS, q_bit, PKS_std, PKS_base_bit, l_rm_PKS, msg_bit);
		
				rlweparams = std::make_shared<FPRLWECryptoParams>( moduliQ, 
						DCRTParamset, N_FP, K_FP, GD_base_bit, l_rm_GD,  FP_std);
				std::cout << "rlweparams is GEN !! " << std::endl;

				m_params =
          std::make_shared<FPRingGSWCryptoParams>(
							lweparams,
							rlweparams,
							DCRTGSWParamset,	// GSW Params
							N_GSW,						// GSW Params
							K_GSW,						// GSW Params
							n_PKS,						// Pre Key Switch params
							KS_std, BK_std, EV_std, PKS_std, // STD
							KS_base_bit, l_rm_KS, // KS info
							EV_base_bit, l_rm_EV, // Evaluation Info
							GD_base_bit, l_rm_GD, // Gadget Info
							PKS_base_bit, l_rm_PKS, // Pre Key Switch Info
							Encode_base_bit, C_base_bit,
							h_FP, h_GSW, h_PKS,
							exp_bit, PARALLEL_FLAG
							); // Encoding Info

	
				std::cout << "m_params is GEN !! " << std::endl;
				break;

			case STD128D:
				N_GSW = 1;
				N_GSW <<= 11;
				K_GSW = 2;
				#ifdef WITH_INTEL_HEXL
				N_FP = 1;
				
				N_FP <<= 8;
				K_FP = 13;
				//N_FP <<= 6;
				//K_FP = 51;
				#else
				N_FP = N_GSW;
				K_FP = K_GSW;
				#endif

				n_PKS = 785;
		
				// Nomi 1
				
				//moduliQ[0] = 35871566856193;
				//moduliQ[1] = 35030827009;



				moduliQ[0] = 286422779035649;
				moduliQ[1] = 279709745153;
				//moduliQ[0] = 144955146241;
				//moduliQ[1] = 141557761;


					

				//moduliQ[0] = 9225589751809;
				//moduliQ[1] = 9009364993;



				// Calculation
				// FP Paramset 
				rootsQ[0] = RootOfUnity<NativeInteger>(N_FP*2,moduliQ[0]);
				rootsQ[1] = RootOfUnity<NativeInteger>(N_FP*2,moduliQ[1]);
				DCRTParamset = std::make_shared<ILDCRTParams<BigInteger>>(N_FP*2,moduliQ,rootsQ);
				ChineseRemainderTransformFTT<NativeVector>::PreCompute(rootsQ, 2*N_FP, moduliQ);
				if(N_FP == N_GSW) {
					DCRTGSWParamset = DCRTParamset;
				} else {
					rootsQ[0] = RootOfUnity<NativeInteger>(N_GSW*2,moduliQ[0]);
					rootsQ[1] = RootOfUnity<NativeInteger>(N_GSW*2,moduliQ[1]);
					DCRTGSWParamset = std::make_shared<ILDCRTParams<BigInteger>>(N_GSW*2,moduliQ,rootsQ);
	
					ChineseRemainderTransformFTT<NativeVector>::PreCompute(rootsQ, 2*N_GSW, moduliQ);
				}
				FP_std		= 3.19;
				KS_std		= 3.19;
				BK_std		= 3.19;
				EV_std		= 3.19;
				PKS_std		= 3.19;
				//Encoding
				Encode_base_bit = 2; // bit 4
				C_base_bit = 4;
				carry_bit = 2;
				//GSW Params Gadget
				GD_base_bit = 11; // 128
				l_rm_GD = 1;

				//KS Key Switching
				KS_base_bit = 15; // 128
				l_rm_KS = 1;


				//EV Evaluation bits
				EV_base_bit = 18; // 128
				l_rm_EV = 1;
	
				//PKS
				PKS_base_bit = 1; // Must 1 !!!
				l_rm_PKS = 0;
				q_bit = 21;
				msg_bit = C_base_bit + carry_bit;
				
				//exp_bit = 10;
				exp_bit = q_bit;


				// sparse h
				h_FP	= 100000000; // this is not works
				h_GSW = 10000000; // this is not works
				h_PKS = 131;
				lweparams = std::make_shared<FPLWECryptoParams>(n_PKS, q_bit, PKS_std, PKS_base_bit, l_rm_PKS, msg_bit);
		
				rlweparams = std::make_shared<FPRLWECryptoParams>( moduliQ, 
						DCRTParamset, N_FP, K_FP, GD_base_bit, l_rm_GD,  FP_std);
				std::cout << "rlweparams is GEN !! " << std::endl;

				m_params =
          std::make_shared<FPRingGSWCryptoParams>(
							lweparams,
							rlweparams,
							DCRTGSWParamset,	// GSW Params
							N_GSW,						// GSW Params
							K_GSW,						// GSW Params
							n_PKS,						// Pre Key Switch params
							KS_std, BK_std, EV_std, PKS_std, // STD
							KS_base_bit, l_rm_KS, // KS info
							EV_base_bit, l_rm_EV, // Evaluation Info
							GD_base_bit, l_rm_GD, // Gadget Info
							PKS_base_bit, l_rm_PKS, // Pre Key Switch Info
							Encode_base_bit, C_base_bit,
							h_FP, h_GSW, h_PKS,
							exp_bit, PARALLEL_FLAG
							); // Encoding Info

	
				std::cout << "m_params is GEN !! " << std::endl;
	
			break;


			case STD128S:
				//uint32_t n_LWE = 4096;
				N_GSW = 1;
				N_GSW <<= 12;
				K_GSW = 1;
				#ifdef WITH_INTEL_HEXL
				N_FP = 1;
				N_FP <<= 8;
				K_FP = 13;
				//N_FP <<= 5;
				//K_FP = 99;
				#else
				N_FP = N_GSW;
				K_FP = K_GSW;
				#endif



				n_PKS = 785;
				
				moduliQ[0] = 74217034874881;
				moduliQ[1] = 144955146241;
				// Calculation
				// FP Paramset 
				rootsQ[0] = RootOfUnity<NativeInteger>(N_FP*2,moduliQ[0]);
				rootsQ[1] = RootOfUnity<NativeInteger>(N_FP*2,moduliQ[1]);
				DCRTParamset = std::make_shared<ILDCRTParams<BigInteger>>(N_FP*2,moduliQ,rootsQ);
				ChineseRemainderTransformFTT<NativeVector>::PreCompute(rootsQ, 2*N_FP, moduliQ);
				if(N_FP == N_GSW) {
					DCRTGSWParamset = DCRTParamset;
				} else {
					rootsQ[0] = RootOfUnity<NativeInteger>(N_GSW*2,moduliQ[0]);
					rootsQ[1] = RootOfUnity<NativeInteger>(N_GSW*2,moduliQ[1]);
					DCRTGSWParamset = std::make_shared<ILDCRTParams<BigInteger>>(N_GSW*2,moduliQ,rootsQ);
	
					ChineseRemainderTransformFTT<NativeVector>::PreCompute(rootsQ, 2*N_GSW, moduliQ);
				}
				FP_std		= 3.19;
				KS_std		= 3.19;
				BK_std		= 3.19;
				EV_std		= 3.19;
				PKS_std		= 3.19;
				//Encoding
				Encode_base_bit = 2; // bit 4
				C_base_bit = 4;
				carry_bit = 2;
				//GSW Params Gadget
				GD_base_bit = 12; // 128
				l_rm_GD = 1;

				//KS Key Switching
				KS_base_bit = 16; // 128
				l_rm_KS = 1;
				//KS_base_bit = 1; // 128
				//l_rm_KS = 30;


				//EV Evaluation bits
				EV_base_bit = 18; // 128
				l_rm_EV = 1;
	
				//PKS
				PKS_base_bit = 1; // Must 1 !!!
				l_rm_PKS = 0;
				q_bit = 21;
				msg_bit = C_base_bit + carry_bit;
				
				//exp_bit = 10;
				exp_bit = q_bit;


				// sparse h
				h_FP	= 100000000; // this is not works
				h_GSW = 10000000; // this is not works
				h_PKS = 131;
				lweparams = std::make_shared<FPLWECryptoParams>(n_PKS, q_bit, PKS_std, PKS_base_bit, l_rm_PKS, msg_bit);
		
				rlweparams = std::make_shared<FPRLWECryptoParams>( moduliQ, 
						DCRTParamset, N_FP, K_FP, GD_base_bit, l_rm_GD,  FP_std);
				std::cout << "rlweparams is GEN !! " << std::endl;

				m_params =
          std::make_shared<FPRingGSWCryptoParams>(
							lweparams,
							rlweparams,
							DCRTGSWParamset,	// GSW Params
							N_GSW,						// GSW Params
							K_GSW,						// GSW Params
							n_PKS,						// Pre Key Switch params
							KS_std, BK_std, EV_std, PKS_std, // STD
							KS_base_bit, l_rm_KS, // KS info
							EV_base_bit, l_rm_EV, // Evaluation Info
							GD_base_bit, l_rm_GD, // Gadget Info
							PKS_base_bit, l_rm_PKS, // Pre Key Switch Info
							Encode_base_bit, C_base_bit,
							h_FP, h_GSW, h_PKS,
							exp_bit, PARALLEL_FLAG
							); // Encoding Info

	
				std::cout << "m_params is GEN !! " << std::endl;
	
			break;
///////////////////////////////////////////////////////////////////////////

			case STD160D:
				//uint32_t n_LWE = 4096;
				N_GSW = 1;
				N_GSW <<= 11;
				K_GSW = 2;
				#ifdef WITH_INTEL_HEXL
				N_FP = 1;
				N_FP <<= 8;
				K_FP = 16;
				//N_FP <<= 6;
				//K_FP = 62;
				#else
				N_FP = N_GSW;
				K_FP = K_GSW;
				#endif
				n_PKS = 1089;

				moduliQ[0] = 286422779035649;
				moduliQ[1] = 279709745153;

				// Calculation
				// FP Paramset 
				rootsQ[0] = RootOfUnity<NativeInteger>(N_FP*2,moduliQ[0]);
				rootsQ[1] = RootOfUnity<NativeInteger>(N_FP*2,moduliQ[1]);
				DCRTParamset = std::make_shared<ILDCRTParams<BigInteger>>(N_FP*2,moduliQ,rootsQ);
				ChineseRemainderTransformFTT<NativeVector>::PreCompute(rootsQ, 2*N_FP, moduliQ);
				if(N_FP == N_GSW) {
					DCRTGSWParamset = DCRTParamset;
				} else {
					rootsQ[0] = RootOfUnity<NativeInteger>(N_GSW*2,moduliQ[0]);
					rootsQ[1] = RootOfUnity<NativeInteger>(N_GSW*2,moduliQ[1]);
					DCRTGSWParamset = std::make_shared<ILDCRTParams<BigInteger>>(N_GSW*2,moduliQ,rootsQ);
	
					ChineseRemainderTransformFTT<NativeVector>::PreCompute(rootsQ, 2*N_GSW, moduliQ);
				}
				FP_std		= 3.19;
				KS_std		= 3.19;
				BK_std		= 3.19;
				EV_std		= 3.19;
				PKS_std		= 3.19;
				//Encoding
				Encode_base_bit = 2; // bit 4
				C_base_bit = 4;
				carry_bit = 2;
				//GSW Params Gadget
				GD_base_bit = 10; // 128
				l_rm_GD = 1;

				//KS Key Switching
				KS_base_bit = 14; // 128
				l_rm_KS = 1;
				//KS_base_bit = 1; // 128
				//l_rm_KS = 30;


				//EV Evaluation bits
				EV_base_bit = 18; // 128
				l_rm_EV = 1;
	
				//PKS
				PKS_base_bit = 1; // Must 1 !!!
				l_rm_PKS = 0;
				q_bit = 21;
				msg_bit = C_base_bit + carry_bit;
				
				//exp_bit = 10;
				exp_bit = q_bit;


				// sparse h
				h_FP	= 100000000; // this is not works
				h_GSW = 10000000; // this is not works
				h_PKS = 131;
				lweparams = std::make_shared<FPLWECryptoParams>(n_PKS, q_bit, PKS_std, PKS_base_bit, l_rm_PKS, msg_bit);
		
				rlweparams = std::make_shared<FPRLWECryptoParams>( moduliQ, 
						DCRTParamset, N_FP, K_FP, GD_base_bit, l_rm_GD,  FP_std);
				std::cout << "rlweparams is GEN !! " << std::endl;

				m_params =
          std::make_shared<FPRingGSWCryptoParams>(
							lweparams,
							rlweparams,
							DCRTGSWParamset,	// GSW Params
							N_GSW,						// GSW Params
							K_GSW,						// GSW Params
							n_PKS,						// Pre Key Switch params
							KS_std, BK_std, EV_std, PKS_std, // STD
							KS_base_bit, l_rm_KS, // KS info
							EV_base_bit, l_rm_EV, // Evaluation Info
							GD_base_bit, l_rm_GD, // Gadget Info
							PKS_base_bit, l_rm_PKS, // Pre Key Switch Info
							Encode_base_bit, C_base_bit,
							h_FP, h_GSW, h_PKS,
							exp_bit, PARALLEL_FLAG
							); // Encoding Info

	
				std::cout << "m_params is GEN !! " << std::endl;
	
			break;


			case STD160S:
				//uint32_t n_LWE = 4096;
				N_GSW = 1;
				N_GSW <<= 12;
				K_GSW = 1;
				#ifdef WITH_INTEL_HEXL
				N_FP = 1;
				N_FP <<= 8;
				K_FP = 15;
				//N_FP <<= 5;
				//K_FP = 120;
				#else
				N_FP = N_GSW;
				K_FP = K_GSW;
				#endif



				n_PKS = 1089;
				
				moduliQ[0] = 74217034874881;
				moduliQ[1] = 144955146241;
				// Calculation
				// FP Paramset 
				rootsQ[0] = RootOfUnity<NativeInteger>(N_FP*2,moduliQ[0]);
				rootsQ[1] = RootOfUnity<NativeInteger>(N_FP*2,moduliQ[1]);
				DCRTParamset = std::make_shared<ILDCRTParams<BigInteger>>(N_FP*2,moduliQ,rootsQ);
				ChineseRemainderTransformFTT<NativeVector>::PreCompute(rootsQ, 2*N_FP, moduliQ);
				if(N_FP == N_GSW) {
					DCRTGSWParamset = DCRTParamset;
				} else {
					rootsQ[0] = RootOfUnity<NativeInteger>(N_GSW*2,moduliQ[0]);
					rootsQ[1] = RootOfUnity<NativeInteger>(N_GSW*2,moduliQ[1]);
					DCRTGSWParamset = std::make_shared<ILDCRTParams<BigInteger>>(N_GSW*2,moduliQ,rootsQ);
	
					ChineseRemainderTransformFTT<NativeVector>::PreCompute(rootsQ, 2*N_GSW, moduliQ);
				}
				FP_std		= 3.19;
				KS_std		= 3.19;
				BK_std		= 3.19;
				EV_std		= 3.19;
				PKS_std		= 3.19;
				//Encoding
				Encode_base_bit = 2; // bit 4
				C_base_bit = 4;
				carry_bit = 2;
				//GSW Params Gadget
				GD_base_bit = 12; // 128
				l_rm_GD = 1;

				//KS Key Switching
				KS_base_bit = 16; // 128
				l_rm_KS = 1;
				//KS_base_bit = 1; // 128
				//l_rm_KS = 30;


				//EV Evaluation bits
				EV_base_bit = 18; // 128
				l_rm_EV = 1;
	
				//PKS
				PKS_base_bit = 1; // Must 1 !!!
				l_rm_PKS = 0;
				q_bit = 21;
				msg_bit = C_base_bit + carry_bit;
				
				//exp_bit = 10;
				exp_bit = q_bit;


				// sparse h
				h_FP	= 100000000; // this is not works
				h_GSW = 10000000; // this is not works
				h_PKS = 131;
				lweparams = std::make_shared<FPLWECryptoParams>(n_PKS, q_bit, PKS_std, PKS_base_bit, l_rm_PKS, msg_bit);
		
				rlweparams = std::make_shared<FPRLWECryptoParams>( moduliQ, 
						DCRTParamset, N_FP, K_FP, GD_base_bit, l_rm_GD,  FP_std);
				std::cout << "rlweparams is GEN !! " << std::endl;

				m_params =
          std::make_shared<FPRingGSWCryptoParams>(
							lweparams,
							rlweparams,
							DCRTGSWParamset,	// GSW Params
							N_GSW,						// GSW Params
							K_GSW,						// GSW Params
							n_PKS,						// Pre Key Switch params
							KS_std, BK_std, EV_std, PKS_std, // STD
							KS_base_bit, l_rm_KS, // KS info
							EV_base_bit, l_rm_EV, // Evaluation Info
							GD_base_bit, l_rm_GD, // Gadget Info
							PKS_base_bit, l_rm_PKS, // Pre Key Switch Info
							Encode_base_bit, C_base_bit,
							h_FP, h_GSW, h_PKS,
							exp_bit, PARALLEL_FLAG
							); // Encoding Info

	
				std::cout << "m_params is GEN !! " << std::endl;
	
			break;

/////////////////////////////////////////////////////////////////////////////////////
			case STD192D:
				std::cout << "STD192D" << std::endl;
				//uint32_t n_LWE = 4096;
				N_GSW = 1;
				N_GSW <<= 11;
				K_GSW = 3;
				#ifdef WITH_INTEL_HEXL
				N_FP = 1;
				N_FP <<= 8;
				K_FP = 19;
				//N_FP <<= 6;
				//K_FP = 73;
				#else
				N_FP = N_GSW;
				K_FP = K_GSW;
				#endif
		
				n_PKS = 1292;
			
				moduliQ[0] = 286422779035649;
				moduliQ[1] = 279709745153;


				// Calculation
				// FP Paramset 
				rootsQ[0] = RootOfUnity<NativeInteger>(N_FP*2,moduliQ[0]);
				rootsQ[1] = RootOfUnity<NativeInteger>(N_FP*2,moduliQ[1]);
				DCRTParamset = std::make_shared<ILDCRTParams<BigInteger>>(N_FP*2,moduliQ,rootsQ);
				ChineseRemainderTransformFTT<NativeVector>::PreCompute(rootsQ, 2*N_FP, moduliQ);
				if(N_FP == N_GSW) {
					DCRTGSWParamset = DCRTParamset;
				} else {
					rootsQ[0] = RootOfUnity<NativeInteger>(N_GSW*2,moduliQ[0]);
					rootsQ[1] = RootOfUnity<NativeInteger>(N_GSW*2,moduliQ[1]);
					DCRTGSWParamset = std::make_shared<ILDCRTParams<BigInteger>>(N_GSW*2,moduliQ,rootsQ);
	
					ChineseRemainderTransformFTT<NativeVector>::PreCompute(rootsQ, 2*N_GSW, moduliQ);
				}
				FP_std		= 3.19;
				KS_std		= 3.19;
				BK_std		= 3.19;
				EV_std		= 3.19;
				PKS_std		= 3.19;
				//Encoding
				Encode_base_bit = 2; // bit 4
				C_base_bit = 4;
				carry_bit = 2;
				//GSW Params Gadget
				GD_base_bit = 10; // 128
				l_rm_GD = 1;

				//KS Key Switching
				KS_base_bit = 13; // 128
				l_rm_KS = 1;


				//EV Evaluation bits
				EV_base_bit = 18; // 128
				l_rm_EV = 1;
	
				//PKS
				PKS_base_bit = 1; // Must 1 !!!
				l_rm_PKS = 0;
				q_bit = 22;
				msg_bit = C_base_bit + carry_bit;
				
				//exp_bit = 10;
				exp_bit = q_bit;


				// sparse h
				h_FP	= 100000000; // this is not works
				h_GSW = 10000000; // this is not works
				h_PKS = 160;
				lweparams = std::make_shared<FPLWECryptoParams>(n_PKS, q_bit, PKS_std, PKS_base_bit, l_rm_PKS, msg_bit);
		
				rlweparams = std::make_shared<FPRLWECryptoParams>( moduliQ, 
						DCRTParamset, N_FP, K_FP, GD_base_bit, l_rm_GD,  FP_std);
				std::cout << "rlweparams is GEN !! " << std::endl;

				m_params =
          std::make_shared<FPRingGSWCryptoParams>(
							lweparams,
							rlweparams,
							DCRTGSWParamset,	// GSW Params
							N_GSW,						// GSW Params
							K_GSW,						// GSW Params
							n_PKS,						// Pre Key Switch params
							KS_std, BK_std, EV_std, PKS_std, // STD
							KS_base_bit, l_rm_KS, // KS info
							EV_base_bit, l_rm_EV, // Evaluation Info
							GD_base_bit, l_rm_GD, // Gadget Info
							PKS_base_bit, l_rm_PKS, // Pre Key Switch Info
							Encode_base_bit, C_base_bit,
							h_FP, h_GSW, h_PKS,
							exp_bit, PARALLEL_FLAG
							); // Encoding Info

	
				std::cout << "m_params is GEN !! " << std::endl;
	
			break;


			case STD192S:
				//uint32_t n_LWE = 4096;
				N_GSW = 1;
				N_GSW <<= 11;
				K_GSW = 3;
				#ifdef WITH_INTEL_HEXL
				N_FP = 1;
				N_FP <<= 8;
				K_FP = 18;
				//N_FP <<= 5;
				//K_FP = 141;
				#else
				N_FP = N_GSW;
				K_FP = K_GSW;
				#endif
		
				n_PKS = 1292;
				
				moduliQ[0] = 74217034874881;
				moduliQ[1] = 144955146241;
				// Calculation
				// FP Paramset 
				rootsQ[0] = RootOfUnity<NativeInteger>(N_FP*2,moduliQ[0]);
				rootsQ[1] = RootOfUnity<NativeInteger>(N_FP*2,moduliQ[1]);
				DCRTParamset = std::make_shared<ILDCRTParams<BigInteger>>(N_FP*2,moduliQ,rootsQ);
				ChineseRemainderTransformFTT<NativeVector>::PreCompute(rootsQ, 2*N_FP, moduliQ);
				if(N_FP == N_GSW) {
					DCRTGSWParamset = DCRTParamset;
				} else {
					rootsQ[0] = RootOfUnity<NativeInteger>(N_GSW*2,moduliQ[0]);
					rootsQ[1] = RootOfUnity<NativeInteger>(N_GSW*2,moduliQ[1]);
					DCRTGSWParamset = std::make_shared<ILDCRTParams<BigInteger>>(N_GSW*2,moduliQ,rootsQ);
	
					ChineseRemainderTransformFTT<NativeVector>::PreCompute(rootsQ, 2*N_GSW, moduliQ);
				}
				FP_std		= 3.19;
				KS_std		= 3.19;
				BK_std		= 3.19;
				EV_std		= 3.19;
				PKS_std		= 3.19;
				//Encoding
				Encode_base_bit = 2; // bit 4
				C_base_bit = 4;
				carry_bit = 2;
				//GSW Params Gadget
				GD_base_bit = 10; // 128
				l_rm_GD = 1;

				//KS Key Switching
				KS_base_bit = 14; // 128
				l_rm_KS = 1;


				//EV Evaluation bits
				EV_base_bit = 18; // 128
				l_rm_EV = 1;
	
				//PKS
				PKS_base_bit = 1; // Must 1 !!!
				l_rm_PKS = 0;
				q_bit = 22;
				msg_bit = C_base_bit + carry_bit;
				
				exp_bit = q_bit;


				// sparse h
				h_FP	= 100000000; // this is not works
				h_GSW = 10000000; // this is not works
				h_PKS = 160;
				lweparams = std::make_shared<FPLWECryptoParams>(n_PKS, q_bit, PKS_std, PKS_base_bit, l_rm_PKS, msg_bit);
		
				rlweparams = std::make_shared<FPRLWECryptoParams>( moduliQ, 
						DCRTParamset, N_FP, K_FP, GD_base_bit, l_rm_GD,  FP_std);
				std::cout << "rlweparams is GEN !! " << std::endl;

				m_params =
          std::make_shared<FPRingGSWCryptoParams>(
							lweparams,
							rlweparams,
							DCRTGSWParamset,	// GSW Params
							N_GSW,						// GSW Params
							K_GSW,						// GSW Params
							n_PKS,						// Pre Key Switch params
							KS_std, BK_std, EV_std, PKS_std, // STD
							KS_base_bit, l_rm_KS, // KS info
							EV_base_bit, l_rm_EV, // Evaluation Info
							GD_base_bit, l_rm_GD, // Gadget Info
							PKS_base_bit, l_rm_PKS, // Pre Key Switch Info
							Encode_base_bit, C_base_bit,
							h_FP, h_GSW, h_PKS,
							exp_bit, PARALLEL_FLAG
							); // Encoding Info

	
				std::cout << "m_params is GEN !! " << std::endl;
	
			break;


			default:
				std::string errMsg = "ERROR: No such parameter set exists afor FPFHE.";
				PALISADE_THROW(config_error, errMsg);
  
		}
  
	}

	//FPLWEPrivateKey FPFHEContext::KeyGen() const {
	//	return m_LWEscheme->KeyGen(m_params->GetLWEParams());
	//}


	//FPLWEPrivateKey FPFHEContext::KeyGenN() const {
	//	return m_LWEscheme->KeyGenN(m_params->GetLWEParams());
	//}

		// Mul
		//
	
	FPRLWEPrivateKey FPFHEContext::RKeyGen() const {
		return m_RLWEscheme -> KeyGen(m_params -> GetRLWEParams(),
				m_params->GetH_FP());
	}

	//FPLWEPrivateKey FPFHEContext::RLWEKey2LWEKey(FPRLWEPrivateKey sk_r) const {	
	//	DCRTPoly sk_poly = sk_r ->GetElement();
	//	return std::make_shared<FPLWEPrivateKeyImpl> (FPLWEPrivateKeyImpl(sk_r->GetElement().GetValues())); 
	
	//}
	

	/*
	FPLWECiphertext FPFHEContext::Encrypt(FPConstLWEPrivateKey sk,
                                     const FPLWEPlaintext &m,
                                     FPFHEOUTPUT output) const {
		//if (output == FRESH) {
			return m_LWEscheme->Encrypt(m_params->GetLWEParams(), sk, m);
		//} else {
		//	auto ct = m_LWEscheme->Encrypt(m_params->GetLWEParams(), sk, m);
		//	return m_RingGSWscheme->Bootstrap(m_params, m_BTKey, ct,1,1);
		//}
	}
	*/


	/*
	FPConstRLWECiphertext FPFHEContext::DivideReduction( 
		FPConstRLWECiphertext ct) const {
		vector<DCRTPoly> ct_a = ct->GetA();
		vector<DCRTPoly> Reduced_a(ct_a.size());
    DCRTPoly ct_b = ct->GetB();

		CTType ct_type = ct->GetType();
		ct_b.SetFormat(Format::COEFFICIENT);
		
		// Level 1 shoud be checked 
		for (uint32_t i = 0; i < ct_a.size(); i ++) {
			ct_a[i].SetFormat(Format::COEFFICIENT);
			Reduced_a[i] = DivideRounding(ct_a[i]);
			Reduced_a[i].SetFormat(Format::EVALUATION);
		}
		DCRTPoly Reduced_b = DivideRounding(ct_b);
		Reduced_b.SetFormat(Format::EVALUATION);
		
		//Sign 
		vector<DCRTPoly> Reduced_a_sign;
		DCRTPoly Reduced_b_sign;
		if (ct_type == INT64) {
			vector<DCRTPoly> ct_a_sign = ct->GetASign();
			Reduced_a_sign.resize(ct_a_sign.size());
			for (uint32_t i = 0; i < ct_a_sign.size(); i ++) {
				ct_a_sign[i].SetFormat(Format::COEFFICIENT);
				Reduced_a_sign[i] = DivideRounding(ct_a_sign[i]);
				Reduced_a_sign[i].SetFormat(Format::EVALUATION);
			}

			DCRTPoly ct_b_sign = ct->GetBSign();
			ct_b_sign.SetFormat(Format::COEFFICIENT);
			Reduced_b_sign = DivideRounding(ct_b_sign);
			Reduced_b_sign.SetFormat(Format::EVALUATION);

		}

		if (ct_type == UINT64) {
			return std::make_shared<FPRLWECiphertextImpl> (FPRLWECiphertextImpl(Reduced_a,Reduced_b,1));
		} else {
			return std::make_shared<FPRLWECiphertextImpl> (FPRLWECiphertextImpl(Reduced_a,Reduced_b, Reduced_a_sign, Reduced_b_sign, 1));
	
		}
	}
	*/






	void FPFHEContext::DivideRounding(DCRTPoly *double_poly, DCRTPoly *out_poly, Format format) const {

    const std::shared_ptr<FPRLWECryptoParams> params = m_params->GetRLWEParams();
		uint32_t Ns = double_poly->GetLength();
		// Invariant to N
		vector<NativeInteger>	QInv = params->GetInv();
		vector<NativeInteger> moduliQ = params->GetModuliQ();
    uint64_t Q1_half = moduliQ[1].ConvertToInt() >> 1; 
	
		Format fm_save = double_poly->GetFormat();
		double_poly->SetFormat(Format::COEFFICIENT);
		out_poly->SetFormat(Format::COEFFICIENT);
		vector<NativePoly> poly = double_poly->GetAllElements();
	
		NativeInteger a0;
		for (uint32_t i = 0; i < Ns; i++) {
			a0 = poly[0][i].ModSub(poly[1][i], moduliQ[0]).ModMul(QInv[0], moduliQ[0]);
      if(poly[1][i].ConvertToInt() >= Q1_half) {
				a0.ModAddEq(1,moduliQ[0]);
      }
      out_poly->GetElementW(0)[i] = a0;
    }
		out_poly->SetFormat(format);
		double_poly->SetFormat(fm_save);
	} 

	DCRTPoly FPFHEContext::DivideRounding(DCRTPoly *in_poly, Format format=Format::COEFFICIENT) const {
    DCRTPoly out_poly = DCRTPoly(in_poly->GetParams(), Format::COEFFICIENT,true);
		
		this->DivideRounding(in_poly, &out_poly, format);
		return out_poly;
  }

	int64_t FPFHEContext::DecryptExp(
			FPConstRLWEPrivateKey sk, 
			const vector<DCRTPoly> &cts) const {
		uint32_t K_FP = m_params->GetRLWEParams()->GetK();

		DCRTPoly exp_result;
		m_RLWEscheme->Decrypt(m_params->GetRLWEParams(), 
			sk, 
			std::make_shared<FPRLWECiphertextImpl>(FPRLWECiphertextImpl(
			cts, cts[K_FP])),
			&exp_result);
		
		exp_result = this->DivideRounding(&exp_result, Format::COEFFICIENT);
		return DecodingExp(exp_result);
	}




	double FPFHEContext::DecryptDouble(
			FPConstRLWEPrivateKey sk, 
			FPConstRLWECiphertext ct,
			const bool DIVIDE
			) const {

		// Handle Frac
		DCRTPoly frac_result;	
		m_RLWEscheme->Decrypt(m_params->GetRLWEParams(), sk, ct, &frac_result);
		if (DIVIDE == true) {
			frac_result = this->DivideRounding(&frac_result, Format::COEFFICIENT);
		}else {
			frac_result.SetFormat(Format::COEFFICIENT);
		}

		// Setting
		NativeInteger Q = m_params -> GetRLWEParams() -> GetModuliQ()[0];
		double m_frac = 0;
		uint64_t tmp_res = 0;

		NativePoly result0 = frac_result.GetElementAtIndex(0);
		result0.SetFormat(Format::COEFFICIENT);	
		uint32_t tmp;
		uint32_t tmp2;

		// 4 ~ 31 are our informations 
		for (uint32_t i = 5; i < 32; i++) {
			//uint32_t save_idx = 31 - i;
			
			if (result0[i] == Q-1) {
				result0[i] = NativeInteger(0);
			}
			
			tmp_res = result0[i].ConvertToInt() / m_params->GetRLWEParams()->GetScaling().ConvertToInt();
			tmp_res = (tmp_res+m_params->GetRLWEParams()->GetBits_over_2()) >> m_params->GetRLWEParams()->GetBits()[1];	
			tmp = (tmp_res % 4) >> 1;
			tmp2 = (tmp_res % 2);
			//std::cout << "i is " << i << ", Val is " << tmp << std::endl;
			m_frac += (double) tmp2;
			m_frac *= 0.5; //pow(4.0,-1);
			m_frac += (double) tmp;
			m_frac *= 0.5; //pow(4.0,-1);



		}	
			
		// Handle Exponent
	
		DCRTPoly exp_result;
		m_RLWEscheme->Decrypt(m_params->GetRLWEParams(), 
				sk, 
				std::make_shared<FPRLWECiphertextImpl>(FPRLWECiphertextImpl(
					ct->GetAExpo(), ct->GetBExpo()
				)),
				&exp_result);
		exp_result = this->DivideRounding(&exp_result, Format::COEFFICIENT);
		int64_t m_exp = DecodingExp(exp_result);
		m_exp -= 511;
		m_exp <<= 1;	

		// Handle Sign
		vector<DCRTPoly> Sign_A = ct->GetASign();
		DCRTPoly Sign_B = ct->GetBSign();
		Sign_B.SetFormat(Format::COEFFICIENT);
		
		Sign_B.GetElementW(0)[0].ModAddEq(
			m_params->GetRLWEParams()->GetModuliQ()[1],
			m_params->GetRLWEParams()->GetModuliQ()[0]);
		
		Sign_B.SetFormat(Format::EVALUATION);	
		uint64_t msg_signed = DecryptUInt64(
				sk,
				std::make_shared<FPRLWECiphertextImpl>(FPRLWECiphertextImpl(
				std::move(Sign_A), std::move(Sign_B)
					)),
				DIVIDE
				);
		
		if (msg_signed == 2) {
			return (m_frac * pow(2.0, m_exp));
		} else if (msg_signed == 0) {
			return - (m_frac * pow(2.0, m_exp));
		} else {
			std::cout <<"Sign Encryption is Invalid!!!!" << std::endl;
			return 0.0;
		}

	}



	float FPFHEContext::DecryptFloat(
			FPConstRLWEPrivateKey sk, 
			FPConstRLWECiphertext ct,
			const bool DIVIDE
			) const {

		// Handle Frac
		DCRTPoly frac_result;	
		m_RLWEscheme->Decrypt(m_params->GetRLWEParams(), sk, ct, &frac_result);
		if (DIVIDE == true) {
			frac_result = this->DivideRounding(&frac_result, Format::COEFFICIENT);
		}else {
			frac_result.SetFormat(Format::COEFFICIENT);
		}

		// Setting
		NativeInteger Q = m_params -> GetRLWEParams() -> GetModuliQ()[0];
		float m_frac = 0;
		uint64_t tmp_res = 0;

		NativePoly result0 = frac_result.GetElementAtIndex(0);
		result0.SetFormat(Format::COEFFICIENT);	
		uint32_t tmp;
		uint32_t tmp2;
		// 4 ~ 31 are our informations 
		for (uint32_t i = 3; i < 16; i++) {
			//uint32_t save_idx = 31 - i;
			
			if (result0[i] == Q-1) {
				result0[i] = NativeInteger(0);
			}
			
			tmp_res = result0[i].ConvertToInt() / m_params->GetRLWEParams()->GetScaling().ConvertToInt();
			tmp_res = (tmp_res+m_params->GetRLWEParams()->GetBits_over_2()) >> m_params->GetRLWEParams()->GetBits()[1];	
			tmp = (tmp_res % 4) >> 1;
			tmp2 = (tmp_res % 2);
			//std::cout << "i is " << i << ", Val is " << tmp << std::endl;
			//if (i != 3) {
			m_frac += (float) tmp2;
			m_frac *= 0.5; //pow(4.0,-1);
			//}
			m_frac += (float) tmp;
			m_frac *= 0.5; //pow(4.0,-1);


		}	
			
		// Handle Exponent
	
		DCRTPoly exp_result;
		m_RLWEscheme->Decrypt(m_params->GetRLWEParams(), 
				sk, 
				std::make_shared<FPRLWECiphertextImpl>(FPRLWECiphertextImpl(
					ct->GetAExpo(), ct->GetBExpo()
				)),
				&exp_result);
		exp_result = this->DivideRounding(&exp_result, Format::COEFFICIENT);
		int64_t m_exp = DecodingExp(exp_result);
		m_exp -= 63;
		m_exp <<= 1;	

		// Handle Sign
		vector<DCRTPoly> Sign_A = ct->GetASign();
		DCRTPoly Sign_B = ct->GetBSign();
		Sign_B.SetFormat(Format::COEFFICIENT);
		
		Sign_B.GetElementW(0)[0].ModAddEq(
			m_params->GetRLWEParams()->GetModuliQ()[1],
			m_params->GetRLWEParams()->GetModuliQ()[0]);
		
		Sign_B.SetFormat(Format::EVALUATION);	
		uint64_t msg_signed = DecryptUInt64(
				sk,
				std::make_shared<FPRLWECiphertextImpl>(FPRLWECiphertextImpl(
				std::move(Sign_A), std::move(Sign_B)
					)),
				DIVIDE
				);
		
		if (msg_signed == 2) {
			return (m_frac * ((float)pow(2.0, m_exp)));
		} else if (msg_signed == 0) {
			return - (m_frac * (float) pow(2.0, m_exp));
		} else {
			std::cout <<"Sign Encryption is Invalid!!!!" << std::endl;
			return 0.0;
		}

	}




	float FPFHEContext::DecryptHalf(
			FPConstRLWEPrivateKey sk, 
			FPConstRLWECiphertext ct,
			const bool DIVIDE
			) const {

		// Handle Frac
		DCRTPoly frac_result;	
		m_RLWEscheme->Decrypt(m_params->GetRLWEParams(), sk, ct, &frac_result);
		if (DIVIDE == true) {
			frac_result = this->DivideRounding(&frac_result, Format::COEFFICIENT);
		}else {
			frac_result.SetFormat(Format::COEFFICIENT);
		}

		// Setting
		NativeInteger Q = m_params -> GetRLWEParams() -> GetModuliQ()[0];
		float m_frac = 0;
		uint64_t tmp_res = 0;

		NativePoly result0 = frac_result.GetElementAtIndex(0);
		result0.SetFormat(Format::COEFFICIENT);	
		uint32_t tmp;
		// 2 ~ 7 are our informations 
		for (uint32_t i = 2; i < 8; i++) {
		
			if (result0[i] == Q-1) {
				result0[i] = NativeInteger(0);
			}
			
			tmp_res = result0[i].ConvertToInt() / m_params->GetRLWEParams()->GetScaling().ConvertToInt();
			tmp_res = (tmp_res+m_params->GetRLWEParams()->GetBits_over_2()) >> m_params->GetRLWEParams()->GetBits()[1];	
			tmp = (tmp_res % 4);
			//std::cout << "i is " << i << ", Val is " << tmp << std::endl;
			m_frac += (float) tmp;
			m_frac *= 0.25; //pow(4.0,-1);
		}	
		//std::cout << "frac is " << m_frac << std::endl;
		// Handle Exponent
	
		DCRTPoly exp_result;
		//std::cout << "ct is " << ct << std::endl:
		m_RLWEscheme->Decrypt(m_params->GetRLWEParams(), 
				sk, 
				std::make_shared<FPRLWECiphertextImpl>(FPRLWECiphertextImpl(
					ct->GetAExpo(), ct->GetBExpo()
				)),
				&exp_result);
		//std::cout << "Decrypt is end " << m_frac << std::endl;
	
		exp_result = this->DivideRounding(&exp_result, Format::COEFFICIENT);
		//std::cout << "Rounding is end " << m_frac << std::endl;
		int64_t m_exp = DecodingExp(exp_result);
		m_exp -= 16;
		m_exp <<= 1;	

		//std::cout << "m_exp is " << m_exp << std::endl;
		// Handle Sign
		vector<DCRTPoly> Sign_A = ct->GetASign();
		DCRTPoly Sign_B = ct->GetBSign();
		Sign_B.SetFormat(Format::COEFFICIENT);
		
		Sign_B.GetElementW(0)[0].ModAddEq(
			m_params->GetRLWEParams()->GetModuliQ()[1],
			m_params->GetRLWEParams()->GetModuliQ()[0]);
		
		Sign_B.SetFormat(Format::EVALUATION);	
		uint64_t msg_signed = DecryptUInt64(
				sk,
				std::make_shared<FPRLWECiphertextImpl>(FPRLWECiphertextImpl(
				std::move(Sign_A), std::move(Sign_B)
					)),
				DIVIDE
				);
		
		//std::cout << "msg_signed is  " << msg_signed << std::endl;
	
		if (msg_signed == 2) {
			return (m_frac * ((float) pow(2.0, m_exp)));
		} else if (msg_signed == 0) {
			return - (m_frac * ((float) pow(2.0, m_exp)));
		} else {
			std::cout <<"Sign Encryption is Invalid!!!!" << std::endl;
			return 0.0;
		}

	}



	FPRLWECiphertext FPFHEContext::EncryptDouble(FPConstRLWEPrivateKey sk, double &m) const {
		// Settings
		uint32_t N = m_params->GetRLWEParams()->GetN();
		vector<NativeInteger> moduliQ = m_params->GetRLWEParams()->GetModuliQ();
		// moduliQ2 is scaling factor	
		NativeVector msg_vector(N);
		msg_vector.SetModulus(moduliQ[0]);
		// Currents base are 2bit
	
		DCRTPoly msg_poly = DCRTPoly(m_params->GetRLWEParams()->GetDCRTParams(), Format::COEFFICIENT,true);
		//std::cout << "Input double m is " << m << std::endl;

		int64_t m_sign;
		double m_frac_double;
		double m_abs;
		int m_exp;
		//double abs_m;
		if (!signbit(m)) {
			//m =(uint64_t) signed_m;
			m_sign = 1;
			m_abs = m;
		} else {
			//m = (uint64_t) (-signed_m);
			m_sign = -1;
			//m_abs = -1.0 * m;
			m_abs = - m;
			//abs_m = -m;
		}
		m_frac_double = frexp(m_abs, &m_exp);


		//std::cout << " sign is " << m_sign << ", frac num is " << m_frac_double << ", exponent is " << m_exp << std::endl;
		if ((m_exp > 1023) || (m_exp < -1022)) {
			PALISADE_THROW(config_error , " double representation overflow (it shoud be in 2^-1024 ~ 2^1023)" );
		} else if (m_exp == 1023) {
			PALISADE_THROW(config_error , "FPFHE 64bit cannot encrypt when exp = 511");
		}
		
		//double recon = m_sign * m_frac_double * pow(2.0, m_exp);
		//std::cout << "reconstructed values are" << recon << std::endl;
		// Make IntFraction
		int32_t tmps = (m_exp % 2);
		if (tmps == -1) {tmps = 1;}
		//std::cout << "tmps are " << tmps << std::endl;
		
		if (tmps == 1) {m_exp += 1;}
		m_exp >>= 1; // Make m_exp base 4, 
		
		uint32_t up;
		uint32_t res_tmp = 0;
		for (uint32_t i = 0; i < 27; i++) {
			uint32_t save_idx = 31 - i;
			
			if ((i == 0) && (tmps == 1)) { // Deviding frac by 2
				m_frac_double *= 2.0;
				up =  (uint32_t) m_frac_double;
				res_tmp += up;
				m_frac_double -= up; 
			} else {
			
				m_frac_double *= 4.0;
				up =  (uint32_t) m_frac_double;
				res_tmp += up;
				m_frac_double -= up;
			}
			msg_poly.GetElementW(0)[save_idx] = res_tmp * moduliQ[1];
			res_tmp = 0;
		}

		DCRTPoly exp_poly = EncodingExp((uint64_t)(m_exp + 511));
		DCRTPoly sign_poly = EncodingSign(m_sign);
		
		auto CT_msg = m_RLWEscheme->Encrypt(m_params->GetRLWEParams(), sk, msg_poly);
		auto CT_sign = m_RLWEscheme->Encrypt(m_params->GetRLWEParams(), sk, sign_poly);
		auto CT_exp = m_RLWEscheme->Encrypt(m_params->GetRLWEParams(), sk, exp_poly);
		
		return std::make_shared<FPRLWECiphertextImpl>(FPRLWECiphertextImpl(
				
				std::move(CT_msg->GetA()), 
				std::move(CT_msg->GetB()),
				std::move(CT_sign->GetA()), 
				std::move(CT_sign->GetB()),
				std::move(CT_exp->GetA()), 
				std::move(CT_exp->GetB())
				));
	}






	FPRLWECiphertext FPFHEContext::EncryptFloat(FPConstRLWEPrivateKey sk, float &m) const {
		// Settings
		uint32_t N = m_params->GetRLWEParams()->GetN();
		vector<NativeInteger> moduliQ = m_params->GetRLWEParams()->GetModuliQ();
		// moduliQ2 is scaling factor	
		NativeVector msg_vector(N);
		msg_vector.SetModulus(moduliQ[0]);
		// Currents base are 2bit
	
		DCRTPoly msg_poly = DCRTPoly(m_params->GetRLWEParams()->GetDCRTParams(), Format::COEFFICIENT,true);
		//std::cout << "Input double m is " << m << std::endl;

		int64_t m_sign;
		float m_frac_float;
		float m_abs;
		int m_exp;
		//double abs_m;
		if (!signbit(m)) {
			//m =(uint64_t) signed_m;
			m_sign = 1;
			m_abs = m;
		} else {
			//m = (uint64_t) (-signed_m);
			m_sign = -1;
			//m_abs = -1.0 * m;
			m_abs = - m;
			//abs_m = -m;
		}
		m_frac_float = frexp(m_abs, &m_exp);

		//std::cout << " sign is " << m_sign << ", frac num is " << m_frac_float << ", exponent is " << m_exp 	<< ", recon is " << m_sign * m_frac_float * pow(2,m_exp) << std::endl;
		if ((m_exp > 128) || (m_exp < -126)) {
			PALISADE_THROW(config_error , " double representation overflow (it shoud be in 2^-126 ~ 2^127)" );
		} 
		//else if (m_exp == 125 && m_exp == 126) {
		//	PALISADE_THROW(config_error , "FPFHE 64bit cannot encrypt when exp = 127");
		//}
		
		//double recon = m_sign * m_frac_double * pow(2.0, m_exp);
		//std::cout << "reconstructed values are" << recon << std::endl;
		// Make IntFraction
		int32_t tmps = (m_exp % 2);
		if (tmps == -1) {tmps = 1;}
		//std::cout << "tmps are " << tmps << std::endl;
		
		if (tmps == 1) {m_exp += 1;}
		m_exp >>= 1; // Make m_exp base 4, 
	
		uint32_t up;
		uint32_t res_tmp = 0;
		for (uint32_t i = 0; i < 13; i++) {
			uint32_t save_idx = 15 - i;
			
			if ((i == 0) && (tmps == 1)) { // Deviding frac by 2
				m_frac_float *= 2.0;
				up =  (uint32_t) m_frac_float;
				res_tmp += up;
				m_frac_float -= up; 
			} else {
			
				m_frac_float *= 4.0;
				up =  (uint32_t) m_frac_float;
				res_tmp += up;
				m_frac_float -= up;
			}
			msg_poly.GetElementW(0)[save_idx] = res_tmp * moduliQ[1];
			res_tmp = 0;
			//std::cout << i << " th frac is "
		}
		//std::cout << "Remain float is " << m_frac_float << std::endl;
		DCRTPoly exp_poly = EncodingExp((uint64_t)(m_exp + 63));
		DCRTPoly sign_poly = EncodingSign(m_sign);
		
		auto CT_msg = m_RLWEscheme->Encrypt(m_params->GetRLWEParams(), sk, msg_poly);
		auto CT_sign = m_RLWEscheme->Encrypt(m_params->GetRLWEParams(), sk, sign_poly);
		auto CT_exp = m_RLWEscheme->Encrypt(m_params->GetRLWEParams(), sk, exp_poly);
		
		return std::make_shared<FPRLWECiphertextImpl>(FPRLWECiphertextImpl(
				
				std::move(CT_msg->GetA()), 
				std::move(CT_msg->GetB()),
				std::move(CT_sign->GetA()), 
				std::move(CT_sign->GetB()),
				std::move(CT_exp->GetA()), 
				std::move(CT_exp->GetB()),
				CTType::FLOAT
				));
	}




	FPConstRLWECiphertext FPFHEContext::EncryptHalf(FPConstRLWEPrivateKey sk, float &m) const {
		// Settings
		uint32_t N = m_params->GetRLWEParams()->GetN();
		vector<NativeInteger> moduliQ = m_params->GetRLWEParams()->GetModuliQ();
		// moduliQ2 is scaling factor	
		NativeVector msg_vector(N);
		msg_vector.SetModulus(moduliQ[0]);
		// Currents base are 2bit
	
		DCRTPoly msg_poly = DCRTPoly(m_params->GetRLWEParams()->GetDCRTParams(), Format::COEFFICIENT,true);
		//std::cout << "Input double m is " << m << std::endl;

		int64_t m_sign;
		float m_frac_half;
		float m_abs;
		int m_exp;
		//double abs_m;
		if (!signbit(m)) {
			//m =(uint64_t) signed_m;
			m_sign = 1;
			m_abs = m;
		} else {
			//m = (uint64_t) (-signed_m);
			m_sign = -1;
			//m_abs = -1.0 * m;
			m_abs = - m;
			//abs_m = -m;
		}
		m_frac_half = frexp(m_abs, &m_exp);


		//std::cout << " sign is " << m_sign << ", frac num is " << m_frac_half << ", exponent is " << m_exp << std::endl;
		if ((m_exp >= 16) || (m_exp < -16)) {
			PALISADE_THROW(config_error , " Half representation overflow (it shoud be in 2^-16 ~ 2^16)" );
		} else if (m_exp == 16) {
			PALISADE_THROW(config_error , "FPFHE 16bit cannot encrypt when exp = 16");
		}
		
		//double recon = m_sign * m_frac_double * pow(2.0, m_exp);
		//std::cout << "reconstructed values are" << recon << std::endl;
		// Make IntFraction
		int32_t tmps = (m_exp % 2);
		if (tmps == -1) {tmps = 1;}
		//std::cout << "tmps are " << tmps << std::endl;
		
		if (tmps == 1) {m_exp += 1;}
		m_exp >>= 1; // Make m_exp base 4, 
		uint32_t up;
		uint32_t res_tmp = 0;
		for (uint32_t i = 0; i < 6; i++) {
			uint32_t save_idx = 7 - i;
			
			if ((i == 0) && (tmps == 1)) { // Deviding frac by 2
				m_frac_half *= 2.0;
				up =  (uint32_t) m_frac_half;
				res_tmp += up;
				m_frac_half -= up; 
			} else {
			
				m_frac_half *= 4.0;
				up =  (uint32_t) m_frac_half;
				res_tmp += up;
				m_frac_half -= up;
			}
			msg_poly.GetElementW(0)[save_idx] = res_tmp * moduliQ[1];
			res_tmp = 0;
		}

		DCRTPoly exp_poly = EncodingExp((uint64_t)(m_exp + 16));
		DCRTPoly sign_poly = EncodingSign(m_sign);
		
		auto CT_msg = m_RLWEscheme->Encrypt(m_params->GetRLWEParams(), sk, msg_poly);
		auto CT_sign = m_RLWEscheme->Encrypt(m_params->GetRLWEParams(), sk, sign_poly);
		auto CT_exp = m_RLWEscheme->Encrypt(m_params->GetRLWEParams(), sk, exp_poly);
		
		return std::make_shared<FPRLWECiphertextImpl>(FPRLWECiphertextImpl(
				std::move(CT_msg->GetA()), 
				std::move(CT_msg->GetB()),
				std::move(CT_sign->GetA()), 
				std::move(CT_sign->GetB()),
				std::move(CT_exp->GetA()), 
				std::move(CT_exp->GetB()),
				CTType::HALF
				));
	}


	DCRTPoly FPFHEContext::EncodingExp( uint64_t m) const {
		//std::cout << "m is " << m << std::endl;
		uint32_t N = m_params->GetRLWEParams()->GetN();
		//uint64_t delta =		
		vector<NativeInteger> moduliQ = m_params->GetRLWEParams()->GetModuliQ();
		NativeInteger Q0 = moduliQ[0];
		NativeInteger Q1 = moduliQ[1];
		// moduliQ2 is scaling factor	
		NativeVector msg_vector(N);
		msg_vector.SetModulus(moduliQ[0]);
		//msg_vector.SetValuesToZero();
		// 64 bit Encoding
		
		const uint32_t scaling_bit = m_params->GetExpBit();
		uint64_t scaling =  m_params->GetRLWEParams()->GetScaling().ConvertToInt();
		scaling <<= scaling_bit;

		DCRTPoly msg_poly = DCRTPoly(m_params->GetRLWEParams()->GetDCRTParams(), Format::COEFFICIENT,true);
		msg_poly.GetElementW(0)[0].ModAddEq(
				Q1.ModMul(scaling,Q0).ModMul(m,Q0), Q0);
		msg_poly.GetElementW(1)[0].ModAddEq(
				Q1.ModMul(scaling,Q1).ModMul(m,Q1), Q1);
		
		return msg_poly;
	}
	




	DCRTPoly FPFHEContext::EncodingUInt64( uint64_t m) const {
		
		uint32_t N = m_params->GetRLWEParams()->GetN();
		vector<NativeInteger> moduliQ = m_params->GetRLWEParams()->GetModuliQ();
		// moduliQ2 is scaling factor	
		NativeVector msg_vector(N);
		msg_vector.SetModulus(moduliQ[0]);
		//msg_vector.SetValuesToZero();
		// 64 bit Encoding
		uint32_t len_encode = m_params->GetLenEnc()[6];
		uint32_t base_encode_bits = m_params->GetBaseEncodeBit();
		uint32_t base_encode = 1 << base_encode_bits;
		
		DCRTPoly msg_poly = DCRTPoly(m_params->GetRLWEParams()->GetDCRTParams(), Format::COEFFICIENT,true);
		
		for (uint32_t i = 0; i < len_encode; i++) {
			//std::cout << i << "'st m is " << m << std::endl;		
			//msg_vector[i] +=  (m % base_encode) *moduliQ[1];
			msg_poly.GetElementW(0)[i] = (m % base_encode) * moduliQ[1];
			//std::cout << i << "'st m mod "<< base_encode<< " is " << msg_vector[i] << std::endl;	
			m = m >> base_encode_bits;
		}
		
		return msg_poly;
	}
	


	FPConstRLWECiphertext FPFHEContext::EncryptUInt64(FPConstRLWEPrivateKey sk,
                                     uint64_t &m) const {
		DCRTPoly msg_poly = EncodingUInt64(m);
		return m_RLWEscheme->Encrypt(m_params->GetRLWEParams(), sk, msg_poly);
	}

	FPConstRLWECiphertext FPFHEContext::EncryptInt64(FPConstRLWEPrivateKey sk,
                                     int64_t &signed_m) const {
		int64_t sign;
		uint64_t m;
		if (signed_m >= 0 ) {
			m =(uint64_t) signed_m;
			sign = 1;
		} else {
			m = (uint64_t) (-signed_m);
			sign = -1;
		}
		
		DCRTPoly msg_poly = EncodingUInt64(m);
		DCRTPoly sign_poly = EncodingSign(sign);
		
		auto CT_msg = m_RLWEscheme->Encrypt(m_params->GetRLWEParams(), sk, msg_poly);
		auto CT_sign = m_RLWEscheme->Encrypt(m_params->GetRLWEParams(), sk, sign_poly);
		


		return std::make_shared<FPRLWECiphertextImpl>(FPRLWECiphertextImpl(
				std::move(CT_msg->GetA()), 
				std::move(CT_msg->GetB()),
				std::move(CT_sign->GetA()), 
				std::move(CT_sign->GetB())));
	}

	DCRTPoly FPFHEContext::EncodingSign( int64_t sign) const {

		if ((sign != 1) && (sign != -1) ) {
			std::cout << "!!!!!!!!!!!!!!  Sign Encoding is error!!!!!!!!!!!!!!" << std::endl;
		}
		uint32_t N = m_params->GetRLWEParams()->GetN();
		vector<NativeInteger> moduliQ = m_params->GetRLWEParams()->GetModuliQ();
		
		NativePoly m_poly1 = NativePoly( m_params->GetRLWEParams()->GetDCRTParams()->GetParams()[0], Format::COEFFICIENT);
		NativePoly m_poly2 = NativePoly( m_params->GetRLWEParams()->GetDCRTParams()->GetParams()[1], Format::COEFFICIENT, true);
		NativeVector msg_vector(N);
		msg_vector.SetModulus(moduliQ[0]);
	
		if (sign == 1) {
			//Set MSG 1
			msg_vector[0] = msg_vector[0].ModAdd(moduliQ[1] , moduliQ[0]);
		} else if (sign == -1) {
			//Set MSG 1
			msg_vector[0] = msg_vector[0].ModSub(moduliQ[1] , moduliQ[0]);
		}
		m_poly1.SetValues(msg_vector, Format::COEFFICIENT);
	
		DCRTPoly msg_poly = DCRTPoly(m_params->GetRLWEParams()->GetDCRTParams(), Format::COEFFICIENT,true);
		msg_poly.SetElementAtIndex(0, m_poly1);
		msg_poly.SetElementAtIndex(1, m_poly2);
	
		return msg_poly;
	}



	/*
	void FPFHEContext::Decrypt(FPConstLWEPrivateKey sk, FPConstLWECiphertext ct,
                            FPLWEPlaintext_double *result) const {
		return m_LWEscheme->Decrypt(m_params->GetLWEParams(), sk, ct, result);
	}
	*/

/* SH
	void FPFHEContext::RDecrypt(FPConstRLWEPrivateKey sk, FPConstRLWECiphertext ct,
                            DCRTPoly *result) const {
		return m_RLWEscheme->Decrypt(m_params->GetRLWEParams(), sk, ct, result);
	}
	DCRTPoly FPFHEContext::MakeDG() const {
		return DCRTPoly(m_params -> GetRLWEParams() -> GetDgg(), m_params-> GetRLWEParams() -> GetDCRTParams(), Format::COEFFICIENT);
	}

*/
	
	uint64_t FPFHEContext::DecodingUInt64(DCRTPoly poly) const {	
		NativeInteger Q = m_params -> GetRLWEParams() -> GetModuliQ()[0];
		uint64_t m = 0;
		uint64_t tmp_res = 0;
		
		uint32_t len_encode = m_params->GetLenEnc()[6];
		uint32_t base_encode_bits = m_params->GetBaseEncodeBit();
		uint32_t base_encode = 1;
						 base_encode <<= base_encode_bits;
	
		NativePoly result0 = poly.GetElementAtIndex(0);
		result0.SetFormat(Format::COEFFICIENT);
	

		for (int32_t i = len_encode-1; i >= 0; i--) {
			m <<= base_encode_bits;
			if (result0[i] == Q-1) {
				result0[i] = NativeInteger(0);
			}
			
			tmp_res = result0[i].ConvertToInt() / m_params->GetRLWEParams()->GetScaling().ConvertToInt();
			tmp_res = (tmp_res+m_params->GetRLWEParams()->GetBits_over_2()) >> m_params->GetRLWEParams()->GetBits()[1];
				
			m += (tmp_res % base_encode);
		}	
			
		return m;
	}

	
	uint64_t FPFHEContext::DecodingExp(DCRTPoly poly) const {	
		// Must SCaling Downs!
		NativeInteger Q = m_params -> GetRLWEParams() -> GetModuliQ()[0];
		const uint32_t scaling_bit = m_params -> GetExpBit();
		uint64_t delta = m_params->GetRLWEParams()->GetScaling().ConvertToInt();
	
		uint64_t scaling = delta;
		scaling <<= scaling_bit;
		uint64_t scaling_over_2 = scaling >> 1;
		
		//std::cout << "delta is" << delta << ", scaling is " << scaling << ", Q is" << Q << std::endl;

		NativePoly result0 = poly.GetElementAtIndex(0);
		result0.SetFormat(Format::COEFFICIENT);
		result0[0].ModAddEq(scaling_over_2, Q);
		

		//std::cout << "m is " << result0[0] << std::endl;
		//uint64_t m = ((result0[0].ConvertToInt() % scaling) <= scaling_over_2 ) ?
		//					( (result0[0].ConvertToInt() >> scaling_bit)  ) : ( (result0[0].ConvertToInt() >> scaling_bit) + 1 )	;
		uint64_t m =  (result0[0].ConvertToInt() / scaling);
		
		//std::cout << "EXP is " << m << std::endl;
		return m;
	}



	int64_t FPFHEContext::DecryptMulPoly(
		FPConstRLWEPrivateKey sk, 
		vector<DCRTPoly> cts,
		const uint32_t mul_max,
		const bool mul_flag,
		const bool DIVIDE) const {

		DCRTPoly result;	
		uint32_t K_FP = m_params->GetRLWEParams()->GetK();
		m_RLWEscheme->Decrypt(m_params->GetRLWEParams(), sk, 
			std::make_shared<FPRLWECiphertextImpl>(FPRLWECiphertextImpl(
			cts, cts[K_FP])),
			&result);
			
		//NativePoly result0;
		if (DIVIDE == true) {
			result = this->DivideRounding(&result, Format::COEFFICIENT);
			//std::cout << "after result is " << result << std::endl;
		} else {
			result.SetFormat(Format::COEFFICIENT);
			//std::cout << "after result is " << result.GetElementAtIndex(0) << std::endl;
		}
		
		if (mul_flag == false) {
			//uint32_t rot_idx = m_params -> GetRLWEParams()->GetN() * 2;
			//rot_idx = (rot_idx - mul_max);
			result.SetFormat(Format::EVALUATION);
			result *= m_params->GetRotMonomial_RLWE(mul_max);
			result.SetFormat(Format::COEFFICIENT);
		
		}

		NativeInteger Q0 = m_params -> GetRLWEParams() -> GetModuliQ()[0];
		NativeInteger Q1 = m_params -> GetRLWEParams() -> GetModuliQ()[1];
		int64_t m = -10000;
		uint64_t tmp_res = 0;
		uint64_t tmp_res2 = 0;
		uint32_t base_encode_bits = m_params->GetBaseEncodeBit();
		uint32_t base_encode = 1;
						 base_encode <<= base_encode_bits;
	
		NativePoly result0 = result.GetElementAtIndex(0);
		result0.SetFormat(Format::COEFFICIENT);

		int32_t signs = 1;
		for (uint32_t i = 0; i <= mul_max; i++) {
		
			// Reverse
			NativeInteger bit_rev = result0[i].ModMul(Q0 - 1,Q0);
			bit_rev.ModAddEq((Q1.ConvertToInt() >> 1), Q0);
			// Orign
			result0[i].ModAddEq((Q1.ConvertToInt() >> 1),Q0);
			

			tmp_res = (uint64_t) (result0[i].ConvertToDouble() / (Q1.ConvertToDouble()));
		
			tmp_res2 = (uint64_t) (bit_rev.ConvertToDouble() / (Q1.ConvertToDouble()));
		

			std::cout << "Mul i is " << i 
				<< ", val is " << tmp_res 
				<< ", rev val is " << tmp_res2 << std::endl;
	


			if (tmp_res == 1 && m == -10000) {
				m = i;
			} else if (tmp_res2 == 1 && m == -10000) {
				m = i;
				signs = -1;
			} else if (tmp_res == 1 && (m != -10000)) {
				m = - 20000;
			} else if (tmp_res2 == 1 && (m != -10000)) {
				m = - 20000;
			}


		}	
		if (mul_flag == true) {
			return m*signs;
		} else {
			return (mul_max - m) * signs;
		}	
	}


	
	uint64_t FPFHEContext::DecryptCompensate(
		FPConstRLWEPrivateKey sk, 
		vector<DCRTPoly> cts,
		const bool DIVIDE) const {

		DCRTPoly result;	
		uint32_t K_FP = m_params->GetRLWEParams()->GetK();
		m_RLWEscheme->Decrypt(m_params->GetRLWEParams(), sk, 
			std::make_shared<FPRLWECiphertextImpl>(FPRLWECiphertextImpl(
			cts, cts[K_FP])),
			&result);
			
		//NativePoly result0;
		if (DIVIDE == true) {
			result = this->DivideRounding(&result, Format::COEFFICIENT);
			//std::cout << "after result is " << result << std::endl;
		} else {
			result.SetFormat(Format::COEFFICIENT);
			//std::cout << "after result is " << result.GetElementAtIndex(0) << std::endl;
		}
		
		NativeInteger Q0 = m_params -> GetRLWEParams() -> GetModuliQ()[0];
		NativeInteger Q1 = m_params -> GetRLWEParams() -> GetModuliQ()[1];
		uint64_t m = 0;
		uint64_t tmp_res = 0;
		
		uint32_t base_encode_bits = m_params->GetBaseEncodeBit();
		uint32_t base_encode = 1;
						 base_encode <<= base_encode_bits;
	
		NativePoly result0 = result.GetElementAtIndex(0);
		result0.SetFormat(Format::COEFFICIENT);

		
		for (int32_t i = 0; i <= 31; i++) {
		
			result0[i].ModAddEq((Q1.ConvertToInt() >> 1),Q0);
			
			tmp_res = (uint64_t) (result0[i].ConvertToDouble() / (Q1.ConvertToDouble()));
			m += tmp_res;
		}
		return m;
	}





	uint64_t FPFHEContext::DecryptUInt64(
		FPConstRLWEPrivateKey sk, 
		FPConstRLWECiphertext ct,
		const bool DIVIDE
		)  const {

		DCRTPoly result;	
		m_RLWEscheme->Decrypt(m_params->GetRLWEParams(), sk, ct, &result);
			//NativePoly result0;
		if (DIVIDE == true) {
			result = this->DivideRounding(&result, Format::COEFFICIENT);
			//std::cout << "after result is " << result << std::endl;
		}else {
			result.SetFormat(Format::COEFFICIENT);
			//std::cout << "after result is " << result.GetElementAtIndex(0) << std::endl;
		}
		return DecodingUInt64(result);			
	}


	int64_t FPFHEContext::DecryptInt64(
			FPConstRLWEPrivateKey sk, 
			FPConstRLWECiphertext ct,
			const bool DIVIDE
			) {

		int64_t msg_unsigned = DecryptUInt64(
				sk,
				std::make_shared<FPRLWECiphertextImpl>(FPRLWECiphertextImpl(
						ct->GetA(), ct->GetB()
					)),
				DIVIDE
				);
		//std::cout << "msg is " << msg_unsigned << std::endl;

		vector<DCRTPoly> Sign_A = ct->GetASign();
		DCRTPoly Sign_B = ct->GetBSign();
		Sign_B.SetFormat(Format::COEFFICIENT);
		//NativePoly Sign_B_poly = Sign_B.GetElementAtIndex(0);
		Sign_B.GetElementW(0)[0].ModAddEq(
			m_params->GetRLWEParams()->GetModuliQ()[1],
			m_params->GetRLWEParams()->GetModuliQ()[0]);
		
		//Sign_B_poly[0] = Sign_B_poly[0].ModAddEq(
		//		m_params->GetRLWEParams()->GetModuliQ()[1],
		//		m_params->GetRLWEParams()->GetModuliQ()[0]);
		//Sign_B.SetElementAtIndex(0, Sign_B_poly);
		Sign_B.SetFormat(Format::EVALUATION);	
		uint64_t msg_signed = DecryptUInt64(
				sk,
				std::make_shared<FPRLWECiphertextImpl>(FPRLWECiphertextImpl(
				std::move(Sign_A), std::move(Sign_B)
					)),
				DIVIDE
				);
		
		if (msg_signed == 2) {
			return msg_unsigned;
		} else if (msg_signed == 0) {
			return -msg_unsigned;
		} else {
			std::cout <<"Sign is Invalid!!!!" << std::endl;
			return 0;
		}

	}


	void FPFHEContext::CheckQ(
			uint32_t bits,
			uint32_t nums,
			vector<NativeInteger> *muls_arr,
			vector<uint32_t> *bits_arr,
			vector<uint32_t> *res_bits_arr
			) const {
		
		NativeInteger Q = PreviousPrime<NativeInteger>(FirstPrime<NativeInteger>(bits, 2048), 2048);
		NativeInteger Q_tmp;
		for (uint32_t i = 0; i < nums; i++) {
			Q = NextPrime(Q,2048);
			Q_tmp = Q - 1;
			while(true) {
				Q_tmp = Q_tmp / 2;
				if (Q_tmp % 2 == 1) {
					break;
				}
			}

			(*bits_arr)[i] = (uint32_t)std::ceil(log(Q.ConvertToDouble()) /	log(static_cast<double>(2)));
			(*res_bits_arr)[i] = (uint32_t)std::ceil(log(Q_tmp.ConvertToDouble()) /	log(static_cast<double>(2)));
			(*muls_arr)[i] = Q_tmp;
		}
		return;
	}

	/*
	std::shared_ptr<FPLWESwitchingKey> FPFHEContext::KeySwitchGen(
    FPConstLWEPrivateKey sk, FPConstLWEPrivateKey skN) const {
		return m_LWEscheme->KeySwitchGen(m_params->GetLWEParams(), sk, skN);
	}
	*/

	
	void FPFHEContext::BTKeyGen( FPConstRLWEPrivateKey sk_r) {
		
		m_BTKey = std::make_shared<FPRingGSWEvalKey> (m_RingGSWscheme->KeyGen(
					m_params, m_LWEscheme
					, m_RLWEscheme, sk_r));
		
		return;
	}

	FPConstRLWECiphertext FPFHEContext::BootstrapOrign(FPConstRLWECiphertext ct1) const {
			return m_RingGSWscheme->BootstrapOrign(m_params, m_LWEscheme, m_BTKey, ct1->GetA(), ct1->GetB());
	}
		 


FPConstRLWECiphertext FPFHEContext::ProductOnly(
		FPConstRLWECiphertext ct1, 
		FPConstRLWECiphertext ct2) const {
	
	std::shared_ptr<vector<DCRTPoly>> Res = m_RingGSWscheme->ProductDCRT(
		m_params, m_BTKey->Evkey,
		ct1->GetA(), ct1->GetB(), ct2->GetA(), ct2->GetB());
	DCRTPoly b = std::move((*Res)[Res->size() - 1]);
		//DCRTPoly b = std::move((*A)[A->size()-1]);
		(*Res).erase((*Res).begin()+((*Res).size()-1));	
		return std::make_shared<FPRLWECiphertextImpl> (FPRLWECiphertextImpl(std::move(*Res), b));
}


FPConstRLWECiphertext FPFHEContext::ProductUOrign(
		FPConstRLWECiphertext ct1, 
		FPConstRLWECiphertext ct2) const {
	

	std::shared_ptr<vector<DCRTPoly>> Res = m_RingGSWscheme->ProductDCRT(
		m_params, m_BTKey->Evkey,
		ct1->GetA(), ct1->GetB(), ct2->GetA(), ct2->GetB());
	
	for (uint32_t i = 0; i < Res->size(); i++ ) {
		(*Res)[i].SetFormat(Format::COEFFICIENT);
		m_RingGSWscheme->DivideRoundingSelf(m_params, &(*Res)[i], Format::COEFFICIENT);
	}
	
	DCRTPoly b = std::move((*Res)[Res->size() - 1]);
	(*Res).erase((*Res).begin()+((*Res).size()-1));
	return 	m_RingGSWscheme->BootstrapOrign(m_params, m_LWEscheme, m_BTKey, *Res, b);

}


	vector<DCRTPoly> FPFHEContext::SGDLevel2(const DCRTPoly &in_poly, uint32_t bits, uint32_t rm_len, Format format) {
		auto res = m_RingGSWscheme->SGDLevel2( m_params, in_poly, bits, rm_len, format);
		return std::move(*res);
	}



FPRLWECiphertext FPFHEContext::ProductDouble(
		FPConstRLWECiphertext ct1, 
		FPConstRLWECiphertext ct2,
		//FPConstRLWEPrivateKey sk,
		std::shared_ptr<FPOFUF> ct_OFUF
		//const FPOFUF &ct_OFUF
		) const {

	CTType ct1_t = ct1->GetType();
	CTType ct2_t = ct2->GetType();	
	if (ct1_t != ct2_t) {
		PALISADE_THROW(config_error, "Current Version does not support cross type operation");
	}
	CTType types = ct1_t;
	uint32_t K_FP = m_params->GetRLWEParams()->GetK();


	// Sign
	std::shared_ptr<vector<DCRTPoly>> Res2 = m_RingGSWscheme->ProductDCRT(
	m_params, m_BTKey->Evkey,
	ct1->GetASign(), ct1->GetBSign(), ct2->GetASign(), ct2->GetBSign());
	for (uint32_t i = 0; i < K_FP+1; i++ ) {
		(*Res2)[i].SetFormat(Format::COEFFICIENT);
		m_RingGSWscheme->DivideRoundingSelf(m_params, &(*Res2)[i], Format::COEFFICIENT);
	}
	
	std::shared_ptr<vector<DCRTPoly>> A_sign =  
	m_RingGSWscheme->BootstrapSign(m_params, m_LWEscheme, (*Res2), m_BTKey);
	DCRTPoly b_sign = std::move((*A_sign)[A_sign->size()-1]);
	A_sign->erase(A_sign->begin()+(A_sign->size()-1));


	// Frac
	std::shared_ptr<vector<DCRTPoly>> Res = m_RingGSWscheme->ProductDCRT(
		m_params, m_BTKey->Evkey,
		ct1->GetA(), ct1->GetB(), ct2->GetA(), ct2->GetB());
	
	for (uint32_t i = 0; i < K_FP+1; i++ ) {
		(*Res)[i].SetFormat(Format::COEFFICIENT);
		m_RingGSWscheme->DivideRoundingSelf(m_params, &(*Res)[i], Format::COEFFICIENT);
	}

	/*
	if (DEBUGSH) {
		std::cout << "Sign is over" << std::endl;
	}*/

	std::shared_ptr<vector<vector<DCRTPoly>>> Frac;
	if (types == CTType::DOUBLE) {
		Frac = m_RingGSWscheme->BootstrapMulDouble64(m_params, m_LWEscheme, (*Res), m_BTKey, 0, 31, 5, types);

	} else if (types == CTType::FLOAT) {
			Frac = m_RingGSWscheme->BootstrapMulDouble64(m_params, m_LWEscheme, (*Res), m_BTKey,0, 15, 3, types);
	
	}


	vector<DCRTPoly> Frac_Results = std::move((*Frac)[Frac->size() - 1]);  
	Frac->erase(Frac->begin() + (Frac->size() - 1));
	
	
	/*
	if (DEBUGSH) { 
		std::cout << "Mul is over" << std::endl;
		
		std::cout << "Frac inner searching..." << std::endl;
		int64_t Frac_inner_result = this ->DecryptMulPoly(sk, (Frac_Results), 32, true);
		std::cout << "After Frac inner result " << Frac_inner_result << std::endl;
		for (uint32_t sh = 0; sh < Frac->size(); sh++) {
			uint64_t tmpss= this ->DecryptUInt64(sk, 
			std::make_shared<FPRLWECiphertextImpl> (FPRLWECiphertextImpl(
			(*Frac)[sh], (*Frac)[sh][K_FP]) ));
			std::cout <<" Zero check, " << sh << "'st zero is " << tmpss << std::endl;			
		}

	}
	*/
	// Expo
	vector<DCRTPoly> A_expo(K_FP+1);
	for (uint32_t i = 0; i < K_FP; i++) {
		A_expo[i] = ct1->GetAExpo()[i] + ct2->GetAExpo()[i];
	}
	A_expo[K_FP] = ct1->GetBExpo() + ct2->GetBExpo();
	

	/*
	if (DEBUGSH) {
		std::cout << "prev Recon2, frac size is " << (*Frac)[0].size() << std::endl;
		std::cout << "Sum expo is "
			<< this->DecryptExp(sk, A_expo) << std::endl;
	}*/
	auto resExp2 = m_RingGSWscheme->BootstrapRecon2(m_params, m_LWEscheme, 
					(*Frac), m_BTKey, 0); 

	/*
	if (DEBUGSH) {
		std::cout << "Recon is over" << std::endl;
		int64_t resExp1dec = this->DecryptExp(sk, (*resExp2)[0]);
		//int64_t resExp2aa = this->DecryptExp(sk, (Exp2_d));
		std::cout << "want to Mul is  " << resExp1dec << std::endl;
		int64_t res25 = this->DecryptExp(sk, (*resExp2)[1]);
		std::cout << "Underflow suspect is  " << res25 << std::endl;

		//std::cout << "Exp res2 is " << resExp2aa << std::endl;
	}	*/

	auto OFUF = m_RingGSWscheme->CheckOFUFMul( m_params, m_LWEscheme,
					(*resExp2)[1], A_expo, types, m_BTKey);
	ct_OFUF->Adding((*OFUF));
	//Minus Bias 
	uint32_t scale_bit = m_params->GetExpBit();
	uint64_t prefix = m_params->GetRLWEParams()->GetScaling().ConvertToInt();
	NativeInteger Q0 = m_params->GetModuliQ()[0];
	NativeInteger Q1 = m_params->GetModuliQ()[1];
	prefix <<= scale_bit;
	if (types == CTType::DOUBLE) {
		prefix *= 511;
	} else if (types == CTType::FLOAT) {
		prefix *= 63;
	}

	A_expo[K_FP].SetFormat(Format::COEFFICIENT);
	A_expo[K_FP].GetElementW(0)[0].ModSubEq(Q1.ModMul(prefix, Q0), Q0);
	A_expo[K_FP].SetFormat(Format::EVALUATION);
	

	/*
	if (DEBUGSH) {
		int64_t res30 = this->DecryptExp(sk, (*OFUF)[0]);
		std::cout << "Overflow is  " << res30<< std::endl;
		int64_t res31 = this->DecryptExp(sk, (*OFUF)[1]);
		std::cout << "Underflow is  " << res31<< std::endl;
	}*/


	// Make Min		
	auto minExp3 =  m_RingGSWscheme->BootstrapReluRev(
					m_params, m_LWEscheme, (A_expo), (*resExp2)[0],  m_BTKey);			
	
	std::shared_ptr<vector<DCRTPoly>> MulPolys;		
	if (types == CTType::DOUBLE) {
		MulPolys =  m_RingGSWscheme->BootstrapMulPolyMake(
		m_params, m_LWEscheme, (*minExp3)[0], m_BTKey, 27, true);
	} else if (types == CTType::FLOAT) {
		MulPolys =  m_RingGSWscheme->BootstrapMulPolyMake(
		m_params, m_LWEscheme, (*minExp3)[0], m_BTKey, 13, true);	
	}
	
	/*
	if (DEBUGSH) {
		std::cout << "Recon Mul Process.." << std::endl;
		int64_t res44 = this ->DecryptMulPoly(sk, (*MulPolys), 27, true);
		std::cout << "Recon mul is  " << res44 << std::endl;
	}*/
	
	// Ending Exp
	for (uint32_t kk = 0; kk < K_FP+1; kk++) {
		A_expo[kk] -= (*minExp3)[0][kk];
	}
	// Expo Make and anding
	DCRTPoly b_expo = std::move(A_expo[K_FP]);
	A_expo.erase(A_expo.begin()+(A_expo.size()-1));


	// Muls
	std::shared_ptr<vector<DCRTPoly>> Muls1;
	if (m_params->GetPARALLEL()) {
		Muls1 = m_RingGSWscheme->ProductDCRTPARALLEL2(m_params, m_BTKey->Evkey, 
			Frac_Results, Frac_Results[K_FP],
			(*MulPolys), (*MulPolys)[K_FP]);

	} else {
		Muls1 = m_RingGSWscheme->ProductDCRT(m_params, m_BTKey->Evkey, 
			Frac_Results, Frac_Results[K_FP],
			(*MulPolys), (*MulPolys)[K_FP]);
	}	
	for (uint32_t i = 0; i < K_FP+1; i++) {
		(*Muls1)[i].SetFormat(Format::COEFFICIENT);
		m_RingGSWscheme->DivideRoundingSelf(m_params, &(*Muls1)[i], Format::COEFFICIENT);
	}		
	
	std::shared_ptr<vector<DCRTPoly>> A;
			if (types == CTType::DOUBLE) {
				A = m_RingGSWscheme->BootstrapOrignPartial(m_params, m_LWEscheme, m_BTKey, (*Muls1), 5, 31);
			
			} else if (types == CTType::FLOAT) {
				A = m_RingGSWscheme->BootstrapOrignPartial(m_params, m_LWEscheme, m_BTKey, (*Muls1), 3, 15);
			}
	/*
	if (DEBUGSH) {
		std::cout << "Final Frac inner searching..." << std::endl;
		int64_t Frac_inner_result2 = this ->DecryptMulPoly(sk, (*A), 32, true);
		std::cout << "Final Frac inner result " << Frac_inner_result2 << std::endl;
	}*/

	DCRTPoly b = std::move((*A)[K_FP]);
	A->erase(A->begin()+(A->size()-1));

		
	// Make Min
	return std::make_shared<FPRLWECiphertextImpl> (FPRLWECiphertextImpl(
		std::move(*A), std::move(b), std::move(*A_sign), std::move(b_sign), 
		std::move(A_expo), std::move(b_expo), types));

}



FPRLWECiphertext FPFHEContext::Product(
		FPConstRLWECiphertext ct1, 
		FPConstRLWECiphertext ct2,
		//FPConstRLWEPrivateKey Rsk,
		//const FPOFUF &ct_OFUF
		std::shared_ptr<FPOFUF> ct_OFUF
		) const {

	CTType ct1_t = ct1->GetType();
	CTType ct2_t = ct2->GetType();	
	if (ct1_t != ct2_t) {
		PALISADE_THROW(config_error, "Current Version does not support cross type operation");
	} 
	
	if(ct1_t == CTType::DOUBLE || ct1_t == CTType::FLOAT) {
		return this->ProductDouble(ct1, ct2, ct_OFUF);
	}


	std::shared_ptr<vector<DCRTPoly>> Res = m_RingGSWscheme->ProductDCRT(
		m_params, m_BTKey->Evkey,
		ct1->GetA(), ct1->GetB(), ct2->GetA(), ct2->GetB());
	
	auto mul_start = std::chrono::system_clock::now();
	for (uint32_t i = 0; i < Res->size(); i++ ) {
		(*Res)[i].SetFormat(Format::COEFFICIENT);
		m_RingGSWscheme->DivideRoundingSelf(m_params, &(*Res)[i], Format::COEFFICIENT);
	}
	
	std::chrono::duration<double> mul_end = std::chrono::system_clock::now() - mul_start;
	std::cout << "Before mul Bootstrapping time is " << mul_end.count() << std::endl     ;	

	std::shared_ptr<vector<DCRTPoly>> A =  
		m_RingGSWscheme->BootstrapMul(m_params, m_LWEscheme, (*Res), m_BTKey);



	DCRTPoly b = std::move((*A)[A->size()-1]);
	A->erase(A->begin()+(A->size()-1));
	
	if (ct1_t == CTType::UINT64) {
		return std::make_shared<FPRLWECiphertextImpl> (FPRLWECiphertextImpl(std::move(*A), b));
	} 
	//else if (ct1_t == CTType::INT64) {
  else {
		std::shared_ptr<vector<DCRTPoly>> Res2 = m_RingGSWscheme->ProductDCRT(
			m_params, m_BTKey->Evkey,
			ct1->GetASign(), ct1->GetBSign(), ct2->GetASign(), ct2->GetBSign());
	
		for (uint32_t i = 0; i < Res->size(); i++ ) {
			(*Res2)[i].SetFormat(Format::COEFFICIENT);
			m_RingGSWscheme->DivideRoundingSelf(m_params, &(*Res2)[i], Format::COEFFICIENT);
		}
	
		std::shared_ptr<vector<DCRTPoly>> A_sign =  
			m_RingGSWscheme->BootstrapSign(m_params, m_LWEscheme, (*Res2), m_BTKey);
		DCRTPoly b_sign = std::move((*A_sign)[A_sign->size()-1]);
		A_sign->erase(A_sign->begin()+(A_sign->size()-1));
	
		return std::make_shared<FPRLWECiphertextImpl> (FPRLWECiphertextImpl(std::move(*A), b, std::move(*A_sign), b_sign));
	}
}
 




		FPRLWECiphertext FPFHEContext::AddCore(
			FPConstRLWECiphertext ct1, 
			FPConstRLWECiphertext ct2,
			//FPConstRLWEPrivateKey sk,
			std::shared_ptr<FPOFUF> ct_OFUF,
			bool adding
			) const {

			CTType types = ct1->GetType();
			uint32_t K_FP = m_params->GetRLWEParams()->GetK();
			vector<DCRTPoly> tmp(K_FP+1);
			DCRTPoly b_tmp;
			std::shared_ptr<vector<DCRTPoly>> resExp;


			// Max
			if (types == CTType::DOUBLE) {
				resExp =  m_RingGSWscheme->BootstrapRelu(
					m_params, m_LWEscheme, 
					ct1->GetAExpo(), ct1->GetBExpo(), ct2->GetAExpo(), ct2->GetBExpo(),  m_BTKey);
			} else if (types == CTType::FLOAT) {
				resExp =  m_RingGSWscheme->BootstrapReluF(
					m_params, m_LWEscheme, 
					ct1->GetAExpo(), ct1->GetBExpo(), ct2->GetAExpo(), ct2->GetBExpo(),  m_BTKey);
			}
			
			// Debug
			/*
			if (DEBUGSH) {
			
				vector<DCRTPoly> tmp_test(K_FP+1);
				for (uint32_t ll=0; ll < K_FP; ll++) {
					tmp_test[ll] = ct2->GetAExpo()[ll] - ct1->GetAExpo()[ll];
				}
				tmp_test[K_FP] = ct2->GetBExpo() - ct1->GetBExpo();
				
				std::cout << "Max are " << std::endl;
				int64_t res0 = this->DecryptExp(sk, (*resExp));
				int64_t resminus = this->DecryptExp(sk, tmp_test);
				std::cout << res0 << std::endl;
				std::cout << "ct2 - ct1 is " << resminus << std::endl;
			}*/
		
			vector<DCRTPoly> Exp1_div(K_FP+1);
			vector<DCRTPoly> Exp2_div(K_FP+1);
			// Min
			for (uint32_t i = 0; i < K_FP+1; i++) {
				if (i == K_FP) {
					Exp1_div[i] = (*resExp)[i] - ct1->GetBExpo() ; 
					Exp2_div[i] = (*resExp)[i] - ct2->GetBExpo() ;
				} else {
					Exp1_div[i] = (*resExp)[i] - ct1->GetAExpo()[i]; 
					Exp2_div[i] = (*resExp)[i] - ct2->GetAExpo()[i];
				}
			}


			std::shared_ptr<vector<DCRTPoly>> min_res1;
			std::shared_ptr<vector<DCRTPoly>> min_res2;

			if (types == CTType::DOUBLE) {
				min_res1 =  m_RingGSWscheme->BootstrapReluRev(
					m_params, m_LWEscheme, Exp1_div, 27,  m_BTKey);
				min_res2 =  m_RingGSWscheme->BootstrapReluRev(
					m_params, m_LWEscheme, Exp2_div, 27,  m_BTKey);
			
			} else if(types == CTType::FLOAT){
				min_res1 =  m_RingGSWscheme->BootstrapReluRevF(
					m_params, m_LWEscheme, Exp1_div, 13,  m_BTKey);		
				min_res2 =  m_RingGSWscheme->BootstrapReluRevF(
					m_params, m_LWEscheme, Exp2_div, 13,  m_BTKey);
			}

			// Debug 
			/*
			if(DEBUGSH) {
				std::cout << "Min operation is end "<< std::endl;
				std::cout << "Before DecryptExp" << std::endl;
				int64_t res1 = this->DecryptExp(sk, (*min_res1));
				int64_t res2 = this->DecryptExp(sk, (*min_res2));
				std::cout << "min res1 is " << res1 << std::endl;
				std::cout << "min res2 is " << res2 << std::endl;
				int64_t resExp1aa = this->DecryptExp(sk, (Exp1_div));
				int64_t resExp2aa = this->DecryptExp(sk, (Exp2_div));
				std::cout << "Exp res1 is " << resExp1aa << std::endl;
				std::cout << "Exp res2 is " << resExp2aa << std::endl;
				
			}*/

			// Make Div, true -> Muls, false -> Divs
			bool  Mul = false;
			std::shared_ptr<vector<DCRTPoly>> div_res1; 
			std::shared_ptr<vector<DCRTPoly>> div_res2;
			if (types == CTType::DOUBLE) {
				div_res1 =  m_RingGSWscheme->BootstrapMulPolyMake0(
					m_params, m_LWEscheme, ct1->GetASign(), ct1->GetBSign(), (*min_res1), m_BTKey, 27,  Mul);
			
				div_res2 =  m_RingGSWscheme->BootstrapMulPolyMake0(
					m_params, m_LWEscheme, ct2->GetASign(), ct2->GetBSign(), (*min_res2), m_BTKey,27,  Mul);
			} else if (types == CTType::FLOAT) {
				div_res1 =  m_RingGSWscheme->BootstrapMulPolyMake0(
					m_params, m_LWEscheme, ct1->GetASign(), ct1->GetBSign(), (*min_res1), m_BTKey, 13,  Mul);
			
				div_res2 =  m_RingGSWscheme->BootstrapMulPolyMake0(
					m_params, m_LWEscheme, ct2->GetASign(), ct2->GetBSign(), (*min_res2), m_BTKey,13,  Mul);
	
			}	
			/*
			if (DEBUGSH) {
				std::cout << "MulPoly is gen" << std::endl;
				
				int64_t res3 = this ->DecryptMulPoly(sk, (*div_res1), 27, Mul);
				int64_t res4 = this ->DecryptMulPoly(sk, (*div_res2), 27, Mul);
				std::cout << "div res3 is " << res3 << std::endl;
				std::cout << "div res4 is " << res4 << std::endl;
				
			}*/
			
			// Multiply

			std::shared_ptr<vector<DCRTPoly>> Muls1;
			std::shared_ptr<vector<DCRTPoly>> Muls2;

			if (m_params->GetPARALLEL()) {
				Muls1 = m_RingGSWscheme->ProductDCRTPARALLEL2(m_params, m_BTKey->Evkey, 
									(*div_res1), (*div_res1)[K_FP],
									ct1->GetA(), ct1->GetB());

				Muls2 = m_RingGSWscheme->ProductDCRTPARALLEL2(m_params, m_BTKey->Evkey, 
									(*div_res2), (*div_res2)[K_FP], 
									ct2->GetA(), ct2->GetB());
			} else {
				Muls1 = m_RingGSWscheme->ProductDCRT(m_params, m_BTKey->Evkey, 
									(*div_res1), (*div_res1)[K_FP],
									ct1->GetA(), ct1->GetB());

				Muls2 = m_RingGSWscheme->ProductDCRT(m_params, m_BTKey->Evkey, 
									(*div_res2), (*div_res2)[K_FP], 
									ct2->GetA(), ct2->GetB());
			}	
			// Adding
			
			if (adding) {
				for (uint32_t i = 0; i < K_FP+1; i++) {
					(*Muls1)[i] +=  (*Muls2)[i];
					(*Muls1)[i].SetFormat(Format::COEFFICIENT);
					m_RingGSWscheme->DivideRoundingSelf(m_params, &(*Muls1)[i], Format::COEFFICIENT);
				}
			} else {
				for (uint32_t i = 0; i < K_FP+1; i++) {
					(*Muls1)[i] -=  (*Muls2)[i];
					(*Muls1)[i].SetFormat(Format::COEFFICIENT);
					m_RingGSWscheme->DivideRoundingSelf(m_params, &(*Muls1)[i], Format::COEFFICIENT);
				}
			}

			/*
			if (DEBUGSH) {
				std::cout << "Adding is over" << std::endl;
			}*/

			//DivideReduction
			std::shared_ptr<vector<vector<DCRTPoly>>> resFrac;
			if (types == CTType::DOUBLE) {
				resFrac = m_RingGSWscheme->BootstrapAddPartialDouble(m_params, m_LWEscheme, (*Muls1), m_BTKey, 3, 32, -1);
			} else {
				resFrac = m_RingGSWscheme->BootstrapAddPartialDouble(m_params, m_LWEscheme, (*Muls1), m_BTKey, 1, 16, -1);
	
			}
			uint32_t Frac_lens = resFrac->size();
			vector<DCRTPoly> A_sign = std::move((*resFrac)[Frac_lens - 1]);  
			DCRTPoly b_sign = std::move(A_sign[A_sign.size()-1]);
			A_sign.erase(A_sign.begin()+(A_sign.size()-1));
			resFrac->erase(resFrac->begin() + (resFrac->size() - 1));
	


			// Sign Checking
			/*
			if (DEBUGSH) {
				std::cout << "sign is over" << std::endl;
				
				int64_t tmpsss= this ->DecryptUInt64(sk, 
					std::make_shared<FPRLWECiphertextImpl> (FPRLWECiphertextImpl(
					A_sign, b_sign)));
				if (tmpsss == 3) {
					tmpsss = -1;
				}
				std::cout <<" Sign is " << (tmpsss) << std::endl;
				
			}*/

			// Exp Sum adding 1 
			uint64_t scaling = m_params->GetRLWEParams()->GetScaling().ConvertToInt();
			scaling <<= m_params->GetExpBit();
			NativeInteger Add1 = m_params->GetModuliQ()[1];
			Add1.ModMulEq(scaling, m_params->GetModuliQ()[0]);
			(*resExp)[K_FP].SetFormat(Format::COEFFICIENT);
			(*resExp)[K_FP].GetElementW(0)[0].ModAddEq(Add1, m_params->GetModuliQ()[0]);
			(*resExp)[K_FP].SetFormat(Format::EVALUATION);
			vector<DCRTPoly> Frac_Results = std::move((*resFrac)[Frac_lens - 2]);  
			resFrac->erase(resFrac->begin() + (resFrac->size() - 1));
			
			/*
			if (DEBUGSH){	
				std::cout << "Frac inner searching..." << std::endl;
				int64_t Frac_inner_result = this ->DecryptMulPoly(sk, (Frac_Results), 32, true);
				std::cout << "After Frac inner result " << Frac_inner_result << std::endl;
				for (uint32_t sh = 0; sh < resFrac->size(); sh++) {
					uint64_t tmpss= this ->DecryptUInt64(sk, 
						std::make_shared<FPRLWECiphertextImpl> (FPRLWECiphertextImpl(
						(*resFrac)[sh], (*resFrac)[sh][K_FP]) ));
					std::cout <<" Zero check, " << sh << "'st zero is " << tmpss << std::endl;			
				}
			}*/

			// real ,...0~ 27 or 0~12, 
			// Homomorphic Counter
			auto resExp2 = m_RingGSWscheme->BootstrapRecon2(m_params, m_LWEscheme, 
					(*resFrac), m_BTKey,0); 
			/*
			if (DEBUGSH) {
				std::cout << "recon2 is over" << std::endl;
				int64_t res15 = this->DecryptExp(sk, (*resExp2)[0]);
				std::cout << "want to up is  " << res15 << std::endl;
				int64_t res25 = this->DecryptExp(sk, (*resExp2)[1]);
				std::cout << "Underflow suspect is  " << res25 << std::endl;
			}*/

			auto OFUF = m_RingGSWscheme->CheckOFUF( m_params, m_LWEscheme,
					(*resExp2)[1], (*resExp), types, m_BTKey);
			ct_OFUF->Adding((*OFUF));
			/*
			if (DEBUGSH) {
				int64_t res30 = this->DecryptExp(sk, (*OFUF)[0]);
				std::cout << "Overflow is  " << res30<< std::endl;
				int64_t res31 = this->DecryptExp(sk, (*OFUF)[1]);
				std::cout << "Underflow is  " << res31<< std::endl;
			}*/



			std::shared_ptr<vector<vector<DCRTPoly>>> minExp3;
			
			if (types == CTType::DOUBLE) {
			// Make Min
				minExp3 =  m_RingGSWscheme->BootstrapReluRev(
				m_params, m_LWEscheme, (*resExp), (*resExp2)[0],  m_BTKey);	
			} else if (types == CTType::FLOAT){ 
				minExp3 =  m_RingGSWscheme->BootstrapReluRevF(
				m_params, m_LWEscheme, (*resExp), (*resExp2)[0],  m_BTKey);	
			}

			/*
			if (DEBUGSH) { 
				std::cout << "Rev relu before input check" << std::endl;
				std::cout << "Res Exp is" << this->DecryptExp(sk, (*resExp)) << std::endl;
				std::cout << "Res Exp3 is" << this->DecryptExp(sk, (*resExp2)[0]) << std::endl;
				
				std::cout << "Rev relu is over" << std::endl;
				int64_t res16 = this->DecryptExp(sk, (*minExp3)[0]);
				std::cout << "min is  " << res16 << std::endl;
			}*/

			for (uint32_t kk = 0; kk < K_FP+1; kk++) {
				(*resExp)[kk] -= (*minExp3)[0][kk];
			}
			/*
			if (DEBUGSH) {
				std::cout << "Prev Recon Mul Process.." << std::endl;
			}*/

			std::shared_ptr<vector<DCRTPoly>> MulPolys;		
			if (types == CTType::DOUBLE) {
				MulPolys =  m_RingGSWscheme->BootstrapMulPolyMake(
					m_params, m_LWEscheme, (*minExp3)[0], m_BTKey, 27, true);
			} else if (types == CTType::FLOAT) {
				MulPolys =  m_RingGSWscheme->BootstrapMulPolyMake(
					m_params, m_LWEscheme, (*minExp3)[0], m_BTKey, 13, true);	
			}
			
			/*
			if (DEBUGSH) {
				std::cout << "Recon Mul Process.." << std::endl;	
				int64_t res44 = this ->DecryptMulPoly(sk, (*MulPolys), 27, true);
				std::cout << "Recon mul is  " << res44 << std::endl;
			}*/

			// Muls
			if (m_params->GetPARALLEL()) {
				Muls1 = m_RingGSWscheme->ProductDCRTPARALLEL2(m_params, m_BTKey->Evkey, 
									Frac_Results, Frac_Results[K_FP],
									(*MulPolys), (*MulPolys)[K_FP]);

			} else {
				Muls1 = m_RingGSWscheme->ProductDCRT(m_params, m_BTKey->Evkey, 
									Frac_Results, Frac_Results[K_FP],
									(*MulPolys), (*MulPolys)[K_FP]);
			}	
			for (uint32_t i = 0; i < K_FP+1; i++) {
				(*Muls1)[i].SetFormat(Format::COEFFICIENT);
				m_RingGSWscheme->DivideRoundingSelf(m_params, &(*Muls1)[i], Format::COEFFICIENT);
			}
		

			std::shared_ptr<vector<DCRTPoly>> A;
			if (types == CTType::DOUBLE) {
				A = m_RingGSWscheme->BootstrapOrignPartial(m_params, m_LWEscheme, m_BTKey, (*Muls1), 5, 31);
			
			} else if (types == CTType::FLOAT) {
				A = m_RingGSWscheme->BootstrapOrignPartial(m_params, m_LWEscheme, m_BTKey, (*Muls1), 3, 15);
			}

			/*
			if (DEBUGSH) {	
				std::cout << "Final Inner Frac search..." << std::endl;
				int64_t Frac_inner_result2 = this ->DecryptMulPoly(sk, (*A), 32, true);
				std::cout << "Final Frac inner result " << Frac_inner_result2 << std::endl;
			}*/


			DCRTPoly b = std::move((*A)[A->size()-1]);
			A->erase(A->begin()+(A->size()-1));
	
			vector<DCRTPoly> A_exp = std::move(*resExp);  
			DCRTPoly b_exp = std::move(A_exp[A_exp.size()-1]);
			A_exp.erase(A_exp.begin()+(A_exp.size()-1));

			return std::make_shared<FPRLWECiphertextImpl> (FPRLWECiphertextImpl(
					std::move(*A), std::move(b), std::move(A_sign), std::move(b_sign), 
					std::move(A_exp), std::move(b_exp), types));

		}

		FPRLWECiphertext FPFHEContext::Add(
			FPConstRLWECiphertext ct1, 
			FPConstRLWECiphertext ct2,
			//FPConstRLWEPrivateKey Rsk,
			std::shared_ptr<FPOFUF> ct_OFUF
			//const FPOFUF &ct_OFUF
		) const {

		CTType ct1_t = ct1->GetType();
		CTType ct2_t = ct2->GetType();	
		if (ct1_t != ct2_t) {
			PALISADE_THROW(config_error, "Current Version does not support cross type operation");
		} 
		// 
		if (ct1_t == CTType::DOUBLE || ct1_t == CTType::FLOAT) {
			return this -> AddCore(ct1, ct2, ct_OFUF, true);
		}

		std::shared_ptr<vector<DCRTPoly>> Muls1;
		std::shared_ptr<vector<DCRTPoly>> Muls2;

		if (m_params->GetPARALLEL()) {
			Muls1 = m_RingGSWscheme->ProductDCRTPARALLEL2(m_params, m_BTKey->Evkey, 
									ct1->GetASign(), ct1->GetBSign(),
									ct1->GetA(), ct1->GetB());

			Muls2 = m_RingGSWscheme->ProductDCRTPARALLEL2(m_params, m_BTKey->Evkey, 
									ct2->GetASign(), ct2->GetBSign(),
									ct2->GetA(), ct2->GetB());
	

		} else {
			Muls1 = m_RingGSWscheme->ProductDCRT(m_params, m_BTKey->Evkey, 
									ct1->GetASign(), ct1->GetBSign(),
									ct1->GetA(), ct1->GetB());

			Muls2 = m_RingGSWscheme->ProductDCRT(m_params, m_BTKey->Evkey, 
									ct2->GetASign(), ct2->GetBSign(),
									ct2->GetA(), ct2->GetB());
		}
		auto start_add = std::chrono::system_clock::now();
	
		for (uint32_t i = 0; i < Muls2->size(); i++) {
			(*Muls1)[i] +=  (*Muls2)[i];
			(*Muls1)[i].SetFormat(Format::COEFFICIENT);
			m_RingGSWscheme->DivideRoundingSelf(m_params, &(*Muls1)[i], Format::COEFFICIENT);
		}
		
		std::chrono::duration<double> end_add = std::chrono::system_clock::now() - start_add;
		std::cout << "Just Adding time is " << end_add.count() << std::endl;

		//DivideReduction
		//auto res = m_RingGSWscheme->BootstrapAdd(m_params, m_LWEscheme, (*Muls1), m_BTKey); 
		auto res = m_RingGSWscheme->BootstrapAddPartial(m_params, m_LWEscheme, (*Muls1), m_BTKey, 0, 31, 0); 
		
		vector<DCRTPoly> A = std::move((*res)[0]);  
		DCRTPoly b = std::move(A[A.size()-1]);
		A.erase(A.begin()+(A.size()-1));
	
		vector<DCRTPoly> A_sign = std::move((*res)[1]);  
		DCRTPoly b_sign = std::move(A_sign[A_sign.size()-1]);
		A_sign.erase(A_sign.begin()+(A_sign.size()-1));
	
		return std::make_shared<FPRLWECiphertextImpl> (FPRLWECiphertextImpl(std::move(A), b, std::move(A_sign), b_sign));
	}

 
	FPRLWECiphertext FPFHEContext::Sub(
			FPConstRLWECiphertext ct1, 
			FPConstRLWECiphertext ct2,
			//FPConstRLWEPrivateKey Rsk,
			std::shared_ptr<FPOFUF> ct_OFUF
			//const FPOFUF &ct_OFUF
			) const {

		CTType ct1_t = ct1->GetType();
		CTType ct2_t = ct2->GetType();	
		if (ct1_t != ct2_t) {
			PALISADE_THROW(config_error, "Current Version does not support cross type operation");
		} 

		if (ct1_t == CTType::DOUBLE || ct1_t == CTType::FLOAT) {
			return this -> AddCore(ct1, ct2 , ct_OFUF, false);
		}



		auto Muls1 = m_RingGSWscheme->ProductDCRT(m_params, m_BTKey->Evkey, 
									ct1->GetASign(), ct1->GetBSign(),
									ct1->GetA(), ct1->GetB());

		auto Muls2 = m_RingGSWscheme->ProductDCRT(m_params, m_BTKey->Evkey, 
									ct2->GetASign(), ct2->GetBSign(),
									ct2->GetA(), ct2->GetB());
		
		
		for (uint32_t i = 0; i < Muls2->size(); i++) {
			(*Muls1)[i] -=  (*Muls2)[i];
			(*Muls1)[i].SetFormat(Format::COEFFICIENT);
			m_RingGSWscheme->DivideRoundingSelf(m_params, &(*Muls1)[i], Format::COEFFICIENT);
		}
		//DivideReduction
		//auto res = m_RingGSWscheme->BootstrapAdd(m_params, m_LWEscheme, (*Muls1), m_BTKey); 
		auto res = m_RingGSWscheme->BootstrapAddPartial(m_params, m_LWEscheme, (*Muls1), m_BTKey, 0, 31,0); 
		
		vector<DCRTPoly> A = std::move((*res)[0]);  
		DCRTPoly b = std::move(A[A.size()-1]);
		A.erase(A.begin()+(A.size()-1));
	
		vector<DCRTPoly> A_sign = std::move((*res)[1]);  
		DCRTPoly b_sign = std::move(A_sign[A_sign.size()-1]);
		A_sign.erase(A_sign.begin()+(A_sign.size()-1));
	
		return std::make_shared<FPRLWECiphertextImpl> (FPRLWECiphertextImpl(std::move(A), b, std::move(A_sign), b_sign));
	}
 


	std::shared_ptr<FPOFUF> FPFHEContext::MakeOverFlowDetectCT() const {
		return std::make_shared<FPOFUF> (FPOFUF());
	}


	void FPFHEContext::CheckOFUF(
			FPConstRLWEPrivateKey Rsk,
			shared_ptr<FPOFUF> ct_OFUF) const {
		int64_t overflow = this->DecryptExp(Rsk, ct_OFUF->GetElements()[0]);
		std::cout << "number of overflow  is    " << overflow << std::endl;

		int64_t underflow = this->DecryptExp(Rsk, ct_OFUF->GetElements()[1]);
		std::cout << "number of underflow is    " << underflow << std::endl;
		return ;
	}
	
	

// TEST Function

/*
	 void FPFHEContext::CheckBS() const {
	 uint64_t M = m_params->GetRLWEParams()->GetN();
	 vector<DCRTPoly> SK = m_BTKey->BK_SK->GetElement();
	 uint64_t Q = (m_params->GetModuliQ()[0].ConvertToInt());
	 uint64_t Q_over_2 = (Q >> 1);
	 for (uint64_t i=0; i < M; i++) {
	 auto tmp_b = (m_BTKey->BSkey);
//auto tmp_c = tmp_b->GetElements();
//std::shared_ptr<FPRingGSWCiphertext> tmp_d = (tmp_b)[0][a0][i];
std::vector<std::vector<std::vector<FPRingGSWCiphertext>>> tmp_d
= tmp_b->GetElements();
uint64_t sizes = tmp_d[0].size();

for (uint64_t j=0; j<sizes; j++) {
DCRTPoly result = (tmp_d[0][0][i].GetElements()[j][1]-(tmp_d[0][0][i].GetElements()[j][0]*SK));

result.SetFormat(Format::COEFFICIENT);
NativePoly result_poly = result.GetElementAtIndex(0);

for (uint64_t k=0; k < m_params->GetGSWN(); k++) {
uint64_t test_val = result_poly[k].ConvertToInt();
if (test_val > Q_over_2) {
test_val = Q - test_val;
} 
if (test_val >= 20) {
//std::cout << "Bootstrapping error!!!!" <<std::endl;
}
}

std::cout <<"i is " << i << ", j is " << j << ", a encrypt res is " << std::endl;
std::cout << result_poly << std::endl;
}

for (uint64_t j=0; j<sizes; j++) {
DCRTPoly result = (tmp_d[0][1][i].GetElements()[j][1]-(tmp_d[0][1][i].GetElements()[j][0]*SK));

result.SetFormat(Format::COEFFICIENT);
NativePoly result_poly = result.GetElementAtIndex(0);

for (uint64_t k=0; k < m_params->GetGSWN(); k++) {
uint64_t test_val = result_poly[k].ConvertToInt();
if (test_val > Q_over_2) {
test_val = Q - test_val;
} 
if (test_val >= 10) {
//std::cout << "Bootstrapping error!!!!" <<std::endl;
}
}
std::cout <<"i is " << i << ", j is " << j << ", b encrypt res is " << std::endl;

std::cout << result_poly << std::endl;

}

}
}
*/

/*
	 void FPFHEContext::CheckKS(FPRLWEPrivateKey sk) const {

	 DCRTPoly SK = sk->GetElement();
	 uint64_t Q = (m_params->GetModuliQ()[0].ConvertToInt());
	 uint64_t Q_over_2 = (Q >> 1);
	 uint64_t M = m_params->GetRLWEParams()->GetN();


	 vector<vector<FPRLWECiphertextImpl>> KSKeyCT = m_BTKey->RKSkey->GetElements();
	 uint32_t size1 = KSKeyCT.size();
	 uint32_t size2 = KSKeyCT[0].size();

	 for (uint32_t i = 0; i < size1; i++) {
	 for (uint32_t j=0; j < size2; j++) {
	 DCRTPoly result =KSKeyCT[i][j].GetB() - (KSKeyCT[i][j].GetA() * SK);

	 result.SetFormat(Format::COEFFICIENT);
	 NativePoly result_poly = result.GetElementAtIndex(0);
	 for (uint32_t k = 0; k < M; k++ ){
	 uint64_t test_val = result_poly[k].ConvertToInt();
	 if (test_val > Q_over_2) {
	 test_val = Q - test_val;
	 } 
	 if (test_val >= 10) {
	 std::cout << "Key Switching error!!!!" <<std::endl;
	 }

	 }
	 }
	 }

	}
	*/

}  // namespace lbcrypto
