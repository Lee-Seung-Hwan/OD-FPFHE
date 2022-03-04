/*************
 * This file is modified from fhew.cpp which is writen by Leo Ducas and Daniele Micciancio.
 * Full paper is listed in : eprint.iacr.org/2022/186.
 * This source is obeyed and applied following copyright law. 
 *
 * If you have any question, please contact us by the email kr3951@hanyang.ac.kr
 * Seunghwan Lee, 2022.03.04
 * *************/


// @file fhew.cpp - FHEW scheme (RingGSW accumulator) implementation
// The scheme is described in https://eprint.iacr.org/2014/816 and in
// Daniele Micciancio and Yuriy Polyakov, "Bootstrapping in FHEW-like
// Cryptosystems", Cryptology ePrint Archive, Report 2020/086,
// https://eprint.iacr.org/2020/086.
//
// Full reference to https://eprint.iacr.org/2014/816:
// @misc{cryptoeprint:2014:816,
//   author = {Leo Ducas and Daniele Micciancio},
//   title = {FHEW: Bootstrapping Homomorphic Encryption in less than a second},
//   howpublished = {Cryptology ePrint Archive, Report 2014/816},
//   year = {2014},
//   note = {\url{https://eprint.iacr.org/2014/816}},
// @author TPOC: contact@palisade-crypto.org
//
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


/******************************/
// This code is
// We modify code fhew.cpp from Daniel Micciancio and Yuriy Polyakov




/******************************/




#include "fpfhew.h"
//#include <ctime>



bool DEBUG = false;


namespace lbcrypto {

	// Encryption as described in Section 5 of https://eprint.iacr.org/2014/816
	// skNTT corresponds to the secret key z
	//
	//  Should be removed.... 
	
	/*
	std::shared_ptr<FPRingGSWCiphertext> FPRingGSWAccumulatorScheme::EncryptAP(
		const std::shared_ptr<FPRingGSWCryptoParams> params, 
		const NativePoly &skNTT,
    const FPLWEPlaintext &m) 
	const {
		NativeInteger Q = params->GetLWEParams()->GetQ();
		int64_t q = params->GetLWEParams()->Getq().ConvertToInt();
		uint32_t N = params->GetLWEParams()->GetN();
		uint32_t digitsG = params->GetDigitsG();
		uint32_t digitsG2 = params->GetDigitsG2();
		const shared_ptr<ILNativeParams> polyParams = params->GetPolyParams();
		auto result = std::make_shared<FPRingGSWCiphertext>(digitsG2, 2);
		DiscreteUniformGeneratorImpl<NativeVector> dug;
		dug.SetModulus(Q);


		// Reduce mod q (dealing with negative number as well)
		int64_t mm = (((m % q) + q) % q) * (2 * N / q);
		int64_t sign = 1;
		if (mm >= N) {
			mm -= N;
			sign = -1;
		}
		// tempA is introduced to minimize the number of NTTs
		std::vector<NativePoly> tempA(digitsG2);
		for (uint32_t i = 0; i < digitsG2; ++i) {
			// populate result[i][0] with uniform random a
			(*result)[i][0] = NativePoly(dug, polyParams, Format::COEFFICIENT);
			tempA[i] = (*result)[i][0];
			// populate result[i][1] with error e
			(*result)[i][1] = NativePoly(params->GetLWEParams()->GetDgg(), polyParams,
                                 Format::COEFFICIENT);
		}

		for (uint32_t i = 0; i < digitsG; ++i) {
			if (sign > 0) {
				// Add G Multiple
				(*result)[2 * i][0][mm].ModAddEq(params->GetGPower()[i], Q);
				// [a,as+e] + X^m*G
				(*result)[2 * i + 1][1][mm].ModAddEq(params->GetGPower()[i], Q);
			} else {
				// Subtract G Multiple
				(*result)[2 * i][0][mm].ModSubEq(params->GetGPower()[i], Q);
				// [a,as+e] - X^m*G
				(*result)[2 * i + 1][1][mm].ModSubEq(params->GetGPower()[i], Q);
			}
		}

		// 3*digitsG2 NTTs are called
		result->SetFormat(Format::EVALUATION);
		for (uint32_t i = 0; i < digitsG2; ++i) {
			tempA[i].SetFormat(Format::EVALUATION);
			(*result)[i][1] += tempA[i] * skNTT;
		}
		
		return result;
	}
	*/

	// Encryption for the GINX variant, as described in "Bootstrapping in FHEW-like
	// Cryptosystems"
	//
	//
	

	std::shared_ptr<FPRingGSWCiphertext> FPRingGSWAccumulatorScheme::EncryptGINX(
    const std::shared_ptr<FPRingGSWCryptoParams> params, 
		const vector<DCRTPoly> &skNTT,
    const FPLWEPlaintext &m) const {
  
		const shared_ptr<ILDCRTParams<BigInteger>> DCRTpolyParams = params->GetGSWPolyParams();
	

		// dim0 : Gadget decom d length of Gadget / adim2: ring nums
   
		uint32_t Gadget_len = params->GetLenGDsave();
		uint32_t GD_Bit			= params->GetBaseGDBit();
		uint32_t GD         = 1;
						 GD <<= GD_Bit;
		uint32_t GD_rm = params->GetLenGDrm();
		uint32_t GSW_K = params->GetGSWK();
		
		auto result = std::make_shared<FPRingGSWCiphertext>(GSW_K+1, GSW_K+1, Gadget_len);
		DCRTPoly::DugType dug;

	
		// tempA is introduced to minimize the number of NTTs
		vector<DCRTPoly> tempA(GSW_K*Gadget_len);

		//Make A And B NTT
		for (uint32_t k_row=0; k_row < GSW_K+1; k_row++) {
			for (uint32_t i = 0; i < (Gadget_len); ++i) {
				for (uint32_t k_col= 0; k_col < GSW_K; k_col++) {
					// Make A
					(*result)[k_row][k_col][i] = DCRTPoly(dug, DCRTpolyParams, Format::COEFFICIENT);
					if (k_col == k_row) {
						tempA[k_col * Gadget_len + i] = (*result)[k_row][k_col][i];
					}
				}
				// Make B
				(*result)[k_row][GSW_K][i] = DCRTPoly(params->GetGSWDgg(), DCRTpolyParams, Format::COEFFICIENT);
			}
		}
		//std::cout << "Uniform is gen" << std::endl; 
		if (m > 0) {
			NativePoly tmp_poly = NativePoly(DCRTpolyParams->GetParams()[0], Format::COEFFICIENT, true);		
			//rm Go
			tmp_poly[0] = 1;
			tmp_poly[0] <<= (GD_Bit * GD_rm);
			DCRTPoly tmp_dcrt_poly = DCRTPoly(DCRTpolyParams,Format::COEFFICIENT);
			tmp_dcrt_poly = tmp_poly;
 
			for (uint32_t i = 0; i < Gadget_len; ++i) {
				for(uint32_t k_col = 0; k_col < GSW_K+1; k_col++) {
					// =========   Orign  ==============
					// Add G Multiple in A part
					(*result)[k_col][k_col][i] += tmp_dcrt_poly;
					// [a,as+e] + G in B part
					//(*result)[2 * i + 1][1] += tmp_dcrt_poly;
				}
				tmp_dcrt_poly *= GD;
			}
		}
		
		// 3*digitsG2 NTTs are called
		result->SetFormat(Format::EVALUATION);
		for (uint32_t i = 0; i < Gadget_len ; ++i) {
			for(uint32_t k_row = 0; k_row < GSW_K+1; k_row++) {
				for(uint32_t k_col = 0; k_col < GSW_K; k_col++) {
					if (k_row == k_col) {
						tempA[k_col * Gadget_len + i].SetFormat(Format::EVALUATION);
						(*result)[k_row][GSW_K][i] += tempA[k_row * Gadget_len + i] * skNTT[k_col];
					} else {
						(*result)[k_row][GSW_K][i] += (*result)[k_row][k_col][i] * skNTT[k_col];
					}
				}
			// Add A.S in B
			//(*result)[i][1] += (tempA[i] *= skNTT);
			}
		}
		return result;
	}


	// wrapper for KeyGen methods


FPRingGSWEvalKey FPRingGSWAccumulatorScheme::KeyGen(
    const std::shared_ptr<FPRingGSWCryptoParams> params,
		const std::shared_ptr<FPLWEEncryptionScheme> lwescheme,
		const std::shared_ptr<FPRLWEEncryptionScheme> rlwescheme,
	const std::shared_ptr<const FPRLWEPrivateKeyImpl> RLWEsk
		) const {
   return KeyGenGINX(params, lwescheme,  rlwescheme,RLWEsk);
}


FPRingGSWEvalKey FPRingGSWAccumulatorScheme::KeyGenGINX(
    const std::shared_ptr<FPRingGSWCryptoParams> params,
	  const std::shared_ptr<FPLWEEncryptionScheme> lwescheme, 
    const std::shared_ptr<FPRLWEEncryptionScheme> rlwescheme, 
		const std::shared_ptr<const FPRLWEPrivateKeyImpl> RLWEsk) const {

	FPRingGSWEvalKey ek; 	

	// EV Key Gen
	ek.Evkey = rlwescheme->EvalKeyGen(params -> GetRLWEParams(), 
			params->GetBaseEVBit(),
			params->GetLenEVsave(),
			params->GetLenEVrm(),
			params->GetEVDgg(), RLWEsk);
	
	std::cout <<"Evaluation Key is Gen" << std::endl;	
	
	// PKS KeyGen 0dim : N, 1dim : Decompose
	
	//LWE KeyGen
	uint32_t h_PKS = params->GetH_PKS();
	std::shared_ptr<FPLWEPrivateKeyImpl> LWE_sk = lwescheme->KeyGen(params->GetLWEParams(), h_PKS);
	
	uint32_t Ns = params->GetRLWEParams()->GetN();
	uint32_t Ks_FP = params->GetRLWEParams()->GetK();
	//std::cout << "Ks_FP is " <<Ks_FP << std::endl;
	uint32_t PKS_len = params->GetLWEParams()->GetPKSLen();
	uint32_t PKS_rm = params->GetLWEParams()->GetPKSLenRm();
	
	vector<NativePoly> RLWE_sk_poly_0(Ks_FP);

	for (uint32_t i = 0; i < Ks_FP; i++) {
		RLWE_sk_poly_0[i] = RLWEsk->GetElement()[i].GetElementAtIndex(0);
		RLWE_sk_poly_0[i].SetFormat(Format::COEFFICIENT);
	}

  int64_t q = (params->GetRLWEParams()->GetModuliQ()[0].ConvertToInt());
	int64_t qHalf = q >> 1;
	uint64_t q_lwe = LWE_sk->GetElement().GetModulus().ConvertToInt();
	uint32_t PKS_base = 1;
					 PKS_base <<= params->GetLWEParams()->GetBasePKSBit();
	ek.PKSkey = std::make_shared<FPLWESwitchingKey>(Ks_FP, Ns, PKS_len- PKS_rm);

	#pragma omp parallel for
	//{
        for (uint32_t ks = 0; ks < Ks_FP; ks++) {
		for (uint32_t i = 0; i < Ns; i++) {
			int64_t s = RLWE_sk_poly_0[ks][i].ConvertToInt();	
			if (s > qHalf) s -= q;
			if (s != -1 && s != 0 && s != 1) {
				PALISADE_THROW(config_error, "LWE sk should be -1, 0, 1");
			}
			for (uint32_t j=0; j<PKS_len; j++) {
				
				if (j < PKS_rm) {
					s *= PKS_base;
					continue;
				} else {
					if (s >= 0) {
						(*ek.PKSkey)[ks][i][j - PKS_rm] = *(lwescheme -> Encrypt_WO_scaling(params -> GetLWEParams(), LWE_sk, s));
					} else {
					// Test
					int64_t test_val = ((int64_t)q_lwe) + s;
					if (test_val <0) {
						PALISADE_THROW(config_error, "PKS msg is error");
					}

						(*ek.PKSkey)[ks][i][j - PKS_rm] = *(lwescheme -> Encrypt_WO_scaling(params -> GetLWEParams(), LWE_sk, q_lwe + s));
					}

					s *= PKS_base;
				}
			}
		}
	}//}

	std::cout <<"Pry KeySwitch is Gen" << std::endl;	
	

	// BK KEy Gen
	uint32_t n_PKS = LWE_sk->GetElement().GetLength();	
	uint32_t K_GSW = params->GetGSWK();
	//std::cout << "GSW_K is "<< K_GSW << std::endl;
	ek.BSkey = std::make_shared<FPRingGSWBTKey>(1, 2, n_PKS);
	// BK secret Generation
	DCRTPoly::TugType tug;
	
	vector<DCRTPoly> sk_GSW_vec(K_GSW);
	for (uint32_t i = 0; i< K_GSW; i++) {
		sk_GSW_vec[i] =  DCRTPoly(tug, params->GetGSWPolyParams(), Format::EVALUATION);
	}

	// handles ternary secrets using signed mod 3 arithmetic; 0 -> {0,0}, 1 ->
	// {1,0}, -1 -> {0,1}
	NativeVector LWE_sk_vector = LWE_sk->GetElement();
	int64_t q_lwe_Half = q_lwe >> 1;

	#pragma omp parallel for
	//{
       for (uint32_t i = 0; i < n_PKS; ++i) {
		int64_t s = LWE_sk_vector[i].ConvertToInt();
		
		if (s > q_lwe_Half) s -= q_lwe;
		switch (s) {
			case 0:
				(*ek.BSkey)[0][0][i] = *(EncryptGINX(params, sk_GSW_vec, 0));
				(*ek.BSkey)[0][1][i] = *(EncryptGINX(params, sk_GSW_vec, 0));
				break;
			case 1:
				(*ek.BSkey)[0][0][i] = *(EncryptGINX(params, sk_GSW_vec, 1));
				(*ek.BSkey)[0][1][i] = *(EncryptGINX(params, sk_GSW_vec, 0));
				break;
			case -1:
				(*ek.BSkey)[0][0][i] = *(EncryptGINX(params, sk_GSW_vec, 0));
				(*ek.BSkey)[0][1][i] = *(EncryptGINX(params, sk_GSW_vec, 1));
				break;
			default:
				std::string errMsg =
         "ERROR: only ternary secret key distributions are supported.";
				PALISADE_THROW(not_implemented_error, errMsg);
		}
	} //}
	std::cout <<"Boostrapping Key is Gen" << std::endl;	
	
	// KS Gen
	uint32_t KS_bit = params->GetBaseKSBit();
	uint32_t KS_len = params -> GetLenKSsave();
	uint32_t KS_rm = params -> GetLenKSrm();
	uint32_t KS         = 1;
						KS <<= KS_bit;
	uint32_t KS_first_mul = 1;
		KS_first_mul <<= (KS_bit * KS_rm);
	uint32_t GSW_N = params->GetGSWN();
	uint32_t GSW_K = params->GetGSWK();
	
	// NTT
	for (uint32_t ks = 0; ks < GSW_K; ks++) {
		sk_GSW_vec[ks].SetFormat(Format::COEFFICIENT); 
	}

	vector<DCRTPoly> sk_poly_vec = RLWEsk->GetElement();
	// sk shouuld be Evaluation

	//DCRT
	// Dim0 is Ks, 1 Ns, 2 len
	std::shared_ptr<vector<vector<vector<FPRLWECiphertextImpl>>>> KSKeyCT
		= std::make_shared<vector<vector<vector<FPRLWECiphertextImpl>>>> (GSW_K);
	//KSKeyCT->resize(GSW_K);
	for (uint32_t i = 0; i < GSW_K; i++){
		(*KSKeyCT)[i].resize(GSW_N);
		for (uint32_t j = 0; j < GSW_N; j ++) {
			(*KSKeyCT)[i][j].resize(KS_len);
		}
	}
	
	//uint32_t 

	#pragma omp parallel for collapse(1)
	//{
       for (uint32_t i = 0; i < GSW_N; i++) {
		for (uint32_t ks= 0; ks < GSW_K; ks++) {
					
			NativePoly tmp_msg = NativePoly(params-> GetRLWEParams()->GetDCRTParams()->GetParams()[0],Format::COEFFICIENT,true);
			tmp_msg[0]= sk_GSW_vec[ks].GetElementAtIndex(0)[i];
		
			DCRTPoly MSG = DCRTPoly(params->GetRLWEParams() -> GetDCRTParams(), Format::COEFFICIENT, true);
			MSG = tmp_msg;
			DCRTPoly::DugType dug; 
			//DCRTPoly::DggType dgg;
			//dgg.SetStd(0.398942280);
			//vector<DCRTPoly> MSG_Poly = MSG.PowersOfBase(KS_bit);
			MSG *= KS_first_mul;
			for (uint32_t j=0; j<KS_len; j++) {
		
				vector<DCRTPoly> as(Ks_FP); 
				DCRTPoly b =  DCRTPoly(params -> GetRLWEParams()-> GetDgg(), params -> GetRLWEParams() -> GetDCRTParams(), Format::COEFFICIENT);
				//DCRTPoly b =  DCRTPoly(dgg, params -> GetRLWEParams() -> GetDCRTParams(), Format::COEFFICIENT);
				//DCRTPoly b =  DCRTPoly(params -> GetRLWEParams() -> GetDCRTParams(), Format::COEFFICIENT, true);
				
				//b+= MSG_Poly[j+KS_rm];
				b += MSG;
				b.SetFormat(Format::EVALUATION);
				for (uint32_t k_fp = 0; k_fp < Ks_FP; k_fp++) {
					as[k_fp] = DCRTPoly(dug, params ->GetRLWEParams() -> GetDCRTParams(), Format::EVALUATION);
					//as[k_fp].SetFormat(Format::EVALUATION);
					b += as[k_fp] * sk_poly_vec[k_fp];
				}
			
				(*KSKeyCT)[ks][i][j] = FPRLWECiphertextImpl(as, b);
				MSG *= KS;
			}
		}
	} //}

	//sk_GSW.SetFormat(Format::EVALUATION);
	//ek.RKSkey = std::make_shared<FPRLWESwitchingKey> (FPRLWESwitchingKey(KSKeyCT));
	ek.RKSkey = std::make_shared<FPRLWESwitchingKey> (FPRLWESwitchingKey((*KSKeyCT)));
	std::cout << "Ring KeySwitching Key Gen" << std::endl;
	
	//ek.BK_SK =  std::make_shared<FPRLWEPrivateKeyImpl> (sk_GSW);
  return ek;
}


	/*
	std::shared_ptr<FPRLWECiphertextImpl> FPRingGSWAccumulatorScheme::EncryptKS(
    const std::shared_ptr<FPRingGSWCryptoParams> params, 
		const vector<DCRTPoly> &skNTT,
    const uint64_t &m) const {
  
		const shared_ptr<ILDCRTParams<BigInteger>> DCRTpolyParams = params->GetGSWPolyParams();
	

		// dim0 : Gadget decom d length of Gadget / adim2: ring nums
   
		uint32_t Gadget_len = params->GetLenGDsave();
		uint32_t GD_Bit			= params->GetBaseGDBit();
		uint32_t GD         = 1;
						 GD <<= GD_Bit;
		uint32_t GD_rm = params->GetLenGDrm();
		uint32_t GSW_K = params->GetGSWK();
		
		auto result = std::make_shared<FPRingGSWCiphertext>(Gadget_len,GSW_K+1, GSW_K+1);
		DCRTPoly::DugType dug;

	
		// tempA is introduced to minimize the number of NTTs
		vector<DCRTPoly> tempA(GSW_K*Gadget_len);

		//Make A And B NTT
		for (uint32_t k_row=0; k_row < GSW_K+1; k_row++) {
			for (uint32_t i = 0; i < (Gadget_len); ++i) {
				for (uint32_t k_col= 0; k_col < GSW_K; k_col++) {
					// Make A
					(*result)[i][k_row][k_col] = DCRTPoly(dug, DCRTpolyParams, Format::COEFFICIENT);
					if (k_col == k_row) {
						tempA[k_col * Gadget_len + i] = (*result)[i][k_row][k_col];
					}
				}
				// Make B
				(*result)[i][k_row][GSW_K] = DCRTPoly(params->GetGSWDgg(), DCRTpolyParams, Format::COEFFICIENT);
			}
		}
		//std::cout << "Uniform is gen" << std::endl; 
		if (m > 0) {
			NativePoly tmp_poly = NativePoly(DCRTpolyParams->GetParams()[0], Format::COEFFICIENT, true);		
			//rm Go
			tmp_poly[0] = 1;
			tmp_poly[0] <<= (GD_Bit * GD_rm);
			DCRTPoly tmp_dcrt_poly = DCRTPoly(DCRTpolyParams,Format::COEFFICIENT);
			tmp_dcrt_poly = tmp_poly;
 
			for (uint32_t i = 0; i < Gadget_len; ++i) {
				for(uint32_t k_col = 0; k_col < GSW_K+1; k_col++) {
					// =========   Orign  ==============
					// Add G Multiple in A part
					(*result)[i][k_col][k_col] += tmp_dcrt_poly;
					// [a,as+e] + G in B part
					//(*result)[2 * i + 1][1] += tmp_dcrt_poly;
				}
				tmp_dcrt_poly *= GD;
			}
		}
		
		// 3*digitsG2 NTTs are called
		result->SetFormat(Format::EVALUATION);
		for (uint32_t i = 0; i < Gadget_len ; ++i) {
			for(uint32_t k_row = 0; k_row < GSW_K+1; k_row++) {
				for(uint32_t k_col = 0; k_col < GSW_K; k_col++) {
					if (k_row == k_col) {
						tempA[k_col * Gadget_len + i].SetFormat(Format::EVALUATION);
						(*result)[i][k_row][GSW_K] += tempA[k_row * Gadget_len + i] * skNTT[k_col];
					} else {
						(*result)[i][k_row][GSW_K] += (*result)[i][k_row][k_col] * skNTT[k_col];
					}
				}
			// Add A.S in B
			//(*result)[i][1] += (tempA[i] *= skNTT);
			}
		}
		return result;
	}
*/






// GINX Accumulation as described in "Bootstrapping in FHEW-like Cryptosystems"
void FPRingGSWAccumulatorScheme::AddToACCGINX(
    const std::shared_ptr<FPRingGSWCryptoParams> params,
    const FPRingGSWCiphertext &input, 
		const uint64_t &a,
    std::shared_ptr<FPRingGSWCiphertext> acc
		) const {
	
	uint32_t K_GSW = params->GetGSWK();
	uint32_t GD_baseBit = params->GetBaseGDBit();
	uint32_t GD_Len = params -> GetLenGDsave();
	uint32_t GD_rm = params -> GetLenGDrm();
	
	// OK

	const shared_ptr<ILDCRTParams<BigInteger>> polyGSWParams = params->GetGSWPolyParams();
	// the reason why we use Index 0 is due to structure of GSW Sample
	
  vector<DCRTPoly> ct = acc->GetElementAtIndex(0,0);	
	vector<vector<DCRTPoly>> DCRT_CT(K_GSW+1);
	// calls K+1  NTTs
	{for (uint32_t i = 0; i < K_GSW+1; i++) {
		ct[i].SetFormat(Format::COEFFICIENT);
		DCRT_CT[i] = std::move(*(this->SGDLevel2(params, ct[i], GD_baseBit, GD_rm, Format::EVALUATION)));
	}}
	/*
	{for (uint32_t i = 0; i < K_GSW+1; i++) {
		for (uint32_t j = 0; j < DCRT_CT[i].size(); j++) {
			DCRT_CT[i][j].SetFormat(Format::EVALUATION);
		}																				
	}}
*/
	//#pragma omp barrier



	// True means that it returns it with EVAL mode.
	//for (uint32_t i = 0; i < K_GSW+1; i++) {
	//DCRT_CT[i] = ct[i].BaseDecompose(GD_baseBit, true);
	//}
	
	//vector<DCRTPoly> dct_0 = ct[0].BaseDecompose(GD_baseBit, true);  // A
  //vector<DCRTPoly> dct_1 = ct[1].BaseDecompose(GD_baseBit, true);  // B
	
	const DCRTPoly &monomial = params->GetMonomial(a);
  // acc = dct * input (matrix product);
  // uses in-place * operators for the last call to dct[i] to gain performance
  // improvement
	// Input is GSW Sample
	// k = 0 -> A, k = 1 -> B
	
	// Output k iterating
	//#pragma omp parallel for
	{for (uint32_t out_k = 0; out_k < K_GSW+1; out_k++) {
		// Zero Index
		DCRTPoly temp1 = DCRTPoly(polyGSWParams, Format::EVALUATION, true );

		//DCRTPoly temp1 = (iter_k < 1) ? (DCRT_CT[iter_k][0] * input[0][][k]) : (DCRT_CT[k][0] *= input[0][k]);	// +A
		
		// Left col & right row index 
		for (uint32_t iter_k = 0; iter_k < K_GSW+1; iter_k++) {
			//temp1 += (k < 1) ? (dct_1[GD_rm] * input[1][k]) : (dct_1[GD_rm] *= input[1][k]);	// +B	
			
			// Gadget Decompositions
			for (uint32_t l = 0; l < GD_Len; l++) {
				
				//temp1 += (DCRT_CT[iter_k][l+GD_rm] * input[l][iter_k][out_k]); 
				temp1 += (DCRT_CT[iter_k][l] * input[l][iter_k][out_k]); 
				/*
				if (k == 0) {
					temp1 += dct_0[GD_rm + l] * input[l*2][k];				// CT_A	Product GSW A
					temp1 += dct_1[GD_rm + l] * input[l*2 + 1][k];		// CT_B	Product GSW A
				} else {
					temp1 += (dct_0[GD_rm + l] *= input[l*2][k]);			// CT_A Product GSW B
					temp1 += (dct_1[GD_rm + l] *= input[l*2 + 1][k]);	// CT_B Product GSW B
				}*/
			}
		}
		(*acc)[0][0][out_k] += (temp1 *= monomial);
	}}	
	//#pragma omp barrier
}	


// GINX Accumulation as described in "Bootstrapping in FHEW-like Cryptosystems"
void FPRingGSWAccumulatorScheme::AddToACCGINXParallel(
    const std::shared_ptr<FPRingGSWCryptoParams> params,
    const FPRingGSWCiphertext &input, 
		const uint64_t &a,
    const vector<DCRTPoly> & ct,
		std::shared_ptr<vector<DCRTPoly>> outs
		) const {
	
	uint32_t K_GSW = params->GetGSWK();
	uint32_t GD_baseBit = params->GetBaseGDBit();
	uint32_t GD_Len = params -> GetLenGDsave();
	uint32_t GD_rm = params -> GetLenGDrm();
	
	// OK

	const shared_ptr<ILDCRTParams<BigInteger>> polyGSWParams = params->GetGSWPolyParams();
	const DCRTPoly &monomial = params->GetMonomial(a);
	vector<vector<DCRTPoly>> DCRT_CT(K_GSW+1);

	
	for (uint32_t i = 0; i < K_GSW+1; i++) {
		DCRT_CT[i] = *(this->SGDLevel2PARALLEL(params, ct[i], GD_baseBit, GD_rm, Format::EVALUATION));
	}
	

	// acc = dct * input (matrix product);
  // uses in-place * operators for the last call to dct[i] to gain performance
  // improvement
	// Input is GSW Sample
	// k = 0 -> A, k = 1 -> B
	
	// Output k iterating
	//
	//
	
	#pragma omp parallel for collapse(1)
	//{
      for (uint32_t out_k = 0; out_k < K_GSW+1; out_k++) {
		//std::cout<< "Inner for idx is " << out_k << ", and thread id is " << omp_get_thread_num() << std::endl;
	
		// Zero Index
		DCRTPoly temp1 = DCRTPoly(polyGSWParams, Format::EVALUATION, true);
		// Left col & right row index 
		for (uint32_t iter_k = 0; iter_k < K_GSW+1; iter_k++) {	
			// Gadget Decompositions
			for (uint32_t l = 0; l < GD_Len; l++) {
				temp1 += (DCRT_CT[iter_k][l] * input[l][iter_k][out_k]); 
			}
		}
		temp1 *= monomial;
		(*outs)[out_k] = std::move(temp1);
	}//}	
	

	/*
	//#pragma omp parallel for collapse(1)
	uint32_t N_GSW = params->GetGSWN();
	NativeInteger Q0 = params->GetModuliQ()[0];
	NativeInteger Q1 = params->GetModuliQ()[1];


	for (uint32_t out_k = 0; out_k < K_GSW+1; out_k++) {
		//std::cout<< "Inner for idx is " << out_k << ", and thread id is " << omp_get_thread_num() << std::endl;

		// Zero Index
		DCRTPoly temp1 = DCRTPoly(polyGSWParams, Format::EVALUATION, true);
		// Left col & right row index 
		#pragma omp parallel for collapse(1)
		{for (uint32_t N_idx = 0 ; N_idx < N_GSW; N_idx++) {
			for (uint32_t iter_k = 0; iter_k < K_GSW+1; iter_k++) {	
				// Gadget Decompositions
				for (uint32_t l = 0; l < GD_Len; l++) {
					temp1.GetElementW(0)[N_idx].ModAddEq( 
							DCRT_CT[iter_k][l].GetElementAtIndex(0)[N_idx] * input[l][iter_k][out_k].GetElementAtIndex(0)[N_idx],  Q0);
					
					temp1.GetElementW(1)[N_idx].ModAddEq( 
							DCRT_CT[iter_k][l].GetElementAtIndex(1)[N_idx] * input[l][iter_k][out_k].GetElementAtIndex(1)[N_idx],  Q1);

				}
			}
		}}
		temp1 *= monomial;
		(*outs)[out_k] = std::move(temp1);
	}	
	*/
}	

std::shared_ptr<vector<vector<DCRTPoly>>> FPRingGSWAccumulatorScheme::Bootstrap_and_KS(
		const std::shared_ptr<FPRingGSWCryptoParams> params, 
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const DCRTPoly &ACC_init,
    const std::shared_ptr<FPRingGSWEvalKey> &EK, const vector<NativeVector> &A, const NativeInteger &b,
		const uint32_t boot_loc_idx, const uint32_t boot_first_bits,
		const uint64_t prefix_add, const uint32_t boot_out_bit,
		const uint32_t out_num) const {

		uint32_t K_GSW = params->GetGSWK();

		auto acc = this -> BootstrapCore(params, LWEscheme,
			ACC_init,EK, A, b, boot_loc_idx, boot_first_bits, prefix_add, boot_out_bit);
  
    // Handle Prefixa
		/*  CHECK !!  */
		//DCRTPoly bNew = std::move((*acc)[0][0][K_GSW]);
		//vector<DCRTPoly> aNew = (*acc)[0][0];
		//aNew.erase(aNew.begin()+K_GSW);
		
		if (params->GetPARALLEL()) {
			#pragma omp parallel for 
			//{
             for (uint32_t k_idx = 0; k_idx < K_GSW+1; k_idx++) {
				(*acc)[0][0][k_idx].SetFormat(Format::COEFFICIENT);
			}//}
			return  this->RingKeySwitchPolyPARALLEL(params, EK->RKSkey, (*acc)[0][0], boot_loc_idx, out_num);
		} else {
			for (uint32_t k_idx = 0; k_idx < K_GSW+1; k_idx++) {
				(*acc)[0][0][k_idx].SetFormat(Format::COEFFICIENT);
			}
			return  this->RingKeySwitchPoly(params, EK->RKSkey, (*acc)[0][0], boot_loc_idx, out_num);
		}
 }

std::shared_ptr<FPLWECiphertextImpl> FPRingGSWAccumulatorScheme::PreKeySwitch(
		const std::shared_ptr<FPRingGSWCryptoParams> params
		,const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme
		,const std::shared_ptr<FPRingGSWEvalKey> &EK
		,const vector<NativeVector> &a, const uint32_t loc_idx, const uint32_t boot_first_bits
		, const uint32_t boot_out_bit
		) const {

		NativeInteger Q0 = params->GetRLWEParams()->GetModuliQ()[0];

		uint32_t N_FP = params->GetRLWEParams()->GetN();
		uint32_t K_FP = params->GetRLWEParams()->GetK();
		uint32_t l_PKS = params->GetLWEParams()->GetPKSLenSave();
		uint32_t l_rm_PKS =  params->GetLWEParams()->GetPKSLenRm();
		uint32_t PKS_Base_bit = (params->GetLWEParams()->GetBasePKSBit());
		uint32_t PKS_Base = 1;
						 PKS_Base <<=	PKS_Base_bit;
		uint32_t n_PKS = (params->GetLWEParams()->Getn());
		uint32_t q_PKS = (params->GetLWEParams()->Getq());
		uint32_t q_PKS_bit =  (params->GetLWEParams()->Getq_bit());
	
		// Test
		uint64_t scaling =  params->GetRLWEParams()->GetScaling().ConvertToInt() 	<< (boot_first_bits + boot_out_bit - q_PKS_bit );
		uint64_t s_over_2 = scaling >> 1;


		// Prepare
		NativeVector a_tot = NativeVector(n_PKS, true);
		a_tot.SetModulus(q_PKS);
		NativeInteger b_tot = NativeInteger(0);
	
		for (uint32_t k = 0; k < K_FP; k++) {
			for (uint32_t i = 0; i < N_FP; i ++) {
				uint32_t real_i_idx;
				uint64_t a_tmp;
				// Deecryption .. sign should be reverse
				if (i <= loc_idx) {
					real_i_idx = loc_idx - i;
					a_tmp = Q0.ModSub(a[k][real_i_idx],Q0).ConvertToInt();
				
				} else {
					real_i_idx = N_FP + loc_idx - i;
					a_tmp = a[k][real_i_idx].ConvertToInt();	
				}
			//Scaling Down
				a_tmp = ( (a_tmp % scaling > s_over_2) ? 
								((a_tmp / scaling) + 1 ) 
							: (a_tmp / scaling) ) % q_PKS;
				for (uint32_t j = 0; j < (l_PKS+l_rm_PKS); j++) {
					uint32_t amps = a_tmp % PKS_Base;
					a_tmp = a_tmp >> PKS_Base_bit;
					if (j < l_rm_PKS) {
						continue;
					}
					if (amps != 0) {
						a_tot += (*EK->PKSkey)[k][i][j-l_rm_PKS].GetA() * amps;
						b_tot = b_tot.ModAddEq((*EK->PKSkey)[k][i][j-l_rm_PKS].GetB() * amps, q_PKS) ;
					
					}
				}
			}
		}
		return std::make_shared<FPLWECiphertextImpl> (FPLWECiphertextImpl(a_tot, b_tot));
	}



std::shared_ptr<FPLWECiphertextImpl> FPRingGSWAccumulatorScheme::PreKeySwitchPARALLEL(
		const std::shared_ptr<FPRingGSWCryptoParams> params
		,const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme
		,const std::shared_ptr<FPRingGSWEvalKey> &EK
		,const vector<NativeVector> &a, const uint32_t loc_idx, const uint32_t boot_first_bits
		, const uint32_t boot_out_bit
		) const {

		NativeInteger Q0 = params->GetRLWEParams()->GetModuliQ()[0];

		uint32_t N_FP = params->GetRLWEParams()->GetN();
		uint32_t K_FP = params->GetRLWEParams()->GetK();
		uint32_t l_PKS = params->GetLWEParams()->GetPKSLenSave();
		uint32_t l_rm_PKS =  params->GetLWEParams()->GetPKSLenRm();
		uint32_t PKS_Base_bit = (params->GetLWEParams()->GetBasePKSBit());
		uint32_t PKS_Base = 1;
						 PKS_Base <<=	PKS_Base_bit;
		uint32_t n_PKS = (params->GetLWEParams()->Getn());
		uint32_t q_PKS = (params->GetLWEParams()->Getq());
		uint32_t q_PKS_bit =  (params->GetLWEParams()->Getq_bit());
	
		// Test
		uint64_t scaling =  params->GetRLWEParams()->GetScaling().ConvertToInt() 	<< (boot_first_bits + boot_out_bit - q_PKS_bit );
		uint64_t s_over_2 = scaling >> 1;


		// Prepare
		NativeVector a_tot = NativeVector(n_PKS, true);
		a_tot.SetModulus(q_PKS);
		NativeInteger b_tot = NativeInteger(0);
	
		for (uint32_t k = 0; k < K_FP; k++) {
			for (uint32_t i = 0; i < N_FP; i ++) {
				uint32_t real_i_idx;
				uint64_t a_tmp;
				// Deecryption .. sign should be reverse
				if (i <= loc_idx) {
					real_i_idx = loc_idx - i;
					a_tmp = Q0.ModSub(a[k][real_i_idx],Q0).ConvertToInt();
				
				} else {
					real_i_idx = N_FP + loc_idx - i;
					a_tmp = a[k][real_i_idx].ConvertToInt();	
				}
			//Scaling Down
				a_tmp = ( (a_tmp % scaling > s_over_2) ? 
								((a_tmp / scaling) + 1 ) 
							: (a_tmp / scaling) ) % q_PKS;
				for (uint32_t j = 0; j < (l_PKS+l_rm_PKS); j++) {
					uint32_t amps = a_tmp % PKS_Base;
					a_tmp = a_tmp >> PKS_Base_bit;
					if (j < l_rm_PKS) {
						continue;
					}
					if (amps != 0) {
						#pragma omp parallel for
						//{
                            for (uint32_t n_idx = 0 ; n_idx < n_PKS+1; n_idx++) {
							if (n_idx == n_PKS)	{
								b_tot.ModAddEq((*EK->PKSkey)[k][i][j-l_rm_PKS].GetB() * amps, q_PKS) ;
							} else {
								a_tot[n_idx].ModAddEq((*EK->PKSkey)[k][i][j-l_rm_PKS].GetA()[n_idx] * amps, q_PKS);
							}
						}//}
					}
				}
			}
		}
		return std::make_shared<FPLWECiphertextImpl> (FPLWECiphertextImpl(a_tot, b_tot));
	}





std::shared_ptr<FPRingGSWCiphertext> FPRingGSWAccumulatorScheme::BootstrapCore(
    const std::shared_ptr<FPRingGSWCryptoParams> params,
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const DCRTPoly &ACC_init,
    const std::shared_ptr<FPRingGSWEvalKey> &EK, const vector<NativeVector> &a, const NativeInteger &b,
		const uint32_t boot_loc_idx, const uint32_t boot_first_bits
		, const uint64_t prefix_add, const uint32_t boot_out_bit
		) const {
	
	auto core_start = std::chrono::system_clock::now();
	
	// a and b shoud be COEFFCIENT status
	if ((EK->BSkey == nullptr) || (EK->RKSkey == nullptr)) {
    std::string errMsg =
        "Bootstrapping keys have not been generated. Please call BTKeyGen "
        "before calling bootstrapping.";
    PALISADE_THROW(config_error, errMsg);
  }
	
  const shared_ptr<ILDCRTParams<BigInteger>> polyParams = params->GetGSWPolyParams();
  //NativeInteger Q0 = params->GetRLWEParams()->GetModuliQ()[0];
  //NativeInteger Q1 = params->GetRLWEParams()->GetModuliQ()[1];


	uint32_t M_GSW = params->GetGSWM();
	//uint32_t N_GSW = params->GetGSWN();
	uint32_t M_GSW_bit = params->GetGSWM_Bit();
	uint32_t K_GSW = params->GetGSWK();

	uint32_t n_PKS = (params->GetLWEParams()->Getn());

	// Make
	uint32_t q_PKS_bit =  (params->GetLWEParams()->Getq_bit());
	uint64_t q_PKS = 1;
					q_PKS <<=  q_PKS_bit;

	int32_t  diff_bt = (q_PKS_bit  - M_GSW_bit);
	int64_t diff_int = 1;
					diff_int <<= diff_bt;

	int32_t q_PKS_test = boot_first_bits + boot_out_bit - q_PKS_bit;
	if (q_PKS_test < 0) {
		PALISADE_THROW(config_error, "qPKS is too large");
	}
	
	// Test
	uint64_t scaling_test =  params->GetRLWEParams()->GetScaling().ConvertToInt();
	scaling_test <<= (boot_first_bits + boot_out_bit - q_PKS_bit );

	uint64_t b_tmp = (b.ConvertToInt()+prefix_add);
					 b_tmp = (b_tmp % scaling_test > (scaling_test >> 1)) ? 
						(	(b_tmp / scaling_test) + 1) : (	b_tmp / scaling_test); 				 
	
					b_tmp = b_tmp % (q_PKS); 

	// Pre Key Switching, We Will Cal B - AS.
	std::shared_ptr<FPLWECiphertextImpl> CT;
	if (params->GetPARALLEL()) {
		CT = PreKeySwitchPARALLEL( params , LWEscheme, EK, a, boot_loc_idx
			, boot_first_bits, boot_out_bit);
	} else {
		CT = PreKeySwitch( params , LWEscheme, EK, a, boot_loc_idx
			, boot_first_bits, boot_out_bit);
	}
	// 2^q moduler cal
	NativeVector	A_boot = CT->GetA();
	NativeInteger B_boot = CT->GetB().ModAdd(NativeInteger(b_tmp),q_PKS);
	
	
	vector<uint64_t>	A_Vals(n_PKS);
	uint64_t B_Vals;
	if (diff_bt < 0) {
		PALISADE_THROW(config_error, "qPKS is greater or equal then M_GSW");
	} else {
		// Exact Reduction
		for (uint32_t i = 0; i < n_PKS; i++) {
			A_Vals[i] = ((A_boot[i].ConvertToInt() % diff_int  > (uint64_t)(diff_int >> 1)) ? 
						( (A_boot[i].ConvertToInt() / diff_int) + 1 ) 
					: (  A_boot[i].ConvertToInt() / diff_int  )
									)  % (M_GSW >> boot_out_bit);
		}
		B_Vals = ((B_boot.ConvertToInt() % diff_int  > (uint64_t)(diff_int >> 1)) ? 
						( (B_boot.ConvertToInt() / diff_int) + 1 ) 
					: (  B_boot.ConvertToInt() / diff_int  )
									)  % (M_GSW >> boot_out_bit);

	}
	if (boot_out_bit > 0) {
		for (uint32_t i = 0; i < n_PKS; i++) {
			A_Vals[i] =  (A_Vals[i] << boot_out_bit);
		}
		B_Vals = (B_Vals << boot_out_bit);
	}
	
	//Reverse
	B_Vals = (M_GSW - B_Vals);
	if (B_Vals == M_GSW) {B_Vals = 0;} 

	
	// Ready To Bootstrap!!!

	std::vector<DCRTPoly> res(K_GSW+1);
	for (uint32_t i = 0; i < K_GSW; i ++) {
		res[i] = DCRTPoly(polyParams, Format::EVALUATION, true);		// A
	}
	res[K_GSW] = (ACC_init * params->GetRotMonomial_GSW(B_Vals));	// B

 // main accumulation computation
  // the following loop is the bottleneck of bootstrapping/binary gate
  // evaluation Mode
  auto acc = std::make_shared<FPRingGSWCiphertext>(1, 1, K_GSW+1);
  //	std::cout << "Res 0 and 1 is made " << b_index <<std::endl;
  

	
	(*acc)[0][0] = std::move(res);
	
	std::chrono::duration<double> core_end = std::chrono::system_clock::now() - core_start;
	if (DEBUG)std::cout << "Bootstrapping Core Prepare time is " <<core_end.count() << std::endl     ;


	//uint32_t real_idx;

	// BsKey [:][0][i] => is 1 if s = 1
	// BsKey [:][1][i] => is 1 if s = -1
	

	// BootstrapParams
	//uint32_t K_GSW = params->GetGSWK();
	uint32_t GD_baseBit = params->GetBaseGDBit();
	uint32_t GD_Len = params -> GetLenGDsave();
	uint32_t GD_rm = params -> GetLenGDrm();
	const shared_ptr<ILDCRTParams<BigInteger>> polyGSWParams = params->GetGSWPolyParams();
	vector<vector<DCRTPoly>> DCRT_CT(K_GSW+1);



	auto starts = std::chrono::system_clock::now();
  // New Code
	//bool New = false;
	//if (New) {
	/*
    
    if (params->GetPARALLEL()) {
        // Has a Bug!!!!1
    vector<DCRTPoly> acc_out0(K_GSW+1);
    vector<DCRTPoly> acc_out1(K_GSW+1);
    vector<DCRTPoly> acc_in(K_GSW+1); // = (*acc)[0][0];


		// Make Zeros
		for (uint32_t kk = 0; kk < K_GSW+1; kk++ ) {
			 acc_out0[kk] =  DCRTPoly(polyGSWParams, Format::EVALUATION);
			 acc_out1[kk] =  DCRTPoly(polyGSWParams, Format::EVALUATION);
			 acc_in[kk] =  DCRTPoly(polyGSWParams, Format::COEFFICIENT);
		}
     for (uint32_t i = 0; i < n_PKS; i++) {
      if (A_Vals[i] == 0) { continue;}

      //vector<DCRTPoly> &acc_in = (*acc)[0][0];  
      // Someting is ...
      #pragma omp parallel for
      //{
          for (uint32_t kk = 0; kk < K_GSW+1; kk++) {
        acc_in[kk] = (*acc)[0][0][kk];
        acc_in[kk].SetFormat(Format::COEFFICIENT);
      }//}

      for (uint32_t idx = 0; idx < K_GSW+1; idx++) {
        DCRT_CT[idx] = std::move(*(this->SGDLevel2PARALLEL(params, acc_in[idx], GD_baseBit, GD_rm, Format::EVALUATION)));
      }
			std::cout << "Before STart" << std::endl;
      //#pragma omp parallel for
      //{
      for (uint32_t mixed_idx = 0; mixed_idx < 2*(K_GSW+1); mixed_idx++) {
        uint32_t two = mixed_idx % 2;
        uint32_t out_k = mixed_idx >> 1;
        const FPRingGSWCiphertext &my_BK_mul = two == 0 ?
          (*EK->BSkey)[0].data()[0].data()[i] :
          (*EK->BSkey)[0].data()[1].data()[i];
        const DCRTPoly &my_monomial = two == 0 ?
          params->GetMonomial(A_Vals[i]) :
          params->GetMonomial(M_GSW - A_Vals[i]);

        DCRTPoly &my_acc_out = two == 0 ?
          acc_out0.data()[out_k] :
          acc_out1.data()[out_k];
				std::
				//my_acc_out =  DCRTPoly(polyGSWParams, Format::EVALUATION);
				uint64_t *op0_Q0 = my_acc_out.GetElementWForHAXL(0);
				uint64_t *op0_Q1 = my_acc_out.GetElementWForHAXL(1);
				uint64_t tmp_vec_address[N_GSW];
	      
				// Zero Index
        for (uint32_t iter_k = 0; iter_k < K_GSW+1; iter_k++) {
					std::cout << "iter idx is " << iter_k << std::endl;
					// Gadget Decompositions
          vector<DCRTPoly> &DCRT_CT_k = DCRT_CT.data()[iter_k];
          const vector<DCRTPoly> my_BK_k = my_BK_mul[iter_k].data()[out_k];
            
					for (uint32_t l = 0; l < GD_Len; l++) {
							uint64_t *op1_Q0 = DCRT_CT_k[l].GetElementWForHAXL(0);
							uint64_t *op1_Q1 = DCRT_CT_k[l].GetElementWForHAXL(1);
							const uint64_t *op2_Q0 = (uint64_t *) (&my_BK_k[l].GetElementAtIndex(0)[0]);	  
							const uint64_t *op2_Q1 = (uint64_t *) (&my_BK_k[l].GetElementAtIndex(1)[0]);	  
							

							if (l == 0 && iter_k == 0) {
								intel::hexl::EltwiseMultMod(op0_Q0, op1_Q0, op2_Q0, N_GSW, Q0.ConvertToInt(),1);
								intel::hexl::EltwiseMultMod(op0_Q1, op1_Q1, op2_Q1, N_GSW, Q1.ConvertToInt(),1);
								//my_acc_out = DCRT_CT_k.data()[l] * my_BK_k.data()[l]; 
						
							} else {
								intel::hexl::EltwiseMultMod(tmp_vec_address, op1_Q0, op2_Q0, N_GSW, Q0.ConvertToInt(),1);
								intel::hexl::EltwiseAddMod(op0_Q0, op0_Q0, tmp_vec_address, N_GSW, Q0.ConvertToInt());
								intel::hexl::EltwiseMultMod(tmp_vec_address, op1_Q1, op2_Q1, N_GSW, Q1.ConvertToInt(),1);
								intel::hexl::EltwiseAddMod(op0_Q1, op0_Q1, tmp_vec_address, N_GSW, Q1.ConvertToInt());
								//my_acc_out += DCRT_CT_k.data()[l] * my_BK_k.data()[l]; 
							}
						}
          }
					const uint64_t *my_monomial_address_Q0 = (uint64_t *) (&my_monomial.GetElementAtIndex(0)[0]);	  
					const uint64_t *my_monomial_address_Q1 = (uint64_t *) (&my_monomial.GetElementAtIndex(1)[0]);	  
					intel::hexl::EltwiseMultMod(op0_Q0, op0_Q0, my_monomial_address_Q0, N_GSW, Q0.ConvertToInt(),1);
					intel::hexl::EltwiseMultMod(op0_Q1, op0_Q1, my_monomial_address_Q1, N_GSW, Q1.ConvertToInt(),1);
					// my_acc_out *= my_monomial;
      }//}

      // Plus
      #pragma omp parallel for
      //{
      for (uint32_t ks = 0; ks < K_GSW+1; ks++) {
        (*acc)[0][0][ks] += acc_out0[ks];
        (*acc)[0][0][ks] += acc_out1[ks];
      }//}

    }
  } else {
		vector<DCRTPoly> acc_in(K_GSW+1); // = (*acc)[0][0];
  	//vector<uint64_t> tmp_vec(N_GSW);
		//uint64_t *tmp_vec_address = tmp_vec.data();
		uint64_t tmp_vec_address[N_GSW];
		
    DCRTPoly Temp = DCRTPoly(polyGSWParams, Format::EVALUATION, true);
		uint64_t *op0_Q0 = Temp.GetElementWForHAXL(0);
		uint64_t *op0_Q1 = Temp.GetElementWForHAXL(1);
	
		//start!!   
		for (uint32_t i = 0; i < n_PKS; i++) {
      if (A_Vals[i] == 0) {continue;}


     {for (uint32_t kk = 0; kk < K_GSW+1; kk++) {
				acc_in[kk] = (*acc)[0][0][kk];
				acc_in[kk].SetFormat(Format::COEFFICIENT);
			}}


      for (uint32_t idx = 0; idx < K_GSW+1; idx++) {
        DCRT_CT[idx] = std::move(*this->SGDLevel2(params, acc_in[idx], GD_baseBit, GD_rm, Format::EVALUATION));
      }

      // + rot
      //uint64_t As = (A_Vals[i] == M_GSW) ?  0 : A_Vals[i];  
      //if (A_Vals[i] == M_GSW) {A_Vals[i] = 0;}
      // handles - (-a) * E(1)  =  + a * Ea
      //athis->AddToACCGINX(params, (*EK->BSkey)[0][0][i], As, acc);
      //uint32_t N_GSW = M_GSW >> 1;

      for (uint32_t out_k = 0; out_k < (K_GSW+1); out_k++) {
        for (uint32_t two = 0; two < 2; two++){
          const FPRingGSWCiphertext &my_BK_mul = two == 0 ?
            (*EK->BSkey)[0].data()[0].data()[i] :
            (*EK->BSkey)[0].data()[1].data()[i];
          const DCRTPoly &my_monomial = two == 0 ?
            params->GetMonomial(A_Vals[i]) :
            params->GetMonomial(M_GSW - A_Vals[i]);
					
					for (uint32_t iter_k = 0; iter_k < K_GSW+1; iter_k++) {
            // Gadget Decompositions
            DCRTPoly *DCRT_CT_k = DCRT_CT.data()[iter_k].data();
            const DCRTPoly *my_BK_k = my_BK_mul[iter_k].data()[out_k].data();

            for (uint32_t l = 0; l < GD_Len; l++) {
							uint64_t *op1_Q0 = DCRT_CT_k[l].GetElementWForHAXL(0);
							uint64_t *op1_Q1 = DCRT_CT_k[l].GetElementWForHAXL(1);
							const uint64_t *op2_Q0 = (uint64_t *) (&my_BK_k[l].GetElementAtIndex(0)[0]);	  
							const uint64_t *op2_Q1 = (uint64_t *) (&my_BK_k[l].GetElementAtIndex(1)[0]);	  
						
							if (l == 0 && iter_k == 0) {
								intel::hexl::EltwiseMultMod(op0_Q0, op1_Q0, op2_Q0, N_GSW, Q0.ConvertToInt(),1);
								intel::hexl::EltwiseMultMod(op0_Q1, op1_Q1, op2_Q1, N_GSW, Q1.ConvertToInt(),1);
								//Temp =  std::move(DCRT_CT_k[l] *= my_BK_k[l]);
							} else {
								intel::hexl::EltwiseMultMod(tmp_vec_address, op1_Q0, op2_Q0, N_GSW, Q0.ConvertToInt(),1);
								intel::hexl::EltwiseAddMod(op0_Q0, op0_Q0, tmp_vec_address, N_GSW, Q0.ConvertToInt());
								intel::hexl::EltwiseMultMod(tmp_vec_address, op1_Q1, op2_Q1, N_GSW, Q1.ConvertToInt(),1);
								intel::hexl::EltwiseAddMod(op0_Q1, op0_Q1, tmp_vec_address, N_GSW, Q1.ConvertToInt());
								//Temp +=  DCRT_CT_k[l] *= my_BK_k[l];
							
							}
            }
          }
					// Making
					const uint64_t *my_monomial_address_Q0 = (uint64_t *) (&my_monomial.GetElementAtIndex(0)[0]);	  
					const uint64_t *my_monomial_address_Q1 = (uint64_t *) (&my_monomial.GetElementAtIndex(1)[0]);	  
									
					uint64_t *acc_out_k_Q0 = (*acc)[0][0][out_k].GetElementWForHAXL(0);
					uint64_t *acc_out_k_Q1 = (*acc)[0][0][out_k].GetElementWForHAXL(1);

					intel::hexl::EltwiseMultMod(op0_Q0, op0_Q0, my_monomial_address_Q0, N_GSW, Q0.ConvertToInt(),1);
					intel::hexl::EltwiseAddMod(acc_out_k_Q0, acc_out_k_Q0, op0_Q0, N_GSW, Q0.ConvertToInt());	
					intel::hexl::EltwiseMultMod(op0_Q1, op0_Q1, my_monomial_address_Q1, N_GSW, Q1.ConvertToInt(),1);
					intel::hexl::EltwiseAddMod(acc_out_k_Q1, acc_out_k_Q1, op0_Q1, N_GSW, Q1.ConvertToInt());
          //Temp *= my_monomial;
          //(*acc)[0][0][out_k] += Temp;
        }
      }
    }
  }
	*/
    //} else { // Orign Code 
	
	if (params->GetPARALLEL()) {
		vector<DCRTPoly> acc_out0(K_GSW+1);
		vector<DCRTPoly> acc_out1(K_GSW+1);
		// Make Zeros
		for (uint32_t ks = 0; ks < K_GSW+1; ks++){
			//#ifdef WITH_INTEL_HEXL
			acc_out0[ks] = DCRTPoly(polyParams, Format::EVALUATION);
			acc_out1[ks] = DCRTPoly(polyParams, Format::EVALUATION);
			//#else
			//acc_out0[ks] = DCRTPoly(polyParams, Format::EVALUATION, true);
			//acc_out1[ks] = DCRTPoly(polyParams, Format::EVALUATION, true);
			//#endif
		}
		vector<DCRTPoly> acc_in(K_GSW+1); // = (*acc)[0][0];
		
		// START!!
		for (uint32_t i = 0; i < n_PKS; i++) {
			if (A_Vals[i] == 0) {continue;}

			/* is Not Working ...
			//#pragma omp parallel for
			{for (uint32_t kk_and_Q = 0; kk_and_Q < (K_GSW+1)*2; kk_and_Q++) {
				uint32_t kk = kk_and_Q >> 1;
				uint32_t Q_idx = kk_and_Q % 2;
				std::cout << "idx is " << kk_and_Q << std::endl;
				acc_in[kk].SetElementAtIndex(Q_idx, (*acc)[0][0][kk].GetElementAtIndex(Q_idx) );
				std::cout << "Save is over" << std::endl;
				acc_in[kk].SetFormatAtIndexW(Format::COEFFICIENT, Q_idx);
				//std::cout << "acc_in is " << acc_in[kk] << std::endl;
			}}
			std::cout << "save is over " << std::endl;
			for (uint32_t kk  = 0; kk < (K_GSW+1); kk++) {
				acc_in[kk].SetOutFormatW(Format::COEFFICIENT);
			}
			*/

			#pragma omp parallel for
			//{
      for (uint32_t kk  = 0; kk  < K_GSW+1; kk++) {
				acc_in[kk] = (*acc)[0][0][kk];
				acc_in[kk].SetFormat(Format::COEFFICIENT);
			}//}
	
			for (uint32_t idx = 0; idx < K_GSW+1; idx++) {
				DCRT_CT[idx] = std::move(*(this->SGDLevel2PARALLEL(params, acc_in[idx], GD_baseBit, GD_rm, Format::EVALUATION)));
			}
	
			//std::cout << "SGDLEvel2 is over" << std::endl;
			#pragma omp parallel for
			//{
      for (uint32_t mixed_idx = 0; mixed_idx < 2*(K_GSW+1); mixed_idx++) {
				uint32_t two = mixed_idx % 2;
				uint32_t out_k = mixed_idx >> 1;
	
				//std::cout<< "out for idx is " << two << ", and thread id is " << omp_get_thread_num() << std::endl;
				if (two == 0) {
					//if (A_Vals[i] == M_GSW) {A_Vals[i] = 0;}
					
					uint64_t As = (A_Vals[i] == M_GSW) ?  0 : A_Vals[i];  
					// handles - (-a) * E(1)  =  + a * E
					//this->AddToACCGINXParallel(params, (*EK->BSkey)[0][0][i], As, acc_in, acc_out0);
					
					// Hardcoding ACCGINX
					const DCRTPoly &monomial = params->GetMonomial(As);	
					//acc_out0[out_k] =  DCRTPoly(polyGSWParams, Format::EVALUATION, true);
					// Zero Index
					for (uint32_t iter_k = 0; iter_k < K_GSW+1; iter_k++) {	
						// Gadget Decompositions
						for (uint32_t l = 0; l < GD_Len; l++) {
							//std::cout << "iter k is " <<iter_k << ", l is " << l << std::endl;

							if (l == 0 && iter_k == 0) {
								acc_out0[out_k] = (DCRT_CT[iter_k][l] * ((*EK->BSkey)[0][0][i])[iter_k][out_k][l]); 
						
							} else {
								acc_out0[out_k] += (DCRT_CT[iter_k][l] * ((*EK->BSkey)[0][0][i])[iter_k][out_k][l]); 
							}
						}
					}
					acc_out0[out_k] *= monomial;

				} else {
					uint64_t As = ( M_GSW - A_Vals[i]);
					if (As == M_GSW) {As = 0;}
					// handles - (-a) * E(-1) =   - a*E
					//this->AddToACCGINXParallel(params, (*EK->BSkey)[0][1][i], As, acc_in, acc_out1);
					// Hardcoding ACCGINX
					const DCRTPoly &monomial = params->GetMonomial(As);	
					//acc_out1[out_k] =  DCRTPoly(polyGSWParams, Format::EVALUATION, true);
					// Zero Index
					for (uint32_t iter_k = 0; iter_k < K_GSW+1; iter_k++) {	
						// Gadget Decompositions
						for (uint32_t l = 0; l < GD_Len; l++) {
							if (l == 0 && iter_k == 0) {
								acc_out1[out_k] = (DCRT_CT[iter_k][l] * ((*EK->BSkey)[0][1][i])[iter_k][out_k][l]); 
							} else {
								acc_out1[out_k] += (DCRT_CT[iter_k][l] * ((*EK->BSkey)[0][1][i])[iter_k][out_k][l]); 
							}
						}
					}
					acc_out1[out_k] *= monomial;
				}
			}//}

			// Plus
			#pragma omp parallel for
			//{
      for (uint32_t ks = 0; ks < K_GSW+1; ks++) {
				(*acc)[0][0][ks] += acc_out0[ks];
				(*acc)[0][0][ks] += acc_out1[ks];
			}//}

		}
	} else {
		// A_Vals should be fliped!! But we do not. 
		vector<DCRTPoly> acc_in(K_GSW+1); //= (*acc)[0][0];	
		
		for (uint32_t i = 0; i < n_PKS; i++) {
			if (A_Vals[i] == 0) {continue;}

			//#pragma omp parallel for
			{for (uint32_t kk = 0; kk < K_GSW+1; kk++) {
				acc_in[kk] = (*acc)[0][0][kk];
				acc_in[kk].SetFormat(Format::COEFFICIENT);
			}}
			for (uint32_t idx = 0; idx < K_GSW+1; idx++) {
				DCRT_CT[idx] = std::move(*(this->SGDLevel2(params, acc_in[idx], GD_baseBit, GD_rm, Format::EVALUATION)));
			}

			// + rot
			uint64_t As = (A_Vals[i] == M_GSW) ?  0 : A_Vals[i];  
			//if (A_Vals[i] == M_GSW) {A_Vals[i] = 0;}
			// handles - (-a) * E(1)  =  + a * Ea
			//this->AddToACCGINX(params, (*EK->BSkey)[0][0][i], As, acc);
			const DCRTPoly &monomial1 = params->GetMonomial(As);	
			// Zero Index
			for (uint32_t out_k = 0; out_k < (K_GSW+1); out_k++) {
				//DCRTPoly Temp; // =  DCRTPoly(polyGSWParams, Format::EVALUATION, true);
				//#ifdef WITH_INTEL_HEXL
					DCRTPoly Temp; // =  DCRTPoly(polyGSWParams, Format::EVALUATION, true);
				//#else
				//	DCRTPoly Temp =  DCRTPoly(polyGSWParams, Format::EVALUATION, true);
				//#endif
	

				for (uint32_t iter_k = 0; iter_k < K_GSW+1; iter_k++) {	
					// Gadget Decompositions
					for (uint32_t l = 0; l < GD_Len; l++) {
						if (l == 0 && iter_k == 0) {
							Temp = (DCRT_CT[iter_k][l] * ((*EK->BSkey)[0][0][i])[iter_k][out_k][l]); 
						} else {
							Temp += (DCRT_CT[iter_k][l] * ((*EK->BSkey)[0][0][i])[iter_k][out_k][l]); 
						}
					}
				}
				Temp *= monomial1;
				(*acc)[0][0][out_k] += Temp;
			}

			// - rot
			// handles - (-a) * E(-1) =   - a*E
			As = ( M_GSW - A_Vals[i]);
			As = ( As == M_GSW) ?  0 : As;  
			//this->AddToACCGINX(params, (*EK->BSkey)[0][1][i], tmp_As, acc);
			
			const DCRTPoly &monomial2 = params->GetMonomial(As);	
			// Zero Index
			for (uint32_t out_k = 0; out_k < (K_GSW+1); out_k++) {
				
				//#ifdef WITH_INTEL_HEXL
					DCRTPoly Temp; // =  DCRTPoly(polyGSWParams, Format::EVALUATION, true);
				//#else
				//	DCRTPoly Temp =  DCRTPoly(polyGSWParams, Format::EVALUATION, true);
				//#endif
				for (uint32_t iter_k = 0; iter_k < K_GSW+1; iter_k++) {	
					// Gadget Decompositions
					for (uint32_t l = 0; l < GD_Len; l++) {
						if (l == 0 && iter_k == 0) {
							if (out_k == K_GSW) {
								DCRT_CT[iter_k][l] *= ((*EK->BSkey)[0][1][i])[iter_k][out_k][l];
								Temp = std::move(DCRT_CT[iter_k][l]);
							} else {
								Temp = (DCRT_CT[iter_k][l] * ((*EK->BSkey)[0][1][i])[iter_k][out_k][l]);
							}
						}	else {
							if (out_k == K_GSW) {
								Temp += (DCRT_CT[iter_k][l] *= ((*EK->BSkey)[0][1][i])[iter_k][out_k][l]);
							} else {
								Temp += (DCRT_CT[iter_k][l] * ((*EK->BSkey)[0][1][i])[iter_k][out_k][l]);
							}
						}
					}
				}
				Temp *= monomial2;
				(*acc)[0][0][out_k] += Temp;
			}
		}	
	}
	//} // NEw

	std::chrono::duration<double> ends = std::chrono::system_clock::now() - starts;
	if (DEBUG) std::cout << "Bootstrapping External Prod time is " << ends.count() << std::endl;
		
	return acc;
}
// Full evaluation as described in "Bootstrapping in FHEW-like
// Cryptosystems"


std::shared_ptr<const FPRLWECiphertextImpl> FPRingGSWAccumulatorScheme::BootstrapOrign(
    const std::shared_ptr<FPRingGSWCryptoParams> params, 
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const std::shared_ptr<FPRingGSWEvalKey> EK,
    //const std::shared_ptr<const FPRLWECiphertextImpl> ct
		const vector<DCRTPoly> &A_poly,
		const DCRTPoly & B_poly
	) const {


	// Level 2 should handle
	uint32_t EN_len = (uint32_t ) std::ceil( (double) 64.0 / (double) params -> GetBaseEncodeBit() );

	uint32_t K_FP		= params->GetRLWEParams()->GetK();
	uint32_t K_GSW	= params->GetGSWK();

	vector<NativeVector> a_poly_vec(K_FP);
	NativePoly tmp;
	for (uint32_t i = 0; i < K_FP; i++) {
		tmp = A_poly[i].GetElementAtIndex(0);
		tmp.SetFormat(Format::COEFFICIENT);
		a_poly_vec[i] = tmp.GetValues();
	}

	NativePoly b_poly = B_poly.GetElementAtIndex(0);
	b_poly.SetFormat(Format::COEFFICIENT);

	uint32_t carry_bit = 2;
	uint32_t out_nums = 1;
	// Padding 1 (it must be zero) + 2bit
	
	uint32_t start_bits = (params->GetRLWEParams()->GetBits()[1]  + 3 + 1);
	uint64_t prefix_add_vals = (params->GetRLWEParams()->GetModuliQ()[1].ConvertToInt() >> 1); 
	
	vector<vector<DCRTPoly>> outs(EN_len);
	
	for (uint32_t idx = 0; idx < EN_len; idx ++) {

		NativeInteger b = b_poly.GetValues()[idx];
		// a and b should handle in BootstrapCore
		
		auto startTimeACC = std::chrono::system_clock::now();
		auto acc = BootstrapCore(params, LWEscheme,
			params->GetACCMul(0),EK, a_poly_vec, b, idx, start_bits, prefix_add_vals, carry_bit);
	
		std::chrono::duration<double> endTimeACC = std::chrono::system_clock::now() - startTimeACC;
		

		auto startTimeKS = std::chrono::system_clock::now();
		// Coeff mode
		//for (uint32_t k_idx = 0; k_idx < K_GSW+1; k_idx++) {
		//	(*acc)[0][0][k_idx].SetFormat(Format::COEFFICIENT);
		//}
		std::shared_ptr<vector<vector<DCRTPoly>>> KS_ct;
		if (params->GetPARALLEL()) { 
			#pragma omp parallel for 
			//{
              for (uint32_t k_idx = 0; k_idx < K_GSW+1; k_idx++) {
				(*acc)[0][0][k_idx].SetFormat(Format::COEFFICIENT);
			}//}
			KS_ct =  this->RingKeySwitchPolyPARALLEL(params, EK->RKSkey, (*acc)[0][0], idx, out_nums);
		}else {
			for (uint32_t k_idx = 0; k_idx < K_GSW+1; k_idx++) {
				(*acc)[0][0][k_idx].SetFormat(Format::COEFFICIENT);
			}
			KS_ct =  this->RingKeySwitchPoly(params, EK->RKSkey, (*acc)[0][0], idx, out_nums);
		}
		outs[idx] = std::move((*KS_ct)[0]);
		std::chrono::duration<double> endTimeKS = std::chrono::system_clock::now() - startTimeKS;

		if (DEBUG) std::cout << "Orign Boots idx is " << idx
			<< ", ACC Time is " << endTimeACC.count()
			<< ", KS  Time is " << endTimeKS.count()
			<< std::endl;
	
	}
	
	// New
	for (uint32_t idx=1; idx < EN_len; idx++) {
		for (uint32_t k_idx = 0; k_idx < K_FP+1; k_idx++) {
			outs[0][k_idx] += (outs[idx][k_idx]);
		}
	}
	// A, B separrate
	DCRTPoly bNew = std::move(outs[0][K_FP]);
	vector<DCRTPoly> aNew = std::move(outs[0]);
	aNew.erase(aNew.begin()+K_FP);

	return std::make_shared<FPRLWECiphertextImpl>(std::move(aNew), std::move(bNew));
	
}



std::shared_ptr<vector<DCRTPoly>> FPRingGSWAccumulatorScheme::BootstrapOrignPartial(
    const std::shared_ptr<FPRingGSWCryptoParams> params, 
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const std::shared_ptr<FPRingGSWEvalKey> EK,
		const vector<DCRTPoly> &CTs,
		const uint32_t from_bit,
		const uint32_t to_bit
	) const {


	// Level 2 should handle
	//uint32_t EN_len = (uint32_t ) std::ceil( (double) 64.0 / (double) params -> GetBaseEncodeBit() );
	uint32_t EN_len = to_bit - from_bit + 1;

	uint32_t K_FP		= params->GetRLWEParams()->GetK();
	uint32_t K_GSW	= params->GetGSWK();
	//uint32_t M_FP   = params->GetRLWEParams()->GetN() * 2;
	vector<NativeVector> a_poly_vec(K_FP);
	NativePoly tmp;
	for (uint32_t i = 0; i < K_FP; i++) {
		tmp = CTs[i].GetElementAtIndex(0);
		tmp.SetFormat(Format::COEFFICIENT);
		a_poly_vec[i] = tmp.GetValues();
	}

	NativePoly b_poly = CTs[K_FP].GetElementAtIndex(0);
	b_poly.SetFormat(Format::COEFFICIENT);

	uint32_t carry_bit = 2;
	uint32_t out_nums = 1;
	// Padding 1 (it must be zero) + 2bit
	
	uint32_t start_bits = (params->GetRLWEParams()->GetBits()[1]  + 3 + 1);
	uint64_t prefix_add_vals = (params->GetRLWEParams()->GetModuliQ()[1].ConvertToInt() >> 1); 
	
	vector<vector<DCRTPoly>> outs(EN_len);
	if (params->GetPARALLEL()) {
		#pragma omp parallel for 	
		//{
            for (uint32_t idx = 0; idx < EN_len; idx++) {
			uint32_t real_idx = idx + from_bit;
			NativeInteger b = b_poly.GetValues()[real_idx];
			// a and b should handle in BootstrapCore
		
			auto startTimeACC = std::chrono::system_clock::now();
			auto acc = BootstrapCore(params, LWEscheme,
				params->GetACCMul(0),EK, a_poly_vec, b, real_idx, start_bits, prefix_add_vals, carry_bit);
	
			std::chrono::duration<double> endTimeACC = std::chrono::system_clock::now() - startTimeACC;
			auto startTimeKS = std::chrono::system_clock::now();
		
			std::shared_ptr<vector<vector<DCRTPoly>>> KS_ct;
			for (uint32_t k_idx = 0; k_idx < K_GSW+1; k_idx++) {
				(*acc)[0][0][k_idx].SetFormat(Format::COEFFICIENT);
			}
			KS_ct =  this->RingKeySwitchPoly(params, EK->RKSkey, (*acc)[0][0], real_idx, out_nums);
			outs[idx] = std::move((*KS_ct)[0]);
			std::chrono::duration<double> endTimeKS = std::chrono::system_clock::now() - startTimeKS;

			if (DEBUG) std::cout << "Orign Boots idx is " << idx
			<< ", ACC Time is " << endTimeACC.count()
			<< ", KS  Time is " << endTimeKS.count()
			<< std::endl;
		}//}
	} else {
		for (uint32_t idx = 0; idx < EN_len; idx++) {
			uint32_t real_idx = idx + from_bit;
			NativeInteger b = b_poly.GetValues()[real_idx];
			// a and b should handle in BootstrapCore
		
			auto startTimeACC = std::chrono::system_clock::now();
			auto acc = BootstrapCore(params, LWEscheme,
				params->GetACCMul(0),EK, a_poly_vec, b, real_idx, start_bits, prefix_add_vals, carry_bit);
	
			std::chrono::duration<double> endTimeACC = std::chrono::system_clock::now() - startTimeACC;
			auto startTimeKS = std::chrono::system_clock::now();
		
			std::shared_ptr<vector<vector<DCRTPoly>>> KS_ct;
			for (uint32_t k_idx = 0; k_idx < K_GSW+1; k_idx++) {
				(*acc)[0][0][k_idx].SetFormat(Format::COEFFICIENT);
			}
			KS_ct =  this->RingKeySwitchPoly(params, EK->RKSkey, (*acc)[0][0], real_idx, out_nums);
		
			outs[idx] = std::move((*KS_ct)[0]);
			std::chrono::duration<double> endTimeKS = std::chrono::system_clock::now() - startTimeKS;

			if (DEBUG) std::cout << "Orign Boots idx is " << idx
				<< ", ACC Time is " << endTimeACC.count()
				<< ", KS  Time is " << endTimeKS.count()
				<< std::endl;
		}
	}

	// New
	//uint32_t rot_idx = M_FP - from_bit;
	//std::cout << "rot is " << rot_idx << std::endl;
	for (uint32_t k_idx = 0; k_idx < K_FP+1; k_idx++) {
		for (uint32_t idx = 1; idx < EN_len; idx++) {
			outs[0][k_idx] += (outs[idx][k_idx]);
		}
		//outs[0][k_idx] *= params->GetRotMonomial_RLWE(from_bit);
	}
	// A, B separrate
	//DCRTPoly bNew = std::move(outs[0][K_FP]);
	//vector<DCRTPoly> aNew = std::move(outs[0]);
	//aNew.erase(aNew.begin()+K_FP);

	return std::make_shared<vector<DCRTPoly>>(std::move(outs[0]));
	
}


  std::shared_ptr<vector<vector<DCRTPoly>>> FPRingGSWAccumulatorScheme::RingKeySwitchPoly(
    const std::shared_ptr<FPRingGSWCryptoParams> params,
    const std::shared_ptr<FPRLWESwitchingKey> K,
    const vector<DCRTPoly> &inputs,
    const uint32_t &orign_idx,
    const uint32_t &out_nums
    ) const {
    auto start = std::chrono::system_clock::now();
  
    //CT should be COEFFICIENT
    if (inputs[0].GetFormat() != Format::COEFFICIENT) {
      PALISADE_THROW(config_error, "RKS input should be COEFF Mode!!");
    } else if ( inputs[inputs.size()-1].GetFormat() != Format::COEFFICIENT) {
      PALISADE_THROW(config_error, "RKS input should be COEFF Mode!!");
    }    

    uint32_t GSW_N = params->GetGSWN();
    uint32_t GSW_K = params->GetGSWK();
    uint32_t FP_K = params->GetRLWEParams()->GetK();
    uint32_t FP_N = params->GetRLWEParams()->GetN();
    uint32_t KS_base_bit = params->GetBaseKSBit();
    uint32_t KS_len = params->GetLenKSsave();
    uint32_t KS_rm = params->GetLenKSrm();
    NativeInteger Q0 = params -> GetModuliQ()[0];
    NativeInteger Q1 = params -> GetModuliQ()[1];
		uint64_t Q0_nat = Q0.ConvertToInt();
		uint64_t Q1_nat = Q1.ConvertToInt();


		#ifdef WITH_INTEL_HEXL
  
		#else
		NativeInteger Mu0 = params->GetRLWEParams()->GetMuFP()[0];
		NativeInteger Mu1 = params->GetRLWEParams()->GetMuFP()[1];
		#endif

		//int64_t Q0_over_2 = Q0.ConvertToInt() >> 1;
		//int64_t a_tmp;
		int32_t a_tmp2;
		bool a_tmp2_positive;
    bool s_positive;
    NativeInteger Tmp_Nat;
    // Make Out
    vector<vector<DCRTPoly>> out_CT(out_nums);
    for (uint32_t i = 0; i < out_nums; i ++) {
      out_CT[i].resize(FP_K+1);
    }    

		vector<std::shared_ptr<vector<DCRTPoly>>> CT_dec(GSW_K);
    vector<std::shared_ptr<vector<vector<int32_t>>>> CT_dec2(GSW_K);
	  //auto Make_start = std::chrono::system_clock::now();
 	
    // K, N, KS
    for (uint32_t i = 0; i < GSW_K; i++) {
      CT_dec2[i] = SGDLevel2ForKS(params, inputs[i], KS_base_bit, KS_rm, Format::COEFFICIENT);
    }    
     // dim0 : N / dim1: K. dim2: KS_len
    const vector<vector<vector<FPRLWECiphertextImpl>>> &KS = K->GetElements();
    for (uint32_t out_idx = 0; out_idx < out_nums; out_idx++) {
      // Setting A
			for (uint32_t out_k = 0; out_k < FP_K; out_k ++) {
        out_CT[out_idx][out_k] =  DCRTPoly(params->GetRLWEParams()->GetDCRTParams(), Format::EVALUATION, true);
      }    
     
      // Setting B 
      out_CT[out_idx][FP_K] = DCRTPoly(params->GetRLWEParams()->GetDCRTParams(), Format::COEFFICIENT,true);
      out_CT[out_idx][FP_K].GetElementW(0)[0] = inputs[GSW_K].GetElementAtIndex(0)[out_idx];
      out_CT[out_idx][FP_K].GetElementW(1)[0] = inputs[GSW_K].GetElementAtIndex(1)[out_idx];
      out_CT[out_idx][FP_K].SetFormat(Format::EVALUATION);
    }    


		//std::chrono::duration<double> Make_end = std::chrono::system_clock::now() - Make_start;
    //if (DEBUG) std::cout << "Making time is " << Make_end.count() << std::endl;
   

    uint32_t real_j_idx;  
    // Makes  
    for(uint32_t in_k_idx = 0; in_k_idx < GSW_K; in_k_idx++) {
      const vector<int32_t> *CT_dec2_k =  (*CT_dec2.data()[in_k_idx]).data();
      const vector<FPRLWECiphertextImpl> *KS_k  = KS.data()[in_k_idx].data();
    
      for (uint32_t out_idx = 0; out_idx < out_nums; out_idx++) {
        DCRTPoly *out_CT_i = out_CT.data()[out_idx].data(); 
        
				for (uint32_t j = 0; j < GSW_N; j++){
          if (j <= out_idx) {
            real_j_idx = out_idx - j; 
            s_positive = true;
          } else {
            real_j_idx = GSW_N + out_idx - j; 
            s_positive = false;
          }    
          const FPRLWECiphertextImpl *KS_k_j = KS_k[real_j_idx].data();
          const int32_t *CT_dec2_k_j = CT_dec2_k[j].data();
    
          for (uint32_t i = 0; i < KS_len; i++) {
            // IDX adjust
            a_tmp2 = CT_dec2_k_j[i];
            //a_tmp2 = (*CT_dec2[in_k_idx])[j][i];
            if (a_tmp2 == 0) { 
              continue;
            } else if (a_tmp2 > 0 ) {
              a_tmp2_positive = true;
            } else {
              a_tmp2_positive = false;
              a_tmp2 *= -1;
            }    
     
            const DCRTPoly *KS_currentA = KS_k_j[i].GetA().data();
            const DCRTPoly &KS_currentB = KS_k_j[i].GetB();
						uint64_t a_tmp2_Q0;
						uint64_t a_tmp2_Q1;
            //Calculation!!
            if (s_positive ^ a_tmp2_positive) {
							a_tmp2_Q0 = a_tmp2;
							a_tmp2_Q1 = a_tmp2;
						} else {
							a_tmp2_Q0 =	Q0_nat - a_tmp2;
							a_tmp2_Q1 =	Q1_nat - a_tmp2;
						} 
            // - a * (s[real_j]) or a * (-s[real_j]), Just Adding
            for (uint32_t out_k = 0; out_k < FP_K; out_k++) {
							//vector Out_CT_Q0
              const NativePoly &KS_A_k_Q0 = KS_currentA[out_k].GetElementAtIndex(0);
              const NativePoly &KS_A_k_Q1 = KS_currentA[out_k].GetElementAtIndex(1);
              NativePoly &out_CT_i_k_Q0 = out_CT_i[out_k].GetElementW(0);
              NativePoly &out_CT_i_k_Q1 = out_CT_i[out_k].GetElementW(1);
							
							#ifdef WITH_INTEL_HEXL
							uint64_t *op1_Q0 = reinterpret_cast<uint64_t  *> (&out_CT_i_k_Q0[0]);
							const uint64_t *op3_Q0 = (uint64_t *) (&KS_A_k_Q0[0]);
							intel::hexl::EltwiseFMAMod(op1_Q0, op3_Q0, a_tmp2_Q0
								, op1_Q0, FP_N, Q0_nat,1);

							uint64_t *op1_Q1 = reinterpret_cast<uint64_t  *> (&out_CT_i_k_Q1[0]);
							const uint64_t *op3_Q1 = (uint64_t *) (&KS_A_k_Q1[0]);
							intel::hexl::EltwiseFMAMod(op1_Q1, op3_Q1, a_tmp2_Q1
								, op1_Q1, FP_N, Q1_nat,1);
							#else
							
              for (uint32_t out_n = 0; out_n < FP_N; out_n++) {
								//out_CT[out_idx][out_k] += (KS[in_k_idx][real_j_idx][i].GetA()[out_k] * a_tmp2) ;
								// A 
								/*	
                Tmp_Nat = KS[in_k_idx][real_j_idx][i].GetA()[out_k].GetElementAtIndex(0)[out_n];    
                Tmp_Nat.MulEq(a_tmp2);
                Tmp_Nat.ModEq(Q0);
                out_CT[out_idx][out_k].GetElementW(0)[out_n].AddEq(Tmp_Nat);
                    
                Tmp_Nat = KS[in_k_idx][real_j_idx][i].GetA()[out_k].GetElementAtIndex(1)[out_n];    
                Tmp_Nat.MulEq(a_tmp2);
                Tmp_Nat.ModEq(Q1);        
                out_CT[out_idx][out_k].GetElementW(1)[out_n].AddEq(Tmp_Nat);
								*/
                  
                Tmp_Nat = KS_A_k_Q0[out_n];
								//if (a_tmp2_Q0 < 2048) {
								//Tmp_Nat.MulEq(a_tmp2_Q0);
								//Tmp_Nat.ModEq(Q0, Mu0);
								//} else {
								Tmp_Nat.ModMulFastEq(a_tmp2_Q0, Q0, Mu0);
								//Tmp_Nat.ModMulEq(a_tmp2_Q0, Q0);
								//}
                out_CT_i_k_Q0[out_n].ModAddFastEq(Tmp_Nat,Q0);
								//out_CT_i_k_Q0[out_n].ModAddEq(Tmp_Nat,Q0);
								//out_CT_i_k_Q0[out_n].AddEq(Tmp_Nat);
								//std::cout <<" is nt>??" << std::endl;
                Tmp_Nat = KS_A_k_Q1[out_n];
								//if (a_tmp2_Q1 < 2048) {
								//	Tmp_Nat.MulEq(a_tmp2_Q1);
								//	Tmp_Nat.ModEq(Q1, Mu1);
								//} else {
								Tmp_Nat.ModMulFastEq(a_tmp2_Q1, Q1, Mu1);
								//Tmp_Nat.ModMulEq(a_tmp2_Q1, Q1);
								//}
								out_CT_i_k_Q1[out_n].ModAddFastEq(Tmp_Nat,Q1);
								//out_CT_i_k_Q1[out_n].ModAddEq(Tmp_Nat,Q1);
								//out_CT_i_k_Q1[out_n].AddEq(Tmp_Nat);
              }
							#endif
            }      
            const NativePoly &KS_B_k_Q0 = KS_currentB.GetElementAtIndex(0);
            const NativePoly &KS_B_k_Q1 = KS_currentB.GetElementAtIndex(1);
            NativePoly &out_CT_i_k_Q0 = out_CT_i[FP_K].GetElementW(0);
            NativePoly &out_CT_i_k_Q1 = out_CT_i[FP_K].GetElementW(1);
            
						
						#ifdef WITH_INTEL_HEXL
						uint64_t *op1_Q0 = reinterpret_cast<uint64_t  *> (&out_CT_i_k_Q0[0]);
						//uint64_t op2		 = (uint64_t) (a_tmp2);
						const uint64_t *op3_Q0 = (uint64_t *) (&KS_B_k_Q0[0]);
						intel::hexl::EltwiseFMAMod(op1_Q0, op3_Q0, a_tmp2_Q0
							, op1_Q0, FP_N, Q0_nat,1);

						uint64_t *op1_Q1 = reinterpret_cast<uint64_t  *> (&out_CT_i_k_Q1[0]);
						const uint64_t *op3_Q1 = (uint64_t *) (&KS_B_k_Q1[0]);
						intel::hexl::EltwiseFMAMod(op1_Q1, op3_Q1, a_tmp2_Q1
							, op1_Q1, FP_N, Q1_nat,1);
						#else
						
						//out_CT[out_idx][FP_K] += (KS[in_k_idx][real_j_idx][i].GetB() * a_tmp2);
            for (uint32_t out_n = 0; out_n < FP_N; out_n++) {
              // B sum
							/*  
              Tmp_Nat = KS[in_k_idx][real_j_idx][i].GetB().GetElementAtIndex(0)[out_n];    
              Tmp_Nat.MulEq(a_tmp2);
              Tmp_Nat.ModEq(Q0);
              out_CT[out_idx][FP_K].GetElementW(0)[out_n].AddEq(Tmp_Nat);
                
              Tmp_Nat = KS[in_k_idx][real_j_idx][i].GetB().GetElementAtIndex(1)[out_n];    
              Tmp_Nat.MulEq(a_tmp2);
              Tmp_Nat.ModEq(Q1);
              out_CT[out_idx][FP_K].GetElementW(1)[out_n].AddEq(Tmp_Nat);
							*/
								
							Tmp_Nat = KS_B_k_Q0[out_n];
							//if (a_tmp2_Q0 < 2048) {
							//	Tmp_Nat.MulEq(a_tmp2_Q0);
							//	Tmp_Nat.ModEq(Q0, Mu1);
							//} else {
							Tmp_Nat.ModMulFastEq(a_tmp2_Q0, Q0, Mu0);
							//Tmp_Nat.ModMulEq(a_tmp2_Q0, Q0);
							//}
							out_CT_i_k_Q0[out_n].ModAddFastEq(Tmp_Nat,Q0);
            	//out_CT_i_k_Q0[out_n].ModAddEq(Tmp_Nat,Q0);
              //out_CT_i_k_Q0[out_n].AddEq(Tmp_Nat);
                    
              Tmp_Nat = KS_B_k_Q1[out_n];
							//if (a_tmp2_Q1 < 2048) {
							//	Tmp_Nat.MulEq(a_tmp2_Q1);
							//	Tmp_Nat.ModEq(Q1, Mu1);
							//} else {
							Tmp_Nat.ModMulFastEq(a_tmp2_Q1, Q1, Mu1);
							//Tmp_Nat.ModMulEq(a_tmp2_Q1, Q1);
							//}
              out_CT_i_k_Q1[out_n].ModAddFastEq(Tmp_Nat,Q1);
              //out_CT_i_k_Q1[out_n].ModAddEq(Tmp_Nat,Q1);
             //out_CT_i_k_Q1[out_n].AddEq(Tmp_Nat);
                
            }
						#endif
           
          }  // out_idx is onver
        } // in GSW_N is over
                
      } // Out idx
      // We must refreshing Q0, Q1 has 
      // Refresing
				
			#ifdef WITH_INTEL_HEXL
				
			#else
			/*
      for (uint32_t out_idx = 0; out_idx < out_nums; out_idx++) {
        DCRTPoly *out_CT_i = out_CT[out_idx].data();
        for (uint32_t out_k = 0; out_k < FP_K; out_k++) {
          NativePoly & out_CT_i_k_Q0 = out_CT_i[out_k].GetElementW(0);
          NativePoly & out_CT_i_k_Q1 = out_CT_i[out_k].GetElementW(1);
          for (uint32_t out_n = 0; out_n < FP_N; out_n++) {
            out_CT_i_k_Q0[out_n].ModEq(Q0, Mu0); 
            out_CT_i_k_Q1[out_n].ModEq(Q1, Mu1); 
                 
          }    
        }        
        NativePoly & out_CT_i_k_Q0 = out_CT_i[FP_K].GetElementW(0);
        NativePoly & out_CT_i_k_Q1 = out_CT_i[FP_K].GetElementW(1);
        for (uint32_t out_n = 0; out_n < FP_N; out_n++) {
          out_CT_i_k_Q0[out_n].ModEq(Q0, Mu0);
          out_CT_i_k_Q1[out_n].ModEq(Q1, Mu1);   
        }    
      } 
			*/
			#endif
				
		} // GSW_K is over 
    // Rotation
    if (orign_idx != 0) { 
      for (uint32_t out_idx = 0; out_idx < out_nums; out_idx++) {
     
        for (uint32_t out_k = 0; out_k < FP_K; out_k++) {
          out_CT[out_idx][out_k] *= params->GetRotMonomial_RLWE(orign_idx);;
        }  
        out_CT[out_idx][FP_K] *= params->GetRotMonomial_RLWE(orign_idx);
      }      
    }          
    std::chrono::duration<double> end = std::chrono::system_clock::now() - start;
    if (DEBUG) std::cout << "Ring Key Switch time is " << end.count() << std::endl;
    return std::make_shared<vector<vector<DCRTPoly>>> (out_CT);
  }



  std::shared_ptr<vector<vector<DCRTPoly>>> FPRingGSWAccumulatorScheme::RingKeySwitchPolyPARALLEL(
    const std::shared_ptr<FPRingGSWCryptoParams> params,
    const std::shared_ptr<FPRLWESwitchingKey> K,
    const vector<DCRTPoly> &inputs,
    const uint32_t &orign_idx,
    const uint32_t &out_nums
    ) const {
    auto start = std::chrono::system_clock::now();

    //CT should be COEFFICIENT
    if (inputs[0].GetFormat() != Format::COEFFICIENT) {
      PALISADE_THROW(config_error, "RKS input should be COEFF Mode!!");
    } else if ( inputs[inputs.size()-1].GetFormat() != Format::COEFFICIENT) {
      PALISADE_THROW(config_error, "RKS input should be COEFF Mode!!");
    }

    uint32_t GSW_N = params->GetGSWN();
    uint32_t GSW_K = params->GetGSWK();
    uint32_t FP_K = params->GetRLWEParams()->GetK();
    uint32_t FP_N = params->GetRLWEParams()->GetN();


    uint32_t KS_base_bit = params->GetBaseKSBit();
    uint32_t KS_len = params->GetLenKSsave();
    uint32_t KS_rm = params->GetLenKSrm();
    NativeInteger Q0 = params -> GetModuliQ()[0];
    NativeInteger Q1 = params -> GetModuliQ()[1];
		uint64_t Q0_nat = Q0.ConvertToInt();
    uint64_t Q1_nat = Q1.ConvertToInt();
		
		#ifdef WITH_INTEL_HEXL
    
		#else
    NativeInteger Mu0 = params -> GetRLWEParams()->GetMuFP()[0];
    NativeInteger Mu1 = params -> GetRLWEParams()->GetMuFP()[1];
    #endif

    vector<vector<DCRTPoly>> out_CT(out_nums);
    for (uint32_t i = 0; i < out_nums; i ++) {
      out_CT[i].resize(FP_K+1);
    }

    vector<std::shared_ptr<vector<vector<int32_t>>>> CT_dec2(GSW_K);


    // dim0 : K / dim1: N. dim2: KS_len
    for (uint32_t i = 0; i < GSW_K; i++) {
      CT_dec2[i] = SGDLevel2ForKSPARALLEL(params, inputs[i], KS_base_bit, KS_rm, Format::COEFFICIENT);
    }
    // dim0 : K / dim1: N. dim2: KS_len
    const vector<vector<vector<FPRLWECiphertextImpl>>> &KS = K->GetElements();

    #pragma omp parallel for collapse(1)  
    //{
      for (uint32_t mix_idx = 0; mix_idx < (out_nums * (FP_K+1)); mix_idx++) {
      uint32_t out_idx = mix_idx % out_nums;
      uint32_t out_k = mix_idx / out_nums;
      DCRTPoly &my_out_CT = out_CT[out_idx][out_k];
      //std::cout << "out_idx is " << out_idx << ", out k is " << out_k << std::endl;
      // Setting A
      //for (uint32_t out_tmp = 0; out_tmp < FP_K; out_tmp ++) {
      //  out_CT[out_idx][out_tmp] =  DCRTPoly(params->GetRLWEParams()->GetDCRTParams(), Format::EVALUATION, true);
      //}

      if (out_k == FP_K) {
        // Setting B 
        my_out_CT = DCRTPoly(params->GetRLWEParams()->GetDCRTParams(), Format::COEFFICIENT,true);
        my_out_CT.GetElementW(0)[0] = inputs[GSW_K].GetElementAtIndex(0)[out_idx];
        my_out_CT.GetElementW(1)[0] = inputs[GSW_K].GetElementAtIndex(1)[out_idx];
        my_out_CT.SetFormat(Format::EVALUATION);
      } else {
        my_out_CT = DCRTPoly(params->GetRLWEParams()->GetDCRTParams(), Format::EVALUATION,true);
        // Setting A
      }
      NativePoly &out_CT_i_k_Q0 = my_out_CT.GetElementW(0);
      NativePoly &out_CT_i_k_Q1 = my_out_CT.GetElementW(1);

      uint32_t real_j_idx;
      int32_t a_tmp2;
      NativeInteger Tmp_Nat;
      bool a_tmp2_positive;
      bool s_positive;
      for(uint32_t in_k_idx = 0; in_k_idx < GSW_K; in_k_idx++) {
        const vector<vector<int32_t>> &CT_dec2_k =  (*CT_dec2[in_k_idx]);
        const vector<vector<FPRLWECiphertextImpl>> &KS_k  = KS[in_k_idx];
        for (uint32_t j = 0; j < GSW_N; j++){
          if (j <= out_idx) {
            real_j_idx = out_idx - j;
            s_positive = true;
          } else {
            real_j_idx = GSW_N + out_idx - j;
            s_positive = false;
          }
          const vector<FPRLWECiphertextImpl> &KS_k_j = KS_k[real_j_idx];
          const vector<int32_t> &CT_dec2_k_j = CT_dec2_k[j];

          for (uint32_t i = 0; i < KS_len; i++) {

            //std::cout << "j is " << j << ", out idx is " << out_idx << "out_k is " << out_k 
            //  <<" KS ien is " << i  <<std::endl;

            // IDX adjust
            a_tmp2 = CT_dec2_k_j[i];
            if (a_tmp2 == 0) {
              continue;
            } else if (a_tmp2 > 0 ) {
              a_tmp2_positive = true;
            } else {
              a_tmp2_positive = false;
              a_tmp2 *= -1;
            }
            //NativePoly &KS_A_k_Q0; 
            const NativePoly &KS_A_k_Q0 = out_k == FP_K ?
              (KS_k_j[i].GetB().GetElementAtIndex(0)) :
              (KS_k_j[i].GetA()[out_k].GetElementAtIndex(0));

            const NativePoly &KS_A_k_Q1 = out_k == FP_K ?
              (KS_k_j[i].GetB().GetElementAtIndex(1)) :
              (KS_k_j[i].GetA()[out_k].GetElementAtIndex(1)) ;


            /*
            std::cout << "out_k is"  << out_k << ", FP_K is" << FP_K 
              << ", KSk_j[i].GetB is" << KS_k_j[i].GetB() << std::endl;
            std::cout << "Made" << std::endl;
            
            std::cout << KS_A_k_Q0 << std::endl;
            */
            //Calculation!!
						uint64_t a_tmp2_Q0;
						uint64_t a_tmp2_Q1;
            if (s_positive ^ a_tmp2_positive) {
							// - a * (s[real_j]) or a * (-s[real_j]), Just Adding
           
							a_tmp2_Q0 = a_tmp2;
							a_tmp2_Q1 = a_tmp2;
						} else {
							a_tmp2_Q0 = Q0_nat - a_tmp2;
							a_tmp2_Q1 = Q1_nat - a_tmp2;
						}
							
						
						#ifdef WITH_INTEL_HEXL
						uint64_t *op1_Q0 = reinterpret_cast<uint64_t  *> (&out_CT_i_k_Q0[0]);
						//uint64_t op2		 = (uint64_t) (a_tmp2);
						const uint64_t *op3_Q0 = (uint64_t *) (&KS_A_k_Q0[0]);
						intel::hexl::EltwiseFMAMod(op1_Q0, op3_Q0, a_tmp2_Q0
							, op1_Q0, FP_N, Q0_nat,1);

						uint64_t *op1_Q1 = reinterpret_cast<uint64_t  *> (&out_CT_i_k_Q1[0]);
						const uint64_t *op3_Q1 = (uint64_t *) (&KS_A_k_Q1[0]);
						intel::hexl::EltwiseFMAMod(op1_Q1, op3_Q1, a_tmp2_Q1
							, op1_Q1, FP_N, Q1_nat,1);
						#else
						
           	for (uint32_t out_n = 0; out_n < FP_N; out_n++) {
							// A
              //out_CT[out_idx][out_k] += (KS[in_k_idx][real_j_idx][i].GetA()[out_k] * mul_a) ;
              /*
              Tmp_Nat = KS[in_k_idx][real_j_idx][i].GetA()[out_k].GetElementAtIndex(0)[out_n];    
              Tmp_Nat.MulEq(a_tmp2);
              Tmp_Nat.ModEq(Q0);
              out_CT[out_idx][out_k].GetElementW(0)[out_n].AddEq(Tmp_Nat);
                    
              Tmp_Nat = KS[in_k_idx][real_j_idx][i].GetA()[out_k].GetElementAtIndex(1)[out_n];    
              Tmp_Nat.MulEq(a_tmp2)
              Tmp_Nat.ModEq(Q1);        
              out_CT[out_idx][out_k].GetElementW(1)[out_n].AddEq(Tmp_Nat);
              */

              Tmp_Nat = KS_A_k_Q0[out_n];
              //Tmp_Nat.MulEq(a_tmp2);
	            //if (a_tmp2_Q0 < 2048) {
								//Tmp_Nat.MulEq(a_tmp2_Q0);
								//Tmp_Nat.ModEq(Q0, Mu0);
							//} else {
								Tmp_Nat.ModMulFastEq(a_tmp2_Q0, Q0, Mu0);
								//Tmp_Nat.ModMulEq(a_tmp2_Q0, Q0);
							//}
							//out_CT_i_k_Q0[out_n].AddEq(Tmp_Nat);
							out_CT_i_k_Q0[out_n].ModAddFastEq(Tmp_Nat,Q0);
         			//out_CT_i_k_Q0[out_n].ModAddEq(Tmp_Nat,Q0);
         			
							Tmp_Nat = KS_A_k_Q1[out_n];
					    //if (a_tmp2_Q0 < 2048) {
								//Tmp_Nat.MulEq(a_tmp2_Q1);
								//Tmp_Nat.ModEq(Q1, Mu1);
							//} else {
								Tmp_Nat.ModMulFastEq(a_tmp2_Q1, Q1, Mu1);
								//Tmp_Nat.ModMulEq(a_tmp2_Q1, Q1);
							//}
							//out_CT_i_k_Q1[out_n].AddEq(Tmp_Nat);
							out_CT_i_k_Q1[out_n].ModAddFastEq(Tmp_Nat,Q1);
          		//out_CT_i_k_Q1[out_n].ModAddEq(Tmp_Nat,Q1);
            }
						#endif
            
					} // KS ending
          // Refresing
					#ifdef WITH_INTEL_HEXL
					#else
					/*
					for (uint32_t out_n = 0; out_n < FP_N; out_n++) {
							out_CT_i_k_Q0[out_n].ModEq(Q0);
							out_CT_i_k_Q1[out_n].ModEq(Q1);
					}
					*/
					#endif

        } // in GSW N is end
      }// in GSW K is end

      // Rotate Index
      if (orign_idx != 0) {
        my_out_CT *= params->GetRotMonomial_RLWE(orign_idx);
      }
    }//}
    std::chrono::duration<double> end = std::chrono::system_clock::now() - start;
    if (DEBUG) std::cout << "Ring Key Switch PARALLEL time is " << end.count() << std::endl;
    return std::make_shared<vector<vector<DCRTPoly>>> (out_CT);
  }



	/*
	std::shared_ptr<vector<vector<DCRTPoly>>> FPRingGSWAccumulatorScheme::RingKeySwitchPolyPARALLEL(
    const std::shared_ptr<FPRingGSWCryptoParams> params,
    const std::shared_ptr<FPRLWESwitchingKey> K,
    const vector<DCRTPoly> &inputs,
    const uint32_t &orign_idx,
    const uint32_t &out_nums
    ) const {
		auto start = std::chrono::system_clock::now();
	
    //CT should be COEFFICIENT
    if (inputs[0].GetFormat() != Format::COEFFICIENT) {
      PALISADE_THROW(config_error, "RKS input should be COEFF Mode!!");
    } else if ( inputs[inputs.size()-1].GetFormat() != Format::COEFFICIENT) {
			PALISADE_THROW(config_error, "RKS input should be COEFF Mode!!");
		}

    uint32_t GSW_N = params->GetGSWN();
    uint32_t GSW_K = params->GetGSWK();
    uint32_t FP_K = params->GetRLWEParams()->GetK();
		uint32_t FP_N = params->GetRLWEParams()->GetN();


    uint32_t KS_base_bit = params->GetBaseKSBit();
    uint32_t KS_len = params->GetLenKSsave();
    uint32_t KS_rm = params->GetLenKSrm();
		NativeInteger Q0 = params -> GetModuliQ()[0];
		NativeInteger Q1 = params -> GetModuliQ()[1];
		
		//int64_t Q0_over_2 = Q0.ConvertToInt() >> 1;
    //int64_t a_tmp;
		//const uint32_t Num_TH = 4;
		// Make Out
    //int32_t a_tmp2;

    vector<vector<DCRTPoly>> out_CT(out_nums);
    for (uint32_t i = 0; i < out_nums; i ++) {
      out_CT[i].resize(FP_K+1);
    }

    vector<std::shared_ptr<vector<vector<int32_t>>>> CT_dec2(GSW_K);
		std::shared_ptr<vector<vector<int32_t>>> tmp_vals;
		for (uint32_t i = 0; i < GSW_K; i++) {

			CT_dec2[i] = SGDLevel2ForKSPARALLEL(params, inputs[i], KS_base_bit, KS_rm, Format::COEFFICIENT);
		}
		
     // dim0 : K / dim1: N. dim2: KS_len
    vector<vector<vector<FPRLWECiphertextImpl>>> KS = K->GetElements();
		
		#pragma omp parallel for collapse(1)	
		{for (uint32_t mix_idx = 0; mix_idx < (out_nums * (FP_K+1)); mix_idx++) {
			uint32_t out_idx = mix_idx % out_nums;
			uint32_t out_k = mix_idx / out_nums;
			//std::cout << "out_idx is " << out_idx << ", out k is " << out_k << std::endl;
			// Setting A
    	//for (uint32_t out_tmp = 0; out_tmp < FP_K; out_tmp ++) {
			//	out_CT[out_idx][out_tmp] =  DCRTPoly(params->GetRLWEParams()->GetDCRTParams(), Format::EVALUATION, true);
      //}
			if (out_k == FP_K) { 
				// Setting B 
				out_CT[out_idx][FP_K] = DCRTPoly(params->GetRLWEParams()->GetDCRTParams(), Format::COEFFICIENT,true);
				out_CT[out_idx][FP_K].GetElementW(0)[0] = inputs[GSW_K].GetElementAtIndex(0)[out_idx];
				out_CT[out_idx][FP_K].GetElementW(1)[0] = inputs[GSW_K].GetElementAtIndex(1)[out_idx];
				out_CT[out_idx][FP_K].SetFormat(Format::EVALUATION);
			} else {
				// Setting A
    		out_CT[out_idx][out_k] =  DCRTPoly(params->GetRLWEParams()->GetDCRTParams(), Format::EVALUATION, true);
			}


     	//#pragma omp parallel for num_threads(Num_TH) 
			//#pragma omp parallel for collapse(1)
			//for (uint32_t out_k = 0; out_k < FP_K+1; out_k++) {
			
			uint32_t real_j_idx;
			int32_t a_tmp2;
			NativeInteger Tmp_Nat;
			bool a_tmp2_positive;
			bool s_positive;
					
			for (uint32_t i = 0; i < KS_len; i++) {
				// Makes
				for (uint32_t j = 0; j < GSW_N; j++){
					// IDX adjust
					for (uint32_t in_k_idx = 0; in_k_idx < GSW_K; in_k_idx++) {
	
						a_tmp2 = (*CT_dec2[in_k_idx])[i][j];
						if (a_tmp2 == 0) {
							continue;
						} else if (a_tmp2 > 0 ) {
							a_tmp2_positive = true;
						} else {
							a_tmp2_positive = false;
							a_tmp2 *= -1;
						}
						if (j <= out_idx) {
								real_j_idx = out_idx - j;
							s_positive = true;
						} else {
							real_j_idx = GSW_N + out_idx - j;
							s_positive = false;
						}
						//Calculation!!
						if (s_positive ^ a_tmp2_positive) {
							// - a * (s[real_j]) or a * (-s[real_j]), Just Adding
							for (uint32_t out_n = 0; out_n < FP_N; out_n++) {
								// A
								if (out_k != FP_K) {
								//out_CT[out_idx][out_k] += (KS[in_k_idx][real_j_idx][i].GetA()[out_k] * mul_a) ;
									Tmp_Nat = KS[in_k_idx][real_j_idx][i].GetA()[out_k].GetElementAtIndex(0)[out_n];    
									Tmp_Nat.MulEq(a_tmp2);
									Tmp_Nat.ModEq(Q0);
									out_CT[out_idx][out_k].GetElementW(0)[out_n].AddEq(Tmp_Nat);
									
									Tmp_Nat = KS[in_k_idx][real_j_idx][i].GetA()[out_k].GetElementAtIndex(1)[out_n];    
									Tmp_Nat.MulEq(a_tmp2);
									Tmp_Nat.ModEq(Q1);				
									out_CT[out_idx][out_k].GetElementW(1)[out_n].AddEq(Tmp_Nat);
							
								} else {
									// B sum
									//out_CT[out_idx][FP_K] += (KS[in_k_idx][real_j_idx][i].GetB() * mul_a);
								
									Tmp_Nat = KS[in_k_idx][real_j_idx][i].GetB().GetElementAtIndex(0)[out_n];    
									Tmp_Nat.MulEq(a_tmp2);
									Tmp_Nat.ModEq(Q0);
									out_CT[out_idx][FP_K].GetElementW(0)[out_n].AddEq(Tmp_Nat);
								
									Tmp_Nat = KS[in_k_idx][real_j_idx][i].GetB().GetElementAtIndex(1)[out_n];    
									Tmp_Nat.MulEq(a_tmp2);
									Tmp_Nat.ModEq(Q1);
									out_CT[out_idx][FP_K].GetElementW(1)[out_n].AddEq(Tmp_Nat);
					
								}
							}
						} else {
							// a * (s[real_j]) or - a * (-s[real_j]), Reverse is needed
							for (uint32_t out_n = 0; out_n < FP_N; out_n++) {
								//out_CT[out_idx][out_k] -= (KS[in_k_idx][real_j_idx][i].GetA()[out_k] * mul_a) ;
								if (out_k != FP_K) {
									Tmp_Nat = Q0;
									Tmp_Nat.SubEq(KS[in_k_idx][real_j_idx][i].GetA()[out_k].GetElementAtIndex(0)[out_n]);    
									Tmp_Nat.MulEq(a_tmp2);
									Tmp_Nat.ModEq(Q0);
									out_CT[out_idx][out_k].GetElementW(0)[out_n].AddEq(Tmp_Nat);
										
									Tmp_Nat = Q1;
									Tmp_Nat.SubEq(KS[in_k_idx][real_j_idx][i].GetA()[out_k].GetElementAtIndex(1)[out_n]);    
									Tmp_Nat.MulEq(a_tmp2);
									Tmp_Nat.ModEq(Q1);				
									out_CT[out_idx][out_k].GetElementW(1)[out_n].AddEq(Tmp_Nat);
						
								} else {
									// B sum
									Tmp_Nat = Q0;
									Tmp_Nat.SubEq(KS[in_k_idx][real_j_idx][i].GetB().GetElementAtIndex(0)[out_n]);    
									Tmp_Nat.MulEq(a_tmp2);
									Tmp_Nat.ModEq(Q0);
									out_CT[out_idx][FP_K].GetElementW(0)[out_n].AddEq(Tmp_Nat);
								
									Tmp_Nat = Q1;
									Tmp_Nat.SubEq(KS[in_k_idx][real_j_idx][i].GetB().GetElementAtIndex(1)[out_n]);    
									Tmp_Nat.MulEq(a_tmp2);
									Tmp_Nat.ModEq(Q1);
									out_CT[out_idx][FP_K].GetElementW(1)[out_n].AddEq(Tmp_Nat);
				
								}
							}
						} // End Calculation 
					} // End GSW_K
				}  // End GSW_N
			} // End KS Len

			for (uint32_t out_n = 0; out_n < FP_N; out_n++) {	
				if (out_k != FP_K) {
					out_CT[out_idx][out_k].GetElementW(0)[out_n].ModEq(Q0);
					out_CT[out_idx][out_k].GetElementW(1)[out_n].ModEq(Q1);
				} else {
					out_CT[out_idx][FP_K].GetElementW(0)[out_n].ModEq(Q0);
					out_CT[out_idx][FP_K].GetElementW(1)[out_n].ModEq(Q1);
				}
			} 
      // Rotate Index
      if (orign_idx != 0) {
			//{for (uint32_t out_k = 0; out_k < FP_K+1; out_k++) {
        if (out_k == FP_K) {	// out B sum
					out_CT[out_idx][FP_K] *= params->GetRotMonomial_RLWE(orign_idx);
				} else {							// out A sum		
					out_CT[out_idx][out_k] *= params->GetRotMonomial_RLWE(orign_idx);;
				}
			}


		}}

		std::chrono::duration<double> end = std::chrono::system_clock::now() - start;
		if (DEBUG) std::cout << "Ring Key Switch PARALLEL time is " << end.count() << std::endl;
 

		return std::make_shared<vector<vector<DCRTPoly>>> (out_CT);
	}
*/



	std::shared_ptr<vector<vector<int32_t>>> FPRingGSWAccumulatorScheme::SGDLevel2ForKS(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const DCRTPoly &in_poly, const uint32_t &bits, 
		const uint32_t &rm_len, const Format format=Format::COEFFICIENT) const {

		// NTT Coeff Check
		if (in_poly.GetFormat() != Format::COEFFICIENT) {
			PALISADE_THROW(config_error, "SGDLevel input should be Coeff Mode!");
    }
		// LEvel is 2

		const uint32_t Ns =  in_poly.GetLength();
    const uint32_t Qbit = params->GetQBit();
    uint32_t total_len = Qbit / bits;
    if (Qbit % bits > 0) {total_len++;}
    
		// Out Poly Making
		//vector<vector<int32_t>> out_poly(total_len - rm_len);
  	vector<vector<int32_t>> out_poly(Ns);
    for (uint32_t i = 0; i < Ns; i++) {
			out_poly[i].resize(total_len - rm_len);
		}
    	
		vector<uint64_t> QOver2Info = params->GetRLWEParams()->GetQOver2Info();

    uint64_t	Delta = params->GetRLWEParams()->GetScaling().ConvertToInt();
    uint64_t	Q0_nat = params->GetModuliQ()[0].ConvertToInt();
    uint64_t  Q1_nat = params->GetModuliQ()[1].ConvertToInt();

		NativeInteger Mu1 = params->GetRLWEParams()->GetMuFP()[1];

		NativeInteger Inv_1 = params->GetRLWEParams()->GetInv()[1];

    uint32_t lower_len = (params->GetRLWEParams()->GetBits()[0]) / bits;
    uint32_t upper_len = total_len - lower_len;
    uint32_t point_bit = bits * lower_len;
    uint32_t shift_len = params->GetRLWEParams()->GetBits()[0] - point_bit;
		uint32_t Base = 1;
						 Base <<= bits;  // bits should be lower than 31

		uint64_t p_mask = 1;
		p_mask <<=  point_bit;
		p_mask -= 1;
    uint64_t p_mask_com = ~p_mask;
		
		//if (NativeInteger::MaxBits)
    NativeInteger::SignedNativeInt gBits =  (NativeInteger::SignedNativeInt) bits ;
    // VARIANT A
    /*
		NativeInteger::SignedNativeInt gBitsMaxBits = NativeInteger::MaxBits() - gBits;
    // Test , Should be Fixed
    
		uint32_t test_bit1 = (uint32_t) std::ceil(log(static_cast<double> 
          ( ((double)Q0_nat) +  ((double)(Q1_nat+ 1)) * ((double)Delta)) ) / log(static_cast<double> (2)));

    if ( (test_bit1 + shift_len) > gBitsMaxBits) {
      PALISADE_THROW(config_error, "SGD is unstable");
    } 
		*/
    // VARIANT B
    
		NativeInteger::SignedNativeInt baseG = NativeInteger::SignedNativeInt(Base);
    NativeInteger::SignedNativeInt gminus1 = Base - 1;
    NativeInteger::SignedNativeInt baseGdiv2 = (Base >> 1)-1;
    

		uint64_t t_down = 0;
    uint64_t  t_up  = 0;
    NativeInteger::SignedNativeInt d_up = 0;
    NativeInteger::SignedNativeInt d_down = 0;
		
		/*
		#ifdef WITH_INTEL_HEXL
		const uint64_t *t_down_vec	 = (uint64_t *)(&in_poly.GetElementAtIndex(0)[0]);
		const uint64_t *t_up_vec     =  (uint64_t *)(&in_poly.GetElementAtIndex(1)[0]);
		uint64_t t_res_vec[Ns];
		intel::hexl::EltwiseReduceMod(t_res_vec, t_down_vec, Ns, Q1_nat, 0,1);
		intel::hexl::EltwiseSubMod(t_res_vec, t_up_vec, t_res_vec, Ns, Q1_nat);
		intel::hexl::EltwiseFMAMod(t_res_vec, t_res_vec,  Inv_1.ConvertToInt(), nullptr , Ns, Q1_nat, 1);
		#endif 
		*/

    for (uint32_t k = 0; k < Ns; k++) {

			t_down = (in_poly.GetElementAtIndex(0)[k]).ConvertToInt(); // 0 ~  Q0-1
			//#ifdef WITH_INTEL_HEXL
			//t_up = t_res_vec[k];
			//#else
			t_up   = (in_poly.GetElementAtIndex(1)[k]).ModSub(in_poly.GetElementAtIndex(0)[k], Q1_nat).ModMulFast(Inv_1, Q1_nat, Mu1).ConvertToInt(); // 0 ~  Q0-1
			//#endif
			// Test If it is Minus
      bool rev = false;
      if ((t_up > QOver2Info[1]) || (  (t_up == QOver2Info[1])   && (t_down > QOver2Info[0])) ) {
				rev = true;

        if ( t_down == 0 ) {
					t_up = Q1_nat - t_up;
        } else {
          t_up    = Q1_nat - t_up - 1;
          t_down  = Q0_nat  - t_down;
        }

      }
      t_down += t_up;
      t_up *= Delta;

      if (rev) {
				d_down  -= (NativeInteger::SignedNativeInt) (t_down & p_mask);
        d_up    -= (NativeInteger::SignedNativeInt) ((t_up << shift_len) + ((t_down & p_mask_com) >> point_bit));
      } else {
        d_down  += (NativeInteger::SignedNativeInt) (t_down & p_mask);
				d_up    += (NativeInteger::SignedNativeInt) ((t_up << shift_len) + ((t_down & p_mask_com) >> point_bit));
      }

      for (uint32_t lo = 0; lo < lower_len; lo++) {

				// This approach gives a slightly better performance
        // VARIANT A
        //NativeInteger::SignedNativeInt r = d_down << gBitsMaxBits;
        //r >>= gBitsMaxBits;
        // VARIANT B
        NativeInteger::SignedNativeInt r = d_down & gminus1;  
        if (r > baseGdiv2) r -= baseG;

        d_down -= r;
        d_down >>= gBits;

        if (lo < rm_len) {
					continue;
        }
        // It should be zeros in values, or It has some intension
				out_poly[k][lo-rm_len] = (int32_t) r;
			}
      d_up += d_down;
        for (uint32_t hi = 0; hi < upper_len; hi++) {
          // This approach gives a slightly better performance
          // VARIANT A
          //NativeInteger::SignedNativeInt r = d_up << gBitsMaxBits;
          //r >>= gBitsMaxBits;

          // VARIANT B
          NativeInteger::SignedNativeInt r = d_up & gminus1;
          if (r > baseGdiv2) r -= baseG;

          d_up -= r;
          d_up >>= gBits;
          if (hi+lower_len < rm_len) {
            continue;
          }
					out_poly[k][hi+lower_len-rm_len] = (int32_t) r;
    
        }
        d_up = 0;
        d_down = 0;
      }
			return std::make_shared<vector<vector<int32_t>>> (out_poly);
    }










	std::shared_ptr<vector<DCRTPoly>> FPRingGSWAccumulatorScheme::SGDLevel2(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const DCRTPoly &in_poly, const uint32_t &bits, 
		const uint32_t &rm_len, const Format format=Format::COEFFICIENT) const {

		// NTT Coeff Check
		if (in_poly.GetFormat() != Format::COEFFICIENT) {
			PALISADE_THROW(config_error, "SGDLevel input should be Coeff Mode!");
    }
			
		const uint32_t Ns =  in_poly.GetLength();
		NativeInteger Q0 = params->GetModuliQ()[0];
		NativeInteger Q1 = params->GetModuliQ()[1];
		// We assme GSW_Q = FP_Q
		NativeInteger Mu1 = params->GetMuGSW()[1];

    const uint32_t Qbit = params->GetQBit();
    uint32_t total_len = Qbit / bits;
    if (Qbit % bits > 0) {total_len++;}
    
		// Out Poly Making
		vector<DCRTPoly> out_poly(total_len - rm_len);
		vector<NativeVector> save_poly_Q0(total_len - rm_len); 
  	vector<NativeVector> save_poly_Q1(total_len - rm_len); 
    for (uint32_t i = 0; i < (total_len - rm_len); i++) {
			//#ifdef WITH_INTEL_HEXL
			out_poly[i] = DCRTPoly(in_poly.GetParams(), Format::COEFFICIENT);
			//#else
			//out_poly[i] = DCRTPoly(in_poly.GetParams(), Format::COEFFICIENT, true);
			//#endif
			save_poly_Q0[i] = NativeVector(Ns, Q0);
			save_poly_Q1[i] = NativeVector(Ns, Q1);
		}
    
	
		//const uint32_t Ns =  in_poly.GetLength();
    vector<uint64_t> QOver2Info = params->GetRLWEParams()->GetQOver2Info();

    uint64_t	Delta = params->GetRLWEParams()->GetScaling().ConvertToInt();
    uint64_t	Q0_nat = Q0.ConvertToInt();
    uint64_t  Q1_nat = Q1.ConvertToInt();
    NativeInteger::SignedNativeInt Q0_sint = Q0_nat;
    NativeInteger::SignedNativeInt Q1_sint = Q1_nat;
    NativeInteger Inv_1 = params->GetRLWEParams()->GetInv()[1];

    uint32_t lower_len = (params->GetRLWEParams()->GetBits()[0]) / bits;
    uint32_t upper_len = total_len - lower_len;
    uint32_t point_bit = bits * lower_len;
    uint32_t shift_len = params->GetRLWEParams()->GetBits()[0] - point_bit;
		uint32_t Base = 1;
						 Base <<= bits;  // bits should be lower than 31

		uint64_t p_mask = 1;
		p_mask <<=  point_bit;
		p_mask -= 1;
    uint64_t p_mask_com = ~p_mask;
		
		//if (NativeInteger::MaxBits)
    NativeInteger::SignedNativeInt gBits =  (NativeInteger::SignedNativeInt) bits ;
    // VARIANT A
    /*
		NativeInteger::SignedNativeInt gBitsMaxBits = NativeInteger::MaxBits() - gBits;
    // Test , Should be Fixed
    
		uint32_t test_bit1 = (uint32_t) std::ceil(log(static_cast<double> 
          ( ((double)Q0_nat) +  ((double)(Q1_nat+ 1)) * ((double)Delta)) ) / log(static_cast<double> (2)));

    if ( (test_bit1 + shift_len) > gBitsMaxBits) {
      PALISADE_THROW(config_error, "SGD is unstable");
    } 
		*/
    // VARIANT B
    
		NativeInteger::SignedNativeInt baseG = NativeInteger::SignedNativeInt(Base);
    NativeInteger::SignedNativeInt gminus1 = Base - 1;
    NativeInteger::SignedNativeInt baseGdiv2 = (Base >> 1)-1;
    

		uint64_t t_down = 0;
    uint64_t  t_up  = 0;
    NativeInteger::SignedNativeInt d_up = 0;
    NativeInteger::SignedNativeInt d_down = 0;

	
		/*
		#ifdef WITH_INTEL_HEXL
		const uint64_t *t_down_vec	 = (uint64_t *)(&in_poly.GetElementAtIndex(0)[0]);
		const uint64_t *t_up_vec     =  (uint64_t *)(&in_poly.GetElementAtIndex(1)[0]);
		uint64_t t_res_vec[Ns];
		intel::hexl::EltwiseReduceMod(t_res_vec, t_down_vec, Ns, Q1_nat, 0,1);
		intel::hexl::EltwiseSubMod(t_res_vec, t_up_vec, t_res_vec, Ns, Q1_nat);
		intel::hexl::EltwiseFMAMod(t_res_vec, t_res_vec,  Inv_1.ConvertToInt(), nullptr , Ns, Q1_nat, 1);
		#endif 
		*/


    for (uint32_t k = 0; k < Ns; k++) {

			t_down = (in_poly.GetElementAtIndex(0)[k]).ConvertToInt(); // 0 ~  Q0-1
			//t_up   = (in_poly.GetElementAtIndex(1)[k]).ModSubFast(in_poly.GetElementAtIndex(0)[k], Q1_nat).ModMulFast(Inv_1, Q1_nat, Mu1).ConvertToInt(); // 0 ~  Q0-1
			//#ifdef WITH_INTEL_HEXL
			//t_up = t_res_vec[k];
			//#else
			t_up   = (in_poly.GetElementAtIndex(1)[k]).ModSub(in_poly.GetElementAtIndex(0)[k], Q1_nat).ModMulFast(Inv_1, Q1_nat, Mu1).ConvertToInt(); // 0 ~  Q0-1
		
			//#endif
			// Test If it is Minus
      bool rev = false;
      if ((t_up > QOver2Info[1]) || (  (t_up == QOver2Info[1])   && (t_down > QOver2Info[0])) ) {
				rev = true;

        if ( t_down == 0 ) {
					t_up = Q1_nat - t_up;
        } else {
          t_up    = Q1_nat - t_up - 1;
          t_down  = Q0_nat  - t_down;
        }

      }
      t_down += t_up;
      t_up *= Delta;

      if (rev) {
				d_down  -= (NativeInteger::SignedNativeInt) (t_down & p_mask);
        d_up    -= (NativeInteger::SignedNativeInt) ((t_up << shift_len) + ((t_down & p_mask_com) >> point_bit));
      } else {
        d_down  += (NativeInteger::SignedNativeInt) (t_down & p_mask);
				d_up    += (NativeInteger::SignedNativeInt) ((t_up << shift_len) + ((t_down & p_mask_com) >> point_bit));
      }

      for (uint32_t lo = 0; lo < lower_len; lo++) {

				// This approach gives a slightly better performance
        // VARIANT A
        //NativeInteger::SignedNativeInt r = d_down << gBitsMaxBits;
        //r >>= gBitsMaxBits;
        // VARIANT B
        NativeInteger::SignedNativeInt r = d_down & gminus1;  
        if (r > baseGdiv2) r -= baseG;

        d_down -= r;
        d_down >>= gBits;

        if (lo < rm_len) {
					continue;
        }
        // It should be zeros in values, or It has some intension
        if (r >= 0) {
          //out_poly[lo-rm_len].GetElementW(0)[k] = NativeInteger(r);
          //out_poly[lo-rm_len].GetElementW(1)[k] = NativeInteger(r);
          save_poly_Q0[lo-rm_len][k] = NativeInteger(r);
          save_poly_Q1[lo-rm_len][k] = NativeInteger(r);
				} else {
					//out_poly[lo-rm_len].GetElementW(0)[k] = NativeInteger(Q0_sint + r);
          //out_poly[lo-rm_len].GetElementW(1)[k] = NativeInteger(Q1_sint + r);
          save_poly_Q0[lo-rm_len][k] = NativeInteger(Q0_sint+r);
          save_poly_Q1[lo-rm_len][k] = NativeInteger(Q1_sint+r);
				
				}
      }
      d_up += d_down;
        for (uint32_t hi = 0; hi < upper_len; hi++) {
          // This approach gives a slightly better performance
          // VARIANT A
          //NativeInteger::SignedNativeInt r = d_up << gBitsMaxBits;
          //r >>= gBitsMaxBits;

          // VARIANT B
          NativeInteger::SignedNativeInt r = d_up & gminus1;
          if (r > baseGdiv2) r -= baseG;

          d_up -= r;
          d_up >>= gBits;
          if (hi+lower_len < rm_len) {
            continue;
          }
          if (r >= 0) {
            //out_poly[hi+lower_len-rm_len].GetElementW(0)[k] = NativeInteger(r);
            //out_poly[hi+lower_len-rm_len].GetElementW(1)[k] = NativeInteger(r);
            save_poly_Q0[hi+lower_len-rm_len][k] = NativeInteger(r);
            save_poly_Q1[hi+lower_len-rm_len][k] = NativeInteger(r);
             
					} else {
						//out_poly[hi+lower_len-rm_len].GetElementW(0)[k] = NativeInteger(Q0_sint+ r);
            //out_poly[hi+lower_len-rm_len].GetElementW(1)[k] = NativeInteger(Q1_sint+ r);
            save_poly_Q0[hi+lower_len-rm_len][k] = NativeInteger(Q0_sint+r);
						save_poly_Q1[hi+lower_len-rm_len][k] = NativeInteger(Q1_sint+r);
					}
        }
        d_up = 0;
        d_down = 0;
      }
			// Move
			
			NativePoly Tmp_Q0; 
			NativePoly Tmp_Q1; 
			for (uint32_t i = 0; i < (total_len - rm_len); i++) {
				//#ifdef WITH_INTEL_HEXL
				Tmp_Q0 = NativePoly(in_poly.GetParams()->GetParams()[0],Format::COEFFICIENT);
				Tmp_Q1 = NativePoly(in_poly.GetParams()->GetParams()[1],Format::COEFFICIENT);
				//#else
				//Tmp_Q0 = NativePoly(in_poly.GetParams()->GetParams()[0],Format::COEFFICIENT,true);
				//Tmp_Q1 = NativePoly(in_poly.GetParams()->GetParams()[1],Format::COEFFICIENT,true);
				//#endif
				Tmp_Q0.SetValues(std::move(save_poly_Q0[i]), Format::COEFFICIENT);
				Tmp_Q1.SetValues(std::move(save_poly_Q1[i]), Format::COEFFICIENT);
				out_poly[i].SetElementAtIndex(0,std::move(Tmp_Q0));
				out_poly[i].SetElementAtIndex(1,std::move(Tmp_Q1));
			}
			
			
      // Format change
      for (uint32_t i  = 0; i < out_poly.size(); i++ ) {
        out_poly[i].SetFormat(format);
      }
			return std::make_shared<vector<DCRTPoly>> (out_poly);
    }


	std::shared_ptr<vector<DCRTPoly>> FPRingGSWAccumulatorScheme::SGDLevel2PARALLEL(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const DCRTPoly &in_poly, const uint32_t &bits, 
		const uint32_t &rm_len, const Format format=Format::COEFFICIENT) const {

		// NTT Coeff Check
		if (in_poly.GetFormat() != Format::COEFFICIENT) {
			PALISADE_THROW(config_error, "SGDLevel input should be Coeff Mode!");
    }
			
    const uint32_t Qbit = params->GetQBit();
    uint32_t total_len = Qbit / bits;
    if (Qbit % bits > 0) {total_len++;}
    const uint32_t Ns =  in_poly.GetLength();
    NativeInteger Q0 = params->GetModuliQ()[0];
		NativeInteger Q1 = params->GetModuliQ()[1];
		NativeInteger Mu1 = params->GetMuGSW()[1];

		uint64_t	Q0_nat = Q0.ConvertToInt();
    uint64_t  Q1_nat = Q1.ConvertToInt();
  
		// Out Poly Making
		vector<DCRTPoly> out_poly(total_len - rm_len);
		vector<NativeVector> save_poly_Q0(total_len - rm_len); 
  	vector<NativeVector> save_poly_Q1(total_len - rm_len); 
    for (uint32_t i = 0; i < (total_len - rm_len); i++) {
			//#ifdef WITH_INTEL_HEXL
			out_poly[i] = DCRTPoly(in_poly.GetParams(), Format::COEFFICIENT);
			//#else
			//out_poly[i] = DCRTPoly(in_poly.GetParams(), Format::COEFFICIENT,true);
			//#endif
			save_poly_Q0[i] = NativeVector(Ns, Q0);
			save_poly_Q1[i] = NativeVector(Ns, Q1);
		}  
	
		vector<uint64_t> QOver2Info = params->GetRLWEParams()->GetQOver2Info();

    uint64_t	Delta = params->GetRLWEParams()->GetScaling().ConvertToInt();
   
		// GSW, FP should be equal
		
		NativeInteger::SignedNativeInt Q0_sint = Q0_nat;
    NativeInteger::SignedNativeInt Q1_sint = Q1_nat;
    NativeInteger Inv_1 = params->GetRLWEParams()->GetInv()[1];

    uint32_t lower_len = (params->GetRLWEParams()->GetBits()[0]) / bits;
    uint32_t upper_len = total_len - lower_len;
    uint32_t point_bit = bits * lower_len;
    uint32_t shift_len = params->GetRLWEParams()->GetBits()[0] - point_bit;
		uint32_t Base = 1;
						 Base <<= bits;  // bits should be lower than 31

		uint64_t p_mask = 1;
		p_mask <<=  point_bit;
		p_mask -= 1;
    uint64_t p_mask_com = ~p_mask;
		
		//if (NativeInteger::MaxBits)
    NativeInteger::SignedNativeInt gBits =  (NativeInteger::SignedNativeInt) bits ;
    // VARIANT A
    /*
		NativeInteger::SignedNativeInt gBitsMaxBits = NativeInteger::MaxBits() - gBits;
    // Test , Should be Fixed
    
		uint32_t test_bit1 = (uint32_t) std::ceil(log(static_cast<double> 
          ( ((double)Q0_nat) +  ((double)(Q1_nat+ 1)) * ((double)Delta)) ) / log(static_cast<double> (2)));

    if ( (test_bit1 + shift_len) > gBitsMaxBits) {
      PALISADE_THROW(config_error, "SGD is unstable");
    } 
		*/
    // VARIANT B
    
		NativeInteger::SignedNativeInt baseG = NativeInteger::SignedNativeInt(Base);
    NativeInteger::SignedNativeInt gminus1 = Base - 1;
    NativeInteger::SignedNativeInt baseGdiv2 = (Base >> 1)-1;
    

		//uint64_t t_down = 0;
    //uint64_t  t_up  = 0;
    //NativeInteger::SignedNativeInt d_up = 0;
    //NativeInteger::SignedNativeInt d_down = 0;
		
		//#pragma omp parallel for
		#pragma omp parallel for
		//{
    for (uint32_t k = 0; k < Ns; k++) {
			
			NativeInteger::SignedNativeInt d_up = 0;
			NativeInteger::SignedNativeInt d_down = 0;
		
			uint64_t t_down = (in_poly.GetElementAtIndex(0)[k]).ConvertToInt(); // 0 ~  Q0-1
			uint64_t t_up   = (in_poly.GetElementAtIndex(1)[k]).ModSub(in_poly.GetElementAtIndex(0)[k], Q1_nat).ModMulFast(Inv_1, Q1_nat, Mu1).ConvertToInt(); // 0 ~  Q0-1

			// Test If it is Minus
      bool rev = false;
      if ((t_up > QOver2Info[1]) || (  (t_up == QOver2Info[1])   && (t_down > QOver2Info[0])) ) {
				rev = true;

        if ( t_down == 0 ) {
					t_up = Q1_nat - t_up;
        } else {
          t_up    = Q1_nat - t_up - 1;
          t_down  = Q0_nat  - t_down;
        }

      }
      t_down += t_up;
      t_up *= Delta;

      if (rev) {
				d_down  -= (NativeInteger::SignedNativeInt) (t_down & p_mask);
        d_up    -= (NativeInteger::SignedNativeInt) ((t_up << shift_len) + ((t_down & p_mask_com) >> point_bit));
      } else {
        d_down  += (NativeInteger::SignedNativeInt) (t_down & p_mask);
				d_up    += (NativeInteger::SignedNativeInt) ((t_up << shift_len) + ((t_down & p_mask_com) >> point_bit));
      }

      for (uint32_t lo = 0; lo < lower_len; lo++) {

				// This approach gives a slightly better performance
        // VARIANT A
        //NativeInteger::SignedNativeInt r = d_down << gBitsMaxBits;
        //r >>= gBitsMaxBits;
        // VARIANT B
        NativeInteger::SignedNativeInt r = d_down & gminus1;  
        if (r > baseGdiv2) r -= baseG;

        d_down -= r;
        d_down >>= gBits;

        if (lo < rm_len) {
					continue;
        }
        // It should be zeros in values, or It has some intension
        if (r >= 0) {
          //out_poly[lo-rm_len].GetElementW(0)[k] = NativeInteger(r);
          //out_poly[lo-rm_len].GetElementW(1)[k] = NativeInteger(r);
          save_poly_Q0[lo-rm_len][k] = NativeInteger(r);
          save_poly_Q1[lo-rm_len][k] = NativeInteger(r);
				} else {
          //out_poly[lo-rm_len].GetElementW(0)[k] = NativeInteger(Q0_sint + r);
          //out_poly[lo-rm_len].GetElementW(1)[k] = NativeInteger(Q1_sint + r);
          save_poly_Q0[lo-rm_len][k] = NativeInteger(Q0_sint + r);
          save_poly_Q1[lo-rm_len][k] = NativeInteger(Q1_sint + r);
				}
      }
      d_up += d_down;
      for (uint32_t hi = 0; hi < upper_len; hi++) {
        // This approach gives a slightly better performance
        // VARIANT A
        //NativeInteger::SignedNativeInt r = d_up << gBitsMaxBits;
        //r >>= gBitsMaxBits;

        // VARIANT B
        NativeInteger::SignedNativeInt r = d_up & gminus1;
        if (r > baseGdiv2) r -= baseG;

        d_up -= r;
        d_up >>= gBits;
        if (hi+lower_len < rm_len) {
					continue;
        }
        if (r >= 0) {
					//out_poly[hi+lower_len-rm_len].GetElementW(0)[k] = NativeInteger(r);
          //out_poly[hi+lower_len-rm_len].GetElementW(1)[k] = NativeInteger(r);
          save_poly_Q0[hi+lower_len-rm_len][k] = NativeInteger(r);
          save_poly_Q1[hi+lower_len-rm_len][k] = NativeInteger(r);
         
				} else {
          //out_poly[hi+lower_len-rm_len].GetElementW(0)[k] = NativeInteger(Q0_sint + r);
          //out_poly[hi+lower_len-rm_len].GetElementW(1)[k] = NativeInteger(Q1_sint + r);
          save_poly_Q0[hi+lower_len-rm_len][k] = NativeInteger(Q0_sint + r);
          save_poly_Q1[hi+lower_len-rm_len][k] = NativeInteger(Q1_sint + r);
				}
      }
      d_up = 0;
      d_down = 0;
    }//}


	
			NativePoly Tmp_Q0; 
			NativePoly Tmp_Q1; 
			for (uint32_t i = 0; i < (total_len - rm_len); i++) {
				//#ifdef WITH_INTEL_HEXL
				Tmp_Q0 = NativePoly(in_poly.GetParams()->GetParams()[0],Format::COEFFICIENT);
				Tmp_Q1 = NativePoly(in_poly.GetParams()->GetParams()[1],Format::COEFFICIENT);
				//#else
				//Tmp_Q0 = NativePoly(in_poly.GetParams()->GetParams()[0],Format::COEFFICIENT,true);
				//Tmp_Q1 = NativePoly(in_poly.GetParams()->GetParams()[1],Format::COEFFICIENT,true);
				//#endif
				Tmp_Q0.SetValues(std::move(save_poly_Q0[i]), Format::COEFFICIENT);
				Tmp_Q1.SetValues(std::move(save_poly_Q1[i]), Format::COEFFICIENT);
				out_poly[i].SetElementAtIndex(0,std::move(Tmp_Q0));
				out_poly[i].SetElementAtIndex(1,std::move(Tmp_Q1));
			}		

    // Format change
    if (format == Format::EVALUATION) { 
			#pragma omp parallel for
			//{
			for (uint32_t i  = 0; i < out_poly.size(); i++ ) {
				out_poly[i].SetFormat(format);
			}//}
		}
		return std::make_shared<vector<DCRTPoly>> (out_poly);
  }




	std::shared_ptr<vector<vector<int32_t>>> FPRingGSWAccumulatorScheme::SGDLevel2ForKSPARALLEL(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const DCRTPoly &in_poly, const uint32_t &bits, 
		const uint32_t &rm_len, const Format format=Format::COEFFICIENT) const {

		// NTT Coeff Check
		if (in_poly.GetFormat() != Format::COEFFICIENT) {
			PALISADE_THROW(config_error, "SGDLevel input should be Coeff Mode!");
    }
		// LEvel is 2

		const uint32_t Ns =  in_poly.GetLength();
    const uint32_t Qbit = params->GetQBit();
    uint32_t total_len = Qbit / bits;
    if (Qbit % bits > 0) {total_len++;}
    
		// Out Poly Making
		//vector<vector<int32_t>> out_poly(total_len - rm_len);
  	vector<vector<int32_t>> out_poly(Ns);
    for (uint32_t i = 0; i < Ns; i++) {
			out_poly[i].resize(total_len - rm_len);
		}
    	
		vector<uint64_t> QOver2Info = params->GetRLWEParams()->GetQOver2Info();

    uint64_t	Delta = params->GetRLWEParams()->GetScaling().ConvertToInt();
    uint64_t	Q0_nat = params->GetModuliQ()[0].ConvertToInt();
    uint64_t  Q1_nat = params->GetModuliQ()[1].ConvertToInt();
		NativeInteger Mu1 = params->GetMuGSW()[1];

		NativeInteger Inv_1 = params->GetRLWEParams()->GetInv()[1];

    uint32_t lower_len = (params->GetRLWEParams()->GetBits()[0]) / bits;
    uint32_t upper_len = total_len - lower_len;
    uint32_t point_bit = bits * lower_len;
    uint32_t shift_len = params->GetRLWEParams()->GetBits()[0] - point_bit;
		uint32_t Base = 1;
						 Base <<= bits;  // bits should be lower than 31

		uint64_t p_mask = 1;
		p_mask <<=  point_bit;
		p_mask -= 1;
    uint64_t p_mask_com = ~p_mask;
		
		//if (NativeInteger::MaxBits)
    NativeInteger::SignedNativeInt gBits =  (NativeInteger::SignedNativeInt) bits ;
    // VARIANT A
    /*
		NativeInteger::SignedNativeInt gBitsMaxBits = NativeInteger::MaxBits() - gBits;
    // Test , Should be Fixed
    
		uint32_t test_bit1 = (uint32_t) std::ceil(log(static_cast<double> 
          ( ((double)Q0_nat) +  ((double)(Q1_nat+ 1)) * ((double)Delta)) ) / log(static_cast<double> (2)));

    if ( (test_bit1 + shift_len) > gBitsMaxBits) {
      PALISADE_THROW(config_error, "SGD is unstable");
    } 
		*/
    // VARIANT B
    
		NativeInteger::SignedNativeInt baseG = NativeInteger::SignedNativeInt(Base);
    NativeInteger::SignedNativeInt gminus1 = Base - 1;
    NativeInteger::SignedNativeInt baseGdiv2 = (Base >> 1)-1;
    

		//uint64_t t_down = 0;
    //uint64_t  t_up  = 0;
		// NativeInteger::SignedNativeInt d_up = 0;
    //NativeInteger::SignedNativeInt d_down = 0;
			
		#pragma omp parallel for
		//{ 
           for (uint32_t k = 0; k < Ns; k++) {
			NativeInteger::SignedNativeInt d_up = 0;
			NativeInteger::SignedNativeInt d_down = 0;
	
			uint64_t t_down = (in_poly.GetElementAtIndex(0)[k]).ConvertToInt(); // 0 ~  Q0-1
			uint64_t t_up   = (in_poly.GetElementAtIndex(1)[k]).ModSub(in_poly.GetElementAtIndex(0)[k], Q1_nat).ModMulFast(Inv_1, Q1_nat, Mu1).ConvertToInt(); // 0 ~  Q0-1

			// Test If it is Minus
      bool rev = false;
      if ((t_up > QOver2Info[1]) || (  (t_up == QOver2Info[1])   && (t_down > QOver2Info[0])) ) {
				rev = true;

        if ( t_down == 0 ) {
					t_up = Q1_nat - t_up;
        } else {
          t_up    = Q1_nat - t_up - 1;
          t_down  = Q0_nat  - t_down;
        }

      }
      t_down += t_up;
      t_up *= Delta;

      if (rev) {
				d_down  -= (NativeInteger::SignedNativeInt) (t_down & p_mask);
        d_up    -= (NativeInteger::SignedNativeInt) ((t_up << shift_len) + ((t_down & p_mask_com) >> point_bit));
      } else {
        d_down  += (NativeInteger::SignedNativeInt) (t_down & p_mask);
				d_up    += (NativeInteger::SignedNativeInt) ((t_up << shift_len) + ((t_down & p_mask_com) >> point_bit));
      }

      for (uint32_t lo = 0; lo < lower_len; lo++) {

				// This approach gives a slightly better performance
        // VARIANT A
        //NativeInteger::SignedNativeInt r = d_down << gBitsMaxBits;
        //r >>= gBitsMaxBits;
        // VARIANT B
        NativeInteger::SignedNativeInt r = d_down & gminus1;  
        if (r > baseGdiv2) r -= baseG;

        d_down -= r;
        d_down >>= gBits;

        if (lo < rm_len) {
					continue;
        }
        // It should be zeros in values, or It has some intension
				out_poly[k][lo-rm_len] = (int32_t) r;
			}
      d_up += d_down;
        for (uint32_t hi = 0; hi < upper_len; hi++) {
          // This approach gives a slightly better performance
          // VARIANT A
          //NativeInteger::SignedNativeInt r = d_up << gBitsMaxBits;
          //r >>= gBitsMaxBits;

          // VARIANT B
          NativeInteger::SignedNativeInt r = d_up & gminus1;
          if (r > baseGdiv2) r -= baseG;

          d_up -= r;
          d_up >>= gBits;
          if (hi+lower_len < rm_len) {
            continue;
          }
					out_poly[k][hi+lower_len-rm_len] = (int32_t) r;
    
        }
        d_up = 0;
        d_down = 0;
      } //}
			return std::make_shared<vector<vector<int32_t>>> (out_poly);
    }










  std::shared_ptr<vector<DCRTPoly>> FPRingGSWAccumulatorScheme::ProductDCRTPARALLEL(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const std::shared_ptr<FPRingEvaluationKey> Evkey, 
		const vector<DCRTPoly> &A1, const DCRTPoly &b1,
    const vector<DCRTPoly> &A2, const DCRTPoly &b2) const {


		auto starts = std::chrono::system_clock::now();
	
    // A1, A2, b1, b2, outA, outb should be Eval mode!!
		//
		if (A1[0].GetFormat() != Format::EVALUATION) {
			PALISADE_THROW(config_error, "ProductDCRT input should be Eval modes!");
		}
		if (A2[0].GetFormat() != Format::EVALUATION) {
			PALISADE_THROW(config_error, "ProductDCRT input should be Eval modes!");
		}
		
    uint32_t Ks_FP = params->GetRLWEParams()->GetK();
    vector<DCRTPoly> out(Ks_FP+1);
    for (uint32_t i = 0; i< Ks_FP+1; i++) {
      out[i] = DCRTPoly(params->GetRLWEParams()->GetDCRTParams(), Format::EVALUATION, true);
		}

		 uint32_t EV_bit       = params->GetBaseEVBit();
    uint32_t EV_len_rm    = params->GetLenEVrm();
    uint32_t EV_len_save  = params->GetLenEVsave();

    // Must FP Params
		// Dim0: Ks_idx, Dim1: GD len, Dim 2: A1,..., Ak, B
    std::vector<std::vector<std::vector<DCRTPoly>>> m_e_key = Evkey->GetElements();





		std::shared_ptr<vector<DCRTPoly>> SGD_vec;
    DCRTPoly A_tmp;
    uint32_t k_idx = 0;
    for (uint32_t k_row = 0; k_row < Ks_FP; k_row++) {
      for (uint32_t k_col = 0; k_col <= k_row; k_col++) {
				
				k_idx = ((k_row * (k_row + 1)) >> 1) +k_col;
				if (k_row == k_col) {
          A_tmp  = A1[k_row];
          A_tmp *= A2[k_col];
        } else {
          A_tmp = A1[k_row];
          A_tmp *= A2[k_col];
          A_tmp += (A1[k_col] * A2[k_row]);
        }
				A_tmp.SetFormat(Format::COEFFICIENT);
				SGD_vec = this->SGDLevel2PARALLEL(params, A_tmp, EV_bit, EV_len_rm, Format::EVALUATION);

				#pragma omp parallel for collapse(1)
				//{
        for (uint32_t out_a_idx = 0; out_a_idx < Ks_FP+1; out_a_idx++) {
					for (uint32_t gd_idx = 0; gd_idx < EV_len_save; gd_idx++) {
						// A1, A2 ,...., Ak, B
						out[out_a_idx] += (m_e_key[k_idx][gd_idx][out_a_idx] * (*SGD_vec)[gd_idx]);
          }
				}//}
      }
    }

    // Collect
    #pragma omp parallel for
		//{
    for (uint32_t k_idx = 0; k_idx < Ks_FP; k_idx++) {
      out[k_idx] += A1[k_idx] * b2;
      out[k_idx] += A2[k_idx] * b1;
    } //}
    out[Ks_FP] += b1 * b2;

		std::chrono::duration<double> ends = std::chrono::system_clock::now() - starts;
		if (DEBUG) std::cout << "Tensor Product PARALLEL 1 time is " << ends.count() << std::endl;

    return std::make_shared<vector<DCRTPoly>>(out);
  };


  std::shared_ptr<vector<DCRTPoly>> FPRingGSWAccumulatorScheme::ProductDCRTPARALLEL2(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const std::shared_ptr<FPRingEvaluationKey> Evkey, 
		const vector<DCRTPoly> &A1, const DCRTPoly &b1,
    const vector<DCRTPoly> &A2, const DCRTPoly &b2) const {

		auto starts = std::chrono::system_clock::now();
	
		// A1, A2, b1, b2, outA, outb should be Eval mode!!
		//
		if (A1[0].GetFormat() != Format::EVALUATION) {
			PALISADE_THROW(config_error, "ProductDCRT input should be Eval modes!");
		}
		if (A2[0].GetFormat() != Format::EVALUATION) {
			PALISADE_THROW(config_error, "ProductDCRT input should be Eval modes!");
		}
		
    uint32_t Ks_FP = params->GetRLWEParams()->GetK();
    vector<DCRTPoly> out(Ks_FP+1);
    for (uint32_t i = 0; i< Ks_FP+1; i++) {
      out[i] = DCRTPoly(params->GetRLWEParams()->GetDCRTParams(), Format::EVALUATION, true);
		}

		uint32_t EV_bit       = params->GetBaseEVBit();
    uint32_t EV_len_rm    = params->GetLenEVrm();
    uint32_t EV_len_save  = params->GetLenEVsave();

    // Must FP Params
		// Dim0: Ks_idx, Dim1: GD len, Dim 2: A1,..., Ak, B
    std::vector<std::vector<std::vector<DCRTPoly>>> m_e_key = Evkey->GetElements();
		uint32_t EV_k_len = m_e_key.size(); 
  


		//std::shared_ptr<vector<DCRTPoly>> SGD_vec;
  	vector<std::shared_ptr<vector<DCRTPoly>>> SGD_vec(EV_k_len);
		//uint32_t k_idx = 0;
		#pragma omp parallel for collapse(1)
		//{
    for (uint32_t k_row = 0; k_row < Ks_FP; k_row++) {
      for (uint32_t k_col = 0; k_col <= k_row; k_col++) {
				DCRTPoly A_tmp;
				uint32_t k_idx = ((k_row * (k_row + 1)) >> 1) +k_col;
				if (k_row == k_col) {
          A_tmp  = A1[k_row];
          A_tmp *= A2[k_col];
        } else {
          A_tmp = A1[k_row];
          A_tmp *= A2[k_col];
          A_tmp += (A1[k_col] * A2[k_row]);
        }
				A_tmp.SetFormat(Format::COEFFICIENT);
				SGD_vec[k_idx] = this->SGDLevel2PARALLEL(params, A_tmp, EV_bit, EV_len_rm, Format::EVALUATION);
			}
		}//}

		//std::shared_ptr<vector<DCRTPoly>> SGD_vec;
    //DCRTPoly A_tmp;
    //uint32_t k_idx = 0;
    #pragma omp parallel for  collapse(1)
		//{
    for (uint32_t out_a_idx = 0; out_a_idx < Ks_FP+1; out_a_idx++) {		
			for (uint32_t k_row = 0; k_row < Ks_FP; k_row++) {
				for (uint32_t k_col = 0; k_col <= k_row; k_col++) {
					uint32_t k_idx = ((k_row * (k_row + 1)) >> 1) +k_col;
					for (uint32_t gd_idx = 0; gd_idx < EV_len_save; gd_idx++) {
						// A1, A2 ,...., Ak, B
						out[out_a_idx] += (m_e_key[k_idx][gd_idx][out_a_idx] * (*(SGD_vec[k_idx]))[gd_idx]);
          }
				}}
      }//}


    // Collect
    #pragma omp parallel for
		//{
    for (uint32_t k_idx = 0; k_idx < Ks_FP; k_idx++) {
      out[k_idx] += A1[k_idx] * b2;
      out[k_idx] += A2[k_idx] * b1;
    }//}
    out[Ks_FP] += b1 * b2;
		
		std::chrono::duration<double> ends = std::chrono::system_clock::now() - starts;
		if (DEBUG) std::cout << "Tensor Product PARALLEL 2 time is " << ends.count() << std::endl;

    return std::make_shared<vector<DCRTPoly>>(out);
  };
	



  std::shared_ptr<vector<DCRTPoly>> FPRingGSWAccumulatorScheme::ProductDCRT(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const std::shared_ptr<FPRingEvaluationKey> Evkey, 
		const vector<DCRTPoly> &A1, const DCRTPoly &b1,
    const vector<DCRTPoly> &A2, const DCRTPoly &b2) const {


		auto starts = std::chrono::system_clock::now();
	
    // A1, A2, b1, b2, outA, outb should be Eval mode!!
		//
		if (A1[0].GetFormat() != Format::EVALUATION) {
			PALISADE_THROW(config_error, "ProductDCRT input should be Eval modes!");
		}
		if (A2[0].GetFormat() != Format::EVALUATION) {
			PALISADE_THROW(config_error, "ProductDCRT input should be Eval modes!");
		}
		
    uint32_t Ks_FP = params->GetRLWEParams()->GetK();
    vector<DCRTPoly> out(Ks_FP+1);
    for (uint32_t i = 0; i< Ks_FP+1; i++) {
      out[i] = DCRTPoly(params->GetRLWEParams()->GetDCRTParams(), Format::EVALUATION, true);
		}

		 uint32_t EV_bit       = params->GetBaseEVBit();
    uint32_t EV_len_rm    = params->GetLenEVrm();
    uint32_t EV_len_save  = params->GetLenEVsave();

    // Must FP Params
		// Dim0: Ks_idx, Dim1: GD len, Dim 2: A1,..., Ak, B
    std::vector<std::vector<std::vector<DCRTPoly>>> m_e_key = Evkey->GetElements();


		std::shared_ptr<vector<DCRTPoly>> SGD_vec;
    DCRTPoly A_tmp;
    uint32_t k_idx = 0;
    for (uint32_t k_row = 0; k_row < Ks_FP; k_row++) {
      for (uint32_t k_col = 0; k_col <= k_row; k_col++) {
				
				k_idx = ((k_row * (k_row + 1)) >> 1) +k_col;
				if (k_row == k_col) {
          A_tmp  = A1[k_row];
          A_tmp *= A2[k_col];
        } else {
          A_tmp = A1[k_row];
          A_tmp *= A2[k_col];
          A_tmp += (A1[k_col] * A2[k_row]);
        }
				A_tmp.SetFormat(Format::COEFFICIENT);
				SGD_vec = this->SGDLevel2(params, A_tmp, EV_bit, EV_len_rm, Format::EVALUATION);

				{for (uint32_t out_a_idx = 0; out_a_idx < Ks_FP+1; out_a_idx++) {
					for (uint32_t gd_idx = 0; gd_idx < EV_len_save; gd_idx++) {
						// A1, A2 ,...., Ak, B
						out[out_a_idx] += (m_e_key[k_idx][gd_idx][out_a_idx] * (*SGD_vec)[gd_idx]);
          }
				}}
      }
    }

    // Collect
		{for (uint32_t k_idx = 0; k_idx < Ks_FP; k_idx++) {
      out[k_idx] += A1[k_idx] * b2;
      out[k_idx] += A2[k_idx] * b1;
    }}
    out[Ks_FP] += b1 * b2;

		std::chrono::duration<double> ends = std::chrono::system_clock::now() - starts;
		if (DEBUG) std::cout << "Tensor Product time is " << ends.count() << std::endl;

    return std::make_shared<vector<DCRTPoly>>(out);
  };

  void FPRingGSWAccumulatorScheme::DivideRounding(
			const std::shared_ptr<FPRingGSWCryptoParams> paramss,
			const DCRTPoly &double_poly, 
			std::shared_ptr<DCRTPoly> out_poly, 
			const Format format) const {

    const std::shared_ptr<FPRLWECryptoParams> params = paramss->GetRLWEParams();
		uint32_t Ns = double_poly.GetLength();
    // Invariant to N
    vector<NativeInteger> QInv = params->GetInv();
    vector<NativeInteger> moduliQ = params->GetModuliQ();
    // FP = GSW assume
		NativeInteger Mu0 = params->GetMuFP()[0];

		uint64_t Q1_half = moduliQ[1].ConvertToInt() >> 1; 

		if (double_poly.GetFormat() != Format::COEFFICIENT) {
			PALISADE_THROW(config_error, "Dividing shoud be Coeff MODE!");
		}

		out_poly->SetFormat(Format::COEFFICIENT);    
		vector<NativePoly> poly = double_poly.GetAllElements();
    NativeInteger a0;
    for (uint32_t i = 0; i < Ns; i++) {
      a0 = poly[0][i].ModSub(poly[1][i], moduliQ[0]).ModMulFast(QInv[0], moduliQ[0], Mu0);
      if(poly[1][i].ConvertToInt() >= Q1_half) {
        a0.ModAddFastEq(1,moduliQ[0]);
      }    
      out_poly->GetElementW(0)[i] = a0;
    }    
    out_poly->SetFormat(format);
  } 





/*
  void FPRingGSWAccumulatorScheme::DivideRounding(
			const std::shared_ptr<FPRingGSWCryptoParams> paramss,
			const DCRTPoly &double_poly, 
			std::shared_ptr<DCRTPoly> out_poly, 
			const Format format) const {

    const std::shared_ptr<FPRLWECryptoParams> params = paramss->GetRLWEParams();
		uint32_t Ns = double_poly.GetLength();
    // Invariant to N
    vector<NativeInteger> QInv = params->GetInv();
    vector<NativeInteger> moduliQ = params->GetModuliQ();
    uint64_t Q1_half = moduliQ[1].ConvertToInt() >> 1; 

		if (double_poly.GetFormat() != Format::COEFFICIENT) {
			PALISADE_THROW(config_error, "Dividing shoud be Coeff MODE!");
		}

		out_poly->SetFormat(Format::COEFFICIENT);    
		vector<NativePoly> poly = double_poly.GetAllElements();
    NativeInteger a0;
    for (uint32_t i = 0; i < Ns; i++) {
      a0 = poly[0][i].ModSub(poly[1][i], moduliQ[0]).ModMul(QInv[0], moduliQ[0]);
      if(poly[1][i].ConvertToInt() >= Q1_half) {
        a0.ModAddEq(1,moduliQ[0]);
      }    
      out_poly->GetElementW(0)[i] = a0;
    }    
    out_poly->SetFormat(format);
  } 
*/

	// Fast
  void FPRingGSWAccumulatorScheme::DivideRoundingSelf(
			const std::shared_ptr<FPRingGSWCryptoParams> paramss,
			DCRTPoly *double_poly, 
			const Format format) const {
		
	  
		const std::shared_ptr<FPRLWECryptoParams> params = paramss->GetRLWEParams();
		uint32_t Ns = double_poly->GetLength();
    // Invariant to N
    vector<NativeInteger> QInv = params->GetInv();
    vector<NativeInteger> moduliQ = params->GetModuliQ();
    uint64_t Q1_half = moduliQ[1].ConvertToInt() >> 1; 
		NativeInteger Mu0 = params->GetMuFP()[0];


		if (double_poly->GetFormat() != Format::COEFFICIENT) {
			double_poly->SetFormat(Format::COEFFICIENT);
			//PALISADE_THROW(config_error, "Dividing shoud be Coeff MODE!");
		}
		vector<NativePoly> poly = double_poly->GetAllElements();
    NativeInteger a0;
    for (uint32_t i = 0; i < Ns; i++) {
      a0 = poly[0][i].ModSub(poly[1][i], moduliQ[0]).ModMulFast(QInv[0], moduliQ[0], Mu0);
      if(poly[1][i].ConvertToInt() >= Q1_half) {
        a0.ModAddFastEq(1,moduliQ[0]);
      }    
      double_poly->GetElementW(0)[i] = a0;
    }    
		double_poly->DropLastElement();
		double_poly->SetFormat(format);
  
	} 


std::shared_ptr<DCRTPoly> FPRingGSWAccumulatorScheme::DivideRounding(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const DCRTPoly &in_poly, 
		const Format format=Format::COEFFICIENT) const {
    
		std::shared_ptr<DCRTPoly> out_poly = std::make_shared<DCRTPoly> (in_poly.GetParams(), Format::COEFFICIENT,true);
		
    this->DivideRounding(params, in_poly, out_poly, format);
    return out_poly;
    }   
	

	std::shared_ptr<vector<NativeVector>>  FPRingGSWAccumulatorScheme::RotandCarryAddA(
  	const std::shared_ptr<FPRingGSWCryptoParams> params,
		const vector<NativePoly> &input,
		vector<DCRTPoly> *CarryInfoA,
		const NativePoly &rot,
		//const uint32_t &idx,
		const bool rot_flag,
		const bool carry_flag
  ) const {
		uint32_t K_FP = input.size();
		vector<NativeVector> out(input.size());
		NativePoly A_Poly_tmp;
	
		// Rotations
		for (uint32_t k = 0; k < K_FP; k++) {
			if (rot_flag) {
				// Copya
				A_Poly_tmp = (input[k] * rot);
			} else {
				A_Poly_tmp = (input[k]);
			}
			A_Poly_tmp.SetFormat(Format::COEFFICIENT);
			// Carry Merging
			if (carry_flag) {
				// Carry and A,B Eval mod
				// Input & Output is COEFF Modesa
				(*CarryInfoA)[k].SetFormat(Format::COEFFICIENT);
				DivideRoundingSelf(params, &((*CarryInfoA)[k]), Format::COEFFICIENT);
				A_Poly_tmp += (*CarryInfoA)[k].GetElementAtIndex(0);
      } 
				out[k] = std::move(A_Poly_tmp.GetValues());
		}
		return std::make_shared<vector<NativeVector>> (out);
	}
	
	std::shared_ptr<vector<vector<DCRTPoly>>>  FPRingGSWAccumulatorScheme::ProductAndBKAndKS(
		const std::shared_ptr<FPRingGSWCryptoParams> params, 
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const DCRTPoly &ACC_init,
    const std::shared_ptr<FPRingGSWEvalKey> &EK, 
		const vector<DCRTPoly> &C1, 
		const vector<DCRTPoly> &C2,
		const uint32_t boot_loc_idx, const uint32_t boot_first_bits,
		const uint64_t prefix_add, const uint32_t boot_out_bit,
		const uint32_t out_num) const {

		uint32_t K_FP = params->GetRLWEParams()->GetK();
		// Product
		std::shared_ptr<vector<DCRTPoly>> mul_res;
		if (params->GetPARALLEL()) {
			mul_res = this->ProductDCRTPARALLEL2(params, EK->Evkey, C1, C1[K_FP], C2, C2[K_FP]); 
		} else {
			mul_res = this->ProductDCRT(params, EK->Evkey, C1, C1[K_FP], C2, C2[K_FP]); 
		}
		// Reduction
		vector<NativeVector> A(K_FP);

		for (uint32_t k = 0; k < K_FP; k++) {
			(*mul_res)[k].SetFormat(Format::COEFFICIENT);
			DivideRoundingSelf(params, &((*mul_res)[k]), Format::COEFFICIENT); 
			A[k] = std::move((*mul_res)[k].GetElementAtIndex(0).GetValues());
		}
		(*mul_res)[K_FP].SetFormat(Format::COEFFICIENT);
		DivideRoundingSelf(params, &((*mul_res)[K_FP]), Format::COEFFICIENT); 
		NativeInteger b = (*mul_res)[K_FP].GetElementAtIndex(0).GetValues()[0];

		return	this->Bootstrap_and_KS(params
            , LWEscheme
            , ACC_init, EK
            , A
            , b
            , boot_loc_idx, boot_first_bits, prefix_add, boot_out_bit, out_num);
	}

	

	
	std::shared_ptr<vector<DCRTPoly>>  FPRingGSWAccumulatorScheme::ProductAndBKAndKSRange(
		const std::shared_ptr<FPRingGSWCryptoParams> params, 
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const DCRTPoly &ACC_init,
    const std::shared_ptr<FPRingGSWEvalKey> &EK, 
		const vector<DCRTPoly> &C1, 
		const vector<DCRTPoly> &C2,
		const uint32_t boot_loc_idx_from, 
		const uint32_t boot_loc_idx_to, 
		const uint32_t boot_first_bits,
		const uint64_t prefix_add, const uint32_t boot_out_bit,
		const uint32_t out_num, bool mul) const {
		if (out_num != 1) {
			PALISADE_THROW(config_error, "PABAKRRev out_num should be 1");
		}
		if (boot_out_bit != 0) {
			PALISADE_THROW(config_error, "PABAKRRev out_but should be 0");
		}

		uint32_t K_FP = params->GetRLWEParams()->GetK();
		uint32_t M_FP = (params->GetRLWEParams()->GetN() * 2);
		vector<NativeVector> A(K_FP);
		vector<DCRTPoly> Outs(K_FP+1);
	
		// Product
		std::shared_ptr<vector<DCRTPoly>> mul_res;
		if (params->GetPARALLEL()) {
			mul_res = this->ProductDCRTPARALLEL2(params, EK->Evkey, C1, C1[K_FP], C2, C2[K_FP]); 
		} else {
			mul_res = this->ProductDCRT(params, EK->Evkey, C1, C1[K_FP], C2, C2[K_FP]); 
		}

		// Reduction
		uint32_t iter_num = boot_loc_idx_to - boot_loc_idx_from + 1;
		for (uint32_t k = 0; k < K_FP; k++) {
			(*mul_res)[k].SetFormat(Format::COEFFICIENT);
			DivideRoundingSelf(params, &((*mul_res)[k]), Format::COEFFICIENT); 
			A[k] = std::move((*mul_res)[k].GetElementAtIndex(0).GetValues());
		}
		(*mul_res)[K_FP].SetFormat(Format::COEFFICIENT);
		DivideRoundingSelf(params, &((*mul_res)[K_FP]), Format::COEFFICIENT); 
		
		if (params->GetPARALLEL() & (iter_num >= 4)) {
			vector<vector<DCRTPoly>> res_save(iter_num);
			for (uint32_t iter_l = 0; iter_l < iter_num; iter_l++) {
				res_save[iter_l].resize(K_FP+1);
			}
			#pragma omp parallel for
			//{
              for (uint32_t iter = 0; iter < iter_num; iter++) {
				NativeInteger b = (*mul_res)[K_FP].GetElementAtIndex(0).GetValues()[iter+boot_loc_idx_from];
				auto tmp_res =	std::move((*this->Bootstrap_and_KS(params
					, LWEscheme
          , ACC_init, EK
          , A
          , b
          , boot_loc_idx_from + iter, boot_first_bits, prefix_add, boot_out_bit, out_num))[0]);
				
				if (mul == false) {  // Make X^-(m)
					// if from =4 , to = 8, vals are X^4, X^5, X^6, X^7, X^8
					// then, we just Adding its...
					// But if not, 4 -> X^-4, 5 -> X^-5 , ..., 8-> X^-8
					// this means, 2*M_FP - 2*(boot_loc_idx_from_iter) mod M_FP
					uint32_t rot_idxs = ((2*M_FP) - 2 * (boot_loc_idx_from + iter)) % M_FP; 
			
					//std::cout << "rot idx is " << rot_idxs << std::endl;
					for (uint32_t k = 0; k < K_FP+1; k++) {	
						/*
						if (iter == 0) {
							Outs[k] = tmp_res[k] * params->GetRotMonomial_RLWE(rot_idxs); 
						} else {
							// it shifted boot_lic_idx_from+iter
							Outs[k] += (tmp_res[k] * params->GetRotMonomial_RLWE(rot_idxs)); 
						}	*/
						res_save[iter][k] = std::move(tmp_res[k]);
						res_save[iter][k] *= params->GetRotMonomial_RLWE(rot_idxs); 
					}
				} else {
					
					/*for (uint32_t k = 0; k < K_FP+1; k++) {
						if (iter == 0) {
							Outs[k] = tmp_res[k]; 
						} else {
							Outs[k] += tmp_res[k]; 
						}
					}*/
					res_save[iter] = std::move(tmp_res);

				}
			}//}
			for (uint32_t iter = 0; iter < iter_num; iter++) {
				if (iter == 0) {
					Outs = std::move(res_save[0]);

				} else {
					for (uint32_t k = 0; k < K_FP+1; k++) {
						Outs[k] += res_save[iter][k]; 
					}
				}
			}
		} else {
			for (uint32_t iter = 0; iter < iter_num; iter++) {
				NativeInteger b = (*mul_res)[K_FP].GetElementAtIndex(0).GetValues()[iter+boot_loc_idx_from];
				auto tmp_res =	std::move((*this->Bootstrap_and_KS(params
            , LWEscheme
            , ACC_init, EK
            , A
            , b
            , boot_loc_idx_from + iter, boot_first_bits, prefix_add, boot_out_bit, out_num))[0]);
				
				if (mul == false) {  // Make X^-(m)
					// if from =4 , to = 8, vals are X^4, X^5, X^6, X^7, X^8
					// then, we just Adding its...
					// But if not, 4 -> X^-4, 5 -> X^-5 , ..., 8-> X^-8
					// this means, 2*M_FP - 2*(boot_loc_idx_from_iter) mod M_FP
					uint32_t rot_idxs = ((2*M_FP) - 2 * (boot_loc_idx_from + iter)) % M_FP; 
			
					//std::cout << "rot idx is " << rot_idxs << std::endl;
					for (uint32_t k = 0; k < K_FP+1; k++) {	
						if (iter == 0) {
							Outs[k] = tmp_res[k] * params->GetRotMonomial_RLWE(rot_idxs); 
						} else {
							// it shifted boot_lic_idx_from+iter
							Outs[k] += (tmp_res[k] * params->GetRotMonomial_RLWE(rot_idxs)); 
						}	
					}
				} else {
					for (uint32_t k = 0; k < K_FP+1; k++) {
						if (iter == 0) {
							Outs[k] = tmp_res[k]; 
						} else {
							Outs[k] += tmp_res[k]; 
						}
					}
				}
			}
		}
		return	std::make_shared<vector<DCRTPoly>> (Outs);
	}





  std::shared_ptr<vector<DCRTPoly>> FPRingGSWAccumulatorScheme::BootstrapIsLower(
    const std::shared_ptr<FPRingGSWCryptoParams> params, 
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<DCRTPoly> & ct,
		const uint64_t &mins,
		const std::shared_ptr<FPRingGSWEvalKey> EK) const {
		
		auto MulPrepare_start = std::chrono::system_clock::now();
    uint32_t K_FP = params->GetRLWEParams()->GetK();
		NativeInteger Q0 = params->GetModuliQ()[0];
    NativeInteger Q1 = params->GetModuliQ()[1];
		uint64_t Q1_Int = Q1.ConvertToInt();
		uint32_t Q1_WO_bits = params->GetRLWEParams()->GetBits()[1];  // <32 bit
		uint32_t start_bits;
		uint64_t prefix_add_vals;
		uint32_t carry_bit;
		uint32_t out_nums;
	
		std::shared_ptr<vector<vector<DCRTPoly>>> BS_out_vec;
		vector<DCRTPoly> BS_out;
		vector<NativeVector> A(K_FP);
		NativeInteger b;
		DCRTPoly DCRT_tmp;
		vector<DCRTPoly> In_Poly(K_FP+1);
	
		vector<vector<DCRTPoly>> Res_Poly(2);
		vector<DCRTPoly> Outs(K_FP+1);
		
		// ct1 - ct 2
		for (uint32_t i = 0; i < K_FP+1; i++) {	
			In_Poly[i] =  ct[i];
		}

		// Substract gives range - 11bit~ 12 bit
		// add 13 bit, bootstrapping, and check whether 13bit is 0 or not
		uint32_t scale_bit = params->GetExpBit();	
		uint64_t prefix = params->GetRLWEParams()->GetScaling().ConvertToInt();
			prefix <<= scale_bit;
		uint64_t bit_check = prefix;	
		bit_check <<= 11;  // 2^10 ? 2^ 11?	
		
		
		// Make

		In_Poly[K_FP].SetFormat(Format::COEFFICIENT);
		In_Poly[K_FP].GetElementW(0)[0].ModAddEq(
				NativeInteger(bit_check).ModMul(Q1, Q0),Q0);
		In_Poly[K_FP].GetElementW(1)[0].ModAddEq(
				NativeInteger(bit_check).ModMul(Q1, Q1),Q1);
		// Sub Values
		//
		In_Poly[K_FP].GetElementW(0)[0].ModSubEq(
				NativeInteger(prefix).ModMul(Q1, Q0).ModMul(mins,Q0),Q0);
		In_Poly[K_FP].GetElementW(1)[0].ModSubEq(
				NativeInteger(prefix).ModMul(Q1, Q1).ModMul(mins,Q1),Q1);
		In_Poly[K_FP].SetFormat(Format::EVALUATION);
		

		////////////////////////////////   STARTS !!!! /////////////////////
		std::chrono::duration<double> MulPrepare_end = std::chrono::system_clock::now() - MulPrepare_start;
		if (DEBUG) {
			std::cout << "bootstrapping Relu prepare time is " << MulPrepare_end.count() << std::endl ;
		}


		// Bootstrap 1
		for (uint32_t i = 0; i < 3; i ++) { 
			
			//std::cout << "i is " << i << std::endl;	
			// Make Polys
			if (i != 2) {
				// Copy A
				for (uint32_t j = 0; j < K_FP; j++) {
					DCRT_tmp = In_Poly[j];
					DCRT_tmp.SetFormat(Format::COEFFICIENT);
					DivideRoundingSelf(params, &DCRT_tmp, Format::COEFFICIENT);
					A[j] = DCRT_tmp.GetElementAtIndex(0).GetValues();
				}
				// Copy B
				DCRT_tmp = In_Poly[K_FP];
				DCRT_tmp.SetFormat(Format::COEFFICIENT);
				DivideRoundingSelf(params, &DCRT_tmp, Format::COEFFICIENT);
				b = DCRT_tmp.GetElementAtIndex(0)[0];
			} else { // i == 2
				for (uint32_t j = 0; j < K_FP; j++) {
					In_Poly[j].SetFormat(Format::COEFFICIENT);
					DivideRoundingSelf(params, &In_Poly[j], Format::COEFFICIENT);
					A[j] = In_Poly[j].GetElementAtIndex(0).GetValues();
				}
				In_Poly[K_FP].SetFormat(Format::COEFFICIENT);
				DivideRoundingSelf(params, &In_Poly[K_FP], Format::COEFFICIENT);
				b = In_Poly[K_FP].GetElementAtIndex(0)[0];
			}
				
			// First Carry, 4+1 (1 bit no information)bit bootstrapping Bootstrapping
      if (i != 2) {
				if (i == 0) {
					start_bits  = scale_bit  + 4 + 1;
					prefix_add_vals = (prefix >> 1);  
				} else if (i == 1) {
					start_bits  = scale_bit  + 8 + 1;
					prefix_add_vals = (prefix << 3);  
				} 
				carry_bit = 1;
				out_nums = 2;
			
				//First Bootstrapping,3 bit + random 1 bit is bootstrapping
				BS_out_vec = this->Bootstrap_and_KS(params
					, LWEscheme
					, params->GetACCRelu(0)
					, EK
					, A
					, b, 0, start_bits, prefix_add_vals, carry_bit, out_nums);
		
				//Second Bootsrapp
				// i = 1, 2, 3 ...
				DCRTPoly ACCs = params -> GetACCRelu(i+1);
				//	ACCs = params -> GetACCRelu(6);
				
				// Orign, 4+1(padding but no infomation) bit bootstrapping msg <= 3bit 
				start_bits  = Q1_WO_bits + 4 + 1;
				carry_bit = 1;
				out_nums = 2;
				prefix_add_vals = ( Q1_Int >> 1);  
				// ================   EVAL MODE SAVING ==================== 
				BS_out_vec = this->ProductAndBKAndKS(params
					, LWEscheme
					, ACCs, EK
					, (*BS_out_vec)[0]
					, (*BS_out_vec)[1]
					, 0, start_bits, prefix_add_vals, carry_bit, out_nums);
			} else if (i == 2) {
				start_bits  = scale_bit  + 12 + 1;
				prefix_add_vals = (prefix << 7);  
				carry_bit = 0;
				out_nums = 1;
			
				//First Bootstrapping,3 bit + random 1 bit is bootstrapping
				BS_out_vec = this->Bootstrap_and_KS(params
					, LWEscheme
					, params->GetACCRelu(19)
					, EK
					, A
					, b, 0, start_bits, prefix_add_vals, carry_bit, out_nums);
				Outs = std::move((*BS_out_vec)[0]);
			}
			// After Handle
			//std::cout << "handle before " << std::endl;
			if (i != 2) {
				
				// remove 4 or 8 bit info
				for (uint32_t j = 0; j < K_FP+1; j++) {
					In_Poly[j]-=(*BS_out_vec)[1][j]; 
				}

				//std::cout << "Minus is end" << std::endl;
				// Saves results 
				Res_Poly[i] = std::move((*BS_out_vec)[0]);
				//std::cout << "Res is moved" << std::endl;
				
			}

		}
		return std::make_shared<vector<DCRTPoly>>(Outs);
	}








  std::shared_ptr<vector<DCRTPoly>> FPRingGSWAccumulatorScheme::BootstrapIsLowerF(
    const std::shared_ptr<FPRingGSWCryptoParams> params, 
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<DCRTPoly> & ct,
		const uint64_t &mins,
		const std::shared_ptr<FPRingGSWEvalKey> EK) const {
	
		auto MulPrepare_start = std::chrono::system_clock::now();
    uint32_t K_FP = params->GetRLWEParams()->GetK();
		NativeInteger Q0 = params->GetModuliQ()[0];
    NativeInteger Q1 = params->GetModuliQ()[1];
		uint64_t Q1_Int = Q1.ConvertToInt();
		uint32_t Q1_WO_bits = params->GetRLWEParams()->GetBits()[1];  // <32 bit
		uint32_t start_bits;
		uint64_t prefix_add_vals;
		uint32_t carry_bit;
		uint32_t out_nums;
	
		std::shared_ptr<vector<vector<DCRTPoly>>> BS_out_vec;
		vector<DCRTPoly> BS_out;
		vector<NativeVector> A(K_FP);
		NativeInteger b;
		DCRTPoly DCRT_tmp;
		vector<DCRTPoly> In_Poly(K_FP+1);
	
		vector<vector<DCRTPoly>> Res_Poly(2);
		vector<DCRTPoly> Outs(K_FP+1);
		
		// ct1 - ct 2
		for (uint32_t i = 0; i < K_FP+1; i++) {	
			In_Poly[i] =  ct[i];
		}

		// Substract gives range - 11bit~ 12 bit
		
		uint32_t scale_bit = params->GetExpBit();	
		uint64_t prefix = params->GetRLWEParams()->GetScaling().ConvertToInt();
			prefix <<= scale_bit;
		uint64_t bit_check = prefix;	
		bit_check <<= 8;  // Check 8 bit	
		
		// Make
		In_Poly[K_FP].SetFormat(Format::COEFFICIENT);
		In_Poly[K_FP].GetElementW(0)[0].ModAddEq(
				NativeInteger(bit_check).ModMul(Q1, Q0),Q0);
		In_Poly[K_FP].GetElementW(1)[0].ModAddEq(
				NativeInteger(bit_check).ModMul(Q1, Q1),Q1);
		
		// Sub Values
		In_Poly[K_FP].GetElementW(0)[0].ModSubEq(
				NativeInteger(prefix).ModMul(Q1, Q0).ModMul(mins,Q0),Q0);
		In_Poly[K_FP].GetElementW(1)[0].ModSubEq(
				NativeInteger(prefix).ModMul(Q1, Q1).ModMul(mins,Q1),Q1);
		In_Poly[K_FP].SetFormat(Format::EVALUATION);
		
		////////////////////////////////   STARTS !!!! /////////////////////
		std::chrono::duration<double> MulPrepare_end = std::chrono::system_clock::now() - MulPrepare_start;
		if (DEBUG) {
			std::cout << "bootstrapping Relu prepare time is " << MulPrepare_end.count() << std::endl ;
		}

		// Bootstrap 1
		for (uint32_t i = 0; i < 3; i ++) { 
			
			//std::cout << "i is " << i << std::endl;	
			// Make Polys
			if (i != 2) {
				// Copy A
				for (uint32_t j = 0; j < K_FP; j++) {
					DCRT_tmp = In_Poly[j];
					DCRT_tmp.SetFormat(Format::COEFFICIENT);
					DivideRoundingSelf(params, &DCRT_tmp, Format::COEFFICIENT);
					A[j] = DCRT_tmp.GetElementAtIndex(0).GetValues();
				}
				// Copy B
				DCRT_tmp = In_Poly[K_FP];
				DCRT_tmp.SetFormat(Format::COEFFICIENT);
				DivideRoundingSelf(params, &DCRT_tmp, Format::COEFFICIENT);
				b = DCRT_tmp.GetElementAtIndex(0)[0];
			} else { // i == 1
				for (uint32_t j = 0; j < K_FP; j++) {
					In_Poly[j].SetFormat(Format::COEFFICIENT);
					DivideRoundingSelf(params, &In_Poly[j], Format::COEFFICIENT);
					A[j] = In_Poly[j].GetElementAtIndex(0).GetValues();
				}
				In_Poly[K_FP].SetFormat(Format::COEFFICIENT);
				DivideRoundingSelf(params, &In_Poly[K_FP], Format::COEFFICIENT);
				b = In_Poly[K_FP].GetElementAtIndex(0)[0];
			}
				
			
			// First Carry, 4+1 (1 bit no information)bit bootstrapping Bootstrapping
      if (i != 2) {

				if (i == 0) {
					start_bits  = scale_bit  + 4 + 1;
					prefix_add_vals = (prefix >> 1);
				} else if (i == 1) {
					start_bits  = scale_bit  + 8 + 1;
					prefix_add_vals = (prefix << 3);
			
				}
				carry_bit = 1;
				out_nums = 2;
			
				//First Bootstrapping,3 bit + random 1 bit is bootstrapping
				BS_out_vec = this->Bootstrap_and_KS(params
					, LWEscheme
					, params->GetACCRelu(0)
					, EK
					, A
					, b, 0, start_bits, prefix_add_vals, carry_bit, out_nums);
				DCRTPoly ACCs = params -> GetACCRelu(1+i);
		 
				//else if (i == 1) {
				//	ACCs = params -> GetACCRelu(2);
				//}
				carry_bit = 1;
				out_nums = 2;
	
				// Orign, 4+1(padding but no infomation) bit bootstrapping msg <= 3bit 
				start_bits  = Q1_WO_bits + 4 + 1;
				prefix_add_vals = ( Q1_Int >> 1);  
				// ================   EVAL MODE SAVING ==================== 
				BS_out_vec = this->ProductAndBKAndKS(params
					, LWEscheme
					, ACCs, EK
					, (*BS_out_vec)[0]
					, (*BS_out_vec)[1]
					, 0, start_bits, prefix_add_vals, carry_bit, out_nums);
			
			} else if (i == 2) {
				start_bits  = scale_bit  + 12 + 1;
				prefix_add_vals = (prefix << 7);  
				carry_bit = 0;
				out_nums = 1;
				//First Bootstrapping,3 bit + random 1 bit is bootstrapping
				Outs= std::move((*this->Bootstrap_and_KS(params
					, LWEscheme
					, params->GetACCRelu(18)
					, EK
					, A
					, b, 0, start_bits, prefix_add_vals, carry_bit, out_nums))[0]);
			}


			// After Handle
			//std::cout << "handle before " << std::endl;
			if (i != 2) {
				
				// remove 4 or 8 bit info
				for (uint32_t j = 0; j < K_FP+1; j++) {
					In_Poly[j]-=(*BS_out_vec)[1][j]; 
				}
				Res_Poly[i] = std::move((*BS_out_vec)[0]);
			} 	
		}
		return std::make_shared<vector<DCRTPoly>>(Outs);
	}







  std::shared_ptr<vector<DCRTPoly>> FPRingGSWAccumulatorScheme::BootstrapReluRevF(
    const std::shared_ptr<FPRingGSWCryptoParams> params, 
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<DCRTPoly> & ct,
		const uint64_t &mins,
		const std::shared_ptr<FPRingGSWEvalKey> EK) const {
		if (mins > 13) {
			PALISADE_THROW(config_error, "Curr version val should be below 13");
		} 

		auto MulPrepare_start = std::chrono::system_clock::now();
    uint32_t K_FP = params->GetRLWEParams()->GetK();
		NativeInteger Q0 = params->GetModuliQ()[0];
    NativeInteger Q1 = params->GetModuliQ()[1];
		uint64_t Q1_Int = Q1.ConvertToInt();
		uint32_t Q1_WO_bits = params->GetRLWEParams()->GetBits()[1];  // <32 bit
		uint32_t start_bits;
		uint64_t prefix_add_vals;
		uint32_t carry_bit;
		uint32_t out_nums;
	
		std::shared_ptr<vector<vector<DCRTPoly>>> BS_out_vec;
		vector<DCRTPoly> BS_out;
		vector<NativeVector> A(K_FP);
		NativeInteger b;
		DCRTPoly DCRT_tmp;
		vector<DCRTPoly> In_Poly(K_FP+1);
	
		vector<vector<DCRTPoly>> Res_Poly(2);
		vector<DCRTPoly> Outs(K_FP+1);
		
		// ct1 - ct 2
		for (uint32_t i = 0; i < K_FP+1; i++) {	
			In_Poly[i] =  ct[i];
		}

		// Substract gives range - 11bit~ 12 bit
		
		uint32_t scale_bit = params->GetExpBit();	
		uint64_t prefix = params->GetRLWEParams()->GetScaling().ConvertToInt();
			prefix <<= scale_bit;
		uint64_t bit_check = prefix;	
		bit_check <<= 7;  // Check 8 bit	
		
		// Make
		In_Poly[K_FP].SetFormat(Format::COEFFICIENT);
		In_Poly[K_FP].GetElementW(0)[0].ModAddEq(
				NativeInteger(bit_check).ModMul(Q1, Q0),Q0);
		In_Poly[K_FP].GetElementW(1)[0].ModAddEq(
				NativeInteger(bit_check).ModMul(Q1, Q1),Q1);
		
		// Sub Values
		In_Poly[K_FP].GetElementW(0)[0].ModSubEq(
				NativeInteger(prefix).ModMul(Q1, Q0).ModMul(mins,Q0),Q0);
		In_Poly[K_FP].GetElementW(1)[0].ModSubEq(
				NativeInteger(prefix).ModMul(Q1, Q1).ModMul(mins,Q1),Q1);
		In_Poly[K_FP].SetFormat(Format::EVALUATION);
		
		////////////////////////////////   STARTS !!!! /////////////////////
		std::chrono::duration<double> MulPrepare_end = std::chrono::system_clock::now() - MulPrepare_start;
		if (DEBUG) {
			std::cout << "bootstrapping Relu prepare time is " << MulPrepare_end.count() << std::endl ;
		}

		// Bootstrap 1
		for (uint32_t i = 0; i < 2; i ++) { 
			
			//std::cout << "i is " << i << std::endl;	
			// Make Polys
			if (i != 1) {
				// Copy A
				for (uint32_t j = 0; j < K_FP; j++) {
					DCRT_tmp = In_Poly[j];
					DCRT_tmp.SetFormat(Format::COEFFICIENT);
					DivideRoundingSelf(params, &DCRT_tmp, Format::COEFFICIENT);
					A[j] = DCRT_tmp.GetElementAtIndex(0).GetValues();
				}
				// Copy B
				DCRT_tmp = In_Poly[K_FP];
				DCRT_tmp.SetFormat(Format::COEFFICIENT);
				DivideRoundingSelf(params, &DCRT_tmp, Format::COEFFICIENT);
				b = DCRT_tmp.GetElementAtIndex(0)[0];
			} else { // i == 1
				for (uint32_t j = 0; j < K_FP; j++) {
					In_Poly[j].SetFormat(Format::COEFFICIENT);
					DivideRoundingSelf(params, &In_Poly[j], Format::COEFFICIENT);
					A[j] = In_Poly[j].GetElementAtIndex(0).GetValues();
				}
				In_Poly[K_FP].SetFormat(Format::COEFFICIENT);
				DivideRoundingSelf(params, &In_Poly[K_FP], Format::COEFFICIENT);
				b = In_Poly[K_FP].GetElementAtIndex(0)[0];
			}
				
			
			// First Carry, 4+1 (1 bit no information)bit bootstrapping Bootstrapping
      if (i != 1) {

				if (i == 0) {
					start_bits  = scale_bit  + 4 + 1;
					prefix_add_vals = (prefix >> 1);
				} 
				//else if(i == 1) {
				//	start_bits  = scale_bit  + 8 + 1;
				//	prefix_add_vals = (prefix << 3);  
				//}
				carry_bit = 1;
				out_nums = 2;
			
				//First Bootstrapping,3 bit + random 1 bit is bootstrapping
				BS_out_vec = this->Bootstrap_and_KS(params
					, LWEscheme
					, params->GetACCRelu(0)
					, EK
					, A
					, b, 0, start_bits, prefix_add_vals, carry_bit, out_nums);
				DCRTPoly ACCs = params -> GetACCRelu(1+i);
		 
				//else if (i == 1) {
				//	ACCs = params -> GetACCRelu(2);
				//}
				carry_bit = 1;
				out_nums = 2;
	
				// Orign, 4+1(padding but no infomation) bit bootstrapping msg <= 3bit 
				start_bits  = Q1_WO_bits + 4 + 1;
				prefix_add_vals = ( Q1_Int >> 1);  
				// ================   EVAL MODE SAVING ==================== 
				BS_out_vec = this->ProductAndBKAndKS(params
					, LWEscheme
					, ACCs, EK
					, (*BS_out_vec)[0]
					, (*BS_out_vec)[1]
					, 0, start_bits, prefix_add_vals, carry_bit, out_nums);
			
			} else if (i == 1) {
				start_bits  = scale_bit  + 8 + 1;
				prefix_add_vals = (prefix << 3);  
				carry_bit = 1;
				out_nums = 2;
				//First Bootstrapping,3 bit + random 1 bit is bootstrapping
				BS_out_vec = this->Bootstrap_and_KS(params
					, LWEscheme
					, params->GetACCRelu(8)
					, EK
					, A
					, b, 0, start_bits, prefix_add_vals, carry_bit, out_nums);
			}


			// After Handle
			//std::cout << "handle before " << std::endl;
			if (i != 1) {
				
				// remove 4 or 8 bit info
				for (uint32_t j = 0; j < K_FP+1; j++) {
					In_Poly[j]-=(*BS_out_vec)[1][j]; 
				}

				//std::cout << "Minus is end" << std::endl;
				// Saves results 
				Res_Poly[i] = std::move((*BS_out_vec)[0]);
				//std::cout << "Res is moved" << std::endl;
				
			} else {// i == 1
				
				// first 4 bit
				start_bits  = Q1_WO_bits + 4 + 1;
				carry_bit = 0;
				out_nums = 1;
				prefix_add_vals = ( Q1_Int >> 1);  
				
				Outs =  std::move((*this->ProductAndBKAndKS(params
        , LWEscheme
				, params->GetACCRelu(4), EK
        , Res_Poly[0]
        , (*BS_out_vec)[0]
        , 0, start_bits, prefix_add_vals, carry_bit, out_nums))[0]);
       	/*
				// Second 4 bit
				start_bits  = Q1_WO_bits + 4 + 1;
				carry_bit = 0;
				out_nums = 1;
				prefix_add_vals = ( Q1_Int >> 1);  
				BS_out =  std::move((*this->ProductAndBKAndKS(params
        , LWEscheme
				, params->GetACCRelu(5), EK
        , Res_Poly[1]
        , (*BS_out_vec)[0]
        , 0, start_bits, prefix_add_vals, carry_bit, out_nums))[0]);
				*/
				// Some Seconds & Third 
				for (uint32_t j = 0; j < K_FP+1; j++) {
					//Outs[j] += BS_out[j];
					Outs[j] += (*BS_out_vec)[1][j];
				}
			}
		}
		
		Outs[K_FP].SetFormat(Format::COEFFICIENT);
		// Add Values

		Outs[K_FP].GetElementW(0)[0].ModAddEq(
				NativeInteger(prefix).ModMul(Q1, Q0).ModMul(mins,Q0),Q0);
		Outs[K_FP].GetElementW(1)[0].ModAddEq(
				NativeInteger(prefix).ModMul(Q1, Q1).ModMul(mins,Q1),Q1);
		
		Outs[K_FP].SetFormat(Format::EVALUATION);
	
		return std::make_shared<vector<DCRTPoly>>(Outs);
	}




  std::shared_ptr<vector<DCRTPoly>> FPRingGSWAccumulatorScheme::BootstrapReluRev(
    const std::shared_ptr<FPRingGSWCryptoParams> params, 
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<DCRTPoly> & ct,
		const uint64_t &mins,
		const std::shared_ptr<FPRingGSWEvalKey> EK) const {
	
		if (mins > 27) {
			PALISADE_THROW(config_error, "Curr version val should be below 27");
		} 

		auto MulPrepare_start = std::chrono::system_clock::now();
    uint32_t K_FP = params->GetRLWEParams()->GetK();
		NativeInteger Q0 = params->GetModuliQ()[0];
    NativeInteger Q1 = params->GetModuliQ()[1];
		uint64_t Q1_Int = Q1.ConvertToInt();
		uint32_t Q1_WO_bits = params->GetRLWEParams()->GetBits()[1];  // <32 bit
		uint32_t start_bits;
		uint64_t prefix_add_vals;
		uint32_t carry_bit;
		uint32_t out_nums;
	
		std::shared_ptr<vector<vector<DCRTPoly>>> BS_out_vec;
		vector<DCRTPoly> BS_out;
		vector<NativeVector> A(K_FP);
		NativeInteger b;
		DCRTPoly DCRT_tmp;
		vector<DCRTPoly> In_Poly(K_FP+1);
	
		vector<vector<DCRTPoly>> Res_Poly(2);
		vector<DCRTPoly> Outs(K_FP+1);
		
		// ct1 - ct 2
		for (uint32_t i = 0; i < K_FP+1; i++) {	
			In_Poly[i] =  ct[i];
		}

		// Substract gives range - 11bit~ 12 bit
		// add 13 bit, bootstrapping, and check whether 13bit is 0 or not
		uint32_t scale_bit = params->GetExpBit();	
		uint64_t prefix = params->GetRLWEParams()->GetScaling().ConvertToInt();
			prefix <<= scale_bit;
		uint64_t bit_check = prefix;	
		bit_check <<= 10;  // 2^10 ? 2^ 11?	
		
		
		// Make

		In_Poly[K_FP].SetFormat(Format::COEFFICIENT);
		In_Poly[K_FP].GetElementW(0)[0].ModAddEq(
				NativeInteger(bit_check).ModMul(Q1, Q0),Q0);
		In_Poly[K_FP].GetElementW(1)[0].ModAddEq(
				NativeInteger(bit_check).ModMul(Q1, Q1),Q1);
		// Sub Values
		//
		In_Poly[K_FP].GetElementW(0)[0].ModSubEq(
				NativeInteger(prefix).ModMul(Q1, Q0).ModMul(mins,Q0),Q0);
		In_Poly[K_FP].GetElementW(1)[0].ModSubEq(
				NativeInteger(prefix).ModMul(Q1, Q1).ModMul(mins,Q1),Q1);
		In_Poly[K_FP].SetFormat(Format::EVALUATION);
		

		////////////////////////////////   STARTS !!!! /////////////////////
		std::chrono::duration<double> MulPrepare_end = std::chrono::system_clock::now() - MulPrepare_start;
		if (DEBUG) {
			std::cout << "bootstrapping Relu prepare time is " << MulPrepare_end.count() << std::endl ;
		}


		// Bootstrap 1
		for (uint32_t i = 0; i < 3; i ++) { 
			
			//std::cout << "i is " << i << std::endl;	
			// Make Polys
			if (i != 2) {
				// Copy A
				for (uint32_t j = 0; j < K_FP; j++) {
					DCRT_tmp = In_Poly[j];
					DCRT_tmp.SetFormat(Format::COEFFICIENT);
					DivideRoundingSelf(params, &DCRT_tmp, Format::COEFFICIENT);
					A[j] = DCRT_tmp.GetElementAtIndex(0).GetValues();
				}
				// Copy B
				DCRT_tmp = In_Poly[K_FP];
				DCRT_tmp.SetFormat(Format::COEFFICIENT);
				DivideRoundingSelf(params, &DCRT_tmp, Format::COEFFICIENT);
				b = DCRT_tmp.GetElementAtIndex(0)[0];
			} else { // i == 2
				for (uint32_t j = 0; j < K_FP; j++) {
					In_Poly[j].SetFormat(Format::COEFFICIENT);
					DivideRoundingSelf(params, &In_Poly[j], Format::COEFFICIENT);
					A[j] = In_Poly[j].GetElementAtIndex(0).GetValues();
				}
				In_Poly[K_FP].SetFormat(Format::COEFFICIENT);
				DivideRoundingSelf(params, &In_Poly[K_FP], Format::COEFFICIENT);
				b = In_Poly[K_FP].GetElementAtIndex(0)[0];
			}
				
			// First Carry, 4+1 (1 bit no information)bit bootstrapping Bootstrapping
      if (i != 2) {
				if (i == 0) {
					start_bits  = scale_bit  + 4 + 1;
					prefix_add_vals = (prefix >> 1);  
				} else if (i == 1) {
					start_bits  = scale_bit  + 8 + 1;
					prefix_add_vals = (prefix << 3);  
				} 
				carry_bit = 1;
				out_nums = 2;
			
				//First Bootstrapping,3 bit + random 1 bit is bootstrapping
				BS_out_vec = this->Bootstrap_and_KS(params
					, LWEscheme
					, params->GetACCRelu(0)
					, EK
					, A
					, b, 0, start_bits, prefix_add_vals, carry_bit, out_nums);
		
				//Second Bootsrapp
				// i = 1, 2, 3 ...
				DCRTPoly ACCs = params -> GetACCRelu(i+1);
				//	ACCs = params -> GetACCRelu(6);
				
				// Orign, 4+1(padding but no infomation) bit bootstrapping msg <= 3bit 
				start_bits  = Q1_WO_bits + 4 + 1;
				carry_bit = 1;
				out_nums = 2;
				prefix_add_vals = ( Q1_Int >> 1);  
				// ================   EVAL MODE SAVING ==================== 
				BS_out_vec = this->ProductAndBKAndKS(params
					, LWEscheme
					, ACCs, EK
					, (*BS_out_vec)[0]
					, (*BS_out_vec)[1]
					, 0, start_bits, prefix_add_vals, carry_bit, out_nums);
			} else if (i == 2) {
				start_bits  = scale_bit  + 12 + 1;
				prefix_add_vals = (prefix << 7);  
				carry_bit = 1;
				out_nums = 2;
			
				//First Bootstrapping,3 bit + random 1 bit is bootstrapping
				BS_out_vec = this->Bootstrap_and_KS(params
					, LWEscheme
					, params->GetACCRelu(6)
					, EK
					, A
					, b, 0, start_bits, prefix_add_vals, carry_bit, out_nums);
			}


			// After Handle
			//std::cout << "handle before " << std::endl;
			if (i != 2) {
				
				// remove 4 or 8 bit info
				for (uint32_t j = 0; j < K_FP+1; j++) {
					In_Poly[j]-=(*BS_out_vec)[1][j]; 
				}

				//std::cout << "Minus is end" << std::endl;
				// Saves results 
				Res_Poly[i] = std::move((*BS_out_vec)[0]);
				//std::cout << "Res is moved" << std::endl;
				
			} else {// i == 2
				
				// first 4 bit
				start_bits  = Q1_WO_bits + 4 + 1;
				carry_bit = 0;
				out_nums = 1;
				prefix_add_vals = ( Q1_Int >> 1);  
				
				Outs =  std::move((*this->ProductAndBKAndKS(params
        , LWEscheme
				, params->GetACCRelu(4), EK
        , Res_Poly[0]
        , (*BS_out_vec)[1]
        , 0, start_bits, prefix_add_vals, carry_bit, out_nums))[0]);
        
				// Second 4 bit
				start_bits  = Q1_WO_bits + 4 + 1;
				carry_bit = 0;
				out_nums = 1;
				prefix_add_vals = ( Q1_Int >> 1);  
				BS_out =  std::move((*this->ProductAndBKAndKS(params
        , LWEscheme
				, params->GetACCRelu(5), EK
        , Res_Poly[1]
        , (*BS_out_vec)[1]
        , 0, start_bits, prefix_add_vals, carry_bit, out_nums))[0]);
        
				// Some Seconds & Third 
				
				for (uint32_t j = 0; j < K_FP+1; j++) {
					Outs[j] += BS_out[j];
					Outs[j] += (*BS_out_vec)[0][j];
				}
				
			}
		}
		
		Outs[K_FP].SetFormat(Format::COEFFICIENT);
		// Add Values
		
		Outs[K_FP].GetElementW(0)[0].ModAddEq(
				NativeInteger(prefix).ModMul(Q1, Q0).ModMul(mins,Q0),Q0);
		Outs[K_FP].GetElementW(1)[0].ModAddEq(
				NativeInteger(prefix).ModMul(Q1, Q1).ModMul(mins,Q1),Q1);
		//Outs[K_FP].SetFormat(Format::EVALUATION);
		// Add vals
		
		/*
		Outs[K_FP].GetElementW(0)[0].ModSubEq(
			NativeInteger(bit_check).ModMul(Q1, Q0),Q0);
		Outs[K_FP].GetElementW(1)[0].ModSubEq(
			NativeInteger(bit_check).ModMul(Q1, Q1),Q1);
		*/

		/*
		Outs[K_FP].GetElementW(0)[0].ModSubEq(
				NativeInteger(prefix).ModMul(Q1, Q0).ModMul(mins,Q0),Q0);
		Outs[K_FP].GetElementW(1)[0].ModSubEq(
				NativeInteger(prefix).ModMul(Q1, Q1).ModMul(mins,Q1),Q1);
		*/



		Outs[K_FP].SetFormat(Format::EVALUATION);
	
		
		
		return std::make_shared<vector<DCRTPoly>>(Outs);
	}







  std::shared_ptr<vector<DCRTPoly>> FPRingGSWAccumulatorScheme::BootstrapMulPolyMake(
    const std::shared_ptr<FPRingGSWCryptoParams> params, 
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<DCRTPoly> & ct_Exp, 
		const std::shared_ptr<FPRingGSWEvalKey> EK,
		const uint32_t max_msg,
		bool mul = true) const {
		auto MulPrepare_start = std::chrono::system_clock::now();
    uint32_t K_FP = params->GetRLWEParams()->GetK();
		NativeInteger Q0 = params->GetModuliQ()[0];
    NativeInteger Q1 = params->GetModuliQ()[1];
		uint64_t Q1_Int = Q1.ConvertToInt();
		uint32_t Q1_WO_bits = params->GetRLWEParams()->GetBits()[1];  // <32 bit
		uint32_t start_bits;
		uint64_t prefix_add_vals;
		uint32_t carry_bit;
		uint32_t out_nums;

		std::shared_ptr<vector<vector<DCRTPoly>>> BS_out_vec;
		vector<DCRTPoly> BS_out(K_FP+1);
		vector<DCRTPoly> Out_tmp(K_FP+1);
		vector<DCRTPoly> BS_out_tmp(K_FP+1);
		vector<NativeVector> A(K_FP);
		NativeInteger b;
		DCRTPoly DCRT_tmp;
		vector<DCRTPoly> In_Poly_Exp(K_FP+1);
	
		// Rotate
		for (uint32_t i = 0; i < K_FP+1; i++) {
			In_Poly_Exp[i] =  ct_Exp[i];
		}
		
		vector<vector<DCRTPoly>> Res_Poly(2);
		vector<DCRTPoly> Outs(K_FP+1);
	

		// Substract gives range - 11bit~ 12 bit
		// add 13 bit, bootstrapping, and check whether 13bit is 0 or not
		uint32_t scale_bit = params->GetExpBit();	
		uint64_t Exp_prefix = params->GetRLWEParams()->GetScaling().ConvertToInt();
			Exp_prefix <<= scale_bit;
		
			////////////////////////////////   STARTS !!!! /////////////////////
		std::chrono::duration<double> MulPrepare_end = std::chrono::system_clock::now() - MulPrepare_start;
		if (DEBUG) std::cout << "bootstrapping Div prepare time is " << MulPrepare_end.count() << std::endl ;
		
		// Bootstrap First
		for (uint32_t j = 0; j < K_FP; j++) {
			DCRT_tmp = ct_Exp[j];
			DCRT_tmp.SetFormat(Format::COEFFICIENT);
			DivideRoundingSelf(params, &DCRT_tmp, Format::COEFFICIENT);
			A[j] = DCRT_tmp.GetElementAtIndex(0).GetValues();
		}
		// Copy B
		DCRT_tmp = ct_Exp[K_FP];
		DCRT_tmp.SetFormat(Format::COEFFICIENT);
		DivideRoundingSelf(params, &DCRT_tmp, Format::COEFFICIENT);
		b = DCRT_tmp.GetElementAtIndex(0)[0];
			
		// First Carry, 4+1 (1 bit no information)bit bootstrapping Bootstrapping
  	start_bits  = scale_bit  + 3 + 1;
		prefix_add_vals = (Exp_prefix >> 1);  
    carry_bit = 2;
    out_nums = 4;
		
		BS_out_vec = this->Bootstrap_and_KS(params
			, LWEscheme
			, params->GetACCDiv(0)
			, EK
			, A
			, b, 0, start_bits, prefix_add_vals, carry_bit, out_nums);

    //First Bootstrapping,3 bit + random 1 bit is bootstrapping
		// Sign * rev
	  start_bits  = Q1_WO_bits  + 2 + 1;
		//prefix_add_vals = (Q1_Int >> 1) + Q1_Int;  
  	prefix_add_vals = (Q1_Int >> 1);  
    carry_bit = 0;
    out_nums = 1;
	  // ================   EVAL MODE SAVING ==================== 
    (*BS_out_vec)[0] = std::move((*this->ProductAndBKAndKS(params
			, LWEscheme
			, params->GetACCDiv(5), EK
  		//, params->GetACCDiv(2), EK
      , (*BS_out_vec)[0]
      , (*BS_out_vec)[3]
      , 0, start_bits, prefix_add_vals, carry_bit, out_nums))[0]);
	
    (*BS_out_vec)[1] = std::move((*this->ProductAndBKAndKS(params
			, LWEscheme
			, params->GetACCDiv(5), EK
  		//, params->GetACCDiv(2), EK
      , (*BS_out_vec)[1]
      , (*BS_out_vec)[3]
      , 0, start_bits, prefix_add_vals, carry_bit, out_nums))[0]);
	  (*BS_out_vec)[2] = std::move((*this->ProductAndBKAndKS(params
			, LWEscheme
			, params->GetACCDiv(5), EK
  		//, params->GetACCDiv(2), EK
      , (*BS_out_vec)[2]
      , (*BS_out_vec)[3]
      , 0, start_bits, prefix_add_vals, carry_bit, out_nums))[0]);
	

		// Should be fix it
		// Make (1 or X) Cipihertext
		for (uint32_t i = 0; i < K_FP+1; i++) {
				if (i == K_FP) {
					BS_out_tmp[i] = (*BS_out_vec)[0][i];
					BS_out_tmp[i].SetFormat(Format::COEFFICIENT);
					BS_out_tmp[i].GetElementW(0)[0].ModSubEq(Q1,Q0);
					BS_out_tmp[i].SetFormat(Format::EVALUATION);
					BS_out[i] = ((*BS_out_vec)[0][i] * params->GetRotMonomial_RLWE(1)) - (BS_out_tmp[i]);
				} else {
					BS_out[i] = ((*BS_out_vec)[0][i] * params->GetRotMonomial_RLWE(1)) - ((*BS_out_vec)[0][i]);
				}
		} 
	

		//std::cout << "ends" << std::endl;
		/*
		// Blind Rotate 1,
		for (uint32_t i = 0; i < K_FP+1; i++) {
				BS_out[i] = In_Poly_Sign[i];
				BS_out[i] *= params->GetRotMonomial_RLWE(1);
				BS_out[i] -= In_Poly_Sign[i];
		} 
		// Get our Sign
	  start_bits  = Q1_WO_bits  + 2 + 1;
		prefix_add_vals = (Q1_Int >> 1) + Q1_Int;  
    carry_bit = 0;
    out_nums = 1;
	  // ================   EVAL MODE SAVING ==================== 
    BS_out = std::move(*this->ProductAndBKAndKSRange(params
			, LWEscheme
			, params->GetACCDiv(2), EK
      , (*BS_out_vec)[0]
      , BS_out
      , 0,1, start_bits, prefix_add_vals, carry_bit, out_nums));
		
		// Adding
		for (uint32_t i = 0; i < K_FP+1; i++) {
			BS_out[i] += In_Poly_Sign[i];
		}
		// BS_out get current vals
		*/


		// Blind Rotation 2,
		for (uint32_t i = 0; i < K_FP+1; i++) {
				BS_out_tmp[i] = BS_out[i];
				BS_out_tmp[i] *= params->GetRotMonomial_RLWE(2);
				BS_out_tmp[i] -= BS_out[i];
		} 

		// Second Construction
		start_bits  = Q1_WO_bits  + 2 + 1;
		prefix_add_vals = (Q1_Int >> 1);  
    carry_bit = 0;
    out_nums = 1;
		BS_out_tmp = std::move(*this->ProductAndBKAndKSRange(params
			, LWEscheme
			, params->GetACCDiv(5), EK
      , BS_out_tmp
      , (*BS_out_vec)[1]
      , 0,3, start_bits, prefix_add_vals, carry_bit, out_nums));
		

	
		// Adding
		for (uint32_t i = 0; i < K_FP+1; i++) {
			BS_out[i] += BS_out_tmp[i];
		}

		// Blind Rotation 3,
		for (uint32_t i = 0; i < K_FP+1; i++) {
				BS_out_tmp[i] = BS_out[i];
				BS_out_tmp[i] *= params->GetRotMonomial_RLWE(4);
				BS_out_tmp[i] -= BS_out[i];
		}

		// Third Construction
		start_bits  = Q1_WO_bits  + 2 + 1;
		prefix_add_vals = (Q1_Int >> 1);  
    carry_bit = 0;
    out_nums = 1;

		BS_out_tmp = std::move(*this->ProductAndBKAndKSRange(params
			, LWEscheme
			, params->GetACCDiv(5), EK
      , BS_out_tmp
      , (*BS_out_vec)[2]
      , 0,7, start_bits, prefix_add_vals, carry_bit, out_nums));
	
		// Adding
		for (uint32_t i = 0; i < K_FP+1; i++) {
			BS_out[i] += BS_out_tmp[i];
		}

		/*
		for (uint32_t i = 0; i < K_FP+1; i++) {
			Out_tmp[i] = BS_out[i];
		}
	*/


		// Remove 3 bit
		start_bits  = scale_bit  + 3 + 1;
		prefix_add_vals = (Exp_prefix >> 1);  
		carry_bit = 1;
    out_nums = 2;
	
		BS_out_vec = this->Bootstrap_and_KS(params
			, LWEscheme
			, params->GetACCDiv(3)
			, EK
			, A
			, b, 0, start_bits, prefix_add_vals, carry_bit, out_nums);
		
		start_bits  = Q1_WO_bits  + 3 + 1;
		prefix_add_vals = (Q1_Int >> 1);  
 		carry_bit = 0;
    out_nums = 1;
	
		BS_out_tmp =  std::move((*this->ProductAndBKAndKS(params
			, LWEscheme
			, params->GetACCDiv(4), EK
      , (*BS_out_vec)[0]
      , (*BS_out_vec)[1]
      , 0, start_bits, prefix_add_vals, carry_bit, out_nums))[0]);


		//std::cout << "Hello" << std::endl;

		// Removes
		
		for (uint32_t i = 0; i < K_FP+1; i++) {
			In_Poly_Exp[i] -=  BS_out_tmp[i];
		}
	
		// Seconds - one
		for (uint32_t j = 0; j < K_FP; j++) {
			DCRT_tmp = In_Poly_Exp[j];
			DCRT_tmp.SetFormat(Format::COEFFICIENT);
			DivideRoundingSelf(params, &DCRT_tmp, Format::COEFFICIENT);
			A[j] = DCRT_tmp.GetElementAtIndex(0).GetValues();
		}
		// Copy B
		DCRT_tmp = In_Poly_Exp[K_FP];
		DCRT_tmp.SetFormat(Format::COEFFICIENT);
		DivideRoundingSelf(params, &DCRT_tmp, Format::COEFFICIENT);
		b = DCRT_tmp.GetElementAtIndex(0)[0];
	


		//std::cout << "Hello2" << std::endl;

		// First Carry, 4+1 (1 bit no information)bit bootstrapping Bootstrapping
  	start_bits  = (scale_bit + 3) + 2 + 1;
		prefix_add_vals = (Exp_prefix << 2);  // -1bit <- 3bit = 2bit  
    carry_bit = 2;
    out_nums = 4;
		
		BS_out_vec = this->Bootstrap_and_KS(params
			, LWEscheme
			, params->GetACCDiv(1)
			, EK
			, A
			, b, 0, start_bits, prefix_add_vals, carry_bit, out_nums);

		//BS rot 3
		for (uint32_t i = 0; i < K_FP+1; i++) {
			if (max_msg >= 8) {
				(*BS_out_vec)[1][i] *= params->GetRotMonomial_RLWE(8); 
				(*BS_out_vec)[0][i] += (*BS_out_vec)[1][i];
			}
			if (max_msg >= 16) {
				(*BS_out_vec)[2][i] *= params->GetRotMonomial_RLWE(16); 
				(*BS_out_vec)[0][i] += (*BS_out_vec)[2][i];
			}
			if (max_msg >= 24) {
				(*BS_out_vec)[3][i] *= params->GetRotMonomial_RLWE(24); 	
				(*BS_out_vec)[0][i] += (*BS_out_vec)[3][i];
			}
			/*
			(*BS_out_vec)[1][i] *= params->GetRotMonomial_RLWE(8); 
			(*BS_out_vec)[2][i] *= params->GetRotMonomial_RLWE(16); 
			(*BS_out_vec)[3][i] *= params->GetRotMonomial_RLWE(24); 	
			// Adding
			(*BS_out_vec)[0][i] += (*BS_out_vec)[1][i];
			(*BS_out_vec)[0][i] += (*BS_out_vec)[2][i];
			(*BS_out_vec)[0][i] += (*BS_out_vec)[3][i];
			*/
		} 
			
		start_bits  = Q1_WO_bits  + 2 + 1;
		prefix_add_vals = (Q1_Int >> 1);  
    carry_bit = 0;
    out_nums = 1;
		
		BS_out = std::move(*this->ProductAndBKAndKSRange(params
			, LWEscheme
			, params->GetACCDiv(5), EK
      , BS_out
      , (*BS_out_vec)[0]
      , 0,max_msg, start_bits, prefix_add_vals, carry_bit, out_nums, mul));
		
		return std::make_shared<vector<DCRTPoly>>(BS_out);
	}




  std::shared_ptr<vector<DCRTPoly>> FPRingGSWAccumulatorScheme::BootstrapMulPolyMake0(
    const std::shared_ptr<FPRingGSWCryptoParams> params, 
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<DCRTPoly> & ct_ASign,
		const DCRTPoly &ct_BSign,
		const vector<DCRTPoly> & ct_Exp, 
		const std::shared_ptr<FPRingGSWEvalKey> EK,
		const uint32_t max_msg,
		bool mul = true) const {
		auto MulPrepare_start = std::chrono::system_clock::now();
    uint32_t K_FP = params->GetRLWEParams()->GetK();
		//NativeInteger Q0 = params->GetModuliQ()[0];
    NativeInteger Q1 = params->GetModuliQ()[1];
		uint64_t Q1_Int = Q1.ConvertToInt();
		uint32_t Q1_WO_bits = params->GetRLWEParams()->GetBits()[1];  // <32 bit
		uint32_t start_bits;
		uint64_t prefix_add_vals;
		uint32_t carry_bit;
		uint32_t out_nums;

		std::shared_ptr<vector<vector<DCRTPoly>>> BS_out_vec;
		vector<DCRTPoly> BS_out(K_FP+1);
		vector<DCRTPoly> Out_tmp(K_FP+1);
		vector<DCRTPoly> BS_out_tmp(K_FP+1);
		vector<NativeVector> A(K_FP);
		NativeInteger b;
		DCRTPoly DCRT_tmp;
		vector<DCRTPoly> In_Poly_Sign(K_FP+1);
		vector<DCRTPoly> In_Poly_Exp(K_FP+1);
	
		vector<vector<DCRTPoly>> Res_Poly(2);
		vector<DCRTPoly> Outs(K_FP+1);
		
		// Rotate
		for (uint32_t i = 0; i < K_FP+1; i++) {
			if (i == K_FP) {
				In_Poly_Sign[i] =  ct_BSign;
				
			} else {
				In_Poly_Sign[i] =  ct_ASign[i];
			}
			In_Poly_Exp[i] =  ct_Exp[i];
		
		}
		
		// Substract gives range - 11bit~ 12 bit
		// add 13 bit, bootstrapping, and check whether 13bit is 0 or not
		uint32_t scale_bit = params->GetExpBit();	
		uint64_t prefix = params->GetRLWEParams()->GetScaling().ConvertToInt();
			prefix <<= scale_bit;
		
			////////////////////////////////   STARTS !!!! /////////////////////
		std::chrono::duration<double> MulPrepare_end = std::chrono::system_clock::now() - MulPrepare_start;
		if (DEBUG) std::cout << "bootstrapping Div prepare time is " << MulPrepare_end.count() << std::endl ;
		
		// Bootstrap First
		for (uint32_t j = 0; j < K_FP; j++) {
			DCRT_tmp = ct_Exp[j];
			DCRT_tmp.SetFormat(Format::COEFFICIENT);
			DivideRoundingSelf(params, &DCRT_tmp, Format::COEFFICIENT);
			A[j] = DCRT_tmp.GetElementAtIndex(0).GetValues();
		}
		// Copy B
		DCRT_tmp = ct_Exp[K_FP];
		DCRT_tmp.SetFormat(Format::COEFFICIENT);
		DivideRoundingSelf(params, &DCRT_tmp, Format::COEFFICIENT);
		b = DCRT_tmp.GetElementAtIndex(0)[0];
			
		// First Carry, 4+1 (1 bit no information)bit bootstrapping Bootstrapping
  	start_bits  = scale_bit  + 3 + 1;
		prefix_add_vals = (prefix >> 1);  
    carry_bit = 2;
    out_nums = 4;
		
		BS_out_vec = this->Bootstrap_and_KS(params
			, LWEscheme
			, params->GetACCDiv(0)
			, EK
			, A
			, b, 0, start_bits, prefix_add_vals, carry_bit, out_nums);

    //First Bootstrapping,3 bit + random 1 bit is bootstrapping
		start_bits  = Q1_WO_bits  + 2 + 1;
		//prefix_add_vals = (Q1_Int >> 1) + Q1_Int;  
  	prefix_add_vals = (Q1_Int >> 1);  
    carry_bit = 0;
    out_nums = 1;
	  // ================   EVAL MODE SAVING ==================== 
    (*BS_out_vec)[0] = std::move((*this->ProductAndBKAndKS(params
			, LWEscheme
			//, params->GetACCDiv(2), EK
  		, params->GetACCDiv(5), EK
      , (*BS_out_vec)[0]
      , (*BS_out_vec)[3]
      , 0, start_bits, prefix_add_vals, carry_bit, out_nums))[0]);
	
    (*BS_out_vec)[1] = std::move((*this->ProductAndBKAndKS(params
			, LWEscheme
			//, params->GetACCDiv(2), EK
  		, params->GetACCDiv(5), EK
      , (*BS_out_vec)[1]
      , (*BS_out_vec)[3]
      , 0, start_bits, prefix_add_vals, carry_bit, out_nums))[0]);
	  (*BS_out_vec)[2] = std::move((*this->ProductAndBKAndKS(params
			, LWEscheme
			//, params->GetACCDiv(2), EK
  		, params->GetACCDiv(5), EK
      , (*BS_out_vec)[2]
      , (*BS_out_vec)[3]
      , 0, start_bits, prefix_add_vals, carry_bit, out_nums))[0]);
	

		// Blind Rotate 1,
		for (uint32_t i = 0; i < K_FP+1; i++) {
				BS_out_tmp[i] = In_Poly_Sign[i];
				BS_out_tmp[i] *= params->GetRotMonomial_RLWE(1);
				BS_out_tmp[i] -= In_Poly_Sign[i];
		} 
	

		// Get our Sign
	  start_bits  = Q1_WO_bits  + 2 + 1;
		prefix_add_vals = (Q1_Int >> 1) + Q1_Int;  
    carry_bit = 0;
    out_nums = 1;
	  // ================   EVAL MODE SAVING ==================== 
    BS_out = std::move(*this->ProductAndBKAndKSRange(params
			, LWEscheme
			, params->GetACCDiv(2), EK
			,	 BS_out_tmp
      , (*BS_out_vec)[0]
			, 0,1, start_bits, prefix_add_vals, carry_bit, out_nums));
		
		// Adding
		for (uint32_t i = 0; i < K_FP+1; i++) {
			BS_out[i] += In_Poly_Sign[i];
		}
		// BS_out get current vals


		// Blind Rotation 2,
		for (uint32_t i = 0; i < K_FP+1; i++) {
				BS_out_tmp[i] = BS_out[i];
				BS_out_tmp[i] *= params->GetRotMonomial_RLWE(2);
				BS_out_tmp[i] -= BS_out[i];
		} 

		// Second Construction
		start_bits  = Q1_WO_bits  + 2 + 1;
		prefix_add_vals = (Q1_Int >> 1) + Q1_Int;  
    carry_bit = 0;
    out_nums = 1;
		BS_out_tmp = std::move(*this->ProductAndBKAndKSRange(params
			, LWEscheme
			, params->GetACCDiv(2), EK
      , BS_out_tmp
      , (*BS_out_vec)[1]
      , 0,3, start_bits, prefix_add_vals, carry_bit, out_nums));
		
		// Adding
		for (uint32_t i = 0; i < K_FP+1; i++) {
			BS_out[i] += BS_out_tmp[i];
		}

		// Blind Rotation 3,
		for (uint32_t i = 0; i < K_FP+1; i++) {
				BS_out_tmp[i] = BS_out[i];
				BS_out_tmp[i] *= params->GetRotMonomial_RLWE(4);
				BS_out_tmp[i] -= BS_out[i];
		}


		// Third Construction
		start_bits  = Q1_WO_bits  + 2 + 1;
		prefix_add_vals = (Q1_Int >> 1) + Q1_Int;  
    carry_bit = 0;
    out_nums = 1;

		BS_out_tmp = std::move(*this->ProductAndBKAndKSRange(params
			, LWEscheme
			, params->GetACCDiv(2), EK
      , BS_out_tmp
      , (*BS_out_vec)[2]
      , 0,7, start_bits, prefix_add_vals, carry_bit, out_nums));
	
		// Adding
		for (uint32_t i = 0; i < K_FP+1; i++) {
			BS_out[i] += BS_out_tmp[i];
		}


		// Remove 3 bit
		start_bits  = scale_bit  + 3 + 1;
		prefix_add_vals = (prefix >> 1);  
		carry_bit = 1;
    out_nums = 2;
	
		BS_out_vec = this->Bootstrap_and_KS(params
			, LWEscheme
			, params->GetACCDiv(3)
			, EK
			, A
			, b, 0, start_bits, prefix_add_vals, carry_bit, out_nums);
		
		start_bits  = Q1_WO_bits  + 3 + 1;
		prefix_add_vals = (Q1_Int >> 1);  
 		carry_bit = 0;
    out_nums = 1;
	
		BS_out_tmp =  std::move((*this->ProductAndBKAndKS(params
			, LWEscheme
			, params->GetACCDiv(4), EK
      , (*BS_out_vec)[0]
      , (*BS_out_vec)[1]
      , 0, start_bits, prefix_add_vals, carry_bit, out_nums))[0]);

		// Removes
		
		for (uint32_t i = 0; i < K_FP+1; i++) {
			In_Poly_Exp[i] -=  BS_out_tmp[i];
		}
	
		// Seconds - one
		for (uint32_t j = 0; j < K_FP; j++) {
			DCRT_tmp = In_Poly_Exp[j];
			DCRT_tmp.SetFormat(Format::COEFFICIENT);
			DivideRoundingSelf(params, &DCRT_tmp, Format::COEFFICIENT);
			A[j] = DCRT_tmp.GetElementAtIndex(0).GetValues();
		}
		// Copy B
		DCRT_tmp = In_Poly_Exp[K_FP];
		DCRT_tmp.SetFormat(Format::COEFFICIENT);
		DivideRoundingSelf(params, &DCRT_tmp, Format::COEFFICIENT);
		b = DCRT_tmp.GetElementAtIndex(0)[0];
	

		// First Carry, 4+1 (1 bit no information)bit bootstrapping Bootstrapping
  	start_bits  = (scale_bit + 3) + 2 + 1;
		prefix_add_vals = (prefix << 2);  // 3bit -> 2bit 
    carry_bit = 2;
    out_nums = 4;
		
		BS_out_vec = this->Bootstrap_and_KS(params
			, LWEscheme
			, params->GetACCDiv(1)
			, EK
			, A
			, b, 0, start_bits, prefix_add_vals, carry_bit, out_nums);

		
		//BS rot 3
		for (uint32_t i = 0; i < K_FP+1; i++) {
			if (max_msg >= 8) {
				(*BS_out_vec)[1][i] *= params->GetRotMonomial_RLWE(8); 
				(*BS_out_vec)[0][i] += (*BS_out_vec)[1][i];
			}
			if (max_msg >= 16) {
				(*BS_out_vec)[2][i] *= params->GetRotMonomial_RLWE(16); 
				(*BS_out_vec)[0][i] += (*BS_out_vec)[2][i];
			}
			if (max_msg >= 24) {
				(*BS_out_vec)[3][i] *= params->GetRotMonomial_RLWE(24); 	
			(*BS_out_vec)[0][i] += (*BS_out_vec)[3][i];
			} 
		}
		start_bits  = Q1_WO_bits  + 2 + 1;
		prefix_add_vals = (Q1_Int >> 1) + Q1_Int;  
    carry_bit = 0;
    out_nums = 1;
		
		BS_out = std::move(*this->ProductAndBKAndKSRange(params
			, LWEscheme
			, params->GetACCDiv(2), EK
      , BS_out
      , (*BS_out_vec)[0]
      , 0,max_msg, start_bits, prefix_add_vals, carry_bit, out_nums, mul));
		

		return std::make_shared<vector<DCRTPoly>>(BS_out);

	}




  std::shared_ptr<vector<vector<DCRTPoly>>> FPRingGSWAccumulatorScheme::BootstrapReluRevF(
    const std::shared_ptr<FPRingGSWCryptoParams> params, 
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<DCRTPoly> & ct_A,
		const vector<DCRTPoly> & ct_B,
		const std::shared_ptr<FPRingGSWEvalKey> EK) const {
	
		auto MulPrepare_start = std::chrono::system_clock::now();
    uint32_t K_FP = params->GetRLWEParams()->GetK();
		NativeInteger Q0 = params->GetModuliQ()[0];
    NativeInteger Q1 = params->GetModuliQ()[1];
		uint64_t Q1_Int = Q1.ConvertToInt();
		uint32_t Q1_WO_bits = params->GetRLWEParams()->GetBits()[1];  // <32 bit
		uint32_t start_bits;
		uint64_t prefix_add_vals;
		uint32_t carry_bit;
		uint32_t out_nums;
	
		std::shared_ptr<vector<vector<DCRTPoly>>> BS_out_vec;
		vector<DCRTPoly> BS_out;
		vector<NativeVector> A(K_FP);
		NativeInteger b;
		DCRTPoly DCRT_tmp;
		vector<DCRTPoly> In_Poly(K_FP+1);
	
		vector<vector<DCRTPoly>> Res_Poly(2);
		vector<DCRTPoly> Outs(K_FP+1);
		vector<vector<DCRTPoly>> Out_Poly(2);


		// ct1 - ct 2
		for (uint32_t i = 0; i < K_FP+1; i++) {	
			In_Poly[i] =  ct_A[i] - ct_B[i];
		}

		// Substract gives range - 11bit~ 12 bit
		// add 13 bit, bootstrapping, and check whether 13bit is 0 or not
		uint32_t scale_bit = params->GetExpBit();	
		uint64_t prefix = params->GetRLWEParams()->GetScaling().ConvertToInt();
			prefix <<= scale_bit;
		uint64_t bit_check = prefix;	
		bit_check <<= 7;  // 2^10 ? 2^ 11?	
		
		// Make

		In_Poly[K_FP].SetFormat(Format::COEFFICIENT);
		In_Poly[K_FP].GetElementW(0)[0].ModAddEq(
				NativeInteger(bit_check).ModMul(Q1, Q0),Q0);
		In_Poly[K_FP].GetElementW(1)[0].ModAddEq(
				NativeInteger(bit_check).ModMul(Q1, Q1),Q1);
		// Sub Values
		In_Poly[K_FP].SetFormat(Format::EVALUATION);
		

		////////////////////////////////   STARTS !!!! /////////////////////
		std::chrono::duration<double> MulPrepare_end = std::chrono::system_clock::now() - MulPrepare_start;
		if (DEBUG) {
			std::cout << "bootstrapping Relu prepare time is " << MulPrepare_end.count() << std::endl ;
		}


		// Bootstrap 1
		for (uint32_t i = 0; i < 2; i ++) { 
			
			//std::cout << "i is " << i << std::endl;	
			// Make Polys
			if (i != 1) {
				// Copy A
				for (uint32_t j = 0; j < K_FP; j++) {
					DCRT_tmp = In_Poly[j];
					DCRT_tmp.SetFormat(Format::COEFFICIENT);
					DivideRoundingSelf(params, &DCRT_tmp, Format::COEFFICIENT);
					A[j] = DCRT_tmp.GetElementAtIndex(0).GetValues();
				}
				// Copy B
				DCRT_tmp = In_Poly[K_FP];
				DCRT_tmp.SetFormat(Format::COEFFICIENT);
				DivideRoundingSelf(params, &DCRT_tmp, Format::COEFFICIENT);
				b = DCRT_tmp.GetElementAtIndex(0)[0];
			} else { // i == 1
				for (uint32_t j = 0; j < K_FP; j++) {
					In_Poly[j].SetFormat(Format::COEFFICIENT);
					DivideRoundingSelf(params, &In_Poly[j], Format::COEFFICIENT);
					A[j] = In_Poly[j].GetElementAtIndex(0).GetValues();
				}
				In_Poly[K_FP].SetFormat(Format::COEFFICIENT);
				DivideRoundingSelf(params, &In_Poly[K_FP], Format::COEFFICIENT);
				b = In_Poly[K_FP].GetElementAtIndex(0)[0];
			}
				

			// First Carry, 4+1 (1 bit no information)bit bootstrapping Bootstrapping
      if (i != 1) {

				if (i == 0) {
					start_bits  = scale_bit  + 4 + 1;
					prefix_add_vals = (prefix >> 1);
				} 
				//else if(i == 1) {
				//	start_bits  = scale_bit  + 8 + 1;
				//	prefix_add_vals = (prefix << 3);  
				//}
				carry_bit = 1;
				out_nums = 2;
			
				//First Bootstrapping,3 bit + random 1 bit is bootstrapping
				BS_out_vec = this->Bootstrap_and_KS(params
					, LWEscheme
					, params->GetACCRelu(0)
					, EK
					, A
					, b, 0, start_bits, prefix_add_vals, carry_bit, out_nums);
				
				DCRTPoly ACCs = params -> GetACCRelu(1+i);
				
				carry_bit = 1;
				out_nums = 2;
	
				// Orign, 4+1(padding but no infomation) bit bootstrapping msg <= 3bit 
				start_bits  = Q1_WO_bits + 4 + 1;
				prefix_add_vals = ( Q1_Int >> 1);  
				// ================   EVAL MODE SAVING ==================== 
				BS_out_vec = this->ProductAndBKAndKS(params
					, LWEscheme
					, ACCs, EK
					, (*BS_out_vec)[0]
					, (*BS_out_vec)[1]
					, 0, start_bits, prefix_add_vals, carry_bit, out_nums);
			
			} else if (i == 1) {
				start_bits  = scale_bit  + 8 + 1;
				prefix_add_vals = (prefix << 3);  
				carry_bit = 1;
				out_nums = 2;
				//First Bootstrapping,3 bit + random 1 bit is bootstrapping
				BS_out_vec = this->Bootstrap_and_KS(params
					, LWEscheme
					, params->GetACCRelu(8)
					, EK
					, A
					, b, 0, start_bits, prefix_add_vals, carry_bit, out_nums);
			}

			// After Handle
			//std::cout << "handle before " << std::endl;
			if (i != 1) {
				// remove 4 or 8 bit info
				for (uint32_t j = 0; j < K_FP+1; j++) {
					In_Poly[j]-=(*BS_out_vec)[1][j]; 
				}

				//std::cout << "Minus is end" << std::endl;
				// Saves results 
				Res_Poly[i] = std::move((*BS_out_vec)[0]);
			
			} else {// i == 1
						
				// first 4 bit
				start_bits  = Q1_WO_bits + 4 + 1;
				carry_bit = 0;
				out_nums = 1;
				prefix_add_vals = ( Q1_Int >> 1);  
				
				Outs =  std::move((*this->ProductAndBKAndKS(params
        , LWEscheme
				, params->GetACCRelu(4), EK
        , Res_Poly[0]
        , (*BS_out_vec)[0]
        , 0, start_bits, prefix_add_vals, carry_bit, out_nums))[0]);
       	
				/*
				// Second 4 bit
				start_bits  = Q1_WO_bits + 4 + 1;
				carry_bit = 0;
				out_nums = 1;
				prefix_add_vals = ( Q1_Int >> 1);  
				BS_out =  std::move((*this->ProductAndBKAndKS(params
        , LWEscheme
				, params->GetACCRelu(5), EK
        , Res_Poly[1]
        , (*BS_out_vec)[0]
        , 0, start_bits, prefix_add_vals, carry_bit, out_nums))[0]);
				*/

				// Some Seconds & Third 
				for (uint32_t j = 0; j < K_FP+1; j++) {
					//Outs[j] += BS_out[j];
					Outs[j] += (*BS_out_vec)[1][j];
			
				}
				
			}
		}

		// Sums mins
		
		for (uint32_t j = 0; j < K_FP+1; j++) {
			Outs[j] +=  ct_B[j];
		}
		
		Out_Poly[0] = std::move(Outs);
		return std::make_shared<vector<vector<DCRTPoly>>>(Out_Poly);
	}






  std::shared_ptr<vector<vector<DCRTPoly>>> FPRingGSWAccumulatorScheme::CheckOFUFMul(
    const std::shared_ptr<FPRingGSWCryptoParams> params, 
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<DCRTPoly> & ct_UF,
		const vector<DCRTPoly> & ct_EXP,
		const CTType types,
		const std::shared_ptr<FPRingGSWEvalKey> EK) const {
	
		
		auto MulPrepare_start = std::chrono::system_clock::now();
    uint32_t K_FP = params->GetRLWEParams()->GetK();
		NativeInteger Q0 = params->GetModuliQ()[0];
    NativeInteger Q1 = params->GetModuliQ()[1];
		uint64_t Q1_Int = Q1.ConvertToInt();
		uint32_t Q1_WO_bits = params->GetRLWEParams()->GetBits()[1];  // <32 bit
		uint32_t start_bits;
		uint64_t prefix_add_vals;
		uint32_t carry_bit;
		uint32_t out_nums;
	
		std::shared_ptr<vector<vector<DCRTPoly>>> BS_out_vec;
		vector<DCRTPoly> BS_out;
		vector<NativeVector> A(K_FP);
		NativeInteger b;
		DCRTPoly DCRT_tmp;
		vector<DCRTPoly> In_Poly(K_FP+1);
	
		vector<vector<DCRTPoly>> Res_Poly(3);
		vector<DCRTPoly> Outs(K_FP+1);
		vector<DCRTPoly> OF(K_FP+1);
		vector<DCRTPoly> UF(K_FP+1);
		vector<DCRTPoly> OF_sus(K_FP+1);
		vector<DCRTPoly> UF_sus(K_FP+1);
		
		vector<vector<DCRTPoly>> Out_Poly(2);


		// ct1 - ct 2
		for (uint32_t i = 0; i < K_FP+1; i++) {	
			In_Poly[i] =  ct_EXP[i];
		}

		// Substract gives range - 11bit~ 12 bit
		// add 13 bit, bootstrapping, and check whether 13bit is 0 or not
		uint32_t scale_bit = params->GetExpBit();	
		uint64_t prefix = params->GetRLWEParams()->GetScaling().ConvertToInt();
			prefix <<= scale_bit;
	

		////////////////////////////////   STARTS !!!! /////////////////////
		std::chrono::duration<double> MulPrepare_end = std::chrono::system_clock::now() - MulPrepare_start;
		if (DEBUG) {
			std::cout << "OFUF prepare time is " << MulPrepare_end.count() << std::endl ;
		}
	

		// Lowerbound test
		//32 bit -> whether it is lower then 63
		//64 bit -> whether it is lower than 511
		//std::shared_ptr<vector<DCRTPoly>> UTMust;
		if (types == CTType::DOUBLE) {
			UF_sus = std::move(*this->BootstrapIsLower(params,  
					LWEscheme,
					ct_EXP,
					512,
					EK));
		} else if (types == CTType::FLOAT) {
			UF_sus = std::move(*this->BootstrapIsLowerF(params,  
					LWEscheme,
					ct_EXP,
					64,
					EK));
		}

		//std::cout << "OKs " << std::endl;
		for (uint32_t j = 0; j < K_FP; j++) {
			UF_sus[j].SetFormat(Format::COEFFICIENT);
			A[j] = UF_sus[j].GetElementAtIndex(0).GetValues();
		}
		UF_sus[K_FP].SetFormat(Format::COEFFICIENT);
		b = UF_sus[K_FP].GetElementAtIndex(0)[0];
					
		start_bits  = Q1_WO_bits + 4 + 1;
		carry_bit = 1;
		out_nums = 2;
		prefix_add_vals = ( Q1_Int >> 1); 
		
		//First Bootstrapping,3 bit + random 1 bit is bootstrapping
		BS_out_vec = this->Bootstrap_and_KS(params
			, LWEscheme
			, params->GetACCRelu(16)
			, EK
			, A
			, b, 0, start_bits, prefix_add_vals, carry_bit, out_nums);
		UF_sus = std::move((*BS_out_vec)[1]);
		UF = std::move((*BS_out_vec)[0]);

		// Upper bound test
	 //Minus Bias 
		uint64_t prefix2 = params->GetRLWEParams()->GetScaling().ConvertToInt();
		prefix2 <<= scale_bit;
		if (types == CTType::DOUBLE) {
			prefix2 *= 511; 
		} else if (types == CTType::FLOAT) {
			prefix2 *= 63;
		}

		In_Poly[K_FP].SetFormat(Format::COEFFICIENT);
		In_Poly[K_FP].GetElementW(0)[0].ModSubEq(Q1.ModMul(prefix2, Q0), Q0); 
		In_Poly[K_FP].SetFormat(Format::EVALUATION);

		uint32_t end_idx=0;
		if (types == CTType::DOUBLE) {
			end_idx = 2;
		} else if(types == CTType::FLOAT) {
			end_idx = 1;
		}

		// Bootstrap 1
		for (uint32_t i = 0; i <= end_idx; i ++) { 
			// Make Polys
			if (i != end_idx) {
				// Copy A
				for (uint32_t j = 0; j < K_FP; j++) {
					DCRT_tmp = In_Poly[j];
					DCRT_tmp.SetFormat(Format::COEFFICIENT);
					DivideRoundingSelf(params, &DCRT_tmp, Format::COEFFICIENT);
					A[j] = DCRT_tmp.GetElementAtIndex(0).GetValues();
				}
				// Copy B
				DCRT_tmp = In_Poly[K_FP];
				DCRT_tmp.SetFormat(Format::COEFFICIENT);
				DivideRoundingSelf(params, &DCRT_tmp, Format::COEFFICIENT);
				b = DCRT_tmp.GetElementAtIndex(0)[0];
			} else { // i == 1 or 2
				for (uint32_t j = 0; j < K_FP; j++) {
					In_Poly[j].SetFormat(Format::COEFFICIENT);
					DivideRoundingSelf(params, &In_Poly[j], Format::COEFFICIENT);
					A[j] = In_Poly[j].GetElementAtIndex(0).GetValues();
				}
				In_Poly[K_FP].SetFormat(Format::COEFFICIENT);
				DivideRoundingSelf(params, &In_Poly[K_FP], Format::COEFFICIENT);
				b = In_Poly[K_FP].GetElementAtIndex(0)[0];
			}
			
			// First Carry, 4+1 (1 bit no information)bit bootstrapping Bootstrapping
      if (i != end_idx) {
				// if end_idx == 1, then i==1 can't come in.
				if (i == 0) {
					start_bits  = scale_bit  + 4 + 1;
					prefix_add_vals = (prefix >> 1);  
				} else if (i == 1) {
					start_bits  = scale_bit  + 8 + 1;
					prefix_add_vals = (prefix << 3);  
				} 
				carry_bit = 1;
				out_nums = 2;
			
				//First Bootstrapping,3 bit + random 1 bit is bootstrapping
				BS_out_vec = this->Bootstrap_and_KS(params
					, LWEscheme
					, params->GetACCRelu(0)
					, EK
					, A
					, b, 0, start_bits, prefix_add_vals, carry_bit, out_nums);
		
				//Second Bootsrapp
				// i = 1, 2, 3 ...
				DCRTPoly ACCs = params -> GetACCRelu(i+11);

				// Orign, 4+1(padding but no infomation) bit bootstrapping msg <= 3bit 
				start_bits  = Q1_WO_bits + 4 + 1;
				carry_bit = 1;
				out_nums = 2;
				prefix_add_vals = ( Q1_Int >> 1);  
				// ================   EVAL MODE SAVING ==================== 
				BS_out_vec = this->ProductAndBKAndKS(params
					, LWEscheme
					, ACCs, EK
					, (*BS_out_vec)[0]
					, (*BS_out_vec)[1]
					, 0, start_bits, prefix_add_vals, carry_bit, out_nums);
			} else if (i == end_idx) {
				DCRTPoly ACCs;

				if (types == CTType::FLOAT) {	
					start_bits  = scale_bit  + 8 + 1; // OF Check
					prefix_add_vals = (prefix << 3);
					ACCs = params->GetACCRelu(13);
				} else if(types == CTType::DOUBLE) {
					start_bits  = scale_bit  + 12 + 1; // OF Check
					prefix_add_vals = (prefix << 7);
					ACCs = params->GetACCRelu(14);
				}

				carry_bit = 1;
				out_nums = 2;
				BS_out_vec = this->Bootstrap_and_KS(params
				, LWEscheme
				, ACCs
				, EK
				, A
				, b, 0, start_bits, prefix_add_vals, carry_bit, out_nums);
				//OF_sus = std::move((*BS_out_vec)[0]);
				OF_sus = std::move((*BS_out_vec)[0]);
				//Res_Poly[end_idx] = std::move((*BS_out_vec)[1]);
			
				start_bits  = Q1_WO_bits + 4 + 1;
				carry_bit = 0;
				out_nums = 1;
				prefix_add_vals = ( Q1_Int >> 1);  
				// ================   EVAL MODE SAVING ==================== 
				
				OF = std::move((*this->ProductAndBKAndKS(params
					, LWEscheme
					, params->GetACCRelu(17), EK
					, OF_sus
					, UF_sus
					, 0, start_bits, prefix_add_vals, carry_bit, out_nums))[0]);
			}	
				
			// After Handle
			//std::cout << "handle before " << std::endl;
			if (i != end_idx) {
				
				// remove 4 or 8 bit info
				for (uint32_t j = 0; j < K_FP+1; j++) {
					In_Poly[j]-=(*BS_out_vec)[1][j]; 
				}

				//std::cout << "Minus is end" << std::endl;
				// Saves results 
				Res_Poly[i] = std::move((*BS_out_vec)[0]);
				//std::cout << "Res is moved" << std::endl;
				
			} else {// i == 2
				// Orign, 4+1(padding but no infomation) bit bootstrapping msg <= 3bit 
				/*
				for (uint32_t j = 0; j < K_FP+1; j++) {
						Res_Poly[0][j] += Res_Poly[1][j];
				}	
				if (types == CTType::DOUBLE) {
					for (uint32_t j = 0; j < K_FP+1; j++) {
						Res_Poly[0][j] += Res_Poly[2][j];
					}	
				}
				// first 4 bit
				start_bits  = Q1_WO_bits + 4 + 1;
				carry_bit = 0;
				out_nums = 1;
				prefix_add_vals = ( Q1_Int >> 1);  
				
				UF =  std::move((*this->ProductAndBKAndKS(params
        , LWEscheme
				, params->GetACCRelu(15), EK
        , Res_Poly[0]
        , ct_UF
        , 0, start_bits, prefix_add_vals, carry_bit, out_nums))[0]);
				//  UF is over
				*/
			}
			
		}
		// Sums mins
		Out_Poly[0] = std::move(OF);
		Out_Poly[1] = std::move(UF);
		return std::make_shared<vector<vector<DCRTPoly>>>(Out_Poly);
	}








  std::shared_ptr<vector<vector<DCRTPoly>>> FPRingGSWAccumulatorScheme::CheckOFUF(
    const std::shared_ptr<FPRingGSWCryptoParams> params, 
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<DCRTPoly> & ct_UF,
		const vector<DCRTPoly> & ct_EXP,
		const CTType types,
		const std::shared_ptr<FPRingGSWEvalKey> EK) const {
	
		
		auto MulPrepare_start = std::chrono::system_clock::now();
    uint32_t K_FP = params->GetRLWEParams()->GetK();
		//NativeInteger Q0 = params->GetModuliQ()[0];
    NativeInteger Q1 = params->GetModuliQ()[1];
		uint64_t Q1_Int = Q1.ConvertToInt();
		uint32_t Q1_WO_bits = params->GetRLWEParams()->GetBits()[1];  // <32 bit
		uint32_t start_bits;
		uint64_t prefix_add_vals;
		uint32_t carry_bit;
		uint32_t out_nums;
	
		std::shared_ptr<vector<vector<DCRTPoly>>> BS_out_vec;
		vector<DCRTPoly> BS_out;
		vector<NativeVector> A(K_FP);
		NativeInteger b;
		DCRTPoly DCRT_tmp;
		vector<DCRTPoly> In_Poly(K_FP+1);
	
		vector<vector<DCRTPoly>> Res_Poly(3);
		vector<DCRTPoly> Outs(K_FP+1);
		vector<DCRTPoly> OF(K_FP+1);
		vector<DCRTPoly> UF(K_FP+1);
		
		vector<vector<DCRTPoly>> Out_Poly(2);


		// ct1 - ct 2
		for (uint32_t i = 0; i < K_FP+1; i++) {	
			In_Poly[i] =  ct_EXP[i];
		}

		// Substract gives range - 11bit~ 12 bit
		// add 13 bit, bootstrapping, and check whether 13bit is 0 or not
		uint32_t scale_bit = params->GetExpBit();	
		uint64_t prefix = params->GetRLWEParams()->GetScaling().ConvertToInt();
			prefix <<= scale_bit;
	

		////////////////////////////////   STARTS !!!! /////////////////////
		std::chrono::duration<double> MulPrepare_end = std::chrono::system_clock::now() - MulPrepare_start;
		if (DEBUG) {
			std::cout << "OFUF prepare time is " << MulPrepare_end.count() << std::endl ;
		}
		
		uint32_t end_idx=0;
		if (types == CTType::DOUBLE) {
			end_idx = 2;
		} else if(types == CTType::FLOAT) {
			end_idx = 1;
		}

		// Bootstrap 1
		for (uint32_t i = 0; i <= end_idx; i ++) { 
			// Make Polys
			if (i != end_idx) {
				// Copy A
				for (uint32_t j = 0; j < K_FP; j++) {
					DCRT_tmp = In_Poly[j];
					DCRT_tmp.SetFormat(Format::COEFFICIENT);
					DivideRoundingSelf(params, &DCRT_tmp, Format::COEFFICIENT);
					A[j] = DCRT_tmp.GetElementAtIndex(0).GetValues();
				}
				// Copy B
				DCRT_tmp = In_Poly[K_FP];
				DCRT_tmp.SetFormat(Format::COEFFICIENT);
				DivideRoundingSelf(params, &DCRT_tmp, Format::COEFFICIENT);
				b = DCRT_tmp.GetElementAtIndex(0)[0];
			} else { // i == 1 or 2
				for (uint32_t j = 0; j < K_FP; j++) {
					In_Poly[j].SetFormat(Format::COEFFICIENT);
					DivideRoundingSelf(params, &In_Poly[j], Format::COEFFICIENT);
					A[j] = In_Poly[j].GetElementAtIndex(0).GetValues();
				}
				In_Poly[K_FP].SetFormat(Format::COEFFICIENT);
				DivideRoundingSelf(params, &In_Poly[K_FP], Format::COEFFICIENT);
				b = In_Poly[K_FP].GetElementAtIndex(0)[0];
			}
			
			// First Carry, 4+1 (1 bit no information)bit bootstrapping Bootstrapping
      if (i != end_idx) {
				// if end_idx == 1, then i==1 can't come in.
				prefix_add_vals =0;
                start_bits  = 0;
                if (i == 0) {
					start_bits  = scale_bit  + 4 + 1;
					prefix_add_vals = (prefix >> 1);  
				} else if (i == 1) {
					start_bits  = scale_bit  + 8 + 1;
					prefix_add_vals = (prefix << 3);  
				} 
				carry_bit = 1;
				out_nums = 2;
			
				//First Bootstrapping,3 bit + random 1 bit is bootstrapping
				BS_out_vec = this->Bootstrap_and_KS(params
					, LWEscheme
					, params->GetACCRelu(0)
					, EK
					, A
					, b, 0, start_bits, prefix_add_vals, carry_bit, out_nums);
		
				//Second Bootsrapp
				// i = 1, 2, 3 ...
				DCRTPoly ACCs = params -> GetACCRelu(i+11);
				//	ACCs = params -> GetACCRelu(6);
				
				// Orign, 4+1(padding but no infomation) bit bootstrapping msg <= 3bit 
				start_bits  = Q1_WO_bits + 4 + 1;
				carry_bit = 1;
				out_nums = 2;
				prefix_add_vals = ( Q1_Int >> 1);  
				// ================   EVAL MODE SAVING ==================== 
				BS_out_vec = this->ProductAndBKAndKS(params
					, LWEscheme
					, ACCs, EK
					, (*BS_out_vec)[0]
					, (*BS_out_vec)[1]
					, 0, start_bits, prefix_add_vals, carry_bit, out_nums);
			} else if (i == end_idx) {
				DCRTPoly ACCs;
                start_bits = 0;
                prefix_add_vals = 0;
				if (types == CTType::FLOAT) {
		
					start_bits  = scale_bit  + 8 + 1; // OF Check
					prefix_add_vals = (prefix << 3);
					ACCs = params->GetACCRelu(13);
				} else if(types == CTType::DOUBLE) {
					start_bits  = scale_bit  + 12 + 1; // OF Check
					prefix_add_vals = (prefix << 7);
					ACCs = params->GetACCRelu(14);
				}
			    out_nums = 2;
                carry_bit = 1;
                
				BS_out_vec = this->Bootstrap_and_KS(params
				, LWEscheme
				, ACCs
				, EK
				, A
				, b, 0, start_bits, prefix_add_vals, carry_bit, out_nums);
				OF = std::move((*BS_out_vec)[0]);
				Res_Poly[end_idx] = std::move((*BS_out_vec)[1]);
			
			}	
		
			// After Handle
			//std::cout << "handle before " << std::endl;
			if (i != end_idx) {
				
				// remove 4 or 8 bit info
				for (uint32_t j = 0; j < K_FP+1; j++) {
					In_Poly[j]-=(*BS_out_vec)[1][j]; 
				}

				//std::cout << "Minus is end" << std::endl;
				// Saves results 
				Res_Poly[i] = std::move((*BS_out_vec)[0]);
				//std::cout << "Res is moved" << std::endl;
				
			} else {// i == 2
				for (uint32_t j = 0; j < K_FP+1; j++) {
						Res_Poly[0][j] += Res_Poly[1][j];
				}	
				if (types == CTType::DOUBLE) {
					for (uint32_t j = 0; j < K_FP+1; j++) {
						Res_Poly[0][j] += Res_Poly[2][j];
					}	
				}
				// first 4 bit
				start_bits  = Q1_WO_bits + 4 + 1;
				carry_bit = 0;
				out_nums = 1;
				prefix_add_vals = ( Q1_Int >> 1);  
				
				UF =  std::move((*this->ProductAndBKAndKS(params
        , LWEscheme
				, params->GetACCRelu(15), EK
        , Res_Poly[0]
        , ct_UF
        , 0, start_bits, prefix_add_vals, carry_bit, out_nums))[0]);
				//  UF is over
			}
		}
		// Sums mins
		Out_Poly[0] = std::move(OF);
		Out_Poly[1] = std::move(UF);
		return std::make_shared<vector<vector<DCRTPoly>>>(Out_Poly);
	}

  std::shared_ptr<vector<vector<DCRTPoly>>> FPRingGSWAccumulatorScheme::BootstrapReluRev(
    const std::shared_ptr<FPRingGSWCryptoParams> params, 
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<DCRTPoly> & ct_A,
		const vector<DCRTPoly> & ct_B,
		const std::shared_ptr<FPRingGSWEvalKey> EK) const {
	
		auto MulPrepare_start = std::chrono::system_clock::now();
    uint32_t K_FP = params->GetRLWEParams()->GetK();
		NativeInteger Q0 = params->GetModuliQ()[0];
    NativeInteger Q1 = params->GetModuliQ()[1];
		uint64_t Q1_Int = Q1.ConvertToInt();
		uint32_t Q1_WO_bits = params->GetRLWEParams()->GetBits()[1];  // <32 bit
		uint32_t start_bits;
		uint64_t prefix_add_vals;
		uint32_t carry_bit;
		uint32_t out_nums;
	
		std::shared_ptr<vector<vector<DCRTPoly>>> BS_out_vec;
		vector<DCRTPoly> BS_out;
		vector<NativeVector> A(K_FP);
		NativeInteger b;
		DCRTPoly DCRT_tmp;
		vector<DCRTPoly> In_Poly(K_FP+1);
	
		vector<vector<DCRTPoly>> Res_Poly(2);
		vector<DCRTPoly> Outs(K_FP+1);
		vector<vector<DCRTPoly>> Out_Poly(2);


		// ct1 - ct 2
		for (uint32_t i = 0; i < K_FP+1; i++) {	
			In_Poly[i] =  ct_A[i] - ct_B[i];
		}

		// Substract gives range - 11bit~ 12 bit
		// add 13 bit, bootstrapping, and check whether 13bit is 0 or not
		uint32_t scale_bit = params->GetExpBit();	
		uint64_t prefix = params->GetRLWEParams()->GetScaling().ConvertToInt();
			prefix <<= scale_bit;
		uint64_t bit_check = prefix;	
		bit_check <<= 10;  // 2^10 ? 2^ 11?	
		
		
		// Make

		In_Poly[K_FP].SetFormat(Format::COEFFICIENT);
		In_Poly[K_FP].GetElementW(0)[0].ModAddEq(
				NativeInteger(bit_check).ModMul(Q1, Q0),Q0);
		In_Poly[K_FP].GetElementW(1)[0].ModAddEq(
				NativeInteger(bit_check).ModMul(Q1, Q1),Q1);
		// Sub Values
		In_Poly[K_FP].SetFormat(Format::EVALUATION);
		

		////////////////////////////////   STARTS !!!! /////////////////////
		std::chrono::duration<double> MulPrepare_end = std::chrono::system_clock::now() - MulPrepare_start;
		if (DEBUG) {
			std::cout << "bootstrapping Relu prepare time is " << MulPrepare_end.count() << std::endl ;
		}


		// Bootstrap 1
		for (uint32_t i = 0; i < 3; i ++) { 
			
			//std::cout << "i is " << i << std::endl;	
			// Make Polys
			if (i != 2) {
				// Copy A
				for (uint32_t j = 0; j < K_FP; j++) {
					DCRT_tmp = In_Poly[j];
					DCRT_tmp.SetFormat(Format::COEFFICIENT);
					DivideRoundingSelf(params, &DCRT_tmp, Format::COEFFICIENT);
					A[j] = DCRT_tmp.GetElementAtIndex(0).GetValues();
				}
				// Copy B
				DCRT_tmp = In_Poly[K_FP];
				DCRT_tmp.SetFormat(Format::COEFFICIENT);
				DivideRoundingSelf(params, &DCRT_tmp, Format::COEFFICIENT);
				b = DCRT_tmp.GetElementAtIndex(0)[0];
			} else { // i == 2
				for (uint32_t j = 0; j < K_FP; j++) {
					In_Poly[j].SetFormat(Format::COEFFICIENT);
					DivideRoundingSelf(params, &In_Poly[j], Format::COEFFICIENT);
					A[j] = In_Poly[j].GetElementAtIndex(0).GetValues();
				}
				In_Poly[K_FP].SetFormat(Format::COEFFICIENT);
				DivideRoundingSelf(params, &In_Poly[K_FP], Format::COEFFICIENT);
				b = In_Poly[K_FP].GetElementAtIndex(0)[0];
			}
			
			// First Carry, 4+1 (1 bit no information)bit bootstrapping Bootstrapping
      if (i != 2) {
				if (i == 0) {
					start_bits  = scale_bit  + 4 + 1;
					prefix_add_vals = (prefix >> 1);  
				} else if (i == 1) {
					start_bits  = scale_bit  + 8 + 1;
					prefix_add_vals = (prefix << 3);  
				} 
				carry_bit = 1;
				out_nums = 2;
			
				//First Bootstrapping,3 bit + random 1 bit is bootstrapping
				BS_out_vec = this->Bootstrap_and_KS(params
					, LWEscheme
					, params->GetACCRelu(0)
					, EK
					, A
					, b, 0, start_bits, prefix_add_vals, carry_bit, out_nums);
		
				//Second Bootsrapp
				// i = 1, 2, 3 ...
				DCRTPoly ACCs = params -> GetACCRelu(i+1);
				//	ACCs = params -> GetACCRelu(6);
				
				// Orign, 4+1(padding but no infomation) bit bootstrapping msg <= 3bit 
				start_bits  = Q1_WO_bits + 4 + 1;
				carry_bit = 1;
				out_nums = 2;
				prefix_add_vals = ( Q1_Int >> 1);  
				// ================   EVAL MODE SAVING ==================== 
				BS_out_vec = this->ProductAndBKAndKS(params
					, LWEscheme
					, ACCs, EK
					, (*BS_out_vec)[0]
					, (*BS_out_vec)[1]
					, 0, start_bits, prefix_add_vals, carry_bit, out_nums);
			} else if (i == 2) {
				start_bits  = scale_bit  + 12 + 1;
				prefix_add_vals = (prefix << 7);  
				carry_bit = 1;
				out_nums = 2;
			
				//First Bootstrapping,3 bit + random 1 bit is bootstrapping
				BS_out_vec = this->Bootstrap_and_KS(params
					, LWEscheme
					, params->GetACCRelu(6)
					, EK
					, A
					, b, 0, start_bits, prefix_add_vals, carry_bit, out_nums);
			}
		
			// After Handle
			//std::cout << "handle before " << std::endl;
			if (i != 2) {
				
				// remove 4 or 8 bit info
				for (uint32_t j = 0; j < K_FP+1; j++) {
					In_Poly[j]-=(*BS_out_vec)[1][j]; 
				}

				//std::cout << "Minus is end" << std::endl;
				// Saves results 
				Res_Poly[i] = std::move((*BS_out_vec)[0]);
				//std::cout << "Res is moved" << std::endl;
				
			} else {// i == 2
						
				// first 4 bit
				start_bits  = Q1_WO_bits + 4 + 1;
				carry_bit = 0;
				out_nums = 1;
				prefix_add_vals = ( Q1_Int >> 1);  
				
				Outs =  std::move((*this->ProductAndBKAndKS(params
        , LWEscheme
				, params->GetACCRelu(4), EK
        , Res_Poly[0]
        , (*BS_out_vec)[1]
        , 0, start_bits, prefix_add_vals, carry_bit, out_nums))[0]);
        
				// Second 4 bit
				start_bits  = Q1_WO_bits + 4 + 1;
				carry_bit = 0;
				out_nums = 1;
				prefix_add_vals = ( Q1_Int >> 1);  
				BS_out =  std::move((*this->ProductAndBKAndKS(params
        , LWEscheme
				, params->GetACCRelu(5), EK
        , Res_Poly[1]
        , (*BS_out_vec)[1]
        , 0, start_bits, prefix_add_vals, carry_bit, out_nums))[0]);
        
				// Some Seconds & Third 
				
				for (uint32_t j = 0; j < K_FP+1; j++) {
					Outs[j] += BS_out[j];
					Outs[j] += (*BS_out_vec)[0][j];
				}
				
			}
		}

		
		// Sums mins
	
		for (uint32_t j = 0; j < K_FP+1; j++) {
			Outs[j] +=  ct_B[j];
		}
		Out_Poly[0] = std::move(Outs);
		return std::make_shared<vector<vector<DCRTPoly>>>(Out_Poly);
	}




  std::shared_ptr<vector<DCRTPoly>> FPRingGSWAccumulatorScheme::BootstrapReluF(
    const std::shared_ptr<FPRingGSWCryptoParams> params, 
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<DCRTPoly> & ct1_A,
		const DCRTPoly &ct1_B,
		const vector<DCRTPoly> & ct2_A,
		const DCRTPoly &ct2_B,
		const std::shared_ptr<FPRingGSWEvalKey> EK) const {
	
		auto MulPrepare_start = std::chrono::system_clock::now();
    uint32_t K_FP = params->GetRLWEParams()->GetK();
		NativeInteger Q0 = params->GetModuliQ()[0];
    NativeInteger Q1 = params->GetModuliQ()[1];
		uint64_t Q1_Int = Q1.ConvertToInt();
		uint32_t Q1_WO_bits = params->GetRLWEParams()->GetBits()[1];  // <32 bit
		uint32_t start_bits;
		uint64_t prefix_add_vals;
		uint32_t carry_bit;
		uint32_t out_nums;
	
		std::shared_ptr<vector<vector<DCRTPoly>>> BS_out_vec;
		vector<DCRTPoly> BS_out;
		vector<NativeVector> A(K_FP);
		NativeInteger b;
		DCRTPoly DCRT_tmp;
		vector<DCRTPoly> In_Poly(K_FP+1);
	
		vector<vector<DCRTPoly>> Res_Poly(2);
		vector<DCRTPoly> Outs(K_FP+1);
		// ct1 - ct 2
		for (uint32_t i = 0; i < K_FP+1; i++) {
			if (i == K_FP) {
				In_Poly[i] =  ct1_B;
				In_Poly[i] -= ct2_B;
			} else {
				In_Poly[i] =  ct1_A[i];
				In_Poly[i] -= ct2_A[i];
			}
		}
		
		// Substract gives range - 11bit~ 12 bit
		// add 13 bit, bootstrapping, and check whether 13bit is 0 or not
		uint32_t scale_bit = params->GetExpBit();	
		uint64_t prefix = params->GetRLWEParams()->GetScaling().ConvertToInt();
			prefix <<= scale_bit;
		uint64_t bit_check = prefix;	
		bit_check <<= 7;  // 2^7 ? 2^ 8?	
		
		In_Poly[K_FP].SetFormat(Format::COEFFICIENT);
		In_Poly[K_FP].GetElementW(0)[0].ModAddEq(
				NativeInteger(bit_check).ModMul(Q1, Q0),Q0);
		In_Poly[K_FP].GetElementW(1)[0].ModAddEq(
				NativeInteger(bit_check).ModMul(Q1, Q1),Q1);
		In_Poly[K_FP].SetFormat(Format::EVALUATION);
		

		////////////////////////////////   STARTS !!!! /////////////////////
		std::chrono::duration<double> MulPrepare_end = std::chrono::system_clock::now() - MulPrepare_start;
		if (DEBUG) std::cout << "bootstrapping Relu prepare time is " << MulPrepare_end.count() << std::endl ;
		
		// Bootstrap 1
		for (uint32_t i = 0; i < 2; i ++) { 
			// Make Polys
			if (i != 1) {
				// Copy A
				for (uint32_t j = 0; j < K_FP; j++) {
					DCRT_tmp = In_Poly[j];
					DCRT_tmp.SetFormat(Format::COEFFICIENT);
					DivideRoundingSelf(params, &DCRT_tmp, Format::COEFFICIENT);
					A[j] = DCRT_tmp.GetElementAtIndex(0).GetValues();
				}
				// Copy B
				DCRT_tmp = In_Poly[K_FP];
				DCRT_tmp.SetFormat(Format::COEFFICIENT);
				DivideRoundingSelf(params, &DCRT_tmp, Format::COEFFICIENT);
				b = DCRT_tmp.GetElementAtIndex(0)[0];
			} else { // i == 2
				for (uint32_t j = 0; j < K_FP; j++) {
					In_Poly[j].SetFormat(Format::COEFFICIENT);
					DivideRoundingSelf(params, &In_Poly[j], Format::COEFFICIENT);
					A[j] = In_Poly[j].GetElementAtIndex(0).GetValues();
				}
				In_Poly[K_FP].SetFormat(Format::COEFFICIENT);
				DivideRoundingSelf(params, &In_Poly[K_FP], Format::COEFFICIENT);
				b = In_Poly[K_FP].GetElementAtIndex(0)[0];
			}
				
			// First Carry, 4+1 (1 bit no information)bit bootstrapping Bootstrapping
      if (i != 1) {

				if (i == 0) {
					start_bits  = scale_bit  + 4 + 1;
					prefix_add_vals = (prefix >> 1);
				} 
				//else if(i == 1) {
				//	start_bits  = scale_bit  + 8 + 1;
				//	prefix_add_vals = (prefix << 3);  
				//}
				carry_bit = 1;
				out_nums = 2;
			
				//First Bootstrapping,3 bit + random 1 bit is bootstrapping
				BS_out_vec = this->Bootstrap_and_KS(params
					, LWEscheme
					, params->GetACCRelu(0)
					, EK
					, A
					, b, 0, start_bits, prefix_add_vals, carry_bit, out_nums);
				DCRTPoly ACCs; 
				if (i == 0) {
					ACCs = params -> GetACCRelu(1);
				} 
				//else if (i == 1) {
				//	ACCs = params -> GetACCRelu(2);
				//}
				carry_bit = 1;
				out_nums = 2;
	
				// Orign, 4+1(padding but no infomation) bit bootstrapping msg <= 3bit 
				start_bits  = Q1_WO_bits + 4 + 1;
				prefix_add_vals = ( Q1_Int >> 1);  
				// ================   EVAL MODE SAVING ==================== 
				BS_out_vec = this->ProductAndBKAndKS(params
					, LWEscheme
					, ACCs, EK
					, (*BS_out_vec)[0]
					, (*BS_out_vec)[1]
					, 0, start_bits, prefix_add_vals, carry_bit, out_nums);
			
			} else if (i == 1) {
				start_bits  = scale_bit  + 8 + 1;
				prefix_add_vals = (prefix << 3);  
				carry_bit = 1;
				out_nums = 2;
				//First Bootstrapping,3 bit + random 1 bit is bootstrapping
				BS_out_vec = this->Bootstrap_and_KS(params
					, LWEscheme
					, params->GetACCRelu(7)
					, EK
					, A
					, b, 0, start_bits, prefix_add_vals, carry_bit, out_nums);
			}
 		
			// After Handle
			//std::cout << "handle before " << std::endl;
			if (i != 1) {
				// remove 4 or 8 bit info
				for (uint32_t j = 0; j < K_FP+1; j++) {
					In_Poly[j]-=(*BS_out_vec)[1][j]; 
				}

				//std::cout << "Minus is end" << std::endl;
				// Saves results 
				Res_Poly[i] = std::move((*BS_out_vec)[0]);
				//std::cout << "Res is moved" << std::endl;
				
			} else {// i == 2
				
				// first 4 bit
				start_bits  = Q1_WO_bits + 4 + 1;
				carry_bit = 0;
				out_nums = 1;
				prefix_add_vals = ( Q1_Int >> 1);  
				
				Outs =  std::move((*this->ProductAndBKAndKS(params
        , LWEscheme
				, params->GetACCRelu(4), EK
        , Res_Poly[0]
        , (*BS_out_vec)[0]
        , 0, start_bits, prefix_add_vals, carry_bit, out_nums))[0]);
				/*
				// Second 4 bit
				start_bits  = Q1_WO_bits + 4 + 1;
				carry_bit = 0;
				out_nums = 1;
				prefix_add_vals = ( Q1_Int >> 1);  
				BS_out =  std::move((*this->ProductAndBKAndKS(params
        , LWEscheme
				, params->GetACCRelu(5), EK
        , Res_Poly[1]
        , (*BS_out_vec)[0]
        , 0, start_bits, prefix_add_vals, carry_bit, out_nums))[0]);
				*/
				// Some Seconds
				
				for (uint32_t j = 0; j < K_FP+1; j++) {
					Outs[j] += (*BS_out_vec)[1][j];
				}
			}
		}

		
		// Sums
		
		for (uint32_t j = 0; j < K_FP+1; j++) {
			if (j == K_FP) {
				Outs[j] +=  ct2_B;
			} else {
				Outs[j] += ct2_A[j]; 
			}
		}	
		return std::make_shared<vector<DCRTPoly>>(Outs);
	}






  std::shared_ptr<vector<DCRTPoly>> FPRingGSWAccumulatorScheme::BootstrapRelu(
    const std::shared_ptr<FPRingGSWCryptoParams> params, 
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<DCRTPoly> & ct1_A,
		const DCRTPoly &ct1_B,
		const vector<DCRTPoly> & ct2_A,
		const DCRTPoly &ct2_B,
		const std::shared_ptr<FPRingGSWEvalKey> EK) const {
	
		auto MulPrepare_start = std::chrono::system_clock::now();
    uint32_t K_FP = params->GetRLWEParams()->GetK();
		NativeInteger Q0 = params->GetModuliQ()[0];
    NativeInteger Q1 = params->GetModuliQ()[1];
		uint64_t Q1_Int = Q1.ConvertToInt();
		uint32_t Q1_WO_bits = params->GetRLWEParams()->GetBits()[1];  // <32 bit
		uint32_t start_bits;
		uint64_t prefix_add_vals;
		uint32_t carry_bit;
		uint32_t out_nums;
	
		std::shared_ptr<vector<vector<DCRTPoly>>> BS_out_vec;
		vector<DCRTPoly> BS_out;
		vector<NativeVector> A(K_FP);
		NativeInteger b;
		DCRTPoly DCRT_tmp;
		vector<DCRTPoly> In_Poly(K_FP+1);
	
		vector<vector<DCRTPoly>> Res_Poly(2);
		vector<DCRTPoly> Outs(K_FP+1);
		// ct1 - ct 2
		for (uint32_t i = 0; i < K_FP+1; i++) {
			if (i == K_FP) {
				In_Poly[i] =  ct1_B;
				In_Poly[i] -= ct2_B;
			} else {
				In_Poly[i] =  ct1_A[i];
				In_Poly[i] -= ct2_A[i];
			}
		}
		
		// Substract gives range - 11bit~ 12 bit
		// add 13 bit, bootstrapping, and check whether 13bit is 0 or not
		uint32_t scale_bit = params->GetExpBit();	
		uint64_t prefix = params->GetRLWEParams()->GetScaling().ConvertToInt();
			prefix <<= scale_bit;
		uint64_t bit_check = prefix;	
		bit_check <<= 10;  // 2^10 ? 2^ 11?	
		
		In_Poly[K_FP].SetFormat(Format::COEFFICIENT);
		In_Poly[K_FP].GetElementW(0)[0].ModAddEq(
				NativeInteger(bit_check).ModMul(Q1, Q0),Q0);
		In_Poly[K_FP].GetElementW(1)[0].ModAddEq(
				NativeInteger(bit_check).ModMul(Q1, Q1),Q1);
		In_Poly[K_FP].SetFormat(Format::EVALUATION);
		

		////////////////////////////////   STARTS !!!! /////////////////////
		std::chrono::duration<double> MulPrepare_end = std::chrono::system_clock::now() - MulPrepare_start;
		if (DEBUG) std::cout << "bootstrapping Relu prepare time is " << MulPrepare_end.count() << std::endl ;
		
		// Bootstrap 1
		for (uint32_t i = 0; i < 3; i ++) { 
			
			//std::cout << "i is " << i << std::endl;	
			// Make Polys
			if (i != 2) {
				// Copy A
				for (uint32_t j = 0; j < K_FP; j++) {
					DCRT_tmp = In_Poly[j];
					DCRT_tmp.SetFormat(Format::COEFFICIENT);
					DivideRoundingSelf(params, &DCRT_tmp, Format::COEFFICIENT);
					A[j] = DCRT_tmp.GetElementAtIndex(0).GetValues();
				}
				// Copy B
				DCRT_tmp = In_Poly[K_FP];
				DCRT_tmp.SetFormat(Format::COEFFICIENT);
				DivideRoundingSelf(params, &DCRT_tmp, Format::COEFFICIENT);
				b = DCRT_tmp.GetElementAtIndex(0)[0];
			} else { // i == 2
				for (uint32_t j = 0; j < K_FP; j++) {
					In_Poly[j].SetFormat(Format::COEFFICIENT);
					DivideRoundingSelf(params, &In_Poly[j], Format::COEFFICIENT);
					A[j] = In_Poly[j].GetElementAtIndex(0).GetValues();
				}
				In_Poly[K_FP].SetFormat(Format::COEFFICIENT);
				DivideRoundingSelf(params, &In_Poly[K_FP], Format::COEFFICIENT);
				b = In_Poly[K_FP].GetElementAtIndex(0)[0];
			}
				
			// First Carry, 4+1 (1 bit no information)bit bootstrapping Bootstrapping
      
			if (i != 2) {
				if (i == 0) {
					start_bits  = scale_bit  + 4 + 1;
					prefix_add_vals = (prefix >> 1);  
				} else if (i == 1) {
					start_bits  = scale_bit  + 8 + 1;
					prefix_add_vals = (prefix << 3);  
				} 
				carry_bit = 1;
				out_nums = 2;
			
				//First Bootstrapping,3 bit + random 1 bit is bootstrapping
				BS_out_vec = this->Bootstrap_and_KS(params
					, LWEscheme
					, params->GetACCRelu(0)
					, EK
					, A
				,	 b, 0, start_bits, prefix_add_vals, carry_bit, out_nums);
		
				//Second Bootsrapp
				// i = 1, 2, 3 ...
				const DCRTPoly& ACCs = params -> GetACCRelu(i+1);
				// Orign, 4+1(padding but no infomation) bit bootstrapping msg <= 3bit 
				start_bits  = Q1_WO_bits + 4 + 1;
				carry_bit = 1;
				out_nums = 2;
				prefix_add_vals = ( Q1_Int >> 1);  
				// ================   EVAL MODE SAVING ==================== 
				BS_out_vec = this->ProductAndBKAndKS(params
					, LWEscheme
					, ACCs, EK
					, (*BS_out_vec)[0]
					, (*BS_out_vec)[1]
					, 0, start_bits, prefix_add_vals, carry_bit, out_nums);
			} else if (i == 2) {
				start_bits  = scale_bit  + 12 + 1;
				prefix_add_vals = (prefix << 7);  		
				carry_bit = 1;
				out_nums = 2;
				//First Bootstrapping,3 bit + random 1 bit is bootstrapping
				BS_out_vec = this->Bootstrap_and_KS(params
					, LWEscheme
					, params->GetACCRelu(3)
					, EK
					, A
				,	 b, 0, start_bits, prefix_add_vals, carry_bit, out_nums);
			}

			// After Handle
			//std::cout << "handle before " << std::endl;
			if (i != 2) {
				
				// remove 4 or 8 bit info
				for (uint32_t j = 0; j < K_FP+1; j++) {
					In_Poly[j]-=(*BS_out_vec)[1][j]; 
				}

				//std::cout << "Minus is end" << std::endl;
				// Saves results 
				Res_Poly[i] = std::move((*BS_out_vec)[0]);
				//std::cout << "Res is moved" << std::endl;
				
			} else {// i == 2
				
				// first 4 bit
				start_bits  = Q1_WO_bits + 4 + 1;
				carry_bit = 0;
				out_nums = 1;
				prefix_add_vals = ( Q1_Int >> 1);  
				
				Outs =  std::move((*this->ProductAndBKAndKS(params
        , LWEscheme
				, params->GetACCRelu(4), EK
        , Res_Poly[0]
        , (*BS_out_vec)[1]
        , 0, start_bits, prefix_add_vals, carry_bit, out_nums))[0]);
        
				// Second 4 bit
				start_bits  = Q1_WO_bits + 4 + 1;
				carry_bit = 0;
				out_nums = 1;
				prefix_add_vals = ( Q1_Int >> 1);  
				BS_out =  std::move((*this->ProductAndBKAndKS(params
        , LWEscheme
				, params->GetACCRelu(5), EK
        , Res_Poly[1]
        , (*BS_out_vec)[1]
        , 0, start_bits, prefix_add_vals, carry_bit, out_nums))[0]);
        
				// Some Seconds & Third 
				
				for (uint32_t j = 0; j < K_FP+1; j++) {
					Outs[j] += BS_out[j];
					Outs[j] += (*BS_out_vec)[0][j];
				}			
			}

		}

		
		// Sums
		
		for (uint32_t j = 0; j < K_FP+1; j++) {
			if (j == K_FP) {
				Outs[j] +=  ct2_B;
			} else {
				Outs[j] += ct2_A[j]; 
			}
		}	
		return std::make_shared<vector<DCRTPoly>>(Outs);
	}




/*
  std::shared_ptr<vector<DCRTPoly>> FPRingGSWAccumulatorScheme::BootstrapRecon(
    const std::shared_ptr<FPRingGSWCryptoParams> params, 
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<DCRTPoly> & ct,
		const std::shared_ptr<FPRingGSWEvalKey> EK,
		const uint32_t from_idx,
		const uint32_t to_idx) const {
	
		auto MulPrepare_start = std::chrono::system_clock::now();
    uint32_t K_FP = params->GetRLWEParams()->GetK();
	  uint32_t M_FP = params->GetRLWEParams()->GetN() * 2;
		
		NativeInteger Q0 = params->GetModuliQ()[0];
    NativeInteger Q1 = params->GetModuliQ()[1];
		uint64_t Q1_Int = Q1.ConvertToInt();
		uint32_t Q1_WO_bits = params->GetRLWEParams()->GetBits()[1];  // <32 bit
		
		uint32_t start_bits;
		uint64_t prefix_add_vals;
		uint32_t carry_bit;
		uint32_t out_nums;
	
		std::shared_ptr<vector<vector<DCRTPoly>>> BS_out_vec;
		vector<DCRTPoly> BS_out;
		vector<NativeVector> A(K_FP);
		NativeVector B;

		NativeInteger b;
		DCRTPoly DCRT_tmp;

		for (uint32_t i = 0; i < K_FP; i++) {
			DCRT_tmp = ct[i];
			DCRT_tmp.SetFormat(Format::COEFFICIENT);
			A[i] = DCRT_tmp.GetElementAtIndex(0).GetValues();
		}
		DCRT_tmp = ct[K_FP];
		DCRT_tmp.SetFormat(Format::COEFFICIENT);
		B = DCRT_tmp.GetElementAtIndex(0).GetValues();
		
		// Substract gives range - 11bit~ 12 bit
		// add 13 bit, bootstrapping, and check whether 13bit is 0 or not
		uint32_t scale_bit = params->GetExpBit();	
		uint64_t prefix = params->GetRLWEParams()->GetScaling().ConvertToInt();
			prefix <<= scale_bit;
		
		////////////////////////////////   STARTS !!!! /////////////////////
		std::chrono::duration<double> MulPrepare_end = std::chrono::system_clock::now() - MulPrepare_start;
		if (DEBUG) std::cout << "bootstrapping Recon prepare time is " << MulPrepare_end.count() << std::endl ;
	
		

		uint32_t boot_nums = to_idx - from_idx + 1;
		vector<DCRTPoly> Adding_Vals(K_FP+1);
		vector<vector<DCRTPoly>> Zero_Or_Not(boot_nums);
		for (uint32_t i = 0; i < boot_nums; i++) {
			Zero_Or_Not[i].resize(K_FP+1);
		}


		for (uint32_t i = from_idx; i < to_idx+1; i ++) {
			uint32_t start_bits_loc  = Q1_WO_bits + 2 + 1;
			uint32_t carry_bit_loc = 1;
			uint32_t out_nums_loc = 2;
			uint64_t prefix_add_vals_loc = ( Q1_Int >> 1);  
			
			uint32_t out_idx_loc = i - from_idx;
	
			NativeInteger b_tmp = B[i];

			Zero_Or_Not[out_idx_loc] = std::move((*this->Bootstrap_and_KS(params
				, LWEscheme
				, params->GetACCReconPoly(0)
				, EK
				, A
				, b_tmp, i, start_bits_loc, prefix_add_vals_loc, carry_bit_loc, out_nums_loc))[0]);
				
			if (i != 0) {
				uint32_t rot_idx = M_FP - i;
				for (uint32_t kk = 0; kk < K_FP+1; kk++) {
					Zero_Or_Not[out_idx_loc][kk] *= params->GetRotMonomial_RLWE(rot_idx); 
				}
			}
			if (i == to_idx) {
				for (uint32_t kk = 0; kk < K_FP+1; kk++) {
					//Adding_Vals[kk] = Zero_Or_Not[out_idx_loc][kk]; 
				}
	
			}
		}
	

		// Adding vals
		vector<NativeVector> A_tmp(K_FP);
		//vector<NativeVector> A_tmp2(K_FP);
		NativeInteger b_tmp;
		//NativeInteger b_tmp2;
		

		uint32_t last_idx = boot_nums - 1;
		//std::cout << "last idx is " << last_idx << std::endl;
		
		std::shared_ptr<vector<DCRTPoly>> mul_res;
		
		for (uint32_t i = 0; i < boot_nums; i++) {
			uint32_t real_idx = last_idx - i;

			// First Product
			std::shared_ptr<vector<DCRTPoly>> mul_res2;
			if (params->GetPARALLEL()) {
				mul_res = this->ProductDCRTPARALLEL2(params, EK->Evkey, 
						Zero_Or_Not[last_idx], Zero_Or_Not[last_idx][K_FP], Zero_Or_Not[real_idx], Zero_Or_Not[real_idx][K_FP]);
			} else {
				mul_res = this->ProductDCRT(params, EK->Evkey, 
						Zero_Or_Not[last_idx], Zero_Or_Not[last_idx][K_FP], Zero_Or_Not[real_idx], Zero_Or_Not[real_idx][K_FP]);
			}
			// Maybe...


			// Adding & reduction
			for (uint32_t k = 0; k < K_FP+1; k++) {
		
				DivideRoundingSelf(params, &((*mul_res)[k]), Format::COEFFICIENT); 
				if (k == K_FP) {
					b_tmp = (*mul_res)[K_FP].GetElementAtIndex(0).GetValues()[0];
				} else {
					A_tmp[k] = std::move((*mul_res)[k].GetElementAtIndex(0).GetValues());
				}
			}
			// Bootstrap
			start_bits  = Q1_WO_bits + 2 + 1;
			carry_bit = 1;
			out_nums = 2;
			prefix_add_vals = ( Q1_Int >> 1);  
		
			auto tmp_vals = this->Bootstrap_and_KS(params
					, LWEscheme
					, params->GetACCReconPoly(1)
					, EK
					, A_tmp
					, b_tmp, 0, start_bits, prefix_add_vals, carry_bit, out_nums);
			
			Zero_Or_Not[last_idx] = std::move((*tmp_vals)[0]);
			
			if (i == 0 ) {
				Adding_Vals = std::move((*tmp_vals)[1]);
			} else {
				for (uint32_t k = 0; k < K_FP+1; k++) {
					Adding_Vals[k] += (*tmp_vals)[1][k];
				}
			}
			

			//if ( i == from_idx) {
			if ( real_idx == from_idx) {
				//std::cout << "Recon Poly is processing ... " << std::endl;
				for (uint32_t k = 0; k < K_FP+1; k++) {
					DivideRoundingSelf(params, &Adding_Vals[k], Format::COEFFICIENT); 
					if (k == K_FP) {
						b_tmp = Adding_Vals[K_FP].GetElementAtIndex(0).GetValues()[0];
					} else {
						A_tmp[k] = std::move(Adding_Vals[k].GetElementAtIndex(0).GetValues());
					}
				}
				start_bits  = Q1_WO_bits + 5 + 1;
				carry_bit = 0;
				out_nums = 1;
				prefix_add_vals = ( Q1_Int >> 1);  
				
				Adding_Vals = std::move((*this->Bootstrap_and_KS(params
					, LWEscheme
					, params->GetACCReconPoly(2)
					, EK
					, A_tmp
					, b_tmp, 0, start_bits, prefix_add_vals, carry_bit, out_nums))[0]);
				break;
			}

		}
	
		return std::make_shared<vector<DCRTPoly>>(Adding_Vals);

	}
*/


  std::shared_ptr<vector<vector<DCRTPoly>>> FPRingGSWAccumulatorScheme::BootstrapRecon2(
    const std::shared_ptr<FPRingGSWCryptoParams> params, 
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<vector<DCRTPoly>> & ct,
		const std::shared_ptr<FPRingGSWEvalKey> EK,
		const uint32_t from_idx) const {
	
		auto MulPrepare_start = std::chrono::system_clock::now();
    uint32_t K_FP = params->GetRLWEParams()->GetK();
	  //uint32_t M_FP = params->GetRLWEParams()->GetN() * 2;
		
		//NativeInteger Q0 = params->GetModuliQ()[0];
    NativeInteger Q1 = params->GetModuliQ()[1];
		uint64_t Q1_Int = Q1.ConvertToInt();
		uint32_t Q1_WO_bits = params->GetRLWEParams()->GetBits()[1];  // <32 bit
		
		uint32_t start_bits;
		uint64_t prefix_add_vals;
		uint32_t carry_bit;
		uint32_t out_nums;
	
		std::shared_ptr<vector<vector<DCRTPoly>>> BS_out_vec;
		vector<DCRTPoly> BS_out;
		vector<NativeVector> A(K_FP);
		NativeVector B;

		NativeInteger b;
		DCRTPoly DCRT_tmp;
	
		// Substract gives range - 11bit~ 12 bit
		// add 13 bit, bootstrapping, and check whether 13bit is 0 or not
		uint32_t scale_bit = params->GetExpBit();	
		uint64_t prefix = params->GetRLWEParams()->GetScaling().ConvertToInt();
			prefix <<= scale_bit;
		
		////////////////////////////////   STARTS !!!! /////////////////////
		std::chrono::duration<double> MulPrepare_end = std::chrono::system_clock::now() - MulPrepare_start;
		if (DEBUG) std::cout << "bootstrapping Recon prepare time is " << MulPrepare_end.count() << std::endl ;
	
		

		uint32_t boot_nums = ct.size();
		vector<DCRTPoly> Adding_Vals(K_FP+1);
		vector<DCRTPoly> Zero_Or_Not_Res(K_FP+1);
		vector<vector<DCRTPoly>> Results(2);

		// Adding vals
		vector<NativeVector> A_tmp(K_FP);
		//vector<NativeVector> A_tmp2(K_FP);
		NativeInteger b_tmp;
		//NativeInteger b_tmp2;
		
		uint32_t last_idx = boot_nums - 1;
		for (uint32_t kk = 0; kk < K_FP+1; kk++) {
			Zero_Or_Not_Res[kk] = ct[last_idx][kk];
		}

		//std::cout << "last idx is " << last_idx << std::endl;
		
		std::shared_ptr<vector<DCRTPoly>> mul_res;
		
		for (uint32_t i = 0; i < boot_nums; i++) {
			
			uint32_t real_idx = last_idx - i;
		
			// First Product
			std::shared_ptr<vector<DCRTPoly>> mul_res2;
			if (params->GetPARALLEL()) {
				mul_res = this->ProductDCRTPARALLEL2(params, EK->Evkey, 
						Zero_Or_Not_Res, Zero_Or_Not_Res[K_FP], ct[real_idx], ct[real_idx][K_FP]);
			} else {
				mul_res = this->ProductDCRT(params, EK->Evkey, 
						Zero_Or_Not_Res, Zero_Or_Not_Res[K_FP], ct[real_idx], ct[real_idx][K_FP]);
			}
			// Adding & reduction
			for (uint32_t k = 0; k < K_FP+1; k++) {
		
				DivideRoundingSelf(params, &((*mul_res)[k]), Format::COEFFICIENT); 
				if (k == K_FP) {
					b_tmp = (*mul_res)[K_FP].GetElementAtIndex(0).GetValues()[0];
				} else {
					A_tmp[k] = std::move((*mul_res)[k].GetElementAtIndex(0).GetValues());
				}
			}
			// Bootstrap
			start_bits  = Q1_WO_bits + 2 + 1;
			carry_bit = 1;
			out_nums = 2;
			prefix_add_vals = ( Q1_Int >> 1);  
		
			auto tmp_vals = this->Bootstrap_and_KS(params
					, LWEscheme
					, params->GetACCReconPoly(1)
					, EK
					, A_tmp
					, b_tmp, 0, start_bits, prefix_add_vals, carry_bit, out_nums);
			
			Zero_Or_Not_Res = std::move((*tmp_vals)[0]);
			
			if (i == 0 ) {
				Adding_Vals = std::move((*tmp_vals)[1]);
			} else {
				for (uint32_t k = 0; k < K_FP+1; k++) {
					Adding_Vals[k] += (*tmp_vals)[1][k];
				}
			}
				
			if ( real_idx == from_idx) {
				//std::cout << "Recon Poly is processing ... " << std::endl;
				for (uint32_t k = 0; k < K_FP+1; k++) {
					DivideRoundingSelf(params, &Adding_Vals[k], Format::COEFFICIENT); 
					if (k == K_FP) {
						b_tmp = Adding_Vals[K_FP].GetElementAtIndex(0).GetValues()[0];
					} else {
						A_tmp[k] = std::move(Adding_Vals[k].GetElementAtIndex(0).GetValues());
					}
				}
				start_bits  = Q1_WO_bits + 5 + 1;
				carry_bit = 0;
				out_nums = 1;
				prefix_add_vals = ( Q1_Int >> 1);  
	
				Adding_Vals = std::move((*this->Bootstrap_and_KS(params
					, LWEscheme
					, params->GetACCReconPoly(2)
					, EK
					, A_tmp
					, b_tmp, 0, start_bits, prefix_add_vals, carry_bit, out_nums))[0]);
			
				break;
			}
		}
		
		// Checking Underflows
	
		/*
		for (uint32_t k = 0; k < K_FP+1; k++) {
			DivideRoundingSelf(params, &Zero_Or_Not_Res[k], Format::COEFFICIENT); 
			if (k == K_FP) {
				b_tmp = Zero_Or_Not_Res[K_FP].GetElementAtIndex(0).GetValues()[0];
			} else {
				A_tmp[k] = std::move(Zero_Or_Not_Res[k].GetElementAtIndex(0).GetValues());
			}
		}

		
		start_bits  = Q1_WO_bits + 5 + 1;
		carry_bit = 0;
		out_nums = 1;
		prefix_add_vals = ( Q1_Int >> 1);  
	
		Zero_Or_Not_Res = std::move((*this->Bootstrap_and_KS(params
			, LWEscheme
			, params->GetACCReconPoly(3)
			, EK
			, A_tmp
			, b_tmp, 0, start_bits, prefix_add_vals, carry_bit, out_nums))[0]);
		*/
		Results[0] = std::move(Adding_Vals);
		Results[1] = std::move(Zero_Or_Not_Res);
		return std::make_shared<vector<vector<DCRTPoly>>> (Results);
	}






  std::shared_ptr<vector<vector<DCRTPoly>>> FPRingGSWAccumulatorScheme::BootstrapMulDouble64(
    const std::shared_ptr<FPRingGSWCryptoParams> params, 
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<DCRTPoly> & ct,
		const std::shared_ptr<FPRingGSWEvalKey> EK,
		const uint32_t from_bits, const uint32_t to_bits,
		const int32_t zero_start_bits,
		const CTType types
  ) const {
	
		auto MulPrepare_start = std::chrono::system_clock::now();
		// Level 2 should handle
		std::shared_ptr<FPRingEvaluationKey> Evkey = EK->Evkey;
		
		uint32_t tot_len_bits = to_bits - from_bits + 1;
		uint32_t Out_len_bits = (int32_t) to_bits - zero_start_bits + 1; 
		if (Out_len_bits < 0) {
            PALISADE_THROW(config_error, "int32_t ...");
        }
        uint32_t end_idx = to_bits * 2 + 1;
		uint32_t start_idx = end_idx - tot_len_bits + 1;

		//std::cout <<"Out_len_bits is " <<Out_len_bits << std::endl;

		uint32_t N_FP = params->GetRLWEParams()->GetN();            // < 32bit 
    uint32_t K_FP = params->GetRLWEParams()->GetK();
		

    NativeInteger Q0 = params->GetModuliQ()[0];
    NativeInteger Q1 = params->GetModuliQ()[1];
    vector<uint64_t> mul_info;
		if (types == CTType::DOUBLE) {
			mul_info = params->GetMulFloatInfo(2); // 6 is 64 bit 
		} else if (types == CTType::FLOAT) {
			mul_info = params->GetMulFloatInfo(1); // 6 is 64 bit 
		}
		//std::cout << "Mul info len is " << mul_info.size() << std::endl;
  	//std::cout << "Iter idx is " << end_idx - start_idx + 2 << std::endl;
    // Merge
    vector<vector<DCRTPoly>> CarryInfo(tot_len_bits);  // Should Eval Mode
    vector<vector<DCRTPoly>> Out_Poly(tot_len_bits);
		//vector<vector<DCRTPoly>> Out_IsZero(tot_len_bits+4 + 1);
		vector<vector<DCRTPoly>> Out_Tot(Out_len_bits+1);


		for (uint32_t k = 0; k < tot_len_bits; k++) {
			CarryInfo[k].resize(K_FP+1);
			Out_Poly[k].resize(K_FP+1);
		}
	
		for (uint32_t k = 0; k < Out_len_bits+1; k++) {
			Out_Tot[k].resize(K_FP+1);
		}
		//vector<DCRTPoly> CarryInfo_B(EN_len); // Should Eva Mode
    
    vector<bool> Carry_occur(tot_len_bits);          // Carry occur checking
    for (uint32_t idx = 0; idx < tot_len_bits; idx ++) {
      Carry_occur[idx] = false;
    }


    // Evaluation Mode
		vector<NativePoly> a_poly_const(K_FP); 
    for (uint32_t k = 0; k < K_FP; k++) {
      a_poly_const[k] = ct[k].GetElementAtIndex(0);
      a_poly_const[k].SetFormat(Format::EVALUATION);    // EVALS!!
    }
    
		NativePoly b_poly_const = ct[K_FP].GetElementAtIndex(0);
		b_poly_const.SetFormat(Format::COEFFICIENT);  // COEFF!!
    
    uint32_t Q1_WO_bits = params->GetRLWEParams()->GetBits()[1];  // <32 bit
    uint64_t Q1_Int = Q1.ConvertToInt();
    uint32_t carry_bit;                           
    uint32_t out_nums;  
    uint32_t start_bits; 
    uint64_t prefix_add_vals; 
    uint32_t In_Out_Index = 0;

    std::shared_ptr<vector<vector<DCRTPoly>>> BS_out_vec;
		vector<DCRTPoly> BS_out;
		vector<NativePoly> A_Poly_tmp(K_FP);
    NativePoly B_Poly_tmp;
		std::shared_ptr<vector<NativeVector>> A;
    vector<NativeVector> A_tmp(K_FP);
    NativeInteger b;
    NativeInteger b_tmp;
		
		DCRTPoly DCRT_tmp;
    NativePoly Native_tmp;
		DCRTPoly ACCs;
	
		////////////////////////////////   STARTS !!!! /////////////////////
		uint32_t rot_idx;
		//uint32_t Test_len = 1;
		std::chrono::duration<double> MulPrepare_end = std::chrono::system_clock::now() - MulPrepare_start;
		if (DEBUG) std::cout << "bootstrapping Mul prepare time is " << MulPrepare_end.count() << std::endl ;

		// result is upto 62 bit. carry info of 62 bit will be updated to 63 bit
		for (uint32_t idx = start_idx; idx < end_idx + 1; idx++) {
   	  uint32_t carry_idx = idx - start_idx;
			int32_t Zero_idx = (int32_t)(idx - start_idx) - zero_start_bits;
			//int32_t out_idx = (int32_t)carry_idx - 4;
			
			rot_idx = (idx==0) ? 0 : (2*N_FP-idx);
			// A Carry & Rotation Handle
			A = RotandCarryAddA( params,
					a_poly_const, 
					&CarryInfo[carry_idx],
					params->GetRotMonomial_RLWE(rot_idx).GetElementAtIndex(0), 
					(idx != 0),				    
					Carry_occur[carry_idx]);

			// Make In & Out Index 0 by multiplying a 
			b = b_poly_const.GetValues()[idx];

      // Carry Merginga
      if (Carry_occur[carry_idx]) {
				CarryInfo[carry_idx][K_FP].SetFormat(Format::COEFFICIENT);
				DivideRoundingSelf(params, &CarryInfo[carry_idx][K_FP], Format::COEFFICIENT);
				b = b.ModAddEq(CarryInfo[carry_idx][K_FP].GetElementAtIndex(0).GetValues()[0] ,params->GetModuliQ()[0]);
      }


			// 4 bit handle
      // 4 Bit is empty 
      if (mul_info[carry_idx] <= 2) {
				//First Bootstrapping,3 bit + random 1 bit is bootstrapping
				start_bits  = Q1_WO_bits  + 2 + 1; // Programmed bit is need 
				prefix_add_vals = (Q1_Int >> 1);  
				if (Zero_idx >= 0) {
					carry_bit = 1;	
					out_nums    = 2;
					ACCs = params->GetACCMul(7);
				} else {
					carry_bit = 0;	
					out_nums    = 1;
					ACCs = params->GetACCMul(6);
	
				}
				BS_out_vec = this->Bootstrap_and_KS(params
				, LWEscheme
				, ACCs
				, EK
				, (*A)
				, b, In_Out_Index, start_bits, prefix_add_vals, carry_bit, out_nums);
        Out_Poly[carry_idx] = std::move((*BS_out_vec)[0]); //BS_out_vec[0] is removed
      
				if (Zero_idx >= 0) {
					Out_Tot[Zero_idx] = std::move((*BS_out_vec)[1]);
				}

			} else {
				// mul_info[carry_idx] >= 3
				if (mul_info[carry_idx] == 3) {
					//First Bootstrapping,3 bit + random 1 bit is bootstrapping
					start_bits  = Q1_WO_bits  + 3 + 1; // Programmed bit is need 
					prefix_add_vals = (Q1_Int >> 1);  
					//const DCRTPoly ACCs;
					if (Zero_idx >= 0) {
						carry_bit = 2;	
						out_nums    = 3;
						ACCs = params->GetACCMul(9);
					} else {
						carry_bit = 1;	
						out_nums    = 2;
						ACCs = params->GetACCMul(8);
	
					}
					BS_out_vec = this->Bootstrap_and_KS(params
					, LWEscheme
					, ACCs
					, EK
					, (*A)
					, b, In_Out_Index, start_bits, prefix_add_vals, carry_bit, out_nums);
					Out_Poly[carry_idx] = std::move((*BS_out_vec)[0]); //BS_out_vec[0] is removed
      
					if (Zero_idx >= 0) {
						Out_Tot[Zero_idx] = std::move((*BS_out_vec)[2]);
					}
					// Continue ...
					BS_out = std::move((*BS_out_vec)[1]);
        
				} else {
					// mul_info[carry_idx] >= 4
			
					//First Bootstrapping,3 bit + random 1 bit is bootstrapping
					start_bits  = Q1_WO_bits  + 3 + 1; // Programmed bit is need 
					carry_bit   = 2;
					out_nums    = 4;
					prefix_add_vals = (Q1_Int >> 1);  

					BS_out_vec = this->Bootstrap_and_KS(params
					, LWEscheme
					, params->GetACCMul(0)
					, EK
					, (*A)
					, b, In_Out_Index, start_bits, prefix_add_vals, carry_bit, out_nums);
        
					// remove 5 bit info
					// Orign, 4+1(padding but no infomation) bit bootstrapping msg <= 3bit 
					start_bits  = Q1_WO_bits + 2 + 1;
					carry_bit = 1;
					prefix_add_vals = (Q1_Int >> 1);  
					if (Zero_idx >= 0) {
						out_nums    = 2;
					}	else {
						out_nums    = 1;
					}
					// ================   EVAL MODE SAVING ==================== 
			

					auto tmp_Res = this->ProductAndBKAndKS(params
						, LWEscheme
						, params->GetACCMul(7),EK
						, (*BS_out_vec)[0]
						, (*BS_out_vec)[2]
						, In_Out_Index, start_bits, prefix_add_vals, carry_bit, out_nums);		
					
					Out_Poly[carry_idx] = std::move((*tmp_Res)[0]); //BS_out_vec[0] is removed
					if (Zero_idx >= 0) {
						Out_Tot[Zero_idx] = std::move((*tmp_Res)[1]);
					}
					// Carry handle 
			
					// First Carry, 4+1 (1 bit no information)bit bootstrapping 
					start_bits  = Q1_WO_bits  + 4 + 1;
					carry_bit = 0;
					out_nums = 1;
					prefix_add_vals = (Q1_Int >> 1);  
					// 0 and 1 Q1Q1 -> Making Carry
					BS_out  = std::move((*this->ProductAndBKAndKS(params
					, LWEscheme
					, params->GetACCMul(3),EK
					, (*BS_out_vec)[1]
					, (*BS_out_vec)[2]
					, In_Out_Index, start_bits, prefix_add_vals, carry_bit, out_nums))[0]);
				
					// Modifiy, 3 is not changed
					auto mf1_start = std::chrono::system_clock::now();
	
					for (uint32_t k = 0; k < K_FP; k++) {
						BS_out[k] -= (*BS_out_vec)[3][k];
						BS_out[k].SetFormat(Format::COEFFICIENT);
						DivideRoundingSelf(params, &BS_out[k], Format::COEFFICIENT);
						A_tmp[k] = std::move(BS_out[k].GetElementAtIndex(0).GetValues());
					}
					BS_out[K_FP] -= (*BS_out_vec)[3][K_FP];
					BS_out[K_FP].SetFormat(Format::COEFFICIENT);
					DivideRoundingSelf(params, &BS_out[K_FP], Format::COEFFICIENT);
					b_tmp = BS_out[K_FP].GetElementAtIndex(0).GetValues()[0];
		    
					std::chrono::duration<double> mf1_end = std::chrono::system_clock::now() - mf1_start;
					if (DEBUG) std::cout << "Modify1 time is " << mf1_end.count() << std::endl;

			 
					// First Carry, 4+1 (1 bit no information)bit bootstrapping 
          start_bits  = Q1_WO_bits + 4 + 1;
          carry_bit = 0;
          out_nums = 1;
          prefix_add_vals = (Q1_Int >> 1);  
          BS_out = std::move((*this->Bootstrap_and_KS(params
            , LWEscheme
            , params->GetACCMul(3),EK
            , A_tmp
            , (b_tmp.ModAdd(Q1,Q0))
            , In_Out_Index, start_bits, prefix_add_vals, carry_bit, out_nums))[0]);
				}
		

				// Carry save
				if (Carry_occur[carry_idx+1]) {
					for (uint32_t k = 0; k < K_FP+1; k++) {
						CarryInfo[carry_idx+1][k] += BS_out[k];
					}
				} else {
					Carry_occur[carry_idx+1] = true;
					for (uint32_t k = 0; k < K_FP+1; k++) {
						CarryInfo[carry_idx+1][k] = BS_out[k];
					}
				}
				
        //std::cout << idx << "'s first carry is saved ..." << std::endl;

				//Second Carry Handle
        // If Carry should calculates
        if (mul_info[carry_idx] > 4) {
          // Subs Carry firsts and Orign
      		// Modifiy
					auto mf2_start = std::chrono::system_clock::now();
			  
					for (uint32_t k = 0; k < K_FP; k++) {
						BS_out[k] -= (*BS_out_vec)[3][k];
						BS_out[k].SetFormat(Format::COEFFICIENT);
						DivideRoundingSelf(params, &BS_out[k], Format::COEFFICIENT);
						(*A)[k] -= (BS_out[k].GetElementAtIndex(0).GetValues()*4);
					}
					BS_out[K_FP] -= (*BS_out_vec)[3][K_FP];
					BS_out[K_FP].SetFormat(Format::COEFFICIENT);
					DivideRoundingSelf(params, &BS_out[K_FP], Format::COEFFICIENT);
					b = b.ModSubEq(BS_out[K_FP].GetElementAtIndex(0).GetValues()[0] * 4, Q0);
        
					for (uint32_t k = 0; k < K_FP; k++) {
						Native_tmp = Out_Poly[carry_idx][k].GetElementAtIndex(0); 
						Native_tmp.SetFormat(Format::COEFFICIENT);
						(*A)[k] -= Native_tmp.GetValues();
					}
					Native_tmp = Out_Poly[carry_idx][K_FP].GetElementAtIndex(0); 
					Native_tmp.SetFormat(Format::COEFFICIENT);
					b = b.ModSubEq( Native_tmp.GetValues()[0], Q0);
	        
					std::chrono::duration<double> mf2_end = std::chrono::system_clock::now() - mf2_start;
					if (DEBUG) std::cout << "Modify2 time is " << mf2_end.count() << std::endl;
          // if 9 bit is zero, then it is okey
          if (mul_info[carry_idx] <= 8) {
            // if my carry bit is 8bit, it is over 
            start_bits = Q1_WO_bits  + 4 + 4 + 1; // 4bit +4bit
            prefix_add_vals = (Q1_Int << 3); // 4 bit - 1   
            carry_bit = 0;
            out_nums = 1;

            BS_out = std::move((*this->Bootstrap_and_KS(params
              , LWEscheme
              , params->GetACCMul(3),EK
              , (*A)
              , b
              , In_Out_Index, start_bits, prefix_add_vals, carry_bit, out_nums))[0]);
          } else {

          // Carry (4bit floated) + 4bit + 1 padding
          // if my carry bit is 8bit, it is over 
            start_bits = Q1_WO_bits  + 4 + 4 + 1; // 4bit +4bit
            prefix_add_vals = (Q1_Int << 3); // 4 bit - 1   
            carry_bit = 1;
            out_nums = 2;

            BS_out_vec = this->Bootstrap_and_KS(params
            , LWEscheme
            , params->GetACCMul(1),EK
            , (*A)
            , b
            , In_Out_Index, start_bits, prefix_add_vals, carry_bit, out_nums);

            //Mul_Divide = Product_RD_Coeff(BS_out_vec[0],BS_out_vec[1]); 
            // Mul Bootstrapping, Handle Sign
            start_bits = Q1_WO_bits  + 4 + 1;
            prefix_add_vals = (Q1_Int >> 1);  
            carry_bit = 0;
            out_nums  = 1;
            
						 BS_out  = std::move((*this->ProductAndBKAndKS(params
							, LWEscheme
							, params->GetACCMul(3),EK
              , (*BS_out_vec)[0]
              , (*BS_out_vec)[1]
              , In_Out_Index, start_bits, prefix_add_vals, carry_bit, out_nums))[0]);
          }
					// Carry save
          if (Carry_occur[carry_idx+2]) {
						for (uint32_t k = 0; k < K_FP+1; k++) {
							CarryInfo[carry_idx+2][k] += BS_out[k];
						} 
		     
					} else {
            Carry_occur[carry_idx+2] = true;
            for (uint32_t k = 0; k < K_FP+1; k++) {
							CarryInfo[carry_idx+2][k] = BS_out[k];
						} 
					
					};
          //std::cout << idx << "'s second carry is saved ..." << std::endl;

          //Third Carry Handle
          if (mul_info[carry_idx] > 8) {

						auto mf3_start = std::chrono::system_clock::now();
	
						// Modifiy, 3 is not changed
						for (uint32_t k = 0; k < K_FP; k++) {
							//BS_out[k] -= (*BS_out_vec)[3][k];
							BS_out[k].SetFormat(Format::COEFFICIENT);
							DivideRoundingSelf(params, &BS_out[k], Format::COEFFICIENT);
							(*A)[k] -= BS_out[k].GetElementAtIndex(0).GetValues() * 256;
						}
						BS_out[K_FP].SetFormat(Format::COEFFICIENT);
						DivideRoundingSelf(params, &BS_out[K_FP], Format::COEFFICIENT);
						b = b.ModSubEq(BS_out[K_FP].GetElementAtIndex(0).GetValues()[0] * 256 , Q0);
          

						std::chrono::duration<double> mf3_end = std::chrono::system_clock::now() - mf3_start;	
						if(DEBUG) std::cout << "Modify3 time is " << mf3_end.count() << std::endl     ;
						//(8bit floated) + 1bit + 1bit padding 
            carry_bit = 0;
            out_nums  = 1;
            start_bits = Q1_WO_bits  + 10; 
            prefix_add_vals = (Q1_Int << 7);    
            BS_out = std::move((*this->Bootstrap_and_KS(params
              , LWEscheme
              , params->GetACCMul(4),EK
              , (*A)
              , b
              , In_Out_Index, start_bits ,prefix_add_vals, carry_bit, out_nums))[0]);
            // FFT & Carry save

            if (Carry_occur[carry_idx+4]) {
            
							for (uint32_t k = 0; k < K_FP+1; k++) {
								CarryInfo[carry_idx+4][k] += BS_out[k];
							} 
		       
						} else {
              Carry_occur[carry_idx+4] = true;
              for (uint32_t k = 0; k < K_FP+1; k++) {
								CarryInfo[carry_idx+4][k] = BS_out[k];
							} 
			    
						};
       
					}
        }
      } 
    }


		//std::cout << "Last" << std::endl;
		//std::cout << "tot_len is " << tot_len_bits << std::endl;
		if (from_bits != 0 ) {
			for (uint32_t k = 0; k < K_FP+1; k++) {
				Out_Poly[0][k] *= params->GetRotMonomial_RLWE(from_bits);
			}
		}
		
		for (uint32_t idx=1; idx < tot_len_bits; idx++) {
			//std::cout << "idx is " << idx << std::endl;
			uint32_t Rot_idx = idx + from_bits; 
			//std::cout << "Rot idx is " << Rot_idx << std::endl;
			for (uint32_t k = 0; k < K_FP+1; k++) {
				Out_Poly[0][k] += (Out_Poly[idx][k] * params->GetRotMonomial_RLWE(Rot_idx));
			}
		}
	
		// Make
		Out_Tot[Out_len_bits] = std::move(Out_Poly[0]);
		
		return std::make_shared<vector<vector<DCRTPoly>>>(std::move(Out_Tot));
  }








  std::shared_ptr<vector<DCRTPoly>> FPRingGSWAccumulatorScheme::BootstrapMul(
    const std::shared_ptr<FPRingGSWCryptoParams> params, 
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<DCRTPoly> & ct,
		const std::shared_ptr<FPRingGSWEvalKey> EK
  ) const {
	
		auto MulPrepare_start = std::chrono::system_clock::now();
		// Level 2 should handle
		std::shared_ptr<FPRingEvaluationKey> Evkey = EK->Evkey;
		
    uint32_t EN_len = (uint32_t ) std::ceil( (double) 64.0 / (double) params -> GetBaseEncodeBit() );
    uint32_t N_FP = params->GetRLWEParams()->GetN();            // < 32bit 
    uint32_t K_FP = params->GetRLWEParams()->GetK();
  
    NativeInteger Q0 = params->GetModuliQ()[0];
    NativeInteger Q1 = params->GetModuliQ()[1];
    vector<uint64_t> mul_info = params->GetMulInfo(6); // 6 is 64 bit 
    
    // Merge
    vector<vector<DCRTPoly>> CarryInfo(EN_len);  // Should Eval Mode
    for (uint32_t k = 0; k < EN_len; k++) {
      CarryInfo[k].resize(K_FP+1);
    }
		//vector<DCRTPoly> CarryInfo_B(EN_len); // Should Eva Mode
    
    vector<bool> Carry_occur(EN_len);          // Carry occur checking
    vector<vector<DCRTPoly>> Out_Poly(EN_len);

    for (uint32_t idx = 0; idx < EN_len; idx ++) {
      Carry_occur[idx] = false;
    }

    // Evaluation Mode
    vector<NativePoly> a_poly_const(K_FP); 
    for (uint32_t k = 0; k < K_FP; k++) {
      a_poly_const[k] = ct[k].GetElementAtIndex(0);
      a_poly_const[k].SetFormat(Format::EVALUATION);    // EVALS!!
    }
    
		NativePoly b_poly_const = ct[K_FP].GetElementAtIndex(0);
		b_poly_const.SetFormat(Format::COEFFICIENT);  // COEFF!!
    
    uint32_t Q1_WO_bits = params->GetRLWEParams()->GetBits()[1];  // <32 bit
    uint64_t Q1_Int = Q1.ConvertToInt();
    uint32_t carry_bit;                           
    uint32_t out_nums;  
    uint32_t start_bits; 
    uint64_t prefix_add_vals; 
    uint32_t In_Out_Index = 0;

    std::shared_ptr<vector<vector<DCRTPoly>>> BS_out_vec;
		vector<DCRTPoly> BS_out;
		vector<NativePoly> A_Poly_tmp(K_FP);
    NativePoly B_Poly_tmp;
		std::shared_ptr<vector<NativeVector>> A;
    vector<NativeVector> A_tmp(K_FP);
    NativeInteger b;
    NativeInteger b_tmp;
		
		DCRTPoly DCRT_tmp;
    NativePoly Native_tmp;
		////////////////////////////////   STARTS !!!! /////////////////////
		uint32_t rot_idx;
		//uint32_t Test_len = 1;
		std::chrono::duration<double> MulPrepare_end = std::chrono::system_clock::now() - MulPrepare_start;
		if (DEBUG) std::cout << "bootstrapping Mul prepare time is " << MulPrepare_end.count() << std::endl ;

		for (uint32_t idx = 0; idx < EN_len; idx ++) {
    //for (uint32_t idx = 0; idx < Test_len; idx ++) {
      //std::cout << idx << "'s processing is on ..." << std::endl;
			   
			rot_idx = (idx==0) ? 0 : (2*N_FP-idx);
			// A Carry & Rotation Handle
			A = RotandCarryAddA( params,
					a_poly_const, 
					&CarryInfo[idx],
					params->GetRotMonomial_RLWE(rot_idx).GetElementAtIndex(0), 
					(idx != 0 ),				    
					Carry_occur[idx]);

			// Make In & Out Index 0 by multiplying a 
			b = b_poly_const.GetValues()[idx];

      // Carry Merging
      if (Carry_occur[idx]) {
				CarryInfo[idx][K_FP].SetFormat(Format::COEFFICIENT);
				DivideRoundingSelf(params, &CarryInfo[idx][K_FP], Format::COEFFICIENT);
				b = b.ModAddEq(CarryInfo[idx][K_FP].GetElementAtIndex(0).GetValues()[0] ,params->GetModuliQ()[0]);
      }

  
      //First Bootstrapping,3 bit + random 1 bit is bootstrapping
      start_bits  = Q1_WO_bits  + 3 + 1; // Programmed bit is need 
      carry_bit   = 2;
      out_nums    = 4;
      prefix_add_vals = (Q1_Int >> 1);  

      BS_out_vec = this->Bootstrap_and_KS(params
      , LWEscheme
      , params->GetACCMul(0)
      , EK
      , (*A)
      , b, In_Out_Index, start_bits, prefix_add_vals, carry_bit, out_nums);

      // 4 Bit is empty 
      if (mul_info[idx] <= 3) {
        Out_Poly[idx] = std::move((*BS_out_vec)[0]); //BS_out_vec[0] is removed
      } else {
        // remove 5 bit info
  
        // Orign, 4+1(padding but no infomation) bit bootstrapping msg <= 3bit 
				start_bits  = Q1_WO_bits + 4 + 1;
        carry_bit = 0;
        out_nums = 1;
        prefix_add_vals = (Q1_Int >> 1);  
        // ================   EVAL MODE SAVING ==================== 
        Out_Poly[idx] = std::move((*this->ProductAndBKAndKS(params
          , LWEscheme
					, params->GetACCMul(2),EK
          , (*BS_out_vec)[0]
          , (*BS_out_vec)[2]
          , In_Out_Index, start_bits, prefix_add_vals, carry_bit, out_nums))[0]);
        
			}

      // Carry handle 
      if ((mul_info[idx] > 2) & (idx < EN_len-1)) {
        if (mul_info[idx] <= 3) {
          BS_out = std::move((*BS_out_vec)[1]);
        } else {
				
					// First Carry, 4+1 (1 bit no information)bit bootstrapping 
          start_bits  = Q1_WO_bits  + 4 + 1;
          carry_bit = 0;
          out_nums = 1;
          prefix_add_vals = (Q1_Int >> 1);  
          // 0 and 1 Q1Q1 -> Making Carry
          BS_out  = std::move((*this->ProductAndBKAndKS(params
            , LWEscheme
						, params->GetACCMul(3),EK
            , (*BS_out_vec)[1]
						, (*BS_out_vec)[2]
            , In_Out_Index, start_bits, prefix_add_vals, carry_bit, out_nums))[0]);
					
					// Modifiy, 3 is not changed
					auto mf1_start = std::chrono::system_clock::now();
	
          for (uint32_t k = 0; k < K_FP; k++) {
						BS_out[k] -= (*BS_out_vec)[3][k];
						BS_out[k].SetFormat(Format::COEFFICIENT);
						DivideRoundingSelf(params, &BS_out[k], Format::COEFFICIENT);
						A_tmp[k] = std::move(BS_out[k].GetElementAtIndex(0).GetValues());
					}
					BS_out[K_FP] -= (*BS_out_vec)[3][K_FP];
					BS_out[K_FP].SetFormat(Format::COEFFICIENT);
					DivideRoundingSelf(params, &BS_out[K_FP], Format::COEFFICIENT);
					b_tmp = BS_out[K_FP].GetElementAtIndex(0).GetValues()[0];
		    
					std::chrono::duration<double> mf1_end = std::chrono::system_clock::now() - mf1_start;
					if (DEBUG) std::cout << "Modify1 time is " << mf1_end.count() << std::endl;


          // First Carry, 4+1 (1 bit no information)bit bootstrapping 
          start_bits  = Q1_WO_bits + 4 + 1;
          carry_bit = 0;
          out_nums = 1;
          prefix_add_vals = (Q1_Int >> 1);  
          BS_out = std::move((*this->Bootstrap_and_KS(params
            , LWEscheme
            , params->GetACCMul(3),EK
            , A_tmp
            , (b_tmp.ModAdd(Q1,Q0))
            , In_Out_Index, start_bits, prefix_add_vals, carry_bit, out_nums))[0]);
        }

        // Carry save
        if (Carry_occur[idx+1]) {
          for (uint32_t k = 0; k < K_FP+1; k++) {
            CarryInfo[idx+1][k] += BS_out[k];
          }
				} else {
          Carry_occur[idx+1] = true;
					for (uint32_t k = 0; k < K_FP+1; k++) {
            CarryInfo[idx+1][k] = BS_out[k];
          }
				};

        //std::cout << idx << "'s first carry is saved ..." << std::endl;

        //Second Carry Handle
        // If Carry should calculates
        if ((mul_info[idx] > 4) & (idx < EN_len-2)) {
          // Subs Carry firsts and Orign
      		// Modifiy
					auto mf2_start = std::chrono::system_clock::now();
	
					for (uint32_t k = 0; k < K_FP; k++) {
						BS_out[k] -= (*BS_out_vec)[3][k];
						BS_out[k].SetFormat(Format::COEFFICIENT);
						DivideRoundingSelf(params, &BS_out[k], Format::COEFFICIENT);
						(*A)[k] -= (BS_out[k].GetElementAtIndex(0).GetValues()*4);
					}
					BS_out[K_FP] -= (*BS_out_vec)[3][K_FP];
					BS_out[K_FP].SetFormat(Format::COEFFICIENT);
					DivideRoundingSelf(params, &BS_out[K_FP], Format::COEFFICIENT);
					b = b.ModSubEq(BS_out[K_FP].GetElementAtIndex(0).GetValues()[0] * 4, Q0);
        
					for (uint32_t k = 0; k < K_FP; k++) {
						Native_tmp = Out_Poly[idx][k].GetElementAtIndex(0); 
						Native_tmp.SetFormat(Format::COEFFICIENT);
						(*A)[k] -= Native_tmp.GetValues();
					}
					Native_tmp = Out_Poly[idx][K_FP].GetElementAtIndex(0); 
					Native_tmp.SetFormat(Format::COEFFICIENT);
					b = b.ModSubEq( Native_tmp.GetValues()[0], Q0);

					std::chrono::duration<double> mf2_end = std::chrono::system_clock::now() - mf2_start;
					if (DEBUG) std::cout << "Modify2 time is " << mf2_end.count() << std::endl;
          // if 9 bit is zero, then it is okey
          if (mul_info[idx] <= 8) {
            // if my carry bit is 8bit, it is over 
            start_bits = Q1_WO_bits  + 4 + 4 + 1; // 4bit +4bit
            prefix_add_vals = (Q1_Int << 3); // 4 bit - 1   
            carry_bit = 0;
            out_nums = 1;

            BS_out = std::move((*this->Bootstrap_and_KS(params
              , LWEscheme
              , params->GetACCMul(3),EK
              , (*A)
              , b
              , In_Out_Index, start_bits, prefix_add_vals, carry_bit, out_nums))[0]);
          } else {

          // Carry (4bit floated) + 4bit + 1 padding
          // if my carry bit is 8bit, it is over 
            start_bits = Q1_WO_bits  + 4 + 4 + 1; // 4bit +4bit
            prefix_add_vals = (Q1_Int << 3); // 4 bit - 1   
            carry_bit = 1;
            out_nums = 2;

            BS_out_vec = this->Bootstrap_and_KS(params
            , LWEscheme
            , params->GetACCMul(1),EK
            , (*A)
            , b
            , In_Out_Index, start_bits, prefix_add_vals, carry_bit, out_nums);

            //Mul_Divide = Product_RD_Coeff(BS_out_vec[0],BS_out_vec[1]); 
            // Mul Bootstrapping, Handle Sign
            start_bits = Q1_WO_bits  + 4 + 1;
            prefix_add_vals = (Q1_Int >> 1);  
            carry_bit = 0;
            out_nums  = 1;
            
						 BS_out  = std::move((*this->ProductAndBKAndKS(params
							, LWEscheme
							, params->GetACCMul(3),EK
              , (*BS_out_vec)[0]
              , (*BS_out_vec)[1]
              , In_Out_Index, start_bits, prefix_add_vals, carry_bit, out_nums))[0]);
          }

          // Carry save
          if (Carry_occur[idx+2]) {
						for (uint32_t k = 0; k < K_FP+1; k++) {
							CarryInfo[idx+2][k] += BS_out[k];
						} 
		     
					} else {
            Carry_occur[idx+2] = true;
            for (uint32_t k = 0; k < K_FP+1; k++) {
							CarryInfo[idx+2][k] = BS_out[k];
						} 
					
					};
          //std::cout << idx << "'s second carry is saved ..." << std::endl;

          //Third Carry Handle
          if ((mul_info[idx] > 8) && ( (idx < EN_len-4))) {

						auto mf3_start = std::chrono::system_clock::now();
	
						// Modifiy, 3 is not changed
						for (uint32_t k = 0; k < K_FP; k++) {
							//BS_out[k] -= (*BS_out_vec)[3][k];
							BS_out[k].SetFormat(Format::COEFFICIENT);
							DivideRoundingSelf(params, &BS_out[k], Format::COEFFICIENT);
							(*A)[k] -= BS_out[k].GetElementAtIndex(0).GetValues() * 256;
						}
						BS_out[K_FP].SetFormat(Format::COEFFICIENT);
						DivideRoundingSelf(params, &BS_out[K_FP], Format::COEFFICIENT);
						b = b.ModSubEq(BS_out[K_FP].GetElementAtIndex(0).GetValues()[0] * 256 , Q0);
          

						std::chrono::duration<double> mf3_end = std::chrono::system_clock::now() - mf3_start;	
						std::cout << "Modify3 time is " << mf3_end.count() << std::endl     ;
						//(8bit floated) + 1bit + 1bit padding 
            carry_bit = 0;
            out_nums  = 1;
            start_bits = Q1_WO_bits  + 10; 
            prefix_add_vals = (Q1_Int << 7);    
            BS_out = std::move((*this->Bootstrap_and_KS(params
              , LWEscheme
              , params->GetACCMul(4),EK
              , (*A)
              , b
              , In_Out_Index, start_bits ,prefix_add_vals, carry_bit, out_nums))[0]);
            // FFT & Carry save

            if (Carry_occur[idx+4]) {
            
							for (uint32_t k = 0; k < K_FP+1; k++) {
								CarryInfo[idx+4][k] += BS_out[k];
							} 
		       
						} else {
              Carry_occur[idx+4] = true;
              for (uint32_t k = 0; k < K_FP+1; k++) {
								CarryInfo[idx+4][k] = BS_out[k];
							} 
			    
						};
            //std::cout << idx << "'s third carry is saved ..." << std::endl;

          }
        }
      } 
    }
    // Collect
		for (uint32_t k = 0; k < K_FP+1; k++) {
			for (uint32_t idx=1; idx < EN_len; idx++) {
				Out_Poly[0][k] += (Out_Poly[idx][k] * params->GetRotMonomial_RLWE(idx));
			}
		}

		return std::make_shared<vector<DCRTPoly>>(Out_Poly[0]);
  }

  std::shared_ptr<vector<DCRTPoly>> FPRingGSWAccumulatorScheme::BootstrapSign(
    const std::shared_ptr<FPRingGSWCryptoParams> params, 
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<DCRTPoly> & ct,
		const std::shared_ptr<FPRingGSWEvalKey> EK
  ) const {
   
		if (ct[0].GetFormat() != Format::COEFFICIENT) {
			PALISADE_THROW(config_error, "BootSign input should be Coefficient");
		}
		//First Bootstrapping,2 bit + random 1 bit is bootstrapping
    NativeInteger Q1 = params->GetModuliQ()[1];
		std::shared_ptr<FPRingEvaluationKey> Evkey = EK->Evkey;
		
		uint32_t Q1_WO_bits = params->GetRLWEParams()->GetBits()[1];  // <32 bit
    uint64_t Q1_Int = Q1.ConvertToInt();
		uint64_t prefix_add_vals = (Q1_Int >> 1) + Q1_Int; // -1, 1 -> 0, 2  
	
		uint32_t K_FP = params->GetRLWEParams()->GetK();  
		uint32_t start_bits  = Q1_WO_bits  + 2 + 1; // Programmed bit is need 
    uint32_t carry_bit   = 0; 
    uint32_t out_nums    = 1; 
		uint32_t loc = 0;
    
		vector<NativeVector> A(K_FP);
		for (uint32_t i = 0; i < K_FP; i++) {
			A[i] = ct[i].GetElementAtIndex(0).GetValues();
		}
		NativeInteger B = ct[K_FP].GetElementAtIndex(0).GetValues()[0];
		vector<DCRTPoly> BS_out = std::move((*this->Bootstrap_and_KS(params
			, LWEscheme
      , params->GetACCMul(5)
      , EK 
      , A
      , B
      , loc, start_bits, prefix_add_vals, carry_bit, out_nums))[0]);
   

		return std::make_shared<vector<DCRTPoly>>(BS_out);
  
	}


	std::shared_ptr<vector<vector<DCRTPoly>>> FPRingGSWAccumulatorScheme::BootstrapAdd(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<DCRTPoly> & ct,
		const std::shared_ptr<FPRingGSWEvalKey> EK
   ) const {
		
		//!!!!!!!!!!!!!!! Should Turn On the Max ,,,
		auto prepare_start = std::chrono::system_clock::now();
	
		// Level 2 should handle
    uint32_t EN_len = (uint32_t ) std::ceil( (double) 64.0 / (double) params -> GetBaseEncodeBit() );

		//NativeInteger Q0 = params->GetModuliQ()[0];
		NativeInteger Q1 = params->GetModuliQ()[1];
		uint32_t K_FP = params->GetRLWEParams()->GetK();
		uint32_t N_FP = params->GetRLWEParams()->GetN();
	
		// Merge
    vector<vector<DCRTPoly>> CarryInfo(EN_len);  // Should Eval Mode
    vector<vector<DCRTPoly>> Out_Poly(EN_len);
    
		for (uint32_t k = 0; k < EN_len; k++) {
      CarryInfo[k].resize(K_FP+1);
			Out_Poly[k].resize(K_FP+1);
    }
	
		// Evaluation Mode
    vector<NativePoly> a_poly_const(K_FP); 
    for (uint32_t k = 0; k < K_FP; k++) {
      a_poly_const[k] = ct[k].GetElementAtIndex(0);
      a_poly_const[k].SetFormat(Format::EVALUATION);    // EVALS!!
    }
		NativePoly b_poly_const = ct[K_FP].GetElementAtIndex(0);
		b_poly_const.SetFormat(Format::COEFFICIENT);  // COEFF!!
    
    uint32_t Q1_WO_bits = params->GetRLWEParams()->GetBits()[1];  // <32 bit
    uint64_t Q1_Int = Q1.ConvertToInt();
    uint32_t carry_bit;                           
    uint32_t out_nums;  
    uint32_t start_bits; 
    uint64_t prefix_add_vals; 
    uint32_t In_Out_Index = 0;

    std::shared_ptr<vector<vector<DCRTPoly>>> BS_out_vec;
		vector<DCRTPoly> BS_out;
    std::shared_ptr<const FPRLWECiphertextImpl> Mul_CT;
    vector<NativePoly> A_Poly_tmp(K_FP);
    NativePoly B_Poly_tmp;
 	
		std::shared_ptr<vector<NativeVector>> A;
    vector<NativeVector> A_tmp(K_FP);
    NativeInteger b;
    NativeInteger b_tmp;
	
		DCRTPoly DCRT_tmp;
    NativePoly Native_tmp;
		uint32_t rot_idx;
		std::chrono::duration<double> prepare_end = std::chrono::system_clock::now() - prepare_start;
		if (DEBUG) std::cout << "Bootstrapping Add Prepare time is " <<prepare_end.count() << std::endl;
		for (uint32_t idx = 0; idx < EN_len; idx ++) {
			// Make In & Out Index 0 by multiplying a 
			//std::cout << "first idx is " << idx << std::endl;
			auto rot_start = std::chrono::system_clock::now();
	
	
			rot_idx = (idx==0) ? 0 : (2*N_FP-idx);
			// A Carry & Rotation Handle
			A = RotandCarryAddA( params,
					a_poly_const, 
					&CarryInfo[idx],
					params->GetRotMonomial_RLWE(rot_idx).GetElementAtIndex(0), 
					(idx != 0),				    
					(idx!=0));
					// Make In & Out Index 0 by multiplying a 
		
			// Handle B
			b = b_poly_const.GetValues()[idx];
      // Carry Merging
      if (idx != 0) {
				CarryInfo[idx][K_FP].SetFormat(Format::COEFFICIENT);
				DivideRoundingSelf(params, &CarryInfo[idx][K_FP], Format::COEFFICIENT);
				b = b.ModAddEq(CarryInfo[idx][K_FP].GetElementAtIndex(0).GetValues()[0] ,params->GetModuliQ()[0]);
      }
			std::chrono::duration<double> rot_end = std::chrono::system_clock::now() - rot_start;
			if (DEBUG) std::cout << "Bootstrapping Add Rotation time is " <<rot_end.count() << std::endl;
			
			//First Bootstrapping,3 bit + random 1 bit is bootstrapping
			start_bits  = Q1_WO_bits  + 4 + 1; // Programmed bit is need 
			carry_bit   = 1;
			out_nums    = 2;
			prefix_add_vals = (Q1_Int >> 1) + 8*Q1_Int;  

			if (idx != EN_len - 1) {
				BS_out_vec = this->Bootstrap_and_KS(params
					, LWEscheme
					, params->GetACCAdd(0)
					, EK
					, (*A)
					, b, In_Out_Index, start_bits, prefix_add_vals, carry_bit, out_nums);
			} else {
				BS_out_vec = this->Bootstrap_and_KS(params
					, LWEscheme
					, params->GetACCAdd(1)
					, EK
					, (*A)
					, b, In_Out_Index, start_bits, prefix_add_vals, carry_bit, out_nums);
			} 
			Out_Poly[idx] = std::move((*BS_out_vec)[0]);

			// Carry save
			if (idx != EN_len - 1) {
				CarryInfo[idx+1] = std::move((*BS_out_vec)[1]);
			} else {
				CarryInfo[0] = std::move((*BS_out_vec)[1]);
			}

		}
		// Tmp Collection
	
		for (uint32_t k = 0; k < K_FP+1; k++) {
			for (uint32_t idx=1; idx < EN_len; idx++) {
				Out_Poly[0][k] += (Out_Poly[idx][k] * params->GetRotMonomial_RLWE(idx));
			}
		}
	
		// Products
		std::shared_ptr<vector<DCRTPoly>> mul_res; 
		if (params->GetPARALLEL()) {
			mul_res = this->ProductDCRTPARALLEL2(params, EK->Evkey, Out_Poly[0], Out_Poly[0][K_FP], CarryInfo[0], CarryInfo[0][K_FP]); 
		} else {
			mul_res = this->ProductDCRT(params, EK->Evkey, Out_Poly[0], Out_Poly[0][K_FP], CarryInfo[0], CarryInfo[0][K_FP]); 
		}
		// Reduction
	
		auto rd_start = std::chrono::system_clock::now();
	

		for (uint32_t k = 0; k < K_FP; k++) {
			(*mul_res)[k].SetFormat(Format::COEFFICIENT);
			DivideRoundingSelf(params, &((*mul_res)[k]), Format::EVALUATION); 
			a_poly_const[k] = std::move((*mul_res)[k].GetElementAtIndex(0));
		}
		(*mul_res)[K_FP].SetFormat(Format::COEFFICIENT);
		DivideRoundingSelf(params, &((*mul_res)[K_FP]), Format::COEFFICIENT); 
		b_poly_const = std::move((*mul_res)[K_FP].GetElementAtIndex(0));
	
		std::chrono::duration<double> rd_end = std::chrono::system_clock::now() - rd_start;
		if (DEBUG) std::cout << "Reduction times is " << rd_end.count() << std::endl;
		

		// Evaluation Mode
		
		//Conpensate Bits One More Times..
		for (uint32_t idx = 0; idx < EN_len; idx ++) {
			//std::cout << "second idx is " << idx << std::endl;
		 
			rot_idx = (idx==0) ? 0 : (2*N_FP-idx);
			// A Carry & Rotation Handle
			A = RotandCarryAddA( params,
					a_poly_const, 
					&CarryInfo[idx],
					params->GetRotMonomial_RLWE(rot_idx).GetElementAtIndex(0), 
					( idx != 0),				    
					(idx!=0));
					// Make In & Out Index 0 by multiplying a 
		
			b = b_poly_const.GetValues()[idx];
      // Carry Merging
      if (idx != 0) {
				CarryInfo[idx][K_FP].SetFormat(Format::COEFFICIENT);
				DivideRoundingSelf(params, &CarryInfo[idx][K_FP], Format::COEFFICIENT);
				b = b.ModAddEq(CarryInfo[idx][K_FP].GetElementAtIndex(0).GetValues()[0] ,params->GetModuliQ()[0]);
      }

			// Orign, 3+1(padding but no infomation) bit bootstrapping msg <= 3bit 
			start_bits  = Q1_WO_bits + 3 + 1;
			carry_bit = 1;
			out_nums = 2;
			prefix_add_vals = (Q1_Int >> 1) + 4*Q1_Int;  
			// ================   EVAL MODE SAVING ==================== 
			BS_out_vec= this->Bootstrap_and_KS(params
				, LWEscheme
				, params->GetACCAdd(3),EK
				, (*A)
				, b
				, In_Out_Index, start_bits, prefix_add_vals, carry_bit, out_nums);
				//  ============================================================
			Out_Poly[idx] = std::move((*BS_out_vec)[0]);
			
			// Carry save
			if (idx < EN_len-1 ) {	
				CarryInfo[idx+1] = std::move((*BS_out_vec)[1]);

			}
		}
		// Collect Again
		for (uint32_t k = 0; k < K_FP+1; k++) {
			for (uint32_t idx=1; idx < EN_len; idx++) {
				Out_Poly[0][k] += (Out_Poly[idx][k] * params->GetRotMonomial_RLWE(idx));
			}
		}
		vector<vector<DCRTPoly>> res(2);
		res[0] = std::move(Out_Poly[0]);
		res[1] = std::move(CarryInfo[0]);

		//if (types == INT64) {
		return std::make_shared<vector<vector<DCRTPoly>>>(res);
	}


	// From, to bits
	std::shared_ptr<vector<vector<DCRTPoly>>> FPRingGSWAccumulatorScheme::BootstrapAddPartial(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<DCRTPoly> & ct,
		const std::shared_ptr<FPRingGSWEvalKey> EK,
		const uint32_t from_bit,
		const uint32_t to_bit,
		//const uint32_t collect_from_bit,
		//const uint32_t collect_to_bit,
		const int32_t muls
   ) const {
		
		//if ((from_bit > collect_from_bit) || (to_bit < collect_to_bit)) {
		//	PALISADE_THROW(config_error, "Adding collectis are error");
		//}
		//!!!!!!!!!!!!!!! Should Turn On the Max ,,,
		auto prepare_start = std::chrono::system_clock::now();
	
		// Level 2 should handle
    //uint32_t EN_len = (uint32_t ) std::ceil( (double) 64.0 / (double) params -> GetBaseEncodeBit() );
		uint32_t EN_len = to_bit - from_bit + 1;

		//NativeInteger Q0 = params->GetModuliQ()[0];
		NativeInteger Q1 = params->GetModuliQ()[1];
		uint32_t K_FP = params->GetRLWEParams()->GetK();
		uint32_t N_FP = params->GetRLWEParams()->GetN();
	
		// Merge
    vector<vector<DCRTPoly>> CarryInfo(EN_len);  // Should Eval Mode
    vector<vector<DCRTPoly>> Out_Poly(EN_len);
    
		for (uint32_t k = 0; k < EN_len; k++) {
      CarryInfo[k].resize(K_FP+1);
			Out_Poly[k].resize(K_FP+1);
    }
	
		// Evaluation Mode
    vector<NativePoly> a_poly_const(K_FP); 
    for (uint32_t k = 0; k < K_FP; k++) {
      a_poly_const[k] = ct[k].GetElementAtIndex(0);
      a_poly_const[k].SetFormat(Format::EVALUATION);    // EVALS!!
    }
		
		NativePoly b_poly_const = ct[K_FP].GetElementAtIndex(0);
		b_poly_const.SetFormat(Format::COEFFICIENT);  // COEFF!!
    

    uint32_t Q1_WO_bits = params->GetRLWEParams()->GetBits()[1];  // <32 bit
    uint64_t Q1_Int = Q1.ConvertToInt();
    uint32_t carry_bit;                           
    uint32_t out_nums;  
    uint32_t start_bits; 
    uint64_t prefix_add_vals; 
    uint32_t In_Out_Index = 0;

    std::shared_ptr<vector<vector<DCRTPoly>>> BS_out_vec;
		vector<DCRTPoly> BS_out;
    std::shared_ptr<const FPRLWECiphertextImpl> Mul_CT;
    vector<NativePoly> A_Poly_tmp(K_FP);
    NativePoly B_Poly_tmp;
 	
		std::shared_ptr<vector<NativeVector>> A;
    vector<NativeVector> A_tmp(K_FP);
    NativeInteger b;
    NativeInteger b_tmp;
	
		DCRTPoly DCRT_tmp;
    NativePoly Native_tmp;
		uint32_t rot_idx;
		std::chrono::duration<double> prepare_end = std::chrono::system_clock::now() - prepare_start;
		if (DEBUG) std::cout << "Bootstrapping Add Prepare time is " <<prepare_end.count() << std::endl;
		for (uint32_t idx = 0; idx < EN_len; idx ++) {
			// Make In & Out Index 0 by multiplying a 
		
			uint32_t loc_idx = from_bit + idx;
			//std::cout << "first idx is " << loc_idx << std::endl;
		
			auto rot_start = std::chrono::system_clock::now();
	
	
			rot_idx = (loc_idx==0) ? 0 : (2*N_FP - loc_idx);
			// A Carry & Rotation Handle
			A = RotandCarryAddA( params,
					a_poly_const, 
					&CarryInfo[idx],
					params->GetRotMonomial_RLWE(rot_idx).GetElementAtIndex(0), 
					( loc_idx != 0),				    
					(idx!=0));
					// Make In & Out Index 0 by multiplying a 
		
			// Handle B
			b = b_poly_const.GetValues()[loc_idx];
      // Carry Merging
      if (idx != 0) {
				CarryInfo[idx][K_FP].SetFormat(Format::COEFFICIENT);
				DivideRoundingSelf(params, &CarryInfo[idx][K_FP], Format::COEFFICIENT);
				b = b.ModAddEq(CarryInfo[idx][K_FP].GetElementAtIndex(0).GetValues()[0] ,params->GetModuliQ()[0]);
      }
			std::chrono::duration<double> rot_end = std::chrono::system_clock::now() - rot_start;
			if (DEBUG) std::cout << "Bootstrapping Add Rotation time is " <<rot_end.count() << std::endl;
			
			//First Bootstrapping,3 bit + random 1 bit is bootstrapping
			start_bits  = Q1_WO_bits  + 4 + 1; // Programmed bit is need 
			carry_bit   = 1;
			out_nums    = 2;
			prefix_add_vals = (Q1_Int >> 1) + 8 * Q1_Int;  

			if (idx != EN_len - 1) {
				BS_out_vec = this->Bootstrap_and_KS(params
					, LWEscheme
					, params->GetACCAdd(0)
					, EK
					, (*A)
					, b, In_Out_Index, start_bits, prefix_add_vals, carry_bit, out_nums);
			} else {
				BS_out_vec = this->Bootstrap_and_KS(params
					, LWEscheme
					, params->GetACCAdd(1)
					, EK
					, (*A)
					, b, In_Out_Index, start_bits, prefix_add_vals, carry_bit, out_nums);
			} 
			Out_Poly[idx] = std::move((*BS_out_vec)[0]);

			// Carry save
			if (idx != EN_len - 1) {
				CarryInfo[idx+1] = std::move((*BS_out_vec)[1]);
			} else {
				CarryInfo[0] = std::move((*BS_out_vec)[1]);
			}

		}
		// Tmp Collection
	
		for (uint32_t k = 0; k < K_FP+1; k++) {
			for (uint32_t idx=1; idx < EN_len; idx++) {
				Out_Poly[0][k] += (Out_Poly[idx][k] * params->GetRotMonomial_RLWE(idx));
			}
		}
	
		// Products
		std::shared_ptr<vector<DCRTPoly>> mul_res; 
		if (params->GetPARALLEL()) {
			mul_res = this->ProductDCRTPARALLEL2(params, EK->Evkey, Out_Poly[0], Out_Poly[0][K_FP], CarryInfo[0], CarryInfo[0][K_FP]); 
		} else {
			mul_res = this->ProductDCRT(params, EK->Evkey, Out_Poly[0], Out_Poly[0][K_FP], CarryInfo[0], CarryInfo[0][K_FP]); 
		}
		// Reduction
	
		auto rd_start = std::chrono::system_clock::now();
	

		for (uint32_t k = 0; k < K_FP; k++) {
			(*mul_res)[k].SetFormat(Format::COEFFICIENT);
			DivideRoundingSelf(params, &((*mul_res)[k]), Format::EVALUATION); 
			a_poly_const[k] = std::move((*mul_res)[k].GetElementAtIndex(0));
		}
		(*mul_res)[K_FP].SetFormat(Format::COEFFICIENT);
		DivideRoundingSelf(params, &((*mul_res)[K_FP]), Format::COEFFICIENT); 
		b_poly_const = std::move((*mul_res)[K_FP].GetElementAtIndex(0));
	

		std::chrono::duration<double> rd_end = std::chrono::system_clock::now() - rd_start;
		if (DEBUG) std::cout << "Reduction times is " << rd_end.count() << std::endl;
		

		// Evaluation Mode
		
		//Conpensate Bits One More Times..
		for (uint32_t idx = 0; idx < EN_len; idx ++) {
	 
			uint32_t loc_idx = idx;
			//std::cout << "second idx is " << loc_idx << std::endl;
		
			rot_idx = (loc_idx==0) ? 0 : (2*N_FP - loc_idx);
			// A Carry & Rotation Handle
			A = RotandCarryAddA( params,
					a_poly_const, 
					&CarryInfo[idx],
					params->GetRotMonomial_RLWE(rot_idx).GetElementAtIndex(0), 
					(loc_idx != 0 ),				    
					(idx!=0));
					// Make In & Out Index 0 by multiplying a 
		
			b = b_poly_const.GetValues()[loc_idx];
      // Carry Merging
      if (idx != 0) {
				CarryInfo[idx][K_FP].SetFormat(Format::COEFFICIENT);
				DivideRoundingSelf(params, &CarryInfo[idx][K_FP], Format::COEFFICIENT);
				b = b.ModAddEq(CarryInfo[idx][K_FP].GetElementAtIndex(0).GetValues()[0] ,params->GetModuliQ()[0]);
      }

			// Orign, 3+1(padding but no infomation) bit bootstrapping msg <= 3bit 
			start_bits  = Q1_WO_bits + 3 + 1;
			carry_bit = 1;
			out_nums = 2;
			prefix_add_vals = (Q1_Int >> 1) + 4*Q1_Int;  
			// ================   EVAL MODE SAVING ==================== 
			BS_out_vec= this->Bootstrap_and_KS(params
				, LWEscheme
				, params->GetACCAdd(3),EK
				, (*A)
				, b
				, In_Out_Index, start_bits, prefix_add_vals, carry_bit, out_nums);
				//  ============================================================
			Out_Poly[idx] = std::move((*BS_out_vec)[0]);
			
			// Carry save
			if (idx < EN_len-1 ) {	
				CarryInfo[idx+1] = std::move((*BS_out_vec)[1]);

			}
		}

		// Collect Again
		// Shifts
		//uint32_t collect_nums = collect_to_bit - collect_from_bit + 1;
		
		uint32_t save_start_bits = (uint32_t)(to_bit + muls - EN_len + 1);
	
		for (uint32_t k = 0; k < K_FP+1; k++) {

			for (uint32_t idx = 0; idx < EN_len; idx++) {
				//uint32_t real_idx = ((N_FP*4) + collect_from_bit + idx) % (2*N_FP);
				uint32_t real_idx = ((N_FP*4) + save_start_bits + idx) % (2*N_FP);
				//std::cout << "real_idx is " << real_idx << std::endl;		
				if ( idx == 0) {
					Out_Poly[0][k] *= params->GetRotMonomial_RLWE(real_idx);
				} else {
					Out_Poly[0][k] += (Out_Poly[idx][k] * params->GetRotMonomial_RLWE(real_idx));
				}
			}
		}
		vector<vector<DCRTPoly>> res(2);
		res[0] = std::move(Out_Poly[0]);
		res[1] = std::move(CarryInfo[0]);
		return std::make_shared<vector<vector<DCRTPoly>>>(res);
	}


	// From, to bits
	std::shared_ptr<vector<vector<DCRTPoly>>> FPRingGSWAccumulatorScheme::BootstrapAddPartialDouble(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<DCRTPoly> & ct,
		const std::shared_ptr<FPRingGSWEvalKey> EK,
		const uint32_t from_bit,
		const uint32_t to_bit,
		//const uint32_t collect_from_bit,
		//const uint32_t collect_to_bit,
		const int32_t muls
   ) const {
		
		//if ((from_bit > collect_from_bit) || (to_bit < collect_to_bit)) {
		//	PALISADE_THROW(config_error, "Adding collectis are error");
		//}
		//!!!!!!!!!!!!!!! Should Turn On the Max ,,,
		auto prepare_start = std::chrono::system_clock::now();
		// 4 ~ 32 after -1 than 3 ~ 31
		uint32_t EN_len = to_bit - from_bit + 1;

		//NativeInteger Q0 = params->GetModuliQ()[0];
		NativeInteger Q1 = params->GetModuliQ()[1];
		uint32_t K_FP = params->GetRLWEParams()->GetK();
		uint32_t N_FP = params->GetRLWEParams()->GetN();
	
		// Merge
    vector<vector<DCRTPoly>> CarryInfo(EN_len);  // Should Eval Mode
    vector<vector<DCRTPoly>> Out_Poly(EN_len);
    vector<vector<DCRTPoly>> Out_Tot(EN_len+2);
    
		for (uint32_t k = 0; k < EN_len; k++) {
      CarryInfo[k].resize(K_FP+1);
			Out_Poly[k].resize(K_FP+1);
    }
		
		for (uint32_t k = 0; k < EN_len+2; k++) {
			Out_Tot[k].resize(K_FP+1);
		}
		

		// Evaluation Mode
    vector<NativePoly> a_poly_const(K_FP); 
    for (uint32_t k = 0; k < K_FP; k++) {
      a_poly_const[k] = ct[k].GetElementAtIndex(0);
      a_poly_const[k].SetFormat(Format::EVALUATION);    // EVALS!!
    }
		
		NativePoly b_poly_const = ct[K_FP].GetElementAtIndex(0);
		b_poly_const.SetFormat(Format::COEFFICIENT);  // COEFF!!
    

    uint32_t Q1_WO_bits = params->GetRLWEParams()->GetBits()[1];  // <32 bit
    uint64_t Q1_Int = Q1.ConvertToInt();
    uint32_t carry_bit;                           
    uint32_t out_nums;  
    uint32_t start_bits; 
    uint64_t prefix_add_vals; 
    uint32_t In_Out_Index = 0;

    std::shared_ptr<vector<vector<DCRTPoly>>> BS_out_vec;
		vector<DCRTPoly> BS_out;
    std::shared_ptr<const FPRLWECiphertextImpl> Mul_CT;
    vector<NativePoly> A_Poly_tmp(K_FP);
    NativePoly B_Poly_tmp;
 	
		std::shared_ptr<vector<NativeVector>> A;
    vector<NativeVector> A_tmp(K_FP);
    NativeInteger b;
    NativeInteger b_tmp;
	
		DCRTPoly DCRT_tmp;
    NativePoly Native_tmp;
		uint32_t rot_idx;
		std::chrono::duration<double> prepare_end = std::chrono::system_clock::now() - prepare_start;
		if (DEBUG) std::cout << "Bootstrapping Add Prepare time is " <<prepare_end.count() << std::endl;
		for (uint32_t idx = 0; idx < EN_len; idx ++) {
			// Make In & Out Index 0 by multiplying a 
		
			uint32_t loc_idx = from_bit + idx;
			//std::cout << "first idx is " << loc_idx << std::endl;
		
			auto rot_start = std::chrono::system_clock::now();
	
	
			rot_idx = (loc_idx==0) ? 0 : (2*N_FP - loc_idx);
			// A Carry & Rotation Handle
			A = RotandCarryAddA( params,
					a_poly_const, 
					&CarryInfo[idx],
					params->GetRotMonomial_RLWE(rot_idx).GetElementAtIndex(0), 
					( loc_idx != 0),	// rot or not		    
					(idx!=0));       // Carry adding or not
					// Make In & Out Index 0 by multiplying a 
		
			// Handle B
			b = b_poly_const.GetValues()[loc_idx];
      // Carry Merging
      if (idx != 0) {
				CarryInfo[idx][K_FP].SetFormat(Format::COEFFICIENT);
				DivideRoundingSelf(params, &CarryInfo[idx][K_FP], Format::COEFFICIENT);
				b = b.ModAddEq(CarryInfo[idx][K_FP].GetElementAtIndex(0).GetValues()[0] ,params->GetModuliQ()[0]);
      }
			std::chrono::duration<double> rot_end = std::chrono::system_clock::now() - rot_start;
			if (DEBUG) std::cout << "Bootstrapping Add Rotation time is " <<rot_end.count() << std::endl;
			
			//First Bootstrapping,3 bit + random 1 bit is bootstrapping
			start_bits  = Q1_WO_bits  + 4 + 1; // Programmed bit is need 
			carry_bit   = 1;
			out_nums    = 2;
			prefix_add_vals = (Q1_Int >> 1) + 8 * Q1_Int;  

			if (idx != EN_len - 1) {
				BS_out_vec = this->Bootstrap_and_KS(params
					, LWEscheme
					, params->GetACCAdd(0)
					, EK
					, (*A)
					, b, In_Out_Index, start_bits, prefix_add_vals, carry_bit, out_nums);
			} else {
				BS_out_vec = this->Bootstrap_and_KS(params
					, LWEscheme
					, params->GetACCAdd(1)
					, EK
					, (*A)
					, b, In_Out_Index, start_bits, prefix_add_vals, carry_bit, out_nums);
			} 
			Out_Poly[idx] = std::move((*BS_out_vec)[0]);

			// Carry save
			if (idx != EN_len - 1) {
				CarryInfo[idx+1] = std::move((*BS_out_vec)[1]);
			} else {
				CarryInfo[0] = std::move((*BS_out_vec)[1]);
			}

		}
		// Tmp Collection
	
		for (uint32_t k = 0; k < K_FP+1; k++) {
			for (uint32_t idx=1; idx < EN_len; idx++) {
				Out_Poly[0][k] += (Out_Poly[idx][k] * params->GetRotMonomial_RLWE(idx));
			}
		}
	
		// Products
		std::shared_ptr<vector<DCRTPoly>> mul_res; 
		if (params->GetPARALLEL()) {
			mul_res = this->ProductDCRTPARALLEL2(params, EK->Evkey, Out_Poly[0], Out_Poly[0][K_FP], CarryInfo[0], CarryInfo[0][K_FP]); 
		} else {
			mul_res = this->ProductDCRT(params, EK->Evkey, Out_Poly[0], Out_Poly[0][K_FP], CarryInfo[0], CarryInfo[0][K_FP]); 
		}
		// Reduction
	
		auto rd_start = std::chrono::system_clock::now();
	
		for (uint32_t k = 0; k < K_FP; k++) {
			(*mul_res)[k].SetFormat(Format::COEFFICIENT);
			DivideRoundingSelf(params, &((*mul_res)[k]), Format::EVALUATION); 
			a_poly_const[k] = std::move((*mul_res)[k].GetElementAtIndex(0));
		}
		(*mul_res)[K_FP].SetFormat(Format::COEFFICIENT);
		DivideRoundingSelf(params, &((*mul_res)[K_FP]), Format::COEFFICIENT); 
		b_poly_const = std::move((*mul_res)[K_FP].GetElementAtIndex(0));
	

		std::chrono::duration<double> rd_end = std::chrono::system_clock::now() - rd_start;
		if (DEBUG) std::cout << "Reduction times is " << rd_end.count() << std::endl;
		

		// Evaluation Mode
		//std::cout << "compensate start" << std::endl;
		//Conpensate Bits One More Times..
		
		for (uint32_t idx = 0; idx < EN_len; idx ++) {
			uint32_t loc_idx = idx;
			//std::cout << "second idx is " << loc_idx << std::endl;
		
			rot_idx = (loc_idx==0) ? 0 : (2*N_FP - loc_idx);
			// A Carry & Rotation Handle
			A = RotandCarryAddA( params,
					a_poly_const, 
					&CarryInfo[idx],
					params->GetRotMonomial_RLWE(rot_idx).GetElementAtIndex(0), 
					(loc_idx != 0 ),				    
					(idx != 0));
					// Make In & Out Index 0 by multiplying a 
		
			b = b_poly_const.GetValues()[loc_idx];
      // Carry Merging
      if (idx != 0) {
				CarryInfo[idx][K_FP].SetFormat(Format::COEFFICIENT);
				DivideRoundingSelf(params, &CarryInfo[idx][K_FP], Format::COEFFICIENT);
				b = b.ModAddEq(CarryInfo[idx][K_FP].GetElementAtIndex(0).GetValues()[0] ,params->GetModuliQ()[0]);
      }

			// Orign, 3+1(padding but no infomation) bit bootstrapping msg <= 3bit 
			start_bits  = Q1_WO_bits + 3 + 1;
			carry_bit = 2;
			out_nums = 3;
			prefix_add_vals = (Q1_Int >> 1) + 4*Q1_Int;  
			// ================   EVAL MODE SAVING ==================== 
			BS_out_vec= this->Bootstrap_and_KS(params
				, LWEscheme
				, params->GetACCAdd(4),EK
				, (*A)
				, b
				, In_Out_Index, start_bits, prefix_add_vals, carry_bit, out_nums);
				//  ============================================================
			Out_Poly[idx] = std::move((*BS_out_vec)[0]);
			
			// Carry save
			Out_Tot[idx] =  std::move((*BS_out_vec)[2]);
			if (idx < EN_len-1 ) {	
				CarryInfo[idx+1] = std::move((*BS_out_vec)[1]);
			}
		}
		
		uint32_t save_start_bits = (uint32_t)(to_bit + muls - EN_len + 1);
		//std::cout <<"save start_bits is " << save_start_bits << std::endl;
		for (uint32_t k = 0; k < K_FP+1; k++) {

			for (uint32_t idx = 0; idx < EN_len; idx++) {
				uint32_t real_idx = ((N_FP*4) + save_start_bits + idx) % (2*N_FP);
				//std::cout << "real_idx is " << real_idx << std::endl;		
				if ( idx == 0) {
					Out_Poly[0][k] *= params->GetRotMonomial_RLWE(real_idx);
				} else {
					Out_Poly[0][k] += (Out_Poly[idx][k] * params->GetRotMonomial_RLWE(real_idx));
				}
			}
		}
		
		Out_Tot[EN_len] = std::move(Out_Poly[0]);
		Out_Tot[EN_len+1] = std::move(CarryInfo[0]);
		return std::make_shared<vector<vector<DCRTPoly>>>(Out_Tot);
	}








};  // namespace lbcrypto
