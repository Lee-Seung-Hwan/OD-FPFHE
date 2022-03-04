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

#include "fprlwe.h"
//#include <ctime>
namespace lbcrypto {

	// Encryption as described in Section 5 of https://eprint.iacr.org/2014/816
	// skNTT corresponds to the secret key z


	// Poly Key Gen
	std::shared_ptr<FPRLWEPrivateKeyImpl> FPRLWEEncryptionScheme::KeyGen(
		
		const std::shared_ptr<FPRLWECryptoParams> params,
		uint32_t h_sparse) const {
		
		DCRTPoly::TugType tug;
		uint32_t K = params->GetK();
		vector<DCRTPoly> SK(K);
		// Sparse not we don't care

		for (uint32_t i = 0; i < K; i++) {
			SK[i] = DCRTPoly(tug, params->GetDCRTParams(), Format::EVALUATION);
		}
			return std::make_shared<FPRLWEPrivateKeyImpl> (SK);

	}

	//std::shared_ptr<FPRLWEPrivateKeyImpl> FPRLWEEncryptionScheme::KeyGenN(
    //const std::shared_ptr<FPLWECryptoParams> params) const {

	//		TernaryUniformGeneratorImpl<NativeVector> tug;
//			return std::make_shared<FPRLWEPrivateKeyImpl>( FPLWEPrivateKeyImpl(tug.GenerateVector(params->GetN(), params->GetQ())));
//		}

	//std::shared_ptr<FPRLWEPrivateKeyImpl> FPRLWEEncryptionScheme::KeyGenPackingN(
  //  const std::shared_ptr<FPLWECryptoParams> params) const {

	//		TernaryUniformGeneratorImpl<NativeVector> tug;
	//		return std::make_shared<FPLWEPrivateKeyImpl>(FPLWEPrivateKeyImpl(tug.GenerateVector(params->GetpackingN(), params->GetQ())));
	//	}


// Divide and Modular reduction


	/* Current not Used!!  */
	
/* SH
	const std::shared_ptr<const FPRLWECiphertextImpl> FPRLWEEncryptionScheme::DivideReduction(
		const std::shared_ptr<FPRLWECryptoParams> params,
		const std::shared_ptr<const FPRLWECiphertextImpl> ct) const {

		vector<uint32_t> Qbits = params-> GetBits();
		vector<NativeInteger> QInv = params->GetInv();
		vector<NativeInteger> moduliQ = params->GetModuliQ();
		uint32_t shifter = Qbits[0] - Qbits[1];
		
		DCRTPoly ct_a = ct->GetA();
		DCRTPoly ct_b = ct->GetB();

		//ct_a 
		vector<NativePoly> a_poly = ct_a.GetAllElements();
		
		a_poly[0].SetFormat(Format::COEFFICIENT);
		a_poly[1].SetFormat(Format::COEFFICIENT);
		
		NativeVector a_vector_0 = a_poly[0].GetValues();
		NativeVector a_vector_1 = a_poly[1].GetValues();
	
		for (uint32_t i = 0; i < a_vector_0.GetLength(); i++) {
			//a_vector_1[i] = (a_vector_1[i] >> shifter);
			a_vector_0[i] = a_vector_0[i].ModMul(QInv[0],moduliQ[0]);
			a_vector_1[i] = a_vector_1[i].ModMul(QInv[1],moduliQ[1]);
			a_vector_0[i] = a_vector_0[i].ModAdd((a_vector_1[i] << shifter), moduliQ[0]);
		}
		
		a_poly[0].SetValues(a_vector_0,Format::COEFFICIENT);
		a_poly[0].SetFormat(Format::EVALUATION);
		//a_avector_0 += a_poly[1];
		
		ct_a.SetElementAtIndex(0,a_poly[0]);
		ct_a.DropLastElement();
		
		ct_a
		
		vector<NativePoly> b_poly = ct_b.GetAllElements();
			
		b_poly[0].SetFormat(Format::COEFFICIENT);
		b_poly[1].SetFormat(Format::COEFFICIENT);
		
		//b_poly[0] *= QInv[0];
		//b_poly[1] *= QInv[1];
		
		NativeVector b_vector_0 = b_poly[0].GetValues();
		NativeVector b_vector_1 = b_poly[1].GetValues();
	
		for (uint32_t i = 0; i < b_vector_0.GetLength(); i++) {
			//a_vector_1[i] = (a_vector_1[i] >> shifter);
			b_vector_0[i] = b_vector_0[i].ModMul(QInv[0],moduliQ[0]);
			b_vector_1[i] = b_vector_1[i].ModMul(QInv[1],moduliQ[1]);
			b_vector_0[i] = b_vector_0[i].ModAdd((b_vector_1[i] >> shifter), moduliQ[0]);
	
			//b_vector_0[i].ModAdd((b_vector_1[i] >> shifter), moduliQ[0]);
			
		}
		
		b_poly[0].SetValues(b_vector_0,Format::COEFFICIENT);
		b_poly[0].SetFormat(Format::EVALUATION);
		//a_avector_0 += a_poly[1];
		
		ct_b.SetElementAtIndex(0,b_poly[0]);
		ct_b.DropLastElement();
		
	return std::make_shared<FPRLWECiphertextImpl> (FPRLWECiphertextImpl(ct_a,ct_b,1));
	
		}
		*/
		


// classical RLWE encryption
// a is a randomly uniform vector of dimension n; with integers mod q
// b = a*s + e + m floor(q/4) is an integer mod q

	const std::shared_ptr<const FPRLWECiphertextImpl> FPRLWEEncryptionScheme::Encrypt(
    const std::shared_ptr<FPRLWECryptoParams> params,
    const std::shared_ptr<const FPRLWEPrivateKeyImpl> sk,
    DCRTPoly &m) const {
		
		// m's Format is COEFFICIENT
		DCRTPoly::DugType dug;
		vector<DCRTPoly> sk_poly = sk->GetElement();	
		// sk shouuld be Evaluation
		
		vector<DCRTPoly> a_vec(sk_poly.size());
		DCRTPoly b =  DCRTPoly(params->GetDgg(), params -> GetDCRTParams(), Format::COEFFICIENT);
		b += m ;
		b.SetFormat(Format::EVALUATION);

		for (uint32_t i = 0; i< sk_poly.size(); i++) {
			a_vec[i] = DCRTPoly(dug, params -> GetDCRTParams(), Format::EVALUATION);
			b += (a_vec[i] * sk_poly[i]);
		}
	
		return std::make_shared<FPRLWECiphertextImpl>(FPRLWECiphertextImpl(a_vec, b));
	
	}



	// classical RLWE decryption
	// m_result = Round(4/q * (b - a*s))
	void FPRLWEEncryptionScheme::Decrypt(
    const std::shared_ptr<FPRLWECryptoParams> params,
    const std::shared_ptr<const FPRLWEPrivateKeyImpl> sk,
    const std::shared_ptr<const FPRLWECiphertextImpl> ct,
    DCRTPoly *result) const {
		// TODO in the future we should add a check to make sure sk parameters match
		// the ct parameters
		//CT should be Evaluation Mode


		// Create local variables to speed up the computations
		vector<DCRTPoly> a = ct->GetA();
		//a.SetFormat(Format::EVALUATION);

		vector<DCRTPoly> s = sk->GetElement();
		//std::cout << "Level is " << ct->GetLevel() << std::endl;
		//std::cout << a[0] << std::endl;
		if (ct->GetLevel() == 1){
			for (uint32_t i = 0; i < a.size(); i++ ){
				s[i].DropLastElement();
				
			}
			std::cout << "Drop Last Elements" << std::endl;
		}
		a[0] *= s[0];
		for (uint32_t i = 1; i < s.size(); i++) {
			a[0] += (a[i] * s[i]);
		}
	
		DCRTPoly b = ct->GetB();
		b -=a[0];
		b.SetFormat(Format::COEFFICIENT);
		*result = b;
	
		// Debug
		
		//std::cout << " ========================= Decrypt Results ========================== " << std::endl;
		//double tmp_res = 0;
		//for (uint32_t i = 0; i < 16; i++) {
		//	tmp_res = (result[0][i]. ConvertToDouble()) / ((double) params->GetScaling());
		//	std::cout << i << "st result is " << tmp_res << std::endl;
		//}
	//std::cout << " ========================= Decrypt Results END ========================== " << std::endl;

		#if defined(FPFHE_DEBUG)
		//double 

		//error = (4.0 * (r.ConvertToDouble() - q.ConvertToInt() / 8)) /
    //                 q.ConvertToDouble() -
    //             static_cast<double>(*result);
		//std::cerr << "error:\t" << error << std::endl;
		#endif
		

		return;
	}
	

	// Switching key as described in Section 3 of https://eprint.iacr.org/2014/816
	
	std::shared_ptr<FPRingEvaluationKey> FPRLWEEncryptionScheme::EvalKeyGen(
    const std::shared_ptr<FPRLWECryptoParams> params,
		uint32_t EV_Bit,
		uint32_t EV_len_save,
		uint32_t EV_len_rm,
		DiscreteGaussianGeneratorImpl<NativeVector> EVDgg,
    const std::shared_ptr<const FPRLWEPrivateKeyImpl> sk) const {

		// 임시 방편
		DCRTPoly::DugType dug;
		vector<DCRTPoly> sk_poly_vec = sk->GetElement();
		DCRTPoly tmp;


		//baseG shoud be bits
		uint32_t K_NP = params->GetK();
		uint32_t len_tot = (K_NP*(K_NP+1) >> 1);
		FPRingEvaluationKey EvalKey(len_tot, EV_len_save, K_NP + 1);
	
		for (uint32_t k_row = 0; k_row < K_NP; k_row ++){
			for (uint32_t k_col = 0; k_col <= k_row; k_col++){
				uint32_t k_idx = ((k_row*(k_row+1)) >> 1) + k_col;
			
				tmp = sk_poly_vec[k_row] * sk_poly_vec[k_col];
				
				// Eval Powers Of Base
				vector<DCRTPoly> mul_2_sk_poly_square = tmp.PowersOfBase(EV_Bit);
	
				for (uint32_t i = 0; i < EV_len_save; i++) {     
					// K=0 Encoding
					EvalKey[k_idx][i][0] = DCRTPoly(dug, params -> GetDCRTParams(), Format::EVALUATION);		
					EvalKey[k_idx][i][K_NP] = EvalKey[k_idx][i][0] * sk_poly_vec[0];
					// Rmain K=1~K_NP setting
					for (uint32_t k_enc = 1; k_enc < K_NP; k_enc++) {
						EvalKey[k_idx][i][k_enc] = DCRTPoly(dug, params -> GetDCRTParams(), Format::EVALUATION);		
						EvalKey[k_idx][i][K_NP] += EvalKey[k_idx][i][k_enc] * sk_poly_vec[k_enc];
					}
					// B setting
					EvalKey[k_idx][i][K_NP] += DCRTPoly(EVDgg, params -> GetDCRTParams(), Format::EVALUATION);
					EvalKey[k_idx][i][K_NP] += mul_2_sk_poly_square[i+EV_len_rm];
				}
			}
		}

		return std::make_shared<FPRingEvaluationKey> (EvalKey);
	}



/* SH
	 const std::shared_ptr<const FPRLWECiphertextImpl> FPRLWEEncryptionScheme::Product(
		const std::shared_ptr<const FPRLWECiphertextImpl> ct1,
		const std::shared_ptr<const FPRLWECiphertextImpl> ct2,
		const std::shared_ptr<FPRLWECryptoParams> r_params,
		uint32_t EV_Bit,
		uint32_t EV_len_rm,
	
		const std::shared_ptr<FPRingEvaluationKey> EvalKey) const {
		
			DCRTPoly ct1_a = ct1->GetA();
			DCRTPoly ct1_b = ct1->GetB();
			
			DCRTPoly ct2_a = ct2->GetA();
			DCRTPoly ct2_b = ct2->GetB();
			
			// Eval transforl

			// a1*a2
			DCRTPoly a_a = ct1_a.Clone();
			a_a *= ct2_a;
			// a1 * b2
			DCRTPoly a_b = ct1_a.Clone();
			a_b *= ct2_b;
			// b1 * a2
			DCRTPoly b_a = ct1_b.Clone();
			b_a *= ct2_a;
			// b1 * b2
			DCRTPoly b_b = ct1_b.Clone();
			b_b *= ct2_b;	

		
			//Relinearlization
			a_a.SetFormat(Format::COEFFICIENT);

			vector<DCRTPoly> a_a_rl = a_a.BaseDecompose(EV_Bit,true);
	
			//a_a.SetFormat(Format::EVALUATION);	
			std::vector<std::vector<DCRTPoly>> m_e_key = EvalKey->GetElements();
		
			Zero Initialization  
			DCRTPoly a_tmp = DCRTPoly(r_params->GetDCRTParams(), Format::EVALUATION, true);
			DCRTPoly b_tmp = DCRTPoly(r_params->GetDCRTParams(), Format::EVALUATION, true);

			DCRTPoly eval_a_tmp;
			DCRTPoly eval_b_tmp;

			for (uint32_t i=0; i < m_e_key.size(); i++) {
			
				eval_a_tmp =(a_a_rl[i+EV_len_rm] * m_e_key[i][0]);
				eval_b_tmp =(a_a_rl[i+EV_len_rm] * m_e_key[i][1]);
			
				a_tmp += eval_a_tmp;
				b_tmp += eval_b_tmp;
			}
				
			b_a += a_b;
			b_a += a_tmp;
			b_b += b_tmp;
			
			// One More Time!!!!
			DCRTPoly a_a_sign;
			DCRTPoly a_b_sign;
			DCRTPoly b_a_sign;
			DCRTPoly b_b_sign;
	
			if (ct1->GetType() == INT64 && ct2->GetType() == INT64) {	
				DCRTPoly ct1_a_sign = ct1->GetASign();
				DCRTPoly ct1_b_sign = ct1->GetBSign();
			
				DCRTPoly ct2_a_sign = ct2->GetASign();
				DCRTPoly ct2_b_sign = ct2->GetBSign();
			
				// a1*a2
				a_a_sign = ct1_a_sign.Clone();
				a_a_sign *= ct2_a_sign;
				// a1 * b2
				a_b_sign = ct1_a_sign.Clone();
				a_b_sign *= ct2_b_sign;
				// b1 * a2
				b_a_sign = ct1_b_sign.Clone();
				b_a_sign *= ct2_a_sign;
				// b1 * b2
				b_b_sign = ct1_b_sign.Clone();
				b_b_sign *= ct2_b_sign;	

		
				Relinearlization
				a_a_sign.SetFormat(Format::COEFFICIENT);

				vector<DCRTPoly> a_a_rl_sign = a_a_sign.BaseDecompose(EV_Bit,true);
	
				//a_a.SetFormat(Format::EVALUATION);	
				//std::vector<std::vector<DCRTPoly>> m_e_key = EvalKey->GetElements();
		
				Zero Initialization  
				DCRTPoly a_tmp_sign = DCRTPoly(r_params->GetDCRTParams(), Format::EVALUATION, true);
				DCRTPoly b_tmp_sign = DCRTPoly(r_params->GetDCRTParams(), Format::EVALUATION, true);

				DCRTPoly eval_a_tmp_sign;
				DCRTPoly eval_b_tmp_sign;

				for (uint32_t i=0; i < m_e_key.size(); i++) {
			
					eval_a_tmp_sign =(a_a_rl_sign[i+EV_len_rm] * m_e_key[i][0]);
					eval_b_tmp_sign =(a_a_rl_sign[i+EV_len_rm] * m_e_key[i][1]);
			
					a_tmp_sign += eval_a_tmp_sign;
					b_tmp_sign += eval_b_tmp_sign;
				}
				
				b_a_sign += a_b_sign;
				b_a_sign += a_tmp_sign;
				b_b_sign += b_tmp_sign;
		
			}
			if (ct1->GetType() == UINT64) { 
				FPRLWECiphertextImpl outputs(b_a, b_b);
				return std::make_shared<const FPRLWECiphertextImpl> (outputs);

			} else { 
				FPRLWECiphertextImpl outputs(b_a, b_b, b_a_sign, b_b_sign);
				return std::make_shared<const FPRLWECiphertextImpl> (outputs);
			}		
		}

	
*/

/* SH
	 const std::shared_ptr<const FPRLWECiphertextImpl> FPRLWEEncryptionScheme::Product(
		const std::shared_ptr<const FPRLWECiphertextImpl> ct1,
		const std::shared_ptr<const FPRLWECiphertextImpl> ct2,
		const std::shared_ptr<FPRLWECryptoParams> r_params,
		uint32_t EV_Bit,
		uint32_t EV_len_rm,
	
		const std::shared_ptr<FPRingEvaluationKey> EvalKey) const {
		
			DCRTPoly ct1_a = ct1->GetA();
			DCRTPoly ct1_b = ct1->GetB();
			
			DCRTPoly ct2_a = ct2->GetA();
			DCRTPoly ct2_b = ct2->GetB();
			
			// Eval transforl

			// a1*a2
			DCRTPoly a_a = ct1_a.Clone();
			a_a *= ct2_a;
			// a1 * b2
			DCRTPoly a_b = ct1_a.Clone();
			a_b *= ct2_b;
			// b1 * a2
			DCRTPoly b_a = ct1_b.Clone();
			b_a *= ct2_a;
			// b1 * b2
			DCRTPoly b_b = ct1_b.Clone();
			b_b *= ct2_b;	

		
			//Relinearlization
			a_a.SetFormat(Format::COEFFICIENT);

			vector<DCRTPoly> a_a_rl = a_a.BaseDecompose(EV_Bit,true);
	
			//a_a.SetFormat(Format::EVALUATION);	
			std::vector<std::vector<DCRTPoly>> m_e_key = EvalKey->GetElements();
		
			Zero Initialization  
			DCRTPoly a_tmp = DCRTPoly(r_params->GetDCRTParams(), Format::EVALUATION, true);
			DCRTPoly b_tmp = DCRTPoly(r_params->GetDCRTParams(), Format::EVALUATION, true);

			DCRTPoly eval_a_tmp;
			DCRTPoly eval_b_tmp;

			for (uint32_t i=0; i < m_e_key.size(); i++) {
			
				eval_a_tmp =(a_a_rl[i+EV_len_rm] * m_e_key[i][0]);
				eval_b_tmp =(a_a_rl[i+EV_len_rm] * m_e_key[i][1]);
			
				a_tmp += eval_a_tmp;
				b_tmp += eval_b_tmp;
			}
				
			b_a += a_b;
			b_a += a_tmp;
			b_b += b_tmp;
			
			// One More Time!!!!
			DCRTPoly a_a_sign;
			DCRTPoly a_b_sign;
			DCRTPoly b_a_sign;
			DCRTPoly b_b_sign;
	
			if (ct1->GetType() == INT64 && ct2->GetType() == INT64) {	
				DCRTPoly ct1_a_sign = ct1->GetASign();
				DCRTPoly ct1_b_sign = ct1->GetBSign();
			
				DCRTPoly ct2_a_sign = ct2->GetASign();
				DCRTPoly ct2_b_sign = ct2->GetBSign();
			
				// a1*a2
				a_a_sign = ct1_a_sign.Clone();
				a_a_sign *= ct2_a_sign;
				// a1 * b2
				a_b_sign = ct1_a_sign.Clone();
				a_b_sign *= ct2_b_sign;
				// b1 * a2
				b_a_sign = ct1_b_sign.Clone();
				b_a_sign *= ct2_a_sign;
				// b1 * b2
				b_b_sign = ct1_b_sign.Clone();
				b_b_sign *= ct2_b_sign;	

		
				Relinearlization
				a_a_sign.SetFormat(Format::COEFFICIENT);

				vector<DCRTPoly> a_a_rl_sign = a_a_sign.BaseDecompose(EV_Bit,true);
	
				//a_a.SetFormat(Format::EVALUATION);	
				//std::vector<std::vector<DCRTPoly>> m_e_key = EvalKey->GetElements();
		
				Zero Initialization  
				DCRTPoly a_tmp_sign = DCRTPoly(r_params->GetDCRTParams(), Format::EVALUATION, true);
				DCRTPoly b_tmp_sign = DCRTPoly(r_params->GetDCRTParams(), Format::EVALUATION, true);

				DCRTPoly eval_a_tmp_sign;
				DCRTPoly eval_b_tmp_sign;

				for (uint32_t i=0; i < m_e_key.size(); i++) {
			
					eval_a_tmp_sign =(a_a_rl_sign[i+EV_len_rm] * m_e_key[i][0]);
					eval_b_tmp_sign =(a_a_rl_sign[i+EV_len_rm] * m_e_key[i][1]);
			
					a_tmp_sign += eval_a_tmp_sign;
					b_tmp_sign += eval_b_tmp_sign;
				}
				
				b_a_sign += a_b_sign;
				b_a_sign += a_tmp_sign;
				b_b_sign += b_tmp_sign;
		
			}
			if (ct1->GetType() == UINT64) { 
				FPRLWECiphertextImpl outputs(b_a, b_b);
				return std::make_shared<const FPRLWECiphertextImpl> (outputs);

			} else { 
				FPRLWECiphertextImpl outputs(b_a, b_b, b_a_sign, b_b_sign);
				return std::make_shared<const FPRLWECiphertextImpl> (outputs);
			}		
		}

		std::shared_ptr<const FPLWECiphertextImpl> FPRLWEEncryptionScheme::SampleExtract(	
			const std::shared_ptr<FPRLWECryptoParams> params,
			const std::shared_ptr<const FPRLWECiphertextImpl> ct,
			uint32_t index_num
			) const {
	
			
//			uint32_t N = params->GetN();
//			NativeInteger Q = params->GetQ();
 // 
		//	NativeVector tmp_a(N);
		//	NativeInteger tmp_b = (ct->GetB())[index_num];
		//	NativePoly a = ct->GetA();
		//	
		//	uint32_t cnt = 0;
		//	for (int32_t i = index_num; i >= 0; i--) {
		//		tmp_a[cnt] =  a[i];
		//		cnt += 1;
		//	}
			
		//for (int32_t i = N-1; i > index_num; i--) {
		//		tmp_a[cnt] = Q - a[i];
			//	cnt += 1;
			//}
			//return std::make_shared<FPLWECiphertextImpl>(FPLWECiphertextImpl(a, b));
			
			//return std::make_shared<const FPLWECiphertextImpl> (FPLWECiphertextImpl(tmp_a,tmp_b));
			
			NativeVector a;
			NativeInteger b;
			return std::make_shared<const FPLWECiphertextImpl> (FPLWECiphertextImpl(a,b));
	}	
	*/


};  // namespace lbcrypto
