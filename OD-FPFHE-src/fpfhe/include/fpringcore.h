/*************
 * This file is modified from ringcore.h which is writen by Leo Ducas and Daniele Micciancio.
 * Full paper is listed in : eprint.iacr.org/2022/186.
 * This source is obeyed and applied following copyright law. 
 *
 * If you have any question, please contact us by the email kr3951@hanyang.ac.kr
 * Seunghwan Lee, 2022.03.04
 * *************/


// @file ringcore.h - Main ring classes for Boolean circuit FHE.
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

#ifndef FPFHE_RINGCORE_H
#define FPFHE_RINGCORE_H

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "fprlwe.h"


#include "lattice/backend.h"
#include "fplwecore.h"

#include "math/backend.h"
#include "math/discretegaussiangenerator.h"
#include "math/nbtheory.h"
#include "math/transfrm.h"
#include "utils/serializable.h"

namespace lbcrypto {

	// enum for all supported binary gates
	enum FPGATE { OR, AND, NOR, NAND, XOR_FAST, XNOR_FAST, XOR, XNOR };

	// Two variants of FHEW are supported based on the bootstrapping technique used:
	// AP and GINX Please see "Bootstrapping in FHEW-like Cryptosystems" for details
	// on both bootstrapping techniques
	enum FPFHEMETHOD { AP, GINX };

	/**
	* @brief Class that stores all parameters for the RingGSW scheme used in
	* bootstrapping
	*/
	class FPRingGSWCryptoParams : public Serializable {
		public:
		FPRingGSWCryptoParams()
      :  m_GD_base_bit(0) {}

  /**
   * Main constructor for RingGSWCryptoParams
   *
   * @param lweparams a shared poiter to an instance of LWECryptoParams
   * @param baseG the gadget base used in the bootstrapping
   * @param baseR the base for the refreshing key
   * @param method bootstrapping method (AP or GINX)
   */
  explicit FPRingGSWCryptoParams(
	
			const std::shared_ptr<FPLWECryptoParams>  lweparams,
  		const std::shared_ptr<FPRLWECryptoParams> rlweparams,
      const std::shared_ptr<ILDCRTParams<BigInteger>> DCRTGSWParamset,
			uint32_t N_GSW,
			uint32_t K_GSW,
			uint32_t n_PKS,
			double KS_std, double BK_std, double EV_std, double PKS_std,
			uint32_t KS_base_bit, uint32_t l_rm_KS,
			uint32_t EV_base_bit, uint32_t l_rm_EV,
			uint32_t GD_base_bit, uint32_t l_rm_GD,
			uint32_t PKS_base_bit, uint32_t l_rm_PKS,
			
			uint32_t Encode_base_bit, uint32_t C_base_bit,
			uint32_t h_FP, uint32_t h_GSW, uint32_t h_PKS,
			uint32_t exp_bit, bool PARALLEL_FLAG
			)

      : 
				m_LWEParams(lweparams),
				m_RLWEParams(rlweparams),
				m_DCRTGSWParams(DCRTGSWParamset),
				m_N_GSW(N_GSW),
				m_K_GSW(K_GSW),
				m_n_PKS(n_PKS),
				m_KS_base_bit(KS_base_bit),
				m_KS_l_rm(l_rm_KS),
				m_EV_base_bit(EV_base_bit),
				m_EV_l_rm(l_rm_EV),
				m_GD_base_bit(GD_base_bit),
				m_GD_l_rm(l_rm_GD),
				m_PKS_base_bit(PKS_base_bit),
				m_PKS_l_rm(l_rm_PKS),
				m_encode_base_bit(Encode_base_bit),
				m_C_base_bit(C_base_bit),
				m_h_FP(h_FP),
				m_h_GSW(h_GSW),
				m_h_PKS(h_PKS),
				m_exp_bit(exp_bit),
				m_PARALLEL(PARALLEL_FLAG)
				{

				m_gsw_dgg.SetStd(BK_std);
				m_ks_dgg.SetStd(KS_std);
				m_ev_dgg.SetStd(EV_std);
				m_pks_dgg.SetStd(PKS_std);
				
    PreCompute();
  }

  /**
   * Performs precomputations based on the supplied parameters
   */
  void PreCompute() {
   
		// TMP NK
		vector<NativeInteger> moduliQ = m_RLWEParams->GetModuliQ();
		/*22.02.03 remove*/ 
		
		/*
		
		vector<NativeInteger> rootsQ_tmp(2);
		rootsQ_tmp[0] = RootOfUnity<NativeInteger>(m_N_GSW*2,moduliQ[0]);
		rootsQ_tmp[1] = RootOfUnity<NativeInteger>(m_N_GSW*2,moduliQ[1]);
		//std::shared_ptr<ILDCRTParams<BigInteger>> 
		m_DCRTGSWNKParams= std::make_shared<ILDCRTParams<BigInteger>>(
				m_N_GSW * 2, moduliQ, rootsQ_tmp);
		*/
		const shared_ptr<FPRLWECryptoParams> rlweparams = m_RLWEParams;
		
		/* Q  Parameter Checking */
		m_moduliQ = moduliQ;
		m_mu_GSW.resize(m_moduliQ.size());
		for (uint32_t i = 0; i < m_moduliQ.size(); i++) {
			m_mu_GSW[i] = m_moduliQ[i].ComputeMu();
		}


		BigInteger Q_total = BigInteger(moduliQ[0]) * BigInteger(moduliQ[1]);
		BigInteger Mod2 = BigInteger(2);
		uint32_t Q_bits = 1;
		/* Check bits*/
		while (true) {
			if (Q_total == 1) {
				break;
			}
			Q_bits +=1;
			if (Q_total.Mod(Mod2) == 1) {
				Q_total = (Q_total - BigInteger(1)).DividedBy(Mod2);
			} else {
				Q_total = (Q_total).DividedBy(Mod2);
			}
		}
		m_Q_bit = Q_bits;
    m_Q_bit2 = m_Q_bit * 2;
	
		rlweparams->SetBitsTotal(Q_bits);
		
		// Calculate carry and overflow information
		// maximum Encryption 64 bits
		for (uint32_t i=0; i < 7; i++){
			uint32_t tmp = (1 << i);
			if (tmp % m_encode_base_bit == 0) {
				m_encode_len.push_back(tmp / m_encode_base_bit); 
			} else {
				m_encode_len.push_back((tmp / m_encode_base_bit)+1); 
			}
		}
	
		//if ((2*m_encode_len[6]) > rlweparams->GetN() ) {
		//	PALISADE_THROW(config_error, "Encoding cannot do 64 bit. square of length of msg is bigger than rlwe degree.");
		//}

	
		///////////////////////////// Encoding Ends/////////////////////

		std::cout << "Q total bits is " << Q_bits << std::endl;
		std::cout << "Encods len bits is " << m_encode_base_bit << std::endl;
	
		m_N_GSW_bit = (uint32_t)std::ceil(log(static_cast<double> (m_N_GSW)) /
                                    log(static_cast<double> (2)));
	

		////////// OTher Bit Info Calculation
	// GSW_N
	m_M_GSW = m_N_GSW * 2;
	m_M_GSW_bit = m_N_GSW_bit + 1;
	
	uint32_t th = 0;
	m_KS_l = m_Q_bit / m_KS_base_bit;
	if (m_Q_bit % m_KS_base_bit > th) {m_KS_l += 1;}
	m_KS_l_save = m_KS_l - m_KS_l_rm;
	

	m_EV_l = m_Q_bit /m_EV_base_bit;
	if (m_Q_bit % m_EV_base_bit > th) {m_EV_l += 1;}
	m_EV_l_save = m_EV_l - m_EV_l_rm;

	m_GD_l = m_Q_bit /m_GD_base_bit;
	if (m_Q_bit % m_GD_base_bit > th) {m_GD_l += 1;}
	m_GD_l_save = m_GD_l - m_GD_l_rm;
	
	m_PKS_l = m_Q_bit /m_PKS_base_bit;
	if (m_Q_bit % m_PKS_base_bit > th) {m_PKS_l += 1;}
	m_PKS_l_save = m_PKS_l - m_PKS_l_rm;
	

	// MSG 
	m_encode_base = 1;
	m_encode_base <<= m_encode_base_bit;
	// CarryUP Bit
	m_C_base = 1;
	m_C_base <<= m_C_base_bit;

	m_carry_mul_infos.resize(7);
	// Carry Calculation
	for (uint32_t j=0; j < 7; j++) {
		vector<uint64_t> carry_info(m_encode_len[j]);
		//vector<uint64_t> carry_info_res(m_encode_len[j]);
		m_carry_mul_infos[j].resize(m_encode_len[j]); 
		// Make zero
		for (uint32_t i = 0; i < m_encode_len[j]; i++ ){
			m_carry_mul_infos[j][i] = 0;
		}			uint32_t tmp = (1 << j); 
		for (uint32_t i = 0; i < m_encode_len[j]; i++) { 
      		if (tmp >= m_encode_base_bit) { 
				carry_info[i] = (1 << m_encode_base_bit) - 1; 
			} else if (tmp == 0){
				break;
			} else { 
				carry_info[i] = (1 << tmp); 
			} 
			tmp -= m_encode_base_bit; 
		}
		CalCarry(carry_info, carry_info, 0, m_carry_mul_infos[j]); 
	}

	m_carry_max = (uint32_t) std::ceil((double) *std::max_element(m_carry_mul_infos[6].begin(), m_carry_mul_infos[6].end()) / ((double) (m_C_base_bit) ));
    std::cout<< "carry max num is " << std::endl 
     << m_carry_max << std::endl
		 << "and vec value 64 bit is " << std::endl 
     << m_carry_mul_infos[6] << std::endl	
		 << "and vec value 32 bit is " << std::endl 
     << m_carry_mul_infos[5] << std::endl	
		 << "and vec value 16 bit is " << std::endl 
     << m_carry_mul_infos[4] << std::endl;	


	
	m_carry_mul_float_infos.resize(3);
	// Carry Calculation
	// 0 -> 16 bit, 1-> 32 bit, 2-> 64 bit
	uint32_t start_bits;
	uint32_t mul_start_bits;
	uint32_t len_bits;
	uint32_t mul_len_bits;
	for (uint32_t j=0; j<3; j++) {
		if (j == 0) {
			start_bits = 2;
			len_bits = 6;
			mul_len_bits = len_bits + 2;
			mul_start_bits = 15 - mul_len_bits + 1;
		} else if (j == 1) {
			start_bits = 3;
			len_bits = 13;
			mul_len_bits = len_bits + 3;
			mul_start_bits =  31 - mul_len_bits + 1;
		} else if (j == 2) {
			start_bits = 5;
			len_bits = 27;
			mul_len_bits = len_bits + 5;
			mul_start_bits =  63 - mul_len_bits + 1;

			}

			vector<uint64_t> carry_info(start_bits + len_bits); 
			m_carry_mul_float_infos[j].resize(mul_len_bits);  // It is enough	
			// Make zero
			for (uint32_t i = 0; i < mul_len_bits; i++ ){
				m_carry_mul_float_infos[j][i] = 0;
			}
			for (uint32_t i = start_bits; i < (start_bits + len_bits); i++) { 
				carry_info[i] = 3; 
			}
			CalCarryFloat(carry_info, carry_info, mul_start_bits, mul_len_bits, 		
					m_carry_mul_float_infos[j]);	
		
	}
		std::cout << "and vec value 64 float bit is " << std::endl 
     << m_carry_mul_float_infos[2] << std::endl	
		 << "and vec value 32 bit is " << std::endl 
     << m_carry_mul_float_infos[1] << std::endl	
		 << "and vec value 16 bit is " << std::endl 
     << m_carry_mul_float_infos[0] << std::endl;	

		// Test Carry 64 one bit
		/*
		vector<uint64_t> carry_info(64); 
		vector<uint64_t> carry_info_res(64);
	
		for (uint32_t i = 0; i < 64; i++) { 
			carry_info[i] = 1; 
		}
		CalCarry_one(carry_info, carry_info, 1, 0, carry_info_res); 
		
		std::cout << "one bit carry vec values is " 
			<< carry_info_res << std::endl; 
		*/

		// q1^2 / 2, q/2 Calculation
		BigInteger Q0 = BigInteger(m_moduliQ[0]);
		BigInteger Q1 = BigInteger(m_moduliQ[1]);

		// Make 3 ..

		MakeBootstrapMulPoly(m_prod_poly_bootstrap);		
		std::cout << "Mul Poly is gen" << std::endl;
		MakeBootstrapMulPolyOne(m_prod_poly_bootstrap_one);		
		std::cout << "Mul PolyOne is gen" << std::endl;
		MakeBootstrapAddPoly(m_add_poly_bootstrap);		
		std::cout << "Add Poly is gen" << std::endl;
		MakeBootstrapAddPolyOne(m_add_poly_bootstrap_one);		
		std::cout << "Add One Poly is gen" << std::endl;
		MakeBootstrapReluPoly(m_relu_poly_bootstrap);
		std::cout << "Relu Poly is gen" << std::endl;
	
		
		MakeBootstrapDivPoly(m_div_poly_bootstrap);		
		std::cout << "Div Poly is gen" << std::endl;
		MakeBootstrapReconPoly(m_recon_poly_bootstrap);		
		std::cout << "Rev Div Boot Poly is End!!" << std::endl;
		
	
		// GINX bootstrapping
		// CHECK!!
		// loop for positive values of m		
  	
		m_monomials.resize(m_N_GSW*2);
		//#pragma omp parallel for 
		{for (uint32_t i = 0; i < m_N_GSW; i++) {
																			
			DCRTPoly aPoly = DCRTPoly(m_DCRTGSWParams, Format::COEFFICIENT);
			NativePoly tmp0 = NativePoly(m_DCRTGSWParams->GetParams()[0], Format::COEFFICIENT, true);
			tmp0[i] = tmp0[i].ModAddEq(NativeInteger(1), m_moduliQ[0]);  // X^m	
			tmp0[0] = tmp0[0].ModSubEq(NativeInteger(1), m_moduliQ[0]);  // - 1
			aPoly = tmp0;		
			aPoly.SetFormat(Format::EVALUATION);   
			m_monomials[i]					 =	aPoly;

		}}

    // loop for negative values of m
		//#pragma omp parallel for 
		{for (uint32_t i = 0; i < m_N_GSW; i++) {

			DCRTPoly aPoly = DCRTPoly(m_DCRTGSWParams, Format::COEFFICIENT);
			NativePoly tmp0 = NativePoly(m_DCRTGSWParams->GetParams()[0], Format::COEFFICIENT, true);

			tmp0[i] = tmp0[i].ModSubEq(NativeInteger(1), m_moduliQ[0]);  // -X^m=X^(-m)
			tmp0[0] = tmp0[0].ModSubEq(NativeInteger(1), m_moduliQ[0]);  // - 1
			aPoly = tmp0; 
			aPoly.SetFormat(Format::EVALUATION);
			m_monomials[i + m_N_GSW] = aPoly;

		}}
		std::cout << "Monomials are Gen!!!!!!" << std::endl;
		////////////////////////////////////////////////////////////////////

		MakeRotMonomial(m_N_GSW, m_DCRTGSWParams, m_rot_monomials_GSW);	
		std::cout << " GSW Rotation are Gen!!!!!!" << std::endl;
		
		MakeRotMonomial(m_RLWEParams->GetN(), m_RLWEParams->GetDCRTParams(),m_rot_monomials_RLWE);
		std::cout << "RLWE Rotation are Gen!!!!!!" << std::endl;
		

		// Orign Bootstrapping Ready
		//uint32_t msg;
		uint32_t cut = m_N_GSW_bit - m_encode_base_bit;
		NativePoly msg_poly0 = NativePoly(m_DCRTGSWParams->GetParams()[0], Format::COEFFICIENT, true);
		NativePoly msg_poly1 = NativePoly(m_DCRTGSWParams->GetParams()[1], Format::COEFFICIENT, true);
		//this parrallel has an error
		// Cal moduliQ with ModMul has an error in parallel processing
		//#pragma omp parallel for 
		{for (uint32_t i=0; i < m_N_GSW; i++) {
			uint32_t msg = i >> (cut);
			NativeInteger Q1_local = m_moduliQ[1];
			NativeInteger Q0_local = m_moduliQ[0];
			// Bootstrapping
			// 
			// msg * Q1 mod Q0
			msg_poly0[i] = Q1_local.ModMul(msg, Q0_local);
					
			// msg * Q1 mod Q1
			msg_poly1[i] = Q1_local.ModMul(msg, Q1_local);
		}}

		m_orign_bootstrap = DCRTPoly(m_DCRTGSWParams, Format::COEFFICIENT, true);
		m_orign_bootstrap.SetElementAtIndex(0, msg_poly0);
		m_orign_bootstrap.SetElementAtIndex(1, msg_poly1);
		m_orign_bootstrap.SetFormat(Format::EVALUATION);
		std::cout << "Bootstrap Orign is Gen!!" << std::endl;

	#if defined(FPFHE_DEBUG)
    
		/*
		std::cerr << "base_g = " << m_GD_base_bit << std::endl;
    std::cerr << "m_Q_bit = " << m_digitsG << std::endl;
    std::cerr << "m_Q_bit2 = " << m_digitsG2 << std::endl;
    //std::cerr << "m_baseR = " << m_baseR << std::endl;
    std::cerr << "m_digitsR = " << m_digitsR << std::endl;
    std::cerr << "m_Gpower = " << m_Gpower << std::endl;
    std::cerr << "n = " << m_LWEParams->Getn() << std::endl;
    std::cerr << "N = " << m_LWEParams->GetN() << std::endl;
    std::cerr << "q = " << m_LWEParams->Getq() << std::endl;
    std::cerr << "Q = " << m_LWEParams->GetQ() << std::endl;
    std::cerr << "baseKS = " << m_LWEParams->GetBaseKS() << std::endl;
    std::cerr << "digitsKS = " << m_LWEParams->GetDigitsKS() << std::endl;
		*/
#endif
		std::cout << "========================== FPRINGCORE PRECOMPUTE IS END =============================" << std::endl;
  }
	/* Carry Calculation*/
	/* It Must be rewrite */
	void CalCarry(vector<uint64_t> max_m1, vector<uint64_t> max_m2, uint32_t start_idx, vector<uint64_t> &out_carry) const{
		uint32_t nums = max_m1.size(); // nums shoud be end idx + 1
		
		// 2 bit + 2 bit Handle
		//m_C_base_bit = 4;
		//m_encode_base_bit = 2;
		//uint32_t Base_Carry = 1 << m_C_base_bit;
		uint32_t Jump = (m_C_base_bit / m_encode_base_bit);
		if (Jump * m_encode_base_bit != m_C_base_bit) {
			PALISADE_THROW(config_error, "carry Bit should be multiples of baseBit.");
		}

		NativePoly P1 = NativePoly(m_DCRTGSWParams->GetParams()[0],Format::COEFFICIENT,true);
		NativePoly P2 = NativePoly(m_DCRTGSWParams->GetParams()[0],
				Format::COEFFICIENT,true);
	
		for (uint32_t i = 0; i< nums; i++) {
			P1[i] = max_m1[i];
			P2[i] = max_m2[i];
		}
		P1.SetFormat(Format::EVALUATION);
		P2.SetFormat(Format::EVALUATION);
		P1 *= P2;
		P1.SetFormat(Format::COEFFICIENT);
		
		
		uint32_t tmp;
		uint32_t carry_idx;
		for (uint32_t i = 0; i < nums; i ++) {
			tmp = P1[i+start_idx].ConvertToInt();
			
			carry_idx = 0;
			while(true) {
				// No carry to give
				if (tmp == 0) {
					break;
				}
				// We will Cut after nums
				if ( (i + start_idx+ (carry_idx)*Jump >= nums)) {
					break;	
				}
				
				if (carry_idx == 0) {
					if (tmp >= 15) {
						out_carry[i] += 4;
					} else {
						out_carry[i] +=	this->GetMSBBitLoc(tmp);
					}
					// carry_idx == 0 Addding
					if ((i + 1 + start_idx) != nums) {
						// Base Carry must 15
						//P1[i+1+start_idx] += (std::min(tmp,Base_Carry) >> m_encode_base_bit);
						P1[i+1+start_idx] += (tmp >> 2) % 4;
					
					}
		
				} else {
					// carry is occurs
					if (tmp >= 15) {
						out_carry[i] += 4;
					} else {
						out_carry[i] +=	this->GetMSBBitLoc(tmp);
					}
					// carry_idx ==0 does not fit
					//P1[i+start_idx+(carry_idx)*Jump] += std::min(tmp,Base_Carry);
					P1[i+start_idx+(carry_idx)*Jump] += tmp % 15;
				}
				// Update
				if (carry_idx == 0) {
					tmp = tmp >> 4;
				} else {
					tmp = tmp >> 4;
				}
				carry_idx += 1;
			}
		}
		return ;
	}

void CalCarryFloat(
		const vector<uint64_t> max_m1, 
		const vector<uint64_t> max_m2, 
		const uint32_t start_idx, 
		const uint32_t frac_len, 
		vector<uint64_t> &out_carry) const{
		
		uint32_t nums = max_m1.size();
		uint32_t nums2 = nums * 2;  // nums2 shoud be end idx + 1
		// 2 bit + 2 bit Handle
		//m_C_base_bit = 4;
		//m_encode_base_bit = 2;
		uint32_t Base_Carry = 1 << m_C_base_bit;
		uint32_t Jump = (m_C_base_bit / m_encode_base_bit);
		if (Jump * m_encode_base_bit != m_C_base_bit) {
			PALISADE_THROW(config_error, "carry Bit should be multiples of baseBit.");
		}

		NativePoly P1 = NativePoly(m_DCRTGSWParams->GetParams()[0],Format::COEFFICIENT,true);
		NativePoly P2 = NativePoly(m_DCRTGSWParams->GetParams()[0],
				Format::COEFFICIENT,true);
	
		//std::cout << "Before Seg test 1" << std::endl;
		for (uint32_t i = 0; i< nums; i++) {
			P1[i] = max_m1[i];
			P2[i] = max_m2[i];
		}
		//std::cout << "Before Seg test 2" << std::endl;
	
		P1.SetFormat(Format::EVALUATION);
		P2.SetFormat(Format::EVALUATION);
		P1 *= P2;
		P1.SetFormat(Format::COEFFICIENT);
		P2.SetFormat(Format::COEFFICIENT);

		std::cout << "Orign val is " << std::endl;
		std::cout << P2 << std::endl;
		std::cout << "Float product is " << std::endl;
		std::cout << P1 << std::endl;
		std::cout << "nums2 are" << nums2 << ", from " << start_idx << " to " << (frac_len + start_idx - 1) << ", total len is "<< frac_len << std::endl;
		

		// After Segments are occur

		for (uint32_t i = 0; i < frac_len; i ++) {
			std::cout << P1[start_idx+i] << ", ";
		}
		std::cout << std::endl;
		uint32_t tmp;
		uint32_t carry_idx;
		
		//std::cout << "Segment Test 1" << std::endl;
		for (uint32_t i = 0; i < frac_len; i ++) {
			tmp = P1[i+start_idx].ConvertToInt();
			//std::cout << "i is " << i << ", tmp is " << tmp << std::endl;
			carry_idx = 0;
			while(true) {
				// No carry to give
				if (tmp == 0) {
					break;
				}

				// We will Cut after nums
				//if ( (i + start_idx+ (carry_idx)*Jump >= nums2)) {
				//	break;	
				//}
				
				if (carry_idx == 0) {
					if (tmp >= 15) {
						out_carry[i] += 4;
					} else {
						out_carry[i] +=	this->GetMSBBitLoc(tmp);
					}

					// carry_idx == 0 Addding
					if ((i + 1 + start_idx) != nums2) {
						//uint32_t three = 3;
						//P1[i+1+start_idx] += std::min(tmp,three);
					
						uint32_t tmp2 = tmp >> 2;
						P1[i+1+start_idx] += tmp2 % 4;
					
					}

				} else {
					// carry is occurs
					if (tmp >= Base_Carry) {
						out_carry[i] += 4;
					} else {
						out_carry[i] +=	this->GetMSBBitLoc(tmp);
					}
					// carry_idx ==0 does not fit
					//uint32_t carry_base = 15;
					//P1[i+start_idx+(carry_idx)*Jump] += std::min(tmp,carry_base);
					P1[i+start_idx+(carry_idx)*Jump] += tmp % 16;
				
				}
				// Update
				if (carry_idx == 0) {
					tmp = tmp >> m_C_base_bit;
				} else {
					tmp = tmp >> m_C_base_bit;
				}
				carry_idx += 1;
			}
		}
		//std::cout << "Segment Test 2" << std::endl;
		//std:: cout << "out carry is " << out_carry << std::endl;
		return ;
	}




	void CalCarry_one(vector<uint64_t> max_m1, vector<uint64_t> max_m2, 
			uint32_t Carry_base_bit,
			uint32_t start_idx, vector<uint64_t> &out_carry) const{
		

		uint32_t nums = max_m1.size();
		uint32_t Base_Carry = 1 << Carry_base_bit;
		
		NativePoly P1 = NativePoly(m_RLWEParams->GetDCRTParams()->GetParams()[0],Format::COEFFICIENT,true);
		NativePoly P2 = NativePoly(m_RLWEParams->GetDCRTParams()->GetParams()[0],
				Format::COEFFICIENT,true);
	
		for (uint32_t i = 0; i< nums; i++) {
			P1[i] = max_m1[i];
			P2[i] = max_m2[i];
		}
		P1.SetFormat(Format::EVALUATION);
		P2.SetFormat(Format::EVALUATION);
		P1 *= P2;
		P1.SetFormat(Format::COEFFICIENT);
		
		uint32_t tmp;
		uint32_t carry_idx;
		for (uint32_t i = 0; i < nums; i ++) {
			tmp = P1[i+start_idx].ConvertToInt();
			
			carry_idx = 0;
			while(true) {
				// No carry to give
				if (tmp == 0) {
					break;
				}
				// We will Cut after nums
				if ( (i + start_idx + (carry_idx) >= nums)) {
					break;	
				}
				
				if (carry_idx == 0) {
					if (tmp >= Base_Carry) {
						out_carry[i] += Carry_base_bit;
					} else {
						out_carry[i] +=	this->GetMSBBitLoc(tmp);
					}
					// carry_idx == 0 Addding
					//P1[i+1+start_idx] += (std::min(tmp,Base_Carry) >> m_encode_base_bit);
			
				} else {
					// carry is occurs
					if (tmp >= Base_Carry) {
						out_carry[i] += Carry_base_bit;
					} else {
						out_carry[i] +=	this->GetMSBBitLoc(tmp);
					}
					// carry_idx ==0 does not fit
					P1[i+start_idx+(carry_idx)] += std::min(tmp,Base_Carry);
				}
				// Update
				tmp = tmp >> Carry_base_bit;
				carry_idx += 1;
			}
		}
		return ;
	}


	uint32_t GetMSBBitLoc(uint32_t tmp ) const {
		uint32_t out = (uint32_t)std::ceil(log(static_cast<double>(tmp)) /
                              log(static_cast<double>(2)));
		if (tmp == (uint32_t) (1<< out)) {
			out += 1;
		}
		return out;
	}


	void MakeBootstrapReluPoly(vector<DCRTPoly> &result) const {
		
		int32_t msg;	
		uint64_t scaling = GetRLWEParams()->GetScaling().ConvertToInt();
		scaling <<= GetExpBit();
		NativeInteger Q1 = m_moduliQ[1];
		NativeInteger Q0 = m_moduliQ[0];

		//NativeInteger scaling_native = NativeInteger(scaling);
		
		result.resize(20); 
		// Out Vector 0
		result[0] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		// 4 & bit reverse
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign
			if (i % 2 == 0) {
				msg = i >> (m_N_GSW_bit - 4); // 4bit		
				// 0*Q1, 1*Q1, 2*Q1, 3*Q1, 0*Q1, 1*Q1, 2*Q1, 3*Q1
				if (msg >= 16) {
					PALISADE_THROW(config_error, "Coding is not correct");
				}
				//msg = msg % 16;
				result[0].GetElementW(0)[i].ModAddEq(Q1 * msg, Q0);
				// ODD: Bootstrapping or Carry
			} else if (i % 2 == 1) {
				result[0].GetElementW(0)[i].ModAddEq(Q1, Q0);
			}
		}
		result[0].SetFormat(Format::EVALUATION);

		
		// Out Vector 1
		result[1] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);

		// 4 bit, Origns and Ups
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign
			if (i % 2 == 0) {
				msg = i >> (m_N_GSW_bit - 4);	// 4bit	
				// 0*Q1, 1*Q1, 2*Q1, 3*Q1, 4*Q1, 5*Q1, 6*Q1, 7*Q1
				result[1].GetElementW(0)[i].ModAddEq(m_moduliQ[1] * msg, m_moduliQ[0]);
		
				// ODD: Bootstrapping or Carry
			} else if (i % 2 == 1) {
				// Adjust Sign
				msg = i >> (m_N_GSW_bit - 4);	// 4bit	
				// 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1
				result[1].GetElementW(0)[i].ModAddEq(
					Q1.ModMul(scaling,Q0).ModMul(msg, Q0), Q0);
				result[1].GetElementW(1)[i].ModAddEq(
						Q1.ModMul(scaling,Q1).ModMul(msg, Q1), Q1);
			}
		}
		result[1].SetFormat(Format::EVALUATION);	
		
		
		// Make
		result[2] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		// 4 bit, Origns and Ups
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign
			if (i % 2 == 0) {
				msg = i >> (m_N_GSW_bit - 4);	// 4bit	
				// 0*Q1, 1*Q1, 2*Q1, 3*Q1, 4*Q1, 5*Q1, 6*Q1, 7*Q1
				result[2].GetElementW(0)[i].ModAddEq(m_moduliQ[1] * msg, m_moduliQ[0]);
		
				// ODD: Bootstrapping or Carry
			} else if (i % 2 == 1) {
				// Adjust Sign
				msg = i >> (m_N_GSW_bit - 4);	// 4bit	
				// 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1
				// 4bit << 
				msg <<= 4;
				result[2].GetElementW(0)[i].ModAddEq(
					Q1.ModMul(scaling,Q0).ModMul(msg, Q0), Q0);
				result[2].GetElementW(1)[i].ModAddEq(
						Q1.ModMul(scaling,Q1).ModMul(msg, Q1), Q1);
				}
		}
		result[2].SetFormat(Format::EVALUATION);	
	

		// Out Vector 2 sign adjust
		result[3] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		// 4bit Orign, 1bit 0 or 1
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign One more bit inserted
			if (i % 2 == 0) {
				msg = i >> (m_N_GSW_bit - 4);	// Get 4bit, 5bit should be 0	
				if (msg >= 16) {
					PALISADE_THROW(config_error, "Coding is not correct");
				}
				if (msg < 4) {
					msg = 0;
				} else {
					msg = msg  - 4; // remove adding effects 2^10
					msg <<= 8;
				}		
				result[3].GetElementW(0)[i].ModAddEq(
					Q1.ModMul(scaling,Q0).ModMul(msg, Q0), Q0);
				result[3].GetElementW(1)[i].ModAddEq(
					Q1.ModMul(scaling,Q1).ModMul(msg, Q1), Q1);
	

				// ODD: Bootstrapping or Carry
			} else if (i % 2 == 1) {

				//msg = i >> (m_N_GSW_bit - 4); // should be 4bit
				//msg = msg >> 3; // MSB 1 bit
				msg = i >> (m_N_GSW_bit - 4); // should be 4bit
				if (msg < 4) {
					msg = 0; 
				} else {
					msg = 1;
				}
				//msg = 0 ;
				result[3].GetElementW(0)[i].ModAddEq(msg * m_moduliQ[1], m_moduliQ[0]);
			}
		}
		result[3].SetFormat(Format::EVALUATION);
		
		// Out Vector 2
		result[4] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		// 1~4bit reconstruction
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign One more bit inserted
			msg = i >> (m_N_GSW_bit - 4); // 4 bits
			// 0Q1 ~ 15Q1
			result[4].GetElementW(0)[i].ModAddEq(
				Q1.ModMul(scaling,Q0).ModMul(msg, Q0), Q0);
			result[4].GetElementW(1)[i].ModAddEq(
				Q1.ModMul(scaling,Q1).ModMul(msg, Q1), Q1);
		}
	
		result[4].SetFormat(Format::EVALUATION);
		

		
		// Out Vector 3 , 4 bit Q1Q1 vals
		result[5] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		
		// Reconstruct 4~7bit
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign One more bit inserted
			msg = i >> (m_N_GSW_bit - 4); // 4 bits
			msg <<= 4;
		
			// 0Q1 ~ 15Q1
			result[5].GetElementW(0)[i].ModAddEq(
				Q1.ModMul(scaling,Q0).ModMul(msg, Q0), Q0);
			result[5].GetElementW(1)[i].ModAddEq(
				Q1.ModMul(scaling,Q1).ModMul(msg, Q1), Q1);
		}
		result[5].SetFormat(Format::EVALUATION);
	
		//////////////// Out Vector 5 /////////////////
		

		// Out Vector 2 sign adjust
		result[6] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		// 4bit Orign, 1bit 0 or 1
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign One more bit inserted
			if (i % 2 == 0) {
				msg = i >> (m_N_GSW_bit - 4);	// Get 4bit, 5bit should be 0	
				if (msg >= 16) {
					PALISADE_THROW(config_error, "Coding is not correct");
				}
				if (msg < 4) {
					// Remove effect 2^10
					msg = 4 - msg;
					msg <<= 8;
				
				} else {
					msg = 0;
				}
			
				result[6].GetElementW(0)[i].ModSubEq(
					Q1.ModMul(scaling,Q0).ModMul(msg, Q0), Q0);
				result[6].GetElementW(1)[i].ModSubEq(
					Q1.ModMul(scaling,Q1).ModMul(msg, Q1), Q1);

				// ODD: Bootstrapping or Carry
			} else if (i % 2 == 1) {

				//msg = i >> (m_N_GSW_bit - 4); // should be 4bit
				//msg = msg >> 3; // MSB 1 bit
				
				msg = i >> (m_N_GSW_bit - 4); // should be 4bit
				if (msg < 4) {
					msg = 1; 
				} else {
					msg = 0;
				}
				
				result[6].GetElementW(0)[i].ModAddEq(msg * m_moduliQ[1], m_moduliQ[0]);
			}
		}
		result[6].SetFormat(Format::EVALUATION);
		


		// Out Vector 7 for Floating Points
		result[7] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		// 1~4bit reconstruction
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign One more bit inserted
			if (i % 2 == 0) {
				msg = i >> (m_N_GSW_bit - 4);	// Get 4bit, 5bit should be 0	
				if (msg >= 16) {
					PALISADE_THROW(config_error, "Coding is not correct");
				}
				if (msg >= 8) {
					msg = 1;
				} else {
					msg = 0;
				}
				result[7].GetElementW(0)[i].ModAddEq(Q1 * msg, Q0);
				result[7].GetElementW(1)[i].ModAddEq(Q1 * msg, Q1);
			} else {

				msg = i >> (m_N_GSW_bit - 4);	// Get 4bit, 5bit should be 0	
				if (msg >= 8) {
					msg = msg - 8;
				} else {
					msg = 0;
				}
				msg <<= 4;
				result[7].GetElementW(0)[i].ModSubEq(
				Q1.ModMul(scaling,Q0).ModMul(msg, Q0), Q0);
				result[7].GetElementW(1)[i].ModSubEq(
				Q1.ModMul(scaling,Q1).ModMul(msg, Q1), Q1);
			}
			// ODD: Bootstrapping or Carry 
		}
		result[7].SetFormat(Format::EVALUATION);
		


		// Out Vector 8 sign adjust
		result[8] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		// 4bit Orign, 1bit 0 or 1
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign One more bit inserted
			if (i % 2 == 0) {
				msg = i >> (m_N_GSW_bit - 4);	// Get 4bit, 5bit should be 0	
				if (msg >= 16) {
					PALISADE_THROW(config_error, "Coding is not correct");
				}
				if (msg < 8) { // minus
					msg = 1;
				} else {
					msg = 0;
				}
				result[8].GetElementW(0)[i].ModAddEq(Q1 * msg, Q0);
				result[8].GetElementW(1)[i].ModAddEq(Q1 * msg, Q1);
			} else {
				msg = i >> (m_N_GSW_bit - 4);	// Get 4bit, 5bit should be 0	
				
				if (msg < 8 ) {
					msg = 8 - msg;
				} else {
					msg = 0;
				}
				msg <<= 4;
				result[8].GetElementW(0)[i].ModSubEq(
				Q1.ModMul(scaling,Q0).ModMul(msg, Q0), Q0);
				result[8].GetElementW(1)[i].ModSubEq(
				Q1.ModMul(scaling,Q1).ModMul(msg, Q1), Q1);
			}
		}
		result[8].SetFormat(Format::EVALUATION);	




		// Out Vector 2 sign adjust
		result[9] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		// 4bit Orign, 1bit 0 or 1
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign One more bit inserte
			msg = i >> (m_N_GSW_bit - 5);	// Get 4bit, 5bit should be 0	
			if (msg >= 12) { // 8 + adding vals 4
				result[9].GetElementW(0)[i].ModAddEq(
					Q1.ModMul(scaling,Q0).ModMul(1, Q0), Q0);
				result[9].GetElementW(1)[i].ModAddEq(
					Q1.ModMul(scaling,Q1).ModMul(1, Q1), Q1);

				// ODD: Bootstrapping or Carry
			}	
		}
		result[9].SetFormat(Format::EVALUATION);
		

		// Out Vector 2 sign adjust
		result[10] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		// 4bit Orign, 1bit 0 or 1
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign One more bit inserte
			msg = i >> (m_N_GSW_bit - 5);	// Get 4bit, 5bit should be 0	
			if (msg >= 24) { // 16 + adding val 8
				result[10].GetElementW(0)[i].ModAddEq(
					Q1.ModMul(scaling,Q0).ModMul(1, Q0), Q0);
				result[10].GetElementW(1)[i].ModAddEq(
					Q1.ModMul(scaling,Q1).ModMul(1, Q1), Q1);

				// ODD: Bootstrapping or Carry
			}	
		}
		result[10].SetFormat(Format::EVALUATION);





	
		// Out Vector 1 OFUF
		result[11] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		// 4 bit, Origns and Ups
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign
			if (i % 2 == 0) {
				msg = i >> (m_N_GSW_bit - 4);	// 4bit	 UF
				if (msg == 0) {
					msg = 1;
				} else {
					msg = 0;
				}
				// 0*Q1, 1*Q1, 2*Q1, 3*Q1, 4*Q1, 5*Q1, 6*Q1, 7*Q1
				result[11].GetElementW(0)[i].ModAddEq(m_moduliQ[1] * msg, m_moduliQ[0]);
		
				// ODD: Bootstrapping or Carry
			} else if (i % 2 == 1) {
				// Adjust Sign
				msg = i >> (m_N_GSW_bit - 4);	// 4bit	
				// 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1
				result[11].GetElementW(0)[i].ModAddEq(
					Q1.ModMul(scaling,Q0).ModMul(msg, Q0), Q0);
				result[11].GetElementW(1)[i].ModAddEq(
						Q1.ModMul(scaling,Q1).ModMul(msg, Q1), Q1);
			}
		}
		result[11].SetFormat(Format::EVALUATION);	
		
		
		// Make
		result[12] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		// 4 bit, Origns and Ups
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign
			if (i % 2 == 0) {
				msg = i >> (m_N_GSW_bit - 4);	// 4bit UF	
				if (msg == 0) {
					msg = 1;
				} else {
					msg = 0;
				}
				// 0*Q1, 1*Q1, 2*Q1, 3*Q1, 4*Q1, 5*Q1, 6*Q1, 7*Q1
				result[12].GetElementW(0)[i].ModAddEq(m_moduliQ[1] * msg, m_moduliQ[0]);
		
				// ODD: Bootstrapping or Carry
			} else if (i % 2 == 1) {
				// Adjust Sign
				msg = i >> (m_N_GSW_bit - 4);	// 4bit	
				// 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1
				// 4bit << 
				msg <<= 4;
				
				result[12].GetElementW(0)[i].ModAddEq(
					Q1.ModMul(scaling,Q0).ModMul(msg, Q0), Q0);
				result[12].GetElementW(1)[i].ModAddEq(
						Q1.ModMul(scaling,Q1).ModMul(msg, Q1), Q1);
				}
		}
		result[12].SetFormat(Format::EVALUATION);	


	
		// Out Vector 1 // 64bit OFUF
		result[13] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		// 4 bit, Origns and Ups
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			if (i % 2 == 0) {
			// Orign
				msg = i >> (m_N_GSW_bit - 4);	// 4bit	 OF
				if (msg >= 8) {
					msg = 1;
				} else {
					msg = 0;
				}
				// 0*Q1, 1*Q1, 2*Q1, 3*Q1, 4*Q1, 5*Q1, 6*Q1, 7*Q1
				result[13].GetElementW(0)[i].ModAddEq(m_moduliQ[1] * msg, m_moduliQ[0]);
				/*
				result[13].GetElementW(0)[i].ModAddEq(
					Q1.ModMul(scaling,Q0).ModMul(msg, Q0), Q0);
				result[13].GetElementW(1)[i].ModAddEq(
					Q1.ModMul(scaling,Q1).ModMul(msg, Q1), Q1);
				*/
			} else {
				msg = i >> (m_N_GSW_bit - 4);	// 4bit	// UF
				if (msg == 0) {
					msg = 1;
				} else {
					msg = 0;
				}
				result[13].GetElementW(0)[i].ModAddEq(m_moduliQ[1] * msg, m_moduliQ[0]);
			}
		}
		result[13].SetFormat(Format::EVALUATION);		
	
		// Out Vector 1 // 32bit OFUF
		result[14] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		// 4 bit, Origns and Ups
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			if (i % 2 == 0) {
				// Orign
				msg = i >> (m_N_GSW_bit - 4);	// 4bit	 OF
				if (msg >= 4) {
					msg = 1;
				} else {
					msg = 0;
				}
				// 0*Q1, 1*Q1, 2*Q1, 3*Q1, 4*Q1, 5*Q1, 6*Q1, 7*Q1
				//result[14].GetElementW(0)[i].ModAddEq(m_moduliQ[1] * msg, m_moduliQ[0]);
				//result[14].GetElementW(0)[i].ModAddEq(
				//	Q1.ModMul(scaling,Q0).ModMul(msg, Q0), Q0);
				//result[14].GetElementW(1)[i].ModAddEq(
				//	Q1.ModMul(scaling,Q1).ModMul(msg, Q1), Q1);
	
				result[14].GetElementW(0)[i].ModAddEq(m_moduliQ[1] * msg, m_moduliQ[0]);
	
			} else {
				msg = i >> (m_N_GSW_bit - 4);	// 4bit	// UF
				if (msg == 0) {
					msg = 1;
				} else {
					msg = 0;
				}
				result[14].GetElementW(0)[i].ModAddEq(m_moduliQ[1] * msg, m_moduliQ[0]);
	
			}
		}
		result[14].SetFormat(Format::EVALUATION);	


		// Out Vector 1 // 32bit OFUF
		result[15] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		// 4 bit, Origns and Ups
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign
			msg = i >> (m_N_GSW_bit - 4);	// 5bit	
			if (msg >= 1) {
				msg = 1;
			} else {
				msg = 0;
			}
			// 0*Q1, 1*Q1, 2*Q1, 3*Q1, 4*Q1, 5*Q1, 6*Q1, 7*Q1
			result[15].GetElementW(0)[i].ModAddEq(
				Q1.ModMul(scaling,Q0).ModMul(msg, Q0), Q0);
			result[15].GetElementW(1)[i].ModAddEq(
				Q1.ModMul(scaling,Q1).ModMul(msg, Q1), Q1);
		}
		result[15].SetFormat(Format::EVALUATION);	


		// Out Vector 1 // 32bit OFUF
		result[16] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		// 4 bit, Origns and Ups
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign
			if (i % 2 == 0) {
				msg = i >> (m_N_GSW_bit - 4);	// 5bit	
				if (msg == 1) {
					msg = 1;
				} else if (msg == 0){
					msg = 0;
				} else {
					msg = 3;
				}
				// 0*Q1, 1*Q1, 2*Q1, 3*Q1, 4*Q1, 5*Q1, 6*Q1, 7*Q1
				result[16].GetElementW(0)[i].ModAddEq(
					Q1.ModMul(scaling,Q0).ModMul(msg, Q0), Q0);
				result[16].GetElementW(1)[i].ModAddEq(
					Q1.ModMul(scaling,Q1).ModMul(msg, Q1), Q1);
			} else {

				msg = i >> (m_N_GSW_bit - 4);	// 5bit	
				if (msg == 1) {
					msg = 0;
				} else if (msg == 0){
					msg = 1;
				} else {
					msg = 3;
				}
				result[16].GetElementW(0)[i].ModAddEq(m_moduliQ[1] * msg, m_moduliQ[0]);
			}
		}
		result[16].SetFormat(Format::EVALUATION);	



		// Out Vector 1 // 32bit OFUF
		result[17] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		// 4 bit, Origns and Ups
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign
			msg = i >> (m_N_GSW_bit - 4);	// 4bit	
			if (msg == 1) {
				msg = 1;
			} else {
				msg = 0;
			}
			// 0*Q1, 1*Q1, 2*Q1, 3*Q1, 4*Q1, 5*Q1, 6*Q1, 7*Q1
			result[17].GetElementW(0)[i].ModAddEq(
				Q1.ModMul(scaling,Q0).ModMul(msg, Q0), Q0);
			result[17].GetElementW(1)[i].ModAddEq(
				Q1.ModMul(scaling,Q1).ModMul(msg, Q1), Q1);
		}	
		result[17].SetFormat(Format::EVALUATION);	


		// Out Vector 8 sign adjust
		result[18] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		// 4bit Orign, 1bit 0 or 1
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign One more bit inserted
			msg = i >> (m_N_GSW_bit - 4);	// Get 4bit, 5bit should be 0	
			if (msg >= 16) {
				PALISADE_THROW(config_error, "Coding is not correct");
			}
			if (msg == 0) { // minus
				msg = 1;
			} else {
				msg = 0;
			}	
			result[18].GetElementW(0)[i].ModAddEq(Q1 * msg, Q0);
			result[18].GetElementW(1)[i].ModAddEq(Q1 * msg, Q1);
			
		}
		result[18].SetFormat(Format::EVALUATION);	



		// Out Vector 2 sign adjust
		result[19] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		// 4bit Orign, 1bit 0 or 1
		for (uint32_t i = 0; i < m_N_GSW; i++) {	
			msg = i >> (m_N_GSW_bit - 4); // should be 4bit
			if (msg < 8) {
				msg = 1; 
			} else {
				msg = 0;
			}
			result[19].GetElementW(0)[i].ModAddEq(msg * m_moduliQ[1], m_moduliQ[0]);	
		}
		result[19].SetFormat(Format::EVALUATION);
		





	}



	void MakeBootstrapDivPoly(vector<DCRTPoly> & result) const {	
		
		int32_t msg;	
		uint64_t scaling = GetRLWEParams()->GetScaling().ConvertToInt();
		scaling <<= GetExpBit();
		NativeInteger Q0 = m_moduliQ[0];
		NativeInteger Q1 = m_moduliQ[1];
		//NativeInteger scaling_native = NativeInteger(scaling);
		//NativeInteger Q1over2 = NativeInteger(m_moduliQ[1].ConvertToInt() >> 1);
		result.resize(6); 
		// Out Vector 0
		result[0] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		//NativeInteger Q1Q1  = m_moduliQ[1].ModMul(m_moduliQ[1], m_moduliQ[0]);
		// 3bit get, 2bit carry & 0,1,2 out, 3 is rev
		for (uint32_t i = 0; i < m_N_GSW; i++) {
		
			// Origna
			if (i % 4 == 0) {
				msg = i >> (m_N_GSW_bit - 3); // 3bit		
				// 0*Q1, 1*Q1, 2*Q1, 3*Q1, 0*Q1, 1*Q1, 2*Q1, 3*Q1
				msg = msg % 2;
				result[0].GetElementW(0)[i].ModAddEq(m_moduliQ[1] * msg, m_moduliQ[0]);
				// ODD: Bootstrapping or Carry
			} else if (i % 4 == 1) {
				msg = i >> (m_N_GSW_bit - 3); // 3bit		
				msg >>= 1;
				msg %=  2;
				result[0].GetElementW(0)[i].ModAddEq(m_moduliQ[1] * msg, m_moduliQ[0]);
			} else if (i % 4 == 2) {
				msg = i >> (m_N_GSW_bit - 3); // 3bit		
				msg >>= 2;
				msg %=  2;
				result[0].GetElementW(0)[i].ModAddEq(m_moduliQ[1] * msg, m_moduliQ[0]);
		
			} else if (i % 4 == 3) {
				// Adjust Flip
				result[0].GetElementW(0)[i].ModAddEq(m_moduliQ[1], m_moduliQ[0]);
			}
		}
		result[0].SetFormat(Format::EVALUATION);

		
				
		// Out Vector 1
		// 2^5 * m^ (-...) 
		// 3bit is adjusted 
		result[1] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);

		// 2 bit, out
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign
			msg = i >> (m_N_GSW_bit - 2);	// 2bit	
			if (i % 4 == 0) {
				if (msg == 0) { // 00
					// 0*Q1, 1*Q1, 2*Q1, 3*Q1, 4*Q1, 5*Q1, 6*Q1, 7*Q1
					result[1].GetElementW(0)[i].ModAddEq(m_moduliQ[1], m_moduliQ[0]);
				}
				// ODD: Bootstrapping or Carry
			} else if (i % 4 == 1) {
				// Adjust Sign
				if (msg == 1) {
					// 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1
					result[1].GetElementW(0)[i].ModAddEq(m_moduliQ[1],  m_moduliQ[0]);
				}
			} else if (i % 4 == 2) {
				if (msg == 2) {
					// 0*Q1, 1*Q1, 2*Q1, 3*Q1, 4*Q1, 5*Q1, 6*Q1, 7*Q1
					result[1].GetElementW(0)[i].ModAddEq(m_moduliQ[1], m_moduliQ[0]);
				}
				// ODD: Bootstrapping or Carry
			} 
			
			//We don't want msg == 0
			else if (i % 4 == 3) {
				if (msg == 3) {
					// 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1
					result[1].GetElementW(0)[i].ModAddEq(m_moduliQ[1],  m_moduliQ[0]);
				}
			}

		}
		result[1].SetFormat(Format::EVALUATION);	
		

		
		// In vals : {-1, 0 ,  1} -> (0, 1,  2). Outs are {-1, 0, 1}
		result[2] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		// Bit reverse
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign One more bit inserted
			msg = i >> (m_N_GSW_bit - 2);	// Get 2bit	
			if (msg == 0) {
				// 0*Q1, 1*Q1, 2*Q1, 3*Q1, 4*Q1, 5*Q1, 6*Q1, 7*Q1
				result[2].GetElementW(0)[i].ModSubEq(m_moduliQ[1], m_moduliQ[0]);
			} else if (msg == 2){
				result[2].GetElementW(0)[i].ModAddEq(m_moduliQ[1], m_moduliQ[0]);
			}
		}
		result[2].SetFormat(Format::EVALUATION);
		

		
		// Out Vector 3
		result[3] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		// 1~4bit reconstruction
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign One more bit inserted
			if (i % 2 == 0 ) {
				msg = i >> (m_N_GSW_bit - 3); // 3 bits

				// 0Q1 ~ 15Q1
				result[3].GetElementW(0)[i].ModAddEq(Q1.ModMul(msg,Q0), Q0);
			
			} else if( i % 2 == 1) {
				result[3].GetElementW(0)[i].ModAddEq(Q1, Q0);
			}
		}
		
		result[3].SetFormat(Format::EVALUATION);
		
		// Out Vector 3
		result[4] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		// 1~4bit reconstruction
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign One more bit inserted
			msg = i >> (m_N_GSW_bit - 3); // 3 bit
			// 0Q1 ~ 15Q1
			result[4].GetElementW(0)[i].ModAddEq(
				Q1.ModMul(scaling,Q0).ModMul(msg, Q0), Q0);
			result[4].GetElementW(1)[i].ModAddEq(
				Q1.ModMul(scaling,Q1).ModMul(msg, Q1), Q1);	
		}
		
		result[4].SetFormat(Format::EVALUATION);
		
		
		//////////////// Out Vector 5 /////////////////
		result[5] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		// Refresh 0, 1 -> 0, 1
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign One more bit inserted
			msg = i >> (m_N_GSW_bit - 2); // 4 bits
			result[5].GetElementW(0)[i].ModAddEq(m_moduliQ[1]*msg, m_moduliQ[0]);
		}
		result[5].SetFormat(Format::EVALUATION);
	
		
	}




	void MakeBootstrapReconPoly(vector<DCRTPoly> & result) const {	
		
		int32_t msg;	
		uint64_t scaling = GetRLWEParams()->GetScaling().ConvertToInt();
		scaling <<= GetExpBit();
		NativeInteger Q0 = m_moduliQ[0];
		NativeInteger Q1 = m_moduliQ[1];
		NativeInteger Q1Q1ModQ0 = Q1.ModMul(Q1,Q0);
		result.resize(4); 
		// Out Vector 0
		result[0] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		//NativeInteger Q1Q1  = m_moduliQ[1].ModMul(m_moduliQ[1], m_moduliQ[0]);
		// 3bit get, 2bit carry & 0,1,2 out, 3 is rev
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			msg = i >> (m_N_GSW_bit - 2); // 2bit		
			if (msg >= 1) {
					msg = 0;
				} else {
					msg = 1;
				};
				result[0].GetElementW(0)[i].ModAddEq(m_moduliQ[1] * msg, m_moduliQ[0]);
			//}
		}
		result[0].SetFormat(Format::EVALUATION);

		
		
		// Out Vector 1
		result[1] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		// 5 bit, Boot Orign
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign
			if (i % 2 == 0) {
				msg = i >> (m_N_GSW_bit - 2);	// 2bit	
				// 0*Q1, 1*Q1, 2*Q1, 3*Q1, 4*Q1, 5*Q1, 6*Q1, 7*Q1
				result[1].GetElementW(0)[i].ModAddEq(m_moduliQ[1] * msg, m_moduliQ[0]);
	
			} else if (i % 2 == 1){
				msg = i >> (m_N_GSW_bit - 2);	// 2bit	
				// 0*Q1, 1*Q1, 2*Q1, 3*Q1, 4*Q1, 5*Q1, 6*Q1, 7*Q1
				result[1].GetElementW(0)[i].ModAddEq(Q1Q1ModQ0.ModMul(msg, Q0), m_moduliQ[0]);
	
			}
		}
		result[1].SetFormat(Format::EVALUATION);	
		

		
		// Out Vector , To Exponent
		result[2] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		// Bit reverse
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign One more bit inserted
			msg = i >> (m_N_GSW_bit - 5);	// Get 5bit	
			result[2].GetElementW(0)[i].ModAddEq(
				Q1.ModMul(scaling,Q0).ModMul(msg, Q0), Q0);
			result[2].GetElementW(1)[i].ModAddEq(
				Q1.ModMul(scaling,Q1).ModMul(msg, Q1), Q1);
		}
		result[2].SetFormat(Format::EVALUATION);
		



		
		// Out Vector , To Exponent
		result[3] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		// Bit reverse
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign One more bit inserted
			msg = i >> (m_N_GSW_bit - 5);	// Get 5bit	
			
			if (msg == 0)  {
				msg = 1;
			} else {
				msg = 0;
			}
			result[3].GetElementW(0)[i].ModAddEq(
				Q1.ModMul(scaling,Q0).ModMul(msg, Q0), Q0);
			result[3].GetElementW(1)[i].ModAddEq(
			Q1.ModMul(scaling,Q1).ModMul(msg, Q1), Q1);
		}
		result[3].SetFormat(Format::EVALUATION);
	
		/*
		// Out Vector 2
		result[3] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		// 1~4bit reconstruction
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign One more bit inserted
			msg = i >> (m_N_GSW_bit - 4); // 4 bits

			// 0Q1 ~ 15Q1
			result[3].GetElementW(0)[i].ModAddEq(scaling_native.ModMul(msg, m_moduliQ[0]), m_moduliQ[0]);
			result[3].GetElementW(1)[i].ModAddEq(scaling_native.ModMul(msg, m_moduliQ[1]), m_moduliQ[1]);
		}
		
		result[3].SetFormat(Format::EVALUATION);
		

		
		// Out Vector 3 , 4 bit Q1Q1 vals
		result[4] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		
		// Reconstruct 4~7bit
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign One more bit inserted
			msg = i >> (m_N_GSW_bit - 4); // 4 bits
			msg <<= 4;
		
			// 0Q1 ~ 15Q1
			result[4].GetElementW(0)[i].ModAddEq(scaling_native.ModMul(msg, m_moduliQ[0]), m_moduliQ[0]);
			result[4].GetElementW(1)[i].ModAddEq(scaling_native.ModMul(msg, m_moduliQ[1]), m_moduliQ[1]);
		}
		result[4].SetFormat(Format::EVALUATION);
	
		//////////////// Out Vector 5 /////////////////
		result[5] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		// Make 8bit ~ 12bit
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign One more bit inserted
			msg = i >> (m_N_GSW_bit - 4); // 4 bits
			msg <<= 8;
			// 0Q1 ~ 15Q1
			result[5].GetElementW(0)[i].ModAddEq(scaling_native.ModMul(msg, m_moduliQ[0]), m_moduliQ[0]);
			result[5].GetElementW(1)[i].ModAddEq(scaling_native.ModMul(msg, m_moduliQ[1]), m_moduliQ[1]);
		}
		result[5].SetFormat(Format::EVALUATION);
	*/

	}


	void MakeBootstrapMulPoly(vector<DCRTPoly> & result) const {
		
		int32_t msg;
		//vector<DCRTPoly> result(6);
		result.resize(20); 
		// Out Vector 0
		result[0] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		NativeInteger Q1Q1  = m_moduliQ[1].ModMul(m_moduliQ[1], m_moduliQ[0]);
		
		// 2bit should be shift
		for (uint32_t i = 0; i < m_N_GSW; i++) {
		
			// Orign
			if (i % 4 == 0) {
				msg = i >> (m_N_GSW_bit - 3); // 3bit		
				// 0*Q1, 1*Q1, 2*Q1, 3*Q1, 0*Q1, 1*Q1, 2*Q1, 3*Q1
				if (msg >= 8) {
					PALISADE_THROW(config_error, "Coding is not correct");
				}
				msg = msg % 4;
				result[0].GetElementW(0)[i].ModAddEq(m_moduliQ[1] * msg, m_moduliQ[0]);
				// ODD: Bootstrapping or Carry
			} else if (i % 4 == 1) {

				msg = i >> (m_N_GSW_bit - 3); // 3bit
				// 0*Q1, 0*Q1, 0*Q1, 0*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1
				msg = msg >> 2;
				if (msg >= 2) {
					PALISADE_THROW(config_error, "Coding is not correct");
				}
				result[0].GetElementW(0)[i].ModAddEq(m_moduliQ[1] * msg, m_moduliQ[0]);
			
			} else if (i % 4 == 2) {
				// 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1
				result[0].GetElementW(0)[i].ModAddEq(m_moduliQ[1], m_moduliQ[0]);
	
			} else if (i % 4 == 3) {
				// 1*Q1Q1, 1*Q1Q1, 1*Q1Q1, 1*Q1Q1, 1*Q1Q1, 1*Q1Q1, 1*Q1Q1
				result[0].GetElementW(0)[i].ModAddEq(Q1Q1, m_moduliQ[0]);
			}
		}
		result[0].SetFormat(Format::EVALUATION);

		
		// Out Vector 1
		result[1] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);

		// 1 bit shift
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign
			if (i % 2 == 0) {
				msg = i >> (m_N_GSW_bit - 4);	// 4bit	
				// 0*Q1, 1*Q1, 2*Q1, 3*Q1, 4*Q1, 5*Q1, 6*Q1, 7*Q1
				result[1].GetElementW(0)[i].ModAddEq(m_moduliQ[1] * msg, m_moduliQ[0]);
		
				// ODD: Bootstrapping or Carry
			} else if (i % 2 == 1) {
				// Adjust Sign
				// 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1
				result[1].GetElementW(0)[i].ModAddEq(m_moduliQ[1], m_moduliQ[0]);
			}
		}
		result[1].SetFormat(Format::EVALUATION);


		// Out Vector 2
		result[2] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);

		// No Shift. 5 bits should be 0
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign One more bit inserted
			msg = i >> (m_N_GSW_bit - 4); // 4 bits
			if (msg >= 16) {
				PALISADE_THROW(config_error, "Coding is not correct");
			}
			//msg = (16 - msg)	
			// 0Q1 ~ 15Q1
			result[2].GetElementW(0)[i].ModAddEq(m_moduliQ[1] * msg, m_moduliQ[0]);
		
		}
		result[2].SetFormat(Format::EVALUATION);
		
		// Out Vector 3 , 4 bit Q1Q1 vals
		result[3] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		// No Shift. 5 bits should be 0
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign One more bit inserted
			msg = i >> (m_N_GSW_bit - 4); // 4 bits
			// 0Q1 ~ 15Q1
			result[3].GetElementW(0)[i].ModAddEq(Q1Q1 * msg, m_moduliQ[0]);
		}
		result[3].SetFormat(Format::EVALUATION);
	
		//////////////// Out Vector 4 /////////////////
		result[4] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		// No Shift. 5 bits should be 0
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign One more bit inserted
			msg = i >> (m_N_GSW_bit - 1); // One Bit
			if (msg >= 2) {
				PALISADE_THROW(config_error, "Coding is not correct");
			}
		
			// 0Q1 ~ 15Q1
			result[4].GetElementW(0)[i].ModAddEq( Q1Q1 * msg, m_moduliQ[0]);	
		}
		result[4].SetFormat(Format::EVALUATION);
	
		
		//////////////// Out Vector 5 /////////////////
		result[5] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		// Sign bit Handle, input is just 0 or 2, 2bit is consider and 3bit is empty
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign One more bit inserted
			msg = i >> (m_N_GSW_bit - 2); // Two Bit
			if (msg >= 4) {
				PALISADE_THROW(config_error, "Coding is not correct");
			}
			if (msg == 0) {
				result[5].GetElementW(0)[i].ModSubEq(m_moduliQ[1] , m_moduliQ[0]);
			} else if (msg == 2) {
				result[5].GetElementW(0)[i].ModAddEq(m_moduliQ[1] , m_moduliQ[0]);
			}
		}
		result[5].SetFormat(Format::EVALUATION);
		
		// Double Mul, 2bit orign
		result[6] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		// Sign bit Handle, input is just 0 or 2, 2bit is consider and 3bit is empty
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign One more bit inserted
			msg = i >> (m_N_GSW_bit - 2); // Two Bit
			result[6].GetElementW(0)[i].ModAddEq(m_moduliQ[1] * msg , m_moduliQ[0]);
		}
		result[6].SetFormat(Format::EVALUATION);
	
		// Double Mul, 2bit orign and Is zero
		result[7] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		// Sign bit Handle, input is just 0 or 2, 2bit is consider and 3bit is empty
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign One more bit inserted
			if (i % 2 == 0) {
				msg = i >> (m_N_GSW_bit - 2); // Two Bit
				result[7].GetElementW(0)[i].ModAddEq(m_moduliQ[1] * msg , m_moduliQ[0]);
			} else if (i % 2 == 1) {
				msg = i >> (m_N_GSW_bit - 2); // Two Bit
				if (msg > 0) {
					msg = 0 ;
				} else if (msg == 0) {
					msg = 1;
				}
				result[7].GetElementW(0)[i].ModAddEq(m_moduliQ[1] * msg, m_moduliQ[0]);
			}
		}
		result[7].SetFormat(Format::EVALUATION);
	

		// Double Mul,  2bit orign, 1bit Carry
		result[8] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		// Sign bit Handle, input is just 0 or 2, 2bit is consider and 3bit is empty
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign
			if (i % 2 == 0) {
				msg = i >> (m_N_GSW_bit - 3); // 3bit		
				// 0*Q1, 1*Q1, 2*Q1, 3*Q1, 0*Q1, 1*Q1, 2*Q1, 3*Q1
				if (msg >= 8) {
					PALISADE_THROW(config_error, "Coding is not correct");
				}
				msg = msg % 4;
				result[8].GetElementW(0)[i].ModAddEq(m_moduliQ[1] * msg, m_moduliQ[0]);
				// ODD: Bootstrapping or Carry
			} else if (i % 2 == 1) {

				msg = i >> (m_N_GSW_bit - 3); // 3bit
				// 0*Q1, 0*Q1, 0*Q1, 0*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1
				msg = msg >> 2;
				if (msg >= 2) {
					PALISADE_THROW(config_error, "Coding is not correct");
				}
				result[8].GetElementW(0)[i].ModAddEq(m_moduliQ[1] * msg, m_moduliQ[0]);
			
			} 
		}
		result[8].SetFormat(Format::EVALUATION);



		// Double Mul,  2bit orign, 1bit Carry, Zero
		result[9] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		// Sign bit Handle, input is just 0 or 2, 2bit is consider and 3bit is empty
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign
            msg = 0;
			if (i % 4 == 0) {
				msg = i >> (m_N_GSW_bit - 3); // 3bit		
				// 0*Q1, 1*Q1, 2*Q1, 3*Q1, 0*Q1, 1*Q1, 2*Q1, 3*Q1
				if (msg >= 8) {
					PALISADE_THROW(config_error, "Coding is not correct");
				}
				msg = msg % 4;
				result[9].GetElementW(0)[i].ModAddEq(m_moduliQ[1] * msg, m_moduliQ[0]);
				// ODD: Bootstrapping or Carry
			} else if (i % 4 == 1) {

				msg = i >> (m_N_GSW_bit - 3); // 3bit
				// 0*Q1, 0*Q1, 0*Q1, 0*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1
				msg = msg >> 2;
				if (msg >= 2) {
					PALISADE_THROW(config_error, "Coding is not correct");
				}
				result[9].GetElementW(0)[i].ModAddEq(m_moduliQ[1] * msg, m_moduliQ[0]);
			
			} else if (i % 4 == 2) {
				msg %= 4;
				if (msg > 0) {
					msg = 0;
				} else {
					msg = 1;
				}
				// 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1
				result[9].GetElementW(0)[i].ModAddEq(m_moduliQ[1] * msg, m_moduliQ[0]);
	
			}		
		}
		result[9].SetFormat(Format::EVALUATION);
	}


	void MakeBootstrapAddPoly(vector<DCRTPoly> &result) const {
		
		int32_t msg;	
		//vector<DCRTPoly> result(4);
		result.resize(6);


		//////////////////////////// Out Vector 0
		result[0] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		NativeInteger Q1Q1  = m_moduliQ[1].ModMul(m_moduliQ[1], m_moduliQ[0]);
		
		// 1bit should be shift
		for (uint32_t i = 0; i < m_N_GSW; i++) {
		
			// Orign 
			if (i % 2 == 0) {
				msg = i >> (m_N_GSW_bit - 4); // 3bit
				// It is Orign values . 0 1 2 3 0 1 2 3 ...
				msg  = msg % 4;
				result[0].GetElementW(0)[i].ModAddEq(m_moduliQ[1] * msg , m_moduliQ[0]);
				// ODD: Bootstrapping or Carry
			} else if (i % 2 == 1) {

				msg = i >> (m_N_GSW_bit - 4); // 3bit
				// -2 -2 -2 -2 -1 -1 -1 -1 0 0 0 0 1 1 1 1 1
				msg = msg >> 2;
				if (msg == 0) {
					result[0].GetElementW(0)[i].ModSubEq(Q1Q1   * 2 , m_moduliQ[0]);
				} else if (msg == 1) {
					result[0].GetElementW(0)[i].ModSubEq(Q1Q1   , m_moduliQ[0]);
				} else if (msg == 2) {
				
				} else if (msg == 3) {
					result[0].GetElementW(0)[i].ModAddEq(Q1Q1		, m_moduliQ[0]);
				} 
			} 		
		}
		result[0].SetFormat(Format::EVALUATION);


		////////////////////////// Out Vector 1, Orign and Sign
		result[1] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);

		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign 
			if (i % 2 == 0) {
				msg = i >> (m_N_GSW_bit - 4); // 4bit
				// It is Orign values . 0 1 2 3 0 1 2 3 ...
				msg  = msg % 4;
				result[1].GetElementW(0)[i].ModAddEq(m_moduliQ[1] * msg , m_moduliQ[0]);
				// ODD: Bootstrapping or Carry
			} else if (i % 2 == 1) {

				msg = i >> (m_N_GSW_bit - 4); // 4bit
				// -1 -1 -1 -1 -1 -1 -1 -1 1 1 1 1 1 1 1 1
				//0=<msg <8 is negative, 8<=msg <15 is positive;
				if (msg < 8) {
					result[1].GetElementW(0)[i].ModSubEq(m_moduliQ[1] , m_moduliQ[0]);
				//} else if (msg == 1) {
				} else {
					result[1].GetElementW(0)[i].ModAddEq(m_moduliQ[1] , m_moduliQ[0]);
				}	
			} 		
		}
		result[1].SetFormat(Format::EVALUATION);


		////////////////////////// Out Vector 2, Orign and zero, Not used
		result[2] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);

		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign 
			if (i % 2 == 0) {
				msg = i >> (m_N_GSW_bit - 4); // 4bit
				// It is Orign values . 0 1 2 3 0 1 2 3 ...
				msg  = msg % 4;
				result[2].GetElementW(0)[i].ModAddEq(m_moduliQ[1] * msg , m_moduliQ[0]);
				// ODD: Bootstrapping or Carry
			} else if (i % 2 == 1) {

				msg = i >> (m_N_GSW_bit - 4); // 4bit
				// 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1
				msg = msg >> 3;
				if (msg == 1) {
					result[2].GetElementW(0)[i].ModAddEq(m_moduliQ[1]   , m_moduliQ[0]);
				}	
			} 		
		}
		result[2].SetFormat(Format::EVALUATION);



		// Out Vector 3 Bit Reverse and Carry 1 Not used
		result[3] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);

		// Out Sign
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign 
			if (i % 2 == 0) {
				msg = i >> (m_N_GSW_bit - 3); // 3bit
				// It is Orign values . 0 1 2 3 0 1 2 3 ...
				msg  = msg % 4;			
				result[3].GetElementW(0)[i].ModAddEq(m_moduliQ[1] * msg , m_moduliQ[0]);	
				// ODD: Bootstrapping or Carry
			} else if (i % 2 == 1) {

				msg = i >> (m_N_GSW_bit - 3); // 4bit
				// -1 -1 -1 -1 0 0 0 0
				msg = msg >> 2;
				if (msg == 0) {
					result[3].GetElementW(0)[i].ModSubEq(Q1Q1 , m_moduliQ[0]);	
				} else if (msg == 1) {

				}	
			}
		}
		result[3].SetFormat(Format::EVALUATION);
		

		// Out Vector 3 Bit Reverse and Carry 1
		result[4] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);


		// Out msg, carry, Is Zero
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign 
			if (i % 4 == 0) {
				msg = i >> (m_N_GSW_bit - 3); // 3bit
				// It is Orign values . 0 1 2 3 0 1 2 3 ...
				msg  = msg % 4;			
				result[4].GetElementW(0)[i].ModAddEq(m_moduliQ[1] * msg , m_moduliQ[0]);	
				// ODD: Bootstrapping or Carry
			} else if (i % 4 == 1) {

				msg = i >> (m_N_GSW_bit - 3); // 3bit
				// -1 -1 -1 -1 0 0 0 0
				msg = msg >> 2;
				if (msg == 0) {
					result[4].GetElementW(0)[i].ModSubEq(Q1Q1 , m_moduliQ[0]);	
				} else if (msg == 1) {

				}	
			} else if (i % 4 == 2) {
				msg = i >> (m_N_GSW_bit - 3); // 3bit
				msg %= 4;
				if (msg > 0 ) {
					msg = 0;
				} else {
					msg = 1;
				}
				result[4].GetElementW(0)[i].ModAddEq(m_moduliQ[1] * msg, m_moduliQ[0]);	
			}
		}
		result[4].SetFormat(Format::EVALUATION);

	}



		//vector<DCRTPoly>  MakeBootstrapMulPolyOne() {
	void MakeBootstrapMulPolyOne(vector<DCRTPoly> &result) const {
		
		int32_t msg;	
		//vector<DCRTPoly> result(6);
		result.resize(6);

		// Out Vector 0
		result[0] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		NativeInteger Q1Q1  = m_moduliQ[1].ModMul(m_moduliQ[1], m_moduliQ[0]);
		
		// 2bit should be shift 3bit infers it
		for (uint32_t i = 0; i < m_N_GSW; i++) {
		
			// Orign
			if (i % 4 == 0) {
				msg = i >> (m_N_GSW_bit - 2); // 2bit		
				// 0*Q1, 1*Q1, 2*Q1, 3*Q1, 0*Q1, 1*Q1, 2*Q1, 3*Q1
				if (msg >= 4) {
					PALISADE_THROW(config_error, "Coding is not correct");
				}
				msg = msg % 2;
				result[0].GetElementW(0)[i].ModAddEq(m_moduliQ[1] * msg, m_moduliQ[0]);
				// ODD: Bootstrapping or Carry
			} else if (i % 4 == 1) {

				msg = i >> (m_N_GSW_bit - 2); // 2bit
				// 0*Q1, 0*Q1, 0*Q1, 0*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1
				msg = msg >> 1;
				if (msg >= 2) {
					PALISADE_THROW(config_error, "Coding is not correct");
				}
				result[0].GetElementW(0)[i].ModAddEq(m_moduliQ[1] * msg, m_moduliQ[0]);
			
			} else if (i % 4 == 2) {
				// 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1
				result[0].GetElementW(0)[i].ModAddEq(m_moduliQ[1], m_moduliQ[0]);
	
			} else if (i % 4 == 3) {
				// 1*Q1Q1, 1*Q1Q1, 1*Q1Q1, 1*Q1Q1, 1*Q1Q1, 1*Q1Q1, 1*Q1Q1
				result[0].GetElementW(0)[i].ModAddEq(Q1Q1, m_moduliQ[0]);
			}
		}
		result[0].SetFormat(Format::EVALUATION);

		// Out Vector 1a
		result[1] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		
		// ...
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign
			if (i % 2 == 0) {
				msg = i >> (m_N_GSW_bit - 2);	// 2bit
				msg = msg % 2;
				// 0*Q1, 1*Q1, 2*Q1, 3*Q1, 4*Q1, 5*Q1, 6*Q1, 7*Q1
				result[1].GetElementW(0)[i].ModAddEq(m_moduliQ[1] * msg, m_moduliQ[0]);
		
				// ODD: Bootstrapping or Carry
			} else if (i % 2 == 1) {
				msg = i >> (m_N_GSW_bit - 2);	// 2bit
				msg = msg >> 1;
				// 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1, 1*Q1
				result[1].GetElementW(0)[i].ModAddEq(Q1Q1 * msg, m_moduliQ[0]);
			}
		}
		result[1].SetFormat(Format::EVALUATION);


		// Out Vector 2
		result[2] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);

		// No Shift. 5 bits should be 0
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign One more bit inserted
			msg = i >> (m_N_GSW_bit - 4); // 4 bits
			if (msg >= 16) {
				PALISADE_THROW(config_error, "Coding is not correct");
			}
			//msg = (16 - msg)	
			// 0Q1 ~ 15Q1
			result[2].GetElementW(0)[i].ModAddEq(m_moduliQ[1] * msg, m_moduliQ[0]);
		
		}
		result[2].SetFormat(Format::EVALUATION);
		
		// Out Vector 3 , 4 bit Q1Q1 vals
		result[3] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		// No Shift. 5 bits should be 0
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign One more bit inserted
			msg = i >> (m_N_GSW_bit - 4); // 4 bits
			// 0Q1 ~ 15Q1
			result[3].GetElementW(0)[i].ModAddEq(Q1Q1 * msg, m_moduliQ[0]);
		}
		result[3].SetFormat(Format::EVALUATION);
	
		//////////////// Out Vector 4 /////////////////
		result[4] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		// No Shift. 5 bits should be 0
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign One more bit inserted
			msg = i >> (m_N_GSW_bit - 1); // One Bit
			if (msg >= 2) {
				PALISADE_THROW(config_error, "Coding is not correct");
			}
		
			// 0Q1 ~ 15Q1
			result[4].GetElementW(0)[i].ModAddEq( Q1Q1 * msg, m_moduliQ[0]);	
		}
		result[4].SetFormat(Format::EVALUATION);
	
		
		//////////////// Out Vector 5 /////////////////
		result[5] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		// Sign bit Handle, input is just 0 or 2, 2bit is consider and 3bit is empty
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign One more bit inserted
			msg = i >> (m_N_GSW_bit - 2); // Two Bit
			if (msg >= 4) {
				PALISADE_THROW(config_error, "Coding is not correct");
			}
			
			if (msg == 0) {
				result[5].GetElementW(0)[i].ModSubEq(m_moduliQ[1] , m_moduliQ[0]);
		
			} else if (msg == 2) {
				result[5].GetElementW(0)[i].ModAddEq(m_moduliQ[1] , m_moduliQ[0]);
	
			}
		
		}
		result[5].SetFormat(Format::EVALUATION);
		//return result;
	}


	//	vector<DCRTPoly>  MakeBootstrapAddPolyOne() {
	void MakeBootstrapAddPolyOne(vector<DCRTPoly> &result) const {
		
		int32_t msg;	
		//vector<DCRTPoly> result(4);
		result.resize(4);


		//////////////////////////// Out Vector 0
		result[0] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);
		NativeInteger Q1Q1  = m_moduliQ[1].ModMul(m_moduliQ[1], m_moduliQ[0]);
		
		// 1bit should be shift
		for (uint32_t i = 0; i < m_N_GSW; i++) {
		
			// Orign 
			if (i % 2 == 0) {
				msg = i >> (m_N_GSW_bit - 4); // 3bit
				// It is Orign values . 0 1 2 3 0 1 2 3 ...
				msg  = msg % 4;
				result[0].GetElementW(0)[i].ModAddEq(m_moduliQ[1] * msg , m_moduliQ[0]);
				// ODD: Bootstrapping or Carry
			} else if (i % 2 == 1) {

				msg = i >> (m_N_GSW_bit - 4); // 3bit
				// -2 -2 -2 -2 -1 -1 -1 -1 0 0 0 0 1 1 1 1 1
				msg = msg >> 2;
				if (msg == 0) {
					result[0].GetElementW(0)[i].ModSubEq(Q1Q1   * 2 , m_moduliQ[0]);
				} else if (msg == 1) {
					result[0].GetElementW(0)[i].ModSubEq(Q1Q1   , m_moduliQ[0]);
				} else if (msg == 2) {
				
				} else if (msg == 3) {
					result[0].GetElementW(0)[i].ModAddEq(Q1Q1		, m_moduliQ[0]);
				} 
			} 		
		}
		result[0].SetFormat(Format::EVALUATION);


		////////////////////////// Out Vector 1, Orign and Sign
		result[1] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);

		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign 
			if (i % 2 == 0) {
				msg = i >> (m_N_GSW_bit - 4); // 4bit
				// It is Orign values . 0 1 2 3 0 1 2 3 ...
				msg  = msg % 4;
				result[1].GetElementW(0)[i].ModAddEq(m_moduliQ[1] * msg , m_moduliQ[0]);
				// ODD: Bootstrapping or Carry
			} else if (i % 2 == 1) {

				msg = i >> (m_N_GSW_bit - 4); // 4bit
				// -1 -1 -1 -1 -1 -1 -1 -1 1 1 1 1 1 1 1 1
				msg = msg >> 3;
				if (msg == 0) {
					result[1].GetElementW(0)[i].ModSubEq(m_moduliQ[1] , m_moduliQ[0]);
				} else if (msg == 1) {
					result[1].GetElementW(0)[i].ModAddEq(m_moduliQ[1] , m_moduliQ[0]);
				}	
			} 		
		}
		result[1].SetFormat(Format::EVALUATION);


		////////////////////////// Out Vector 2, Orign and zero & 1
		result[2] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);

		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign 
			if (i % 2 == 0) {
				msg = i >> (m_N_GSW_bit - 4); // 4bit
				// It is Orign values . 0 1 2 3 0 1 2 3 ...
				msg  = msg % 4;
				result[2].GetElementW(0)[i].ModAddEq(m_moduliQ[1] * msg , m_moduliQ[0]);
				// ODD: Bootstrapping or Carry
			} else if (i % 2 == 1) {

				msg = i >> (m_N_GSW_bit - 4); // 4bit
				// 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1
				msg = msg >> 3;
				if (msg == 1) {
					result[2].GetElementW(0)[i].ModAddEq(m_moduliQ[1]   , m_moduliQ[0]);
				}	
			} 		
		}
		result[2].SetFormat(Format::EVALUATION);



		// Out Vector 3 Bit Reverse and Carry 1
		result[3] = DCRTPoly(m_DCRTGSWParams,Format::COEFFICIENT,true);

		// Out Sign
		for (uint32_t i = 0; i < m_N_GSW; i++) {
			// Orign 
			if (i % 2 == 0) {
				msg = i >> (m_N_GSW_bit - 3); // 3bit
				// It is Orign values . 0 1 2 3 0 1 2 3 ...
				msg  = msg % 4;			
				result[3].GetElementW(0)[i].ModAddEq(m_moduliQ[1] * msg , m_moduliQ[0]);	
				// ODD: Bootstrapping or Carry
			} else if (i % 2 == 1) {

				msg = i >> (m_N_GSW_bit - 3); // 4bit
				// -1 -1 -1 -1 0 0 0 0
				msg = msg >> 2;
				if (msg == 0) {
					result[3].GetElementW(0)[i].ModSubEq(Q1Q1 , m_moduliQ[0]);	
				} else if (msg == 1) {

				}	
			}
		}
		result[3].SetFormat(Format::EVALUATION);
		//return result;
	}



	//vector<DCRTPoly> MakeRotMonomial(uint32_t Ns,
	void MakeRotMonomial(uint32_t Ns,
		const shared_ptr<ILDCRTParams<BigInteger>> DCRTParams,
		vector<DCRTPoly> &rot_monomials ) const {			
		//vector<DCRTPoly> rot_monomials(Ns * 2);
		rot_monomials.resize(Ns * 2);
		// loop for positive values of m		
  	
		//#pragma omp parallel for 
		{for (uint32_t i = 0; i < Ns; i++) {
			//rot_monomials[i] = DCRTPoly(DCRTParams, Format::COEFFICIENT, true);
			//Local Varlabie
			NativeInteger Q0_loc = m_moduliQ[0]; 
      DCRTPoly aPoly = DCRTPoly(DCRTParams, Format::COEFFICIENT, true);
      NativePoly tmp0 = NativePoly(DCRTParams->GetParams()[0], Format::COEFFICIENT, true);
			

			tmp0[i] = tmp0[i].ModAddEq(NativeInteger(1), Q0_loc);  // X^m
			aPoly = tmp0;
			aPoly.SetFormat(Format::EVALUATION);
			rot_monomials[i]	=	std::move(aPoly);
    }}

    // loop for negative values of m
		//#pragma omp parallel for 
		{for (uint32_t i = 0; i < Ns; i++) {
      //rot_monomials[i+Ns] = DCRTPoly(DCRTParams, Format::COEFFICIENT, true);
			NativeInteger Q0_loc = m_moduliQ[0]; 
			DCRTPoly aPoly = DCRTPoly(DCRTParams, Format::COEFFICIENT, true);
      NativePoly tmp0 = NativePoly(DCRTParams->GetParams()[0], Format::COEFFICIENT, true);
			
			tmp0[i] = tmp0[i].ModSubEq(NativeInteger(1), Q0_loc);  // X^-m
			aPoly = tmp0;
			aPoly.SetFormat(Format::EVALUATION);
			rot_monomials[i+Ns]	=	std::move(aPoly);
    }}
		
		//return rot_monomials;

	}

	//Get Parameters 
  //const std::vector<NativeInteger>& GetDigitsR() const { return m_digitsR; }
	const std::shared_ptr<FPLWECryptoParams> GetLWEParams() const {
		return m_LWEParams;
	}
	const std::shared_ptr<FPRLWECryptoParams> GetRLWEParams() const {
		return m_RLWEParams;
	}
	const shared_ptr<ILDCRTParams<BigInteger>> GetGSWPolyParams() const {
    return m_DCRTGSWParams;
  }
	const shared_ptr<ILDCRTParams<BigInteger>> GetGSWNKPolyParams() const {
    return m_DCRTGSWNKParams;
  }
	const uint32_t GetGSWN() const{ return m_N_GSW;}
	const uint32_t GetGSWN_Bit() const { return m_N_GSW_bit; }
	const uint32_t GetGSWM() const{ return m_M_GSW;}
	const uint32_t GetGSWM_Bit() const { return m_M_GSW_bit; }
	const uint32_t GetGSWK() const { return m_K_GSW; }
	const uint32_t GetGSWNK() const { return m_K_GSW * m_N_GSW; }
	
	const DiscreteGaussianGeneratorImpl<NativeVector> &GetGSWDgg() const {
		return m_gsw_dgg;
  }
  const DiscreteGaussianGeneratorImpl<NativeVector> &GetKSDgg() const {
		return m_ks_dgg;
  }
	const DiscreteGaussianGeneratorImpl<NativeVector> &GetEVDgg() const {
		return m_ev_dgg;
  }
 	const DiscreteGaussianGeneratorImpl<NativeVector> &GetPKSDgg() const {
		return m_pks_dgg;
  }
 

	uint32_t GetBaseKSBit() const { return m_KS_base_bit; }
 	uint32_t GetLenKS() const { return m_KS_l; }
	uint32_t GetLenKSsave() const { return m_KS_l_save; }
	uint32_t GetLenKSrm() const { return m_KS_l_rm; }
 
	uint32_t GetBaseEVBit() const { return m_EV_base_bit; }
 	uint32_t GetLenEV() const { return m_EV_l; }
	uint32_t GetLenEVsave() const { return m_EV_l_save; }
	uint32_t GetLenEVrm() const { return m_EV_l_rm; }
 
	uint32_t GetBaseGDBit() const { return m_GD_base_bit; }
 	uint32_t GetLenGD() const { return m_GD_l; }
	uint32_t GetLenGDsave() const { return m_GD_l_save; }
	uint32_t GetLenGDrm() const { return m_GD_l_rm; }
 

	uint32_t GetBaseEncodeBit() const { return m_encode_base_bit; }
 	uint32_t GetBaseEncode() const { return m_encode_base; }
 	
	uint32_t GetBaseCarryBit() const { return m_C_base_bit; }
 	uint32_t GetBaseCarry() const { return m_C_base; }
 	
	vector<NativeInteger> GetModuliQ() const { return m_moduliQ;} 
  vector<NativeInteger> GetMuGSW() const { return m_mu_GSW;} 
  
	uint32_t GetQBit() const { return m_Q_bit; }
  uint32_t GetQ2Bit() const { return m_Q_bit2; }
	vector<uint64_t> GetMulInfo(uint32_t i) const {return m_carry_mul_infos[i];}
	vector<uint64_t> GetMulFloatInfo(uint32_t i) const {return m_carry_mul_float_infos[i];}
	vector<uint32_t> GetLenEnc() const { return m_encode_len;}
 
	const DCRTPoly& GetMonomial(uint32_t i) const { return m_monomials[i]; }
  const DCRTPoly& GetRotMonomial_GSW(uint32_t i) const { return m_rot_monomials_GSW[i]; }
  const DCRTPoly& GetRotMonomial_RLWE(uint32_t i) const { return m_rot_monomials_RLWE[i]; }


	// Bootstrap
	const DCRTPoly& GetACCMul(uint32_t i) const {return  m_prod_poly_bootstrap[i];}
	const DCRTPoly& GetACCMulOne(uint32_t i) const {return  m_prod_poly_bootstrap_one[i];}
	const DCRTPoly& GetACCAdd(uint32_t i) const {return  m_add_poly_bootstrap[i];}
	const DCRTPoly& GetACCAddOne(uint32_t i) const {return  m_add_poly_bootstrap_one[i];}
	const DCRTPoly& GetACCOrign() const { return m_orign_bootstrap;}
	const DCRTPoly& GetACCRelu(uint32_t i) const { return m_relu_poly_bootstrap[i];}
	const DCRTPoly& GetACCDiv(uint32_t i) const { return m_div_poly_bootstrap[i];}
	const DCRTPoly& GetACCReconPoly(uint32_t i) const { return m_recon_poly_bootstrap[i];}
	const DCRTPoly& GetPrefixBaseQ1Q1over2() const {return m_Q1Q1_base_over_2;}
  const DCRTPoly& GetPrefixBaseQ1over2() const {return m_Q1_base_over_2;}
  const bool &GetPARALLEL() const {return m_PARALLEL;}
	const uint32_t &GetExpBit() const {return m_exp_bit;}
	FPFHEMETHOD GetMethod() const { return m_method; }

	uint32_t GetH_FP() const { return m_h_FP;}
	uint32_t GetH_GSW() const { return m_h_GSW;}
	uint32_t GetH_PKS() const { return m_h_PKS;}

  //bool operator==(const FPRingGSWCryptoParams& other) const {
  //  return 
			//*m_LWEParams == *other.m_LWEParams && m_baseR == other.m_baseR &&
  //         m_GD_base_bit == other.m_baseGBit && m_method == other.m_method;
  //}
  //bool operator!=(const FPRingGSWCryptoParams& other) const {
  //  return !(*this == other);
  //}

  template <class Archive>
  void save(Archive& ar, std::uint32_t const version) const {
    //ar(::cereal::make_nvp("params", m_LWEParams));
    //ar(::cereal::make_nvp("bR", m_baseR));
    ar(::cereal::make_nvp("bG", m_GD_base_bit));
    ar(::cereal::make_nvp("method", m_method));
  }

  template <class Archive>
  void load(Archive& ar, std::uint32_t const version) {
    if (version > SerializedVersion()) {
      PALISADE_THROW(deserialize_error,
                     "serialized object version " + std::to_string(version) +
                         " is from a later version of the library");
    }
    //ar(::cereal::make_nvp("params", m_LWEParams));
    //ar(::cereal::make_nvp("bR", m_baseR));
    ar(::cereal::make_nvp("bG", m_GD_base_bit));
    ar(::cereal::make_nvp("method", m_method));

    this->PreCompute();
  }

  std::string SerializedObjectName() const { return "FPRingGSWCryptoParams"; }
  static uint32_t SerializedVersion() { return 1; }

 private:
	std::shared_ptr<FPLWECryptoParams> m_LWEParams;
  
	// Ring_params
	std::shared_ptr<FPRLWECryptoParams> m_RLWEParams;
  // Parameters for polynomials in RingGSW/RingLWE
  shared_ptr<ILDCRTParams<BigInteger>> m_DCRTGSWParams;
	shared_ptr<ILDCRTParams<BigInteger>> m_DCRTGSWNKParams;
	

	// GSW_N
  uint32_t m_N_GSW;
	uint32_t m_N_GSW_bit;
	uint32_t m_K_GSW;
	uint32_t m_M_GSW;
	uint32_t m_M_GSW_bit;
	uint32_t m_n_PKS;
	uint32_t m_n_PKS_bit;

	DiscreteGaussianGeneratorImpl<NativeVector> m_gsw_dgg;
	DiscreteGaussianGeneratorImpl<NativeVector> m_ks_dgg;
	DiscreteGaussianGeneratorImpl<NativeVector> m_ev_dgg;
	DiscreteGaussianGeneratorImpl<NativeVector> m_pks_dgg;


	uint32_t m_KS_base_bit;
	uint32_t m_KS_l;
	uint32_t m_KS_l_rm;
	uint32_t m_KS_l_save;
	
	uint32_t m_EV_base_bit;
	uint32_t m_EV_l;
	uint32_t m_EV_l_rm;
	uint32_t m_EV_l_save;

	// gadget base used in bootstrapping
	uint32_t m_GD_base_bit;
	uint32_t m_GD_l;
	uint32_t m_GD_l_rm;
	uint32_t m_GD_l_save;

	// gadget base used in bootstrapping
	uint32_t m_PKS_base_bit;
	uint32_t m_PKS_l;
	uint32_t m_PKS_l_rm;
	uint32_t m_PKS_l_save;

	// MSG 
	uint32_t m_encode_base_bit;
	uint32_t m_encode_base;
	// CarryUP Bit
	uint32_t m_C_base_bit;
	uint32_t m_C_base;

	vector<NativeInteger> m_moduliQ; 

  // number of digits in decomposing integers mod Q
  uint32_t m_Q_bit;
  // twice the number of digits in decomposing integers mod Q
  uint32_t m_Q_bit2;
  
	// Encoding len;
	vector<uint32_t> m_encode_len;
	
	// Carry mul Info
	vector<vector<uint64_t>> m_carry_mul_infos;
	// Carry float Info
	vector<vector<uint64_t>> m_carry_mul_float_infos;
	uint32_t m_carry_max;
	// Q over 2
	vector<NativeInteger> m_Q1Q1_prefix_0;
	vector<NativeInteger> m_Q1Q1_prefix_1;
	vector<NativeInteger> m_Q1_prefix_0;
	vector<NativeInteger> m_Q1_prefix_1;
	DCRTPoly m_Q1Q1_base_over_2;
	DCRTPoly m_Q1_base_over_2;



	// Boostrapp Polys
	// Mul ACC
	vector<DCRTPoly> m_prod_poly_bootstrap;
	vector<DCRTPoly> m_prod_poly_bootstrap_one;
	vector<DCRTPoly> m_add_poly_bootstrap;
	vector<DCRTPoly> m_add_poly_bootstrap_one;
	vector<DCRTPoly> m_relu_poly_bootstrap;
	vector<DCRTPoly> m_div_poly_bootstrap;
	vector<DCRTPoly> m_recon_poly_bootstrap;
	// Orign recover ACC
	DCRTPoly m_orign_bootstrap;

  
	// A vector of powers of baseG
  std::vector<NativeInteger> m_Gpower;


  // Precomputed polynomials in Format::EVALUATION representation for X^m - 1
  // (used only for GINX bootstrapping)
  std::vector<DCRTPoly> m_monomials;

	// Precomputed X^m
	std::vector<DCRTPoly> m_rot_monomials_GSW;
	std::vector<DCRTPoly> m_rot_monomials_RLWE;

  // Bootstrapping method (AP or GINX)
  FPFHEMETHOD m_method;


	uint32_t m_h_FP;
	uint32_t m_h_GSW;
	uint32_t m_h_PKS;
	uint32_t m_exp_bit;
	bool m_PARALLEL;
	vector<NativeInteger> m_mu_GSW;

};

/**
 * @brief Class that stores a RingGSW ciphertext; a two-dimensional vector of
 * ring elements
 */

class FPRingGSWCiphertext : public Serializable {
 public:
  FPRingGSWCiphertext() {}

  FPRingGSWCiphertext(uint32_t gadget_dim, uint32_t rowSize, uint32_t colSize) {
  	m_elements.resize(gadget_dim);
		for (uint32_t j=0; j < gadget_dim; j++) {
			m_elements[j].resize(rowSize); 
			for (uint32_t i = 0; i < rowSize; i++) {
				m_elements[j][i].resize(colSize);
			}
		}
  }
  explicit FPRingGSWCiphertext(
      const std::vector<std::vector<std::vector<DCRTPoly>>>& elements)
      : m_elements(elements) {}

  explicit FPRingGSWCiphertext(const FPRingGSWCiphertext& rhs) {
    this->m_elements = rhs.m_elements;
  }

  explicit FPRingGSWCiphertext(const FPRingGSWCiphertext&& rhs) {
    this->m_elements = std::move(rhs.m_elements);
  }

  const FPRingGSWCiphertext& operator=(const FPRingGSWCiphertext& rhs) {
    this->m_elements = rhs.m_elements;
    return *this;
  }

  const FPRingGSWCiphertext& operator=(const FPRingGSWCiphertext&& rhs) {
    this->m_elements = rhs.m_elements;
    return *this;
  }

  const std::vector<std::vector<std::vector<DCRTPoly>>>& GetElements() const {
    return m_elements;
  }

	
  const DCRTPoly& GetElementAtIndex(uint32_t i, uint32_t j, uint32_t k)
      const { return (m_elements[i][j][k]);}

  const vector<DCRTPoly>& GetElementAtIndex(uint32_t i, uint32_t j)
      const { return (m_elements[i][j]);}



  void SetElements(const std::vector<std::vector<std::vector<DCRTPoly>>>& elements) {
    m_elements = elements;
  }

  /**
   * Switches between COEFFICIENT and Format::EVALUATION polynomial
   * representations using NTT
   */
  void SetFormat(const Format format) {
    for (uint32_t i = 0; i < m_elements.size(); i++) {
      // column size is assume to be the same
      for (uint32_t j = 0; j < m_elements[0].size(); j++) {
				for (uint32_t k = 0; k < m_elements[0][0].size(); k++) {
					m_elements[i][j][k].SetFormat(format);
				}
			}
		}
  }
		
	std::vector<std::vector<DCRTPoly>>& operator[](uint32_t i) { return m_elements[i]; }

  const std::vector<std::vector<DCRTPoly>>& operator[](usint i) const {
    return m_elements[i];
  }

  bool operator==(const FPRingGSWCiphertext& other) const {
    return m_elements == other.m_elements;
  }

  bool operator!=(const FPRingGSWCiphertext& other) const {
    return !(*this == other);
  }

  template <class Archive>
  void save(Archive& ar, std::uint32_t const version) const {
    ar(::cereal::make_nvp("elements", m_elements));
  }

  template <class Archive>
  void load(Archive& ar, std::uint32_t const version) {
    if (version > SerializedVersion()) {
      PALISADE_THROW(deserialize_error,
                     "serialized object version " + std::to_string(version) +
                         " is from a later version of the library");
    }
    ar(::cereal::make_nvp("elements", m_elements));
  }

  std::string SerializedObjectName() const { return "FPRingGSWCiphertext"; }
  static uint32_t SerializedVersion() { return 1; }

 private:
	std::vector<std::vector<std::vector<DCRTPoly>>> m_elements;
};

/**
 * @brief Class that stores the refreshing key (used in bootstrapping)
 * A three-dimensional vector of RingGSW ciphertexts
 */






class FPRingGSWBTKey : public Serializable {
 public:
  FPRingGSWBTKey() {}

	/*dim1:  dim2: module dim e.g., a,b is 2, dim3: degree of N */
  explicit FPRingGSWBTKey(uint32_t dim1, uint32_t dim2, uint32_t dim3) {
    m_key.resize(dim1);
    for (uint32_t i = 0; i < dim1; i++) {
      m_key[i].resize(dim2);
      for (uint32_t j = 0; j < dim2; j++) {
				m_key[i][j].resize(dim3);
			}
    }
  }

  explicit FPRingGSWBTKey(
      const std::vector<std::vector<std::vector<FPRingGSWCiphertext>>>& key)
      : m_key(key) {}

  explicit FPRingGSWBTKey(const FPRingGSWBTKey& rhs) { this->m_key = rhs.m_key; }

  explicit FPRingGSWBTKey(const FPRingGSWBTKey&& rhs) {
    this->m_key = std::move(rhs.m_key);
  }

  const FPRingGSWBTKey& operator=(const FPRingGSWBTKey& rhs) {
    this->m_key = rhs.m_key;
    return *this;
  }

  const FPRingGSWBTKey& operator=(const FPRingGSWBTKey&& rhs) {
    this->m_key = std::move(rhs.m_key);
    return *this;
  }

  const std::vector<std::vector<std::vector<FPRingGSWCiphertext>>>& GetElements()
      const {
    return m_key;
  }



  void SetElements(
      const std::vector<std::vector<std::vector<FPRingGSWCiphertext>>>& key) {
    m_key = key;
  }

  std::vector<std::vector<FPRingGSWCiphertext>>& operator[](uint32_t i) {
    return m_key[i];
  }

  const std::vector<std::vector<FPRingGSWCiphertext>>& operator[](usint i) const {
    return m_key[i];
  }

  bool operator==(const FPRingGSWBTKey& other) const {
    return m_key == other.m_key;
  }

  bool operator!=(const FPRingGSWBTKey& other) const { return !(*this == other); }

  template <class Archive>
  void save(Archive& ar, std::uint32_t const version) const {
    ar(::cereal::make_nvp("key", m_key));
  }

  template <class Archive>
  void load(Archive& ar, std::uint32_t const version) {
    if (version > SerializedVersion()) {
      PALISADE_THROW(deserialize_error,
                     "serialized object version " + std::to_string(version) +
                         " is from a later version of the library");
    }
    ar(::cereal::make_nvp("key", m_key));
  }

  std::string SerializedObjectName() const { return "FPRingGSWBTKey"; }
  static uint32_t SerializedVersion() { return 1; }

 private:
  std::vector<std::vector<std::vector<FPRingGSWCiphertext>>> m_key;
};





// The struct for storing bootstrapping keys
typedef struct {
  // refreshing key
  std::shared_ptr<FPRingGSWBTKey> BSkey;
  // switching key
  std::shared_ptr<FPLWESwitchingKey> PKSkey;

	// Mul key
	std::shared_ptr<FPRingEvaluationKey> Evkey;
	
	//////////////////////////////////
	// Ring Switching Key
	std::shared_ptr<FPRLWESwitchingKey> RKSkey;

	//////// TEST Values
	std::shared_ptr<const FPRLWEPrivateKeyImpl> BK_SK;

} FPRingGSWEvalKey;




}  // namespace lbcrypto

#endif
