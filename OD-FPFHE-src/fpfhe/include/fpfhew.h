/*************
 * This file is modified from fhew.h which is writen by Leo Ducas and Daniele Micciancio.
 * Full paper is listed in : eprint.iacr.org/2022/186.
 * This source is obeyed and applied following copyright law. 
 *
 * If you have any question, please contact us by the email kr3951@hanyang.ac.kr
 * Seunghwan Lee, 2022.03.04
 * *************/


// @file fhew.h - FHEW scheme header file
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

#ifndef FHE_FP_FHEW_H
#define FHE_FP_FHEW_H

#include "fplwe.h"
#include "fpringcore.h"

namespace lbcrypto {

/**
 * @brief Ring GSW accumulator schemes described in
 * https://eprint.iacr.org/2014/816 and "Bootstrapping in FHEW-like
 * Cryptosystems"
 */
class FPRingGSWAccumulatorScheme {
 public:
  FPRingGSWAccumulatorScheme() {}

  /**
   * Generates a refreshing key
   *
   * @param params a shared pointer to RingGSW scheme parameters
   * @param lwescheme a shared pointer to additive LWE scheme
   * @param LWEsk a shared pointer to the secret key of the underlying additive
   * LWE scheme
   * @return a shared pointer to the refreshing key
   */
  FPRingGSWEvalKey KeyGen(
    const std::shared_ptr<FPRingGSWCryptoParams> params,
    const std::shared_ptr<FPLWEEncryptionScheme> lwescheme,
		const std::shared_ptr<FPRLWEEncryptionScheme> rlwescheme,
		const std::shared_ptr<const FPRLWEPrivateKeyImpl> RLWEsk) const;

  /**
   * Evaluates a binary gate (calls bootstrapping as a subroutine)
   *
   * @param params a shared pointer to RingGSW scheme parameters
   * @param gate the gate; can be AND, OR, NAND, NOR, XOR, or XOR
   * @param &EK a shared pointer to the bootstrapping keys
   * @param ct1 first ciphertext
   * @param ct2 second ciphertext
   * @param lwescheme a shared pointer to additive LWE scheme
   * @return a shared pointer to the resulting ciphertext
   */

	/*
  std::shared_ptr<FPLWECiphertextImpl> EvalBinGate(
      const std::shared_ptr<FPRingGSWCryptoParams> params, const FPGATE gate,
      const FPRingGSWEvalKey &EK,
      const std::shared_ptr<const FPLWECiphertextImpl> ct1,
      const std::shared_ptr<const FPLWECiphertextImpl> ct2,
      const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme) const;
*/
  /**
   * Evaluates NOT gate
   *
   * @param params a shared pointer to RingGSW scheme parameters
   * @param ct1 the input ciphertext
   * @return a shared pointer to the resulting ciphertext
   */
	/*
  std::shared_ptr<FPLWECiphertextImpl> EvalNOT(
      const std::shared_ptr<FPRingGSWCryptoParams> params,
      const std::shared_ptr<const FPLWECiphertextImpl> ct1) const;
	*/
  /**
   * Bootstraps a fresh ciphertext
   *
   * @param params a shared pointer to RingGSW scheme parameters
   * @param &EK a shared pointer to the bootstrapping keys
   * @param ct1 input ciphertext
   * @param lwescheme a shared pointer to additive LWE scheme
   * @return a shared pointer to the resulting ciphertext
   */
  std::shared_ptr<const FPRLWECiphertextImpl> BootstrapOrign(
      const std::shared_ptr<FPRingGSWCryptoParams> params,
			const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
      const std::shared_ptr<FPRingGSWEvalKey> EK,
      const vector<DCRTPoly> &A_poly,
			const DCRTPoly & B_poly
	) const ;
 
	std::shared_ptr<vector<DCRTPoly>> BootstrapOrignPartial(
      const std::shared_ptr<FPRingGSWCryptoParams> params,
			const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
      const std::shared_ptr<FPRingGSWEvalKey> EK,
      const vector<DCRTPoly> &CTs,
			const uint32_t from_bit, const uint32_t to_bit
	) const ;


	vector<std::shared_ptr<const FPRLWECiphertextImpl>> RingKeySwitch(
     const std::shared_ptr<FPRingGSWCryptoParams> params,
     const std::shared_ptr<FPRLWESwitchingKey> K,
     const std::shared_ptr<const FPRLWECiphertextImpl> ct,
     uint32_t idx,
		 uint32_t out_nums
     ) ;
 
  std::shared_ptr<vector<vector<DCRTPoly>>> RingKeySwitchPoly(
     const std::shared_ptr<FPRingGSWCryptoParams> params,
     const std::shared_ptr<FPRLWESwitchingKey> K,
     const vector<DCRTPoly> &input,
     const uint32_t &idx,
		 const uint32_t &out_nums
     ) const ;
 
  std::shared_ptr<vector<vector<DCRTPoly>>> RingKeySwitchPolyPARALLEL(
     const std::shared_ptr<FPRingGSWCryptoParams> params,
     const std::shared_ptr<FPRLWESwitchingKey> K,
     const vector<DCRTPoly> &input,
     const uint32_t &idx,
		 const uint32_t &out_nums
     ) const ;
 
  std::shared_ptr<vector<vector<DCRTPoly>>> RingKeySwitchPolyPARALLEL2(
     const std::shared_ptr<FPRingGSWCryptoParams> params,
     const std::shared_ptr<FPRLWESwitchingKey> K,
     const vector<DCRTPoly> &input,
     const uint32_t &idx,
		 const uint32_t &out_nums
     ) const ;


 std::shared_ptr<vector<vector<DCRTPoly>>> Bootstrap_and_KS( 
    const std::shared_ptr<FPRingGSWCryptoParams> params, 
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const DCRTPoly &ACC_init, 
    const std::shared_ptr<FPRingGSWEvalKey> &EK, const vector<NativeVector> &A, const NativeInteger &b, 
    const uint32_t boot_loc_idx, const uint32_t boot_first_bits,
    const uint64_t prefix_add, const uint32_t boot_out_bit,
    const uint32_t out_num) const ;




 std::shared_ptr<vector<vector<DCRTPoly>>>  ProductAndBKAndKS(
    const std::shared_ptr<FPRingGSWCryptoParams> params, 
    const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const DCRTPoly &ACC_init,
    const std::shared_ptr<FPRingGSWEvalKey> &EK, 
    const vector<DCRTPoly> &C1, 
    const vector<DCRTPoly> &C2, 
    const uint32_t boot_loc_idx, const uint32_t boot_first_bits,
    const uint64_t prefix_add, const uint32_t boot_out_bit,
    const uint32_t out_num) const ;

 
 std::shared_ptr<vector<DCRTPoly>>  ProductAndBKAndKSRange(
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
    const uint32_t out_num, bool mul = true) const ;


 std::shared_ptr<FPLWECiphertextImpl> PreKeySwitch(
		const std::shared_ptr<FPRingGSWCryptoParams> l_params
		,const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme
		,const std::shared_ptr<FPRingGSWEvalKey> &EK
		,const vector<NativeVector> &a, const uint32_t loc_idx, const uint32_t boot_first_bits
		, const uint32_t boot_out_bit) const;

 std::shared_ptr<FPLWECiphertextImpl> PreKeySwitchPARALLEL(
		const std::shared_ptr<FPRingGSWCryptoParams> l_params
		,const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme
		,const std::shared_ptr<FPRingGSWEvalKey> &EK
		,const vector<NativeVector> &a, const uint32_t loc_idx, const uint32_t boot_first_bits
		, const uint32_t boot_out_bit) const;



 std::shared_ptr<vector<DCRTPoly>> SGDLevel2(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const DCRTPoly &in_poly, 
		const uint32_t &bits, 
		const uint32_t &rm_len, 
		const Format format) const ;

 std::shared_ptr<vector<DCRTPoly>> SGDLevel2PARALLEL(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const DCRTPoly &in_poly, 
		const uint32_t &bits, 
		const uint32_t &rm_len, 
		const Format format) const ;


 std::shared_ptr<vector<vector<int32_t>>> SGDLevel2ForKS(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const DCRTPoly &in_poly, 
		const uint32_t &bits, 
		const uint32_t &rm_len, 
		const Format format) const ;

 std::shared_ptr<vector<vector<int32_t>>> SGDLevel2ForKSPARALLEL(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const DCRTPoly &in_poly, 
		const uint32_t &bits, 
		const uint32_t &rm_len, 
		const Format format) const ;


/*
 std::shared_ptr<vector<vector<int32_t>>> SGDLevel2CompactPARALLEL(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const DCRTPoly &in_poly, 
		const uint32_t &bits, 
		const uint32_t &rm_len, 
		const Format format) const ;
*/


 std::shared_ptr<DCRTPoly> DivideRounding(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const DCRTPoly &poly,
		const Format format) const;
 
 void DivideRounding(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const DCRTPoly &poly,
		std::shared_ptr<DCRTPoly> out,
		const Format format) const;

 
 void DivideRoundingSelf(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		DCRTPoly *poly,
		const Format format) const;

std::shared_ptr<vector<NativeVector>> RotandCarryAddA(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const vector<NativePoly> &input,
		vector<DCRTPoly> *CarryInfoA,
		const NativePoly &rot,
		//const uint32_t &idx,
		const bool rot_flag,
		const bool carry_flag) const ;



 std::shared_ptr<vector<DCRTPoly>> BootstrapRelu(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<DCRTPoly> & ct1_A,
		const DCRTPoly &ct1_B,
		const vector<DCRTPoly> & ct2_A,
		const DCRTPoly &ct2_B,
		const std::shared_ptr<FPRingGSWEvalKey> EK) const;


 std::shared_ptr<vector<DCRTPoly>> BootstrapReluF(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<DCRTPoly> & ct1_A,
		const DCRTPoly &ct1_B,
		const vector<DCRTPoly> & ct2_A,
		const DCRTPoly &ct2_B,
		const std::shared_ptr<FPRingGSWEvalKey> EK) const;





 std::shared_ptr<vector<DCRTPoly>> BootstrapRecon(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<DCRTPoly> & ct,
		const std::shared_ptr<FPRingGSWEvalKey> EK,
		const uint32_t from_idx,
		const uint32_t to_idx) const;


 std::shared_ptr<vector<vector<DCRTPoly>>> BootstrapRecon2(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<vector<DCRTPoly>> & ct,
		const std::shared_ptr<FPRingGSWEvalKey> EK,
		const uint32_t from_idx) const;




 
 std::shared_ptr<vector<DCRTPoly>> BootstrapMulPolyMake0(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<DCRTPoly> & ct_ASign,
		const DCRTPoly &ct_Bsign,
		const vector<DCRTPoly> & ct_Exp,
		const std::shared_ptr<FPRingGSWEvalKey> EK,
		const uint32_t max_msg,
		bool mul) const;


 
 std::shared_ptr<vector<DCRTPoly>> BootstrapMulPolyMake(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<DCRTPoly> & ct_Exp,
		const std::shared_ptr<FPRingGSWEvalKey> EK,
		const uint32_t max_msg,
		bool mul) const;



 std::shared_ptr<vector<DCRTPoly>> BootstrapIsLower(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<DCRTPoly> & ct,
		const uint64_t &mins,
		const std::shared_ptr<FPRingGSWEvalKey> EK) const;





 std::shared_ptr<vector<DCRTPoly>> BootstrapIsLowerF(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<DCRTPoly> & ct,
		const uint64_t &mins,
		const std::shared_ptr<FPRingGSWEvalKey> EK) const;




 std::shared_ptr<vector<DCRTPoly>> BootstrapIsLoweR(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<DCRTPoly> & ct,
		const uint64_t &mins,
		const std::shared_ptr<FPRingGSWEvalKey> EK) const;






 std::shared_ptr<vector<DCRTPoly>> BootstrapReluRev(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<DCRTPoly> & ct,
		const uint64_t &mins,
		const std::shared_ptr<FPRingGSWEvalKey> EK) const;

 std::shared_ptr<vector<DCRTPoly>> BootstrapReluRevF(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<DCRTPoly> & ct,
		const uint64_t &mins,
		const std::shared_ptr<FPRingGSWEvalKey> EK) const;




 std::shared_ptr<vector<vector<DCRTPoly>>> CheckOFUF(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<DCRTPoly> & ct_UP,
		const vector<DCRTPoly> & ct_EXP,
		const CTType types,
		const std::shared_ptr<FPRingGSWEvalKey> EK) const;


 std::shared_ptr<vector<vector<DCRTPoly>>> CheckOFUFMul(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<DCRTPoly> & ct_UP,
		const vector<DCRTPoly> & ct_EXP,
		const CTType types,
		const std::shared_ptr<FPRingGSWEvalKey> EK) const;



 std::shared_ptr<vector<vector<DCRTPoly>>> BootstrapReluRev(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<DCRTPoly> & ct_A,
		const vector<DCRTPoly> & ct_B,
		const std::shared_ptr<FPRingGSWEvalKey> EK) const;



 std::shared_ptr<vector<vector<DCRTPoly>>> BootstrapReluRevF(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<DCRTPoly> & ct_A,
		const vector<DCRTPoly> & ct_B,
		const std::shared_ptr<FPRingGSWEvalKey> EK) const;





 std::shared_ptr<vector<DCRTPoly>> BootstrapMul(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<DCRTPoly> & ct,
		const std::shared_ptr<FPRingGSWEvalKey> EK) const;


 std::shared_ptr<vector<vector<DCRTPoly>>> BootstrapMulDouble64(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<DCRTPoly> & ct,
		const std::shared_ptr<FPRingGSWEvalKey> EK,
		const uint32_t from_bits, const uint32_t to_bits, 
		const int32_t zero_start_bits,
		const CTType types
		) const;

/*
 std::shared_ptr<vector<vector<DCRTPoly>>> BootstrapMulFloat(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<DCRTPoly> & ct,
		const std::shared_ptr<FPRingGSWEvalKey> EK,
		const uint32_t from_bits, const uint32_t to_bits, 
		const int32_t zero_start_bits) const;

*/

 std::shared_ptr<vector<DCRTPoly>> BootstrapMulOne(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<DCRTPoly> & ct,
		const std::shared_ptr<FPRingGSWEvalKey> EK) const;


std::shared_ptr<vector<vector<DCRTPoly>>> BootstrapAdd(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<DCRTPoly> & ct,
		const std::shared_ptr<FPRingGSWEvalKey> EK) const;


std::shared_ptr<vector<vector<DCRTPoly>>> BootstrapAddPartial(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<DCRTPoly> & ct,
		const std::shared_ptr<FPRingGSWEvalKey> EK,
		const uint32_t from_bit, 
		const uint32_t to_bit,
		const int32_t muls) const;


std::shared_ptr<vector<vector<DCRTPoly>>> BootstrapAddPartialDouble(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<DCRTPoly> & ct,
		const std::shared_ptr<FPRingGSWEvalKey> EK,
		const uint32_t from_bit, 
		const uint32_t to_bit,
		const int32_t muls) const;





 std::shared_ptr<vector<DCRTPoly>> BootstrapSign(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
		const vector<DCRTPoly> & ct,
		const std::shared_ptr<FPRingGSWEvalKey> EK) const;


 /**
   * Core bootstrapping operation
   *
   * @param params a shared pointer to RingGSW scheme parameters
   * @param &EK a shared pointer to the bootstrapping keys
   * @param gate the gate; can be AND, OR, NAND, NOR, XOR, or XOR
   * @param &a first part of the input LWE ciphertext
   * @param &b second part of the input LWE ciphertext
   * @param lwescheme a shared pointer to additive LWE scheme
   * @return the output RingLWE accumulator
   */
  std::shared_ptr<FPRingGSWCiphertext> BootstrapCore(
      const std::shared_ptr<FPRingGSWCryptoParams> params, 
			const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme,
			const DCRTPoly &ACC_init,
      const std::shared_ptr<FPRingGSWEvalKey> &EK, const vector<NativeVector> &a, const NativeInteger &b,
			const uint32_t boot_loc_idx, const uint32_t boot_start_bits
			, const uint64_t prefix_add
			,const uint32_t boot_out_bit
			//,const std::shared_ptr<FPLWEEncryptionScheme> LWEscheme
			) const;

 
	std::shared_ptr<vector<DCRTPoly>> ProductDCRT(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const std::shared_ptr<FPRingEvaluationKey> Evkey,
		const vector<DCRTPoly> &A1, const DCRTPoly &b1,
		const vector<DCRTPoly> &A2, const DCRTPoly &b2) const;

	std::shared_ptr<vector<DCRTPoly>> ProductDCRTPARALLEL(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const std::shared_ptr<FPRingEvaluationKey> Evkey,
		const vector<DCRTPoly> &A1, const DCRTPoly &b1,
		const vector<DCRTPoly> &A2, const DCRTPoly &b2) const;

	std::shared_ptr<vector<DCRTPoly>> ProductDCRTPARALLEL2(
		const std::shared_ptr<FPRingGSWCryptoParams> params,
		const std::shared_ptr<FPRingEvaluationKey> Evkey,
		const vector<DCRTPoly> &A1, const DCRTPoly &b1,
		const vector<DCRTPoly> &A2, const DCRTPoly &b2) const;





 private:
  /**
   * Generates a refreshing key - GINX variant
   *
   * @param params a shared pointer to RingGSW scheme parameters
   * @param lwescheme a shared pointer to additive LWE scheme
   * @param LWEsk a shared pointer to the secret key of the underlying additive
   * LWE scheme
   * @return a shared pointer to the refreshing key
   */
  FPRingGSWEvalKey KeyGenGINX(
      const std::shared_ptr<FPRingGSWCryptoParams> params,
			const std::shared_ptr<FPLWEEncryptionScheme> lwescheme,
			const std::shared_ptr<FPRLWEEncryptionScheme> rlwescheme,
			const std::shared_ptr<const FPRLWEPrivateKeyImpl> RLWEsk) const;

  /**
   * Generates a refreshing key - AP variant
   *
   * @param params a shared pointer to RingGSW scheme parameters
   * @param lwescheme a shared pointer to additive LWE scheme
   * @param LWEsk a shared pointer to the secret key of the underlying additive
   * LWE scheme
   * @return a shared pointer to the refreshing key
   */
  FPRingGSWEvalKey KeyGenAP(
      const std::shared_ptr<FPRingGSWCryptoParams> params,
      const std::shared_ptr<FPLWEEncryptionScheme> lwescheme,
      const std::shared_ptr<const FPLWEPrivateKeyImpl> LWEsk) const;

  /**
   * Internal RingGSW encryption used in generating the refreshing key - AP
   * variant
   *
   * @param params a shared pointer to RingGSW scheme parameters
   * @param skFFT secret key polynomial in the EVALUATION representation
   * @param m plaintext (corresponds to a lookup entry for the LWE scheme secret
   * key)
   * @return a shared pointer to the resulting ciphertext
   */
  std::shared_ptr<FPRingGSWCiphertext> EncryptAP(
      const std::shared_ptr<FPRingGSWCryptoParams> params,
      const NativePoly &skFFT, const FPLWEPlaintext &m) const;

  /**
   * Internal RingGSW encryption used in generating the refreshing key - GINX
   * variant
   *
   * @param params a shared pointer to RingGSW scheme parameters
   * @param skFFT secret key polynomial in the EVALUATION representation
   * @param m plaintext (corresponds to a lookup entry for the LWE scheme secret
   * key)
   * @return a shared pointer to the resulting ciphertext
   */
  std::shared_ptr<FPRingGSWCiphertext> EncryptGINX(
      const std::shared_ptr<FPRingGSWCryptoParams> params,
      const vector<DCRTPoly> &skFFT, const FPLWEPlaintext &m) const;

  std::shared_ptr<FPRingGSWCiphertext> EncryptKS(
      const std::shared_ptr<FPRingGSWCryptoParams> params,
      const vector<DCRTPoly> &skFFT, const uint64_t &m) const;



  /**
   * Main accumulator function used in bootstrapping - AP variant
   *
   * @param params a shared pointer to RingGSW scheme parameters
   * @param &input input ciphertext
   * @param acc previous value of the accumulator
   */
  void AddToACCAP(const std::shared_ptr<FPRingGSWCryptoParams> params,
                  const FPRingGSWCiphertext &input,
                  std::shared_ptr<FPRingGSWCiphertext> acc) const;

  /**
   * Main accumulator function used in bootstrapping - GINX variant
   *
   * @param params a shared pointer to RingGSW scheme parameters
   * @param &input input ciphertext
   * @param &a integer a in each step of GINX accumulation
   * @param acc previous value of the accumulator
   */
  void AddToACCGINX(const std::shared_ptr<FPRingGSWCryptoParams> params,
                    const FPRingGSWCiphertext &input, const uint64_t &a,
                    std::shared_ptr<FPRingGSWCiphertext> acc) const;
  
	void AddToACCGINXParallel(const std::shared_ptr<FPRingGSWCryptoParams> params,
                    const FPRingGSWCiphertext &input, const uint64_t &a,
                    const vector<DCRTPoly> & ct,
										std::shared_ptr<vector<DCRTPoly>> outs) const;


  /**
   * Takes an RLWE ciphertext input and outputs a vector of its digits, i.e., an
   * RLWE' ciphertext
   *
   * @param params a shared pointer to RingGSW scheme parameters
   * @param &input input RLWE ciphertext
   * @param *output input RLWE ciphertext
   */
  inline void SignedDigitDecompose(
      const std::shared_ptr<FPRingGSWCryptoParams> params,
      const std::vector<NativePoly> &input,
      std::vector<NativePoly> *output) const;
	};
	
}  // namespace lbcrypto

#endif
