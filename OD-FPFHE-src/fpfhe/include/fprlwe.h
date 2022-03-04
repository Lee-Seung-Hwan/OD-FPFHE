/*************
 * This file is modified from fplwe.h which is writen by Leo Ducas and Daniele Micciancio.
 * Full paper is listed in : eprint.iacr.org/2022/186.
 * This source is obeyed and applied following copyright law. 
 *
 * If you have any question, please contact us by the email kr3951@hanyang.ac.kr
 * Seunghwan Lee, 2022.03.04
 * *************/



// @file fplwe.h - Main ring classes for Boolean circuit FHE.
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

#ifndef FPFHE_RLWE_MECRO_H
#define FPFHE_RLWE_MECRO_H

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "lattice/backend.h"
//#include "fpringcore.h"

#include "fplwecore.h"
#include "math/backend.h"
#include "math/discretegaussiangenerator.h"
#include "math/nbtheory.h"
#include "math/transfrm.h"
#include "utils/serializable.h"



namespace lbcrypto {
	
	typedef NativeVector FPRLWEPlaintext; 

	/**
	 * @brief Class that stores all parameters for the RingGSW scheme used in
   * bootstrapping
   */
	

	// Evaluation Key!!
	// SH
	class FPRingEvaluationKey : public Serializable {
		public:
		FPRingEvaluationKey() {}

		explicit FPRingEvaluationKey(uint32_t k_time_k_over_2, uint32_t decompose_dim , uint32_t module_dim) {
			
			m_e_key.resize(k_time_k_over_2);
			for (uint32_t i = 0; i < k_time_k_over_2; i++) {
				m_e_key[i].resize(decompose_dim);
				for (uint32_t j = 0; j < decompose_dim; j++) {
					m_e_key[i][j].resize(module_dim);
				}
			}		
		}
		explicit FPRingEvaluationKey(
			const std::vector<std::vector<std::vector<DCRTPoly>>>& key)
      : m_e_key(key) {}

		explicit FPRingEvaluationKey(const FPRingEvaluationKey& rhs) { this->m_e_key = rhs.m_e_key; }

		explicit FPRingEvaluationKey(const FPRingEvaluationKey&& rhs) {
			this->m_e_key = std::move(rhs.m_e_key);
		}

		const FPRingEvaluationKey& operator=(const FPRingEvaluationKey& rhs) {
			this->m_e_key = rhs.m_e_key;
			return *this;
		}

		const FPRingEvaluationKey& operator=(const FPRingEvaluationKey&& rhs) {
			this->m_e_key = std::move(rhs.m_e_key);
			return *this;
		}

		const std::vector<std::vector<std::vector<DCRTPoly>>>& GetElements()
      const {
			return m_e_key;
		}

		void SetElements(
      const std::vector<std::vector<std::vector<DCRTPoly>>>& key) {
			m_e_key = key;
		}

		std::vector<std::vector<DCRTPoly>>& operator[](uint32_t i) {
			return m_e_key[i];
		}

		const std::vector<std::vector<DCRTPoly>>& operator[](usint i) const {
			return m_e_key[i];
		}

		bool operator==(const FPRingEvaluationKey& other) const {
			return m_e_key == other.m_e_key;
		}

		bool operator!=(const FPRingEvaluationKey& other) const { return !(*this == other); }

		template <class Archive>
		void save(Archive& ar, std::uint32_t const version) const {
			ar(::cereal::make_nvp("key", m_e_key));
		}

		template <class Archive>
		void load(Archive& ar, std::uint32_t const version) {
			if (version > SerializedVersion()) {
      PALISADE_THROW(deserialize_error,
                     "serialized object version " + std::to_string(version) +
                         " is from a later version of the library");
			}
			ar(::cereal::make_nvp("key", m_e_key));
		}

		std::string SerializedObjectName() const { return "FPRingEvaluationKey"; }
		static uint32_t SerializedVersion() { return 1; }

		private:
		std::vector<std::vector<std::vector<DCRTPoly>>> m_e_key;
	};



	class FPRLWECryptoParams : public Serializable {
		public:
			FPRLWECryptoParams()
				: m_baseGD(0), m_digitsG(0), m_digitsG2(0), m_r_std(0.0) {}
	
  /**
   * Main constructor for RingCryptoParams
   *
   * @param lweparams a shared poiter to an instance of LWECryptoParams
   * @param baseG the gadget base used in the bootstrapping
   * @param baseR the base for the refreshing key
   * @param method bootstrapping method (AP or GINX)
   */	
		explicit FPRLWECryptoParams(vector<NativeInteger> moduliQ,
															 const std::shared_ptr<ILDCRTParams<BigInteger>> DCRTParamset,
															 uint32_t degree_N,
															 uint32_t K_FP,
                               uint32_t baseGD, uint32_t l_rm_GD, double r_std)
				:
				  m_baseGD(baseGD),
					//m_l_rm_GD(l_rm_GD),
					m_N(degree_N),
					m_K(K_FP),
					m_r_std(r_std),
					m_moduliQ(moduliQ),
					m_DCRTpolyParams(DCRTParamset)
					{
						
						//if (!IsPowerOfTwo(baseG)) {
						//	PALISADE_THROW(config_error, "Gadget base should be a power of two.");
						//}
						m_r_dgg.SetStd(r_std);
						PreCompute();
					}
	
  /**
   * Performs precomputations based on the supplied parameters
   */
			void PreCompute() {
				
				// Precompute mu
				m_mu_FP.resize(m_moduliQ.size());
				for (uint32_t i = 0; i < m_mu_FP.size(); i++) {
					m_mu_FP[i] = m_moduliQ[i].ComputeMu();
				
				}

				// ??
				//NativeInteger Q = m_moduliQ[0];
				m_digitsG.push_back((uint32_t)std::ceil(log(m_moduliQ[0].ConvertToDouble()) /
                                    log(static_cast<double>(m_baseGD))));
				m_digitsG.push_back((uint32_t)std::ceil(log(m_moduliQ[1].ConvertToDouble()) /
                                    log(static_cast<double>(m_baseGD))));
			
				m_digitsG2.push_back(m_digitsG[0] * 2);
				m_digitsG2.push_back(m_digitsG[1] * 2);
				
				// Computes baseG^i
				
				vector<NativeInteger> vTemp;

				vTemp.push_back(NativeInteger(1));
				vTemp.push_back(NativeInteger(1));

				for (uint32_t i = 0; i < (m_digitsG[0]+m_digitsG[1]); i++) {
					m_Gpower.push_back(vTemp);
					vTemp[0] = vTemp[0].ModMul(NativeInteger(m_baseGD), m_moduliQ[0]);
					vTemp[1] = vTemp[1].ModMul(NativeInteger(m_baseGD), m_moduliQ[1]);
				}

				// each bits saved
				uint32_t bits_Q0 = 0;
				uint32_t bits_Q1 = 0;
				NativeInteger scaling_0 = m_moduliQ[0] - 1;
				NativeInteger scaling_1 = m_moduliQ[1] - 1;
				uint64_t get_bits_over_2 = 1;
				
				while(true) {
					scaling_0 = scaling_0 / 2;
					bits_Q0 += 1;
					if (scaling_0  % 2 == 1) {
						break;
					}
				}

				while(true) {
					scaling_1 = scaling_1 / 2;
					bits_Q1 += 1;
					get_bits_over_2 *= 2;
					
					if (scaling_1  % 2 == 1) {
						break;
					}
				}
				if (scaling_1 != scaling_0) {
					std::string errMSG = "The scaling factor Q and P does not fit";
					PALISADE_THROW(config_error, errMSG);
				} else {
					std::cout << "Delta is " <<  scaling_1 << std::endl;
				}
				m_scaling_r = scaling_0;
				m_bits_Q.push_back(bits_Q0);
				m_bits_Q.push_back(bits_Q1);
				m_getbits_over_2 = get_bits_over_2 / 2;


				// Inv
				NativeInteger Q0_basis			= NativeInteger(m_moduliQ[1]);
				NativeInteger Q1_basis			= NativeInteger(m_moduliQ[0]);

				m_QInv.push_back(Q0_basis.ModInverse(m_moduliQ[0])); // q1^-1 mod q0
				m_QInv.push_back(Q1_basis.ModInverse(m_moduliQ[1])); // q0^-1 mod q1
				

				if (m_QInv[0].ModMul(m_moduliQ[1],m_moduliQ[0]).ConvertToInt() != 1) {
					std::string errMSG = "The QInv is not correct";
					PALISADE_THROW(config_error, errMSG);
				}
				if (m_QInv[1].ModMul(m_moduliQ[0],m_moduliQ[1]).ConvertToInt() != 1) {
					std::string errMSG = "The QInv is not correct";
					PALISADE_THROW(config_error, errMSG);
				}
				// Make Q/2
				BigInteger Q_over_2 = (BigInteger(m_moduliQ[0]) * BigInteger(m_moduliQ[1])).DivideAndRound(BigInteger(2));
				m_Q_over_2.resize(2);
				m_Q_over_2[0] = Q_over_2.Mod(BigInteger(m_moduliQ[0])).ConvertToInt();
				m_Q_over_2[1] = Q_over_2.Mod(BigInteger(m_moduliQ[1])).ConvertToInt();


				// double_inv
				m_QInv_double.push_back(1./static_cast<double>(m_moduliQ[0].ConvertToInt()));
				m_QInv_double.push_back(1./static_cast<double>(m_moduliQ[1].ConvertToInt()));
			
				/*
				const BigInteger BarrettBase128Bit("340282366920938463463374607431768211456");
				const BigInteger TwoPower64("18446744073709551616");
				m_modqBarrettMu.resize(m_moduliQ.size());
				for (uint32_t i = 0; i < m_moduliQ.size(); i++) {
					BigInteger mu = BarrettBase128Bit / BigInteger(m_moduliQ[i]);
					uint64_t val[2];
					val[0] = (mu % TwoPower64).ConvertToInt();
					val[1] = mu.RShift(64).ConvertToInt();
					memcpy(&m_modqBarrettMu[i], val, sizeof(DoubleNativeInt));
				}
				*/
				//usint qMSB = m_moduliQ[0].GetMSB();
				//usint sizeQMSB = GetMSB64(sizeQ);
				//m_tQHatInvModqDivqModt.resize(sizeQ);
				//m_tQHatInvModqDivqModtPrecon.resize(sizeQ);
				//m_tQHatInvModqDivqFrac.resize(sizeQ);



				// Computes polynomials X^m - 1 that are needed in the accumulator for the
				// GINX bootstrapping
				//if (m_method == GINX) {
				// loop for positive values of m
				//
			/*
				for (uint32_t i = 0; i < N; i++) {
					NativePoly aPoly = NativePoly(m_polyParams, Format::COEFFICIENT, true);
					aPoly[i].ModAddEq(NativeInteger(1), Q);  // X^m
					aPoly[0].ModSubEq(NativeInteger(1), Q);  // -1
					aPoly.SetFormat(Format::EVALUATION);
					m_monomials.push_back(aPoly);
				}

      // loop for negative values of m
				for (uint32_t i = 0; i < N; i++) {
					NativePoly aPoly = NativePoly(m_polyParams, Format::COEFFICIENT, true);
					aPoly[i].ModSubEq(NativeInteger(1), Q);  // -X^m
					aPoly[0].ModSubEq(NativeInteger(1), Q);  // -1
					aPoly.SetFormat(Format::EVALUATION);
					m_monomials.push_back(aPoly);
				}
			*/
#if defined(FPFHE_DEBUG)
				/*
				std::cerr << "base_g = " << m_baseG << std::endl;
				std::cerr << "m_digitsG = " << m_digitsG << std::endl;
				std::cerr << "m_digitsG2 = " << m_digitsG2 << std::endl;
				std::cerr << "m_baseR = " << m_baseR << std::endl;
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
			}
			
			const std::shared_ptr<ILDCRTParams<BigInteger>> GetDCRTParams() const {
				return m_DCRTpolyParams;
			}
			uint32_t GetBaseG() const { return m_baseGD; }

			vector<uint32_t> GetDigitsG() const { return m_digitsG; }
			
			const DiscreteGaussianGeneratorImpl<NativeVector> &GetDgg() const {
				return m_r_dgg;
			}

			vector<uint32_t> GetDigitsG2() const { return m_digitsG2; }

			uint32_t GetN() const {return m_N;}
			uint32_t GetK() const {return m_K;}
			uint32_t GetNK() const {return m_K * m_N;}



			uint64_t GetBits_over_2() const {return m_getbits_over_2;}
				
			uint32_t GetBitsTotal() const { return m_bits_Q_total;} 
			void SetBitsTotal(uint32_t bits) { m_bits_Q_total = bits;} 
				
			NativeInteger GetScaling() const {return m_scaling_r;}
			
			vector<NativeInteger> GetInv() const {return m_QInv;} 
			vector<double> GetInvDouble() const {return m_QInv_double;} 
			vector<uint32_t> GetBits() const {return m_bits_Q;}
			vector<uint64_t> GetQOver2Info() const {return m_Q_over_2;}
			NativeInteger GetQ() const {return m_q;}
			vector<NativeInteger> GetModuliQ() const {return m_moduliQ;}
			vector<NativeInteger> GetMuFP() const {return m_mu_FP;}


			//uint32_t GetBaseR() const { return m_baseR; }

			const std::vector<NativeInteger>& GetDigitsR() const { return m_digitsR; }

			//const shared_ptr<ILNativeParams> GetRPolyParams() const {
			//	return m_polyParams;
			//}

			const std::vector<vector<NativeInteger>>& GetGPower() const { return m_Gpower; }
			const std::vector<NativeInteger>& GetGateConst() const { return m_gateConst; }
			const NativePoly& GetMonomial(uint32_t i) const { return m_monomials[i]; }


			//bool operator==(const FPRLWECryptoParams& other) const {
			//	return (*m_LWEParams == *other.m_LWEParams && m_baseR == other.m_baseR &&
      //     m_baseGD == other.m_baseGD;
			//}

			//bool operator!=(const FPRLWECryptoParams& other) const {
			//	return !(*this == other);
			//}

			template <class Archive>
			void save(Archive& ar, std::uint32_t const version) const {
				//ar(::cereal::make_nvp("params", m_LWEParams));
				//ar(::cereal::make_nvp("bR", m_baseR));
				ar(::cereal::make_nvp("bG", m_baseGD));
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
				ar(::cereal::make_nvp("bG", m_baseGD));

				this->PreCompute();
			}

			std::string SerializedObjectName() const { return "FPRingGSWCryptoParams"; }
			static uint32_t SerializedVersion() { return 1; }
		

		private:
			// gadget base used in bootstrapping
			uint32_t m_baseGD;
			// number of digits in decomposing integers mod Q
			vector<uint32_t> m_digitsG;
			// twice the number of digits in decomposing integers mod Q
			vector<uint32_t> m_digitsG2;
			// base used for the refreshing key (used only for AP bootstrapping)
			//uint32_t m_l_rm_GD;
			// N
			uint32_t m_N;
			uint32_t m_K;
			// Q
			NativeInteger m_q;
			NativeInteger m_scaling_r;
			uint64_t m_getbits_over_2;
			// standard deviation
			double m_r_std;
			// Discrete Gaussian Generator
			DiscreteGaussianGeneratorImpl<NativeVector> m_r_dgg;
			// powers of m_baseR (used only for AP bootstrapping)
			std::vector<NativeInteger> m_digitsR;
			// A vector of powers of baseG
			std::vector<std::vector<NativeInteger>> m_Gpower;
			// Parameters for polynomials in RingGSW/RingLWE
			shared_ptr<ILNativeParams> m_polyParams;
			
			// 
			// moduliQ
			vector<NativeInteger> m_moduliQ;
			vector<double> m_QInv_double;
			vector<NativeInteger> m_QInv;
			vector<DoubleNativeInt> m_modqBarrettMu;
			vector<uint32_t> m_bits_Q;
			vector<uint64_t> m_Q_over_2;
			uint32_t m_bits_Q_total;
			// DCRTPolyParam
			const shared_ptr<ILDCRTParams<BigInteger>> m_DCRTpolyParams;
			// Constants used in evaluating binary gates
			std::vector<NativeInteger> m_gateConst;
			// Precomputed polynomials in Format::EVALUATION representation for X^m - 1
			// (used only for GINX bootstrapping)
			std::vector<NativePoly> m_monomials;
			// Bootstrapping method (AP or GINX)
			vector<NativeInteger> m_mu_FP;					
	};

/* RING LWE Private Key manegement!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

	class FPRLWEPrivateKeyImpl : public Serializable {
		public:
			FPRLWEPrivateKeyImpl() {}

			explicit FPRLWEPrivateKeyImpl(const vector<DCRTPoly> &s) : m_s(s) {}

			explicit FPRLWEPrivateKeyImpl(const FPRLWEPrivateKeyImpl &rhs) {
				this->m_s = rhs.m_s;
			}

			explicit FPRLWEPrivateKeyImpl(const FPRLWEPrivateKeyImpl &&rhs) {
				this->m_s = std::move(rhs.m_s);
			}

			const FPRLWEPrivateKeyImpl &operator=(const FPRLWEPrivateKeyImpl &rhs) {
				this->m_s = rhs.m_s;
				return *this;
			}

			const FPRLWEPrivateKeyImpl &operator=(const FPRLWEPrivateKeyImpl &&rhs) {
				this->m_s = std::move(rhs.m_s);
				return *this;
			}

			const vector<DCRTPoly> &GetElement() const { return m_s; }

			void SetElement(const vector<DCRTPoly> &s) { m_s = s; }

			bool operator==(const FPRLWEPrivateKeyImpl &other) const {
				return m_s == other.m_s;
			}

			bool operator!=(const FPRLWEPrivateKeyImpl &other) const {
				return !(*this == other);
			}

			template <class Archive>
			void save(Archive &ar, std::uint32_t const version) const {
				ar(::cereal::make_nvp("s", m_s));
			}

			template <class Archive>
			void load(Archive &ar, std::uint32_t const version) {
				if (version > SerializedVersion()) {
					PALISADE_THROW(deserialize_error,
                     "serialized object version " + std::to_string(version) +
                         " is from a later version of the library");
				}

				ar(::cereal::make_nvp("s", m_s));
			}

			std::string SerializedObjectName() const { return "FPLWEPrivateKey"; }
			static uint32_t SerializedVersion() { return 1; }

			private:
			vector<DCRTPoly> m_s;
 };

/*   RLWE Cipheretext!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         */

enum CTType {
	UINT64,
	UINT32,
	UINT16,
	INT64,
	INT32,
	INT16,
	DOUBLE,
	FLOAT,
	HALF
};


/**
 * @brief Class that stores a RLWE scheme ciphertext; composed of a Polynomial "a"
 * and polynomial "b"
 */
class FPRLWECiphertextImpl : public Serializable {
 public:

  FPRLWECiphertextImpl() {}

	// Unsigned
  explicit FPRLWECiphertextImpl(const vector<DCRTPoly> &a, const DCRTPoly &b)
      : m_a(a), m_b(b), m_type(UINT64), m_level(2) {}
  
	explicit FPRLWECiphertextImpl(const vector<DCRTPoly> &a, const DCRTPoly &b, uint32_t level)
      : m_a(a), m_b(b),  m_type(UINT64), m_level(level) {}
	
	// signed
  explicit FPRLWECiphertextImpl(
		const vector<DCRTPoly> &a, const DCRTPoly &b,
		const vector<DCRTPoly> &a_s, const DCRTPoly &b_s)
      : m_a(a), m_b(b), m_a_sign(a_s), m_b_sign(b_s), m_type(INT64), m_level(2) {}

	explicit FPRLWECiphertextImpl(
		const vector<DCRTPoly> &a, const DCRTPoly &b,
		const vector<DCRTPoly> &a_s, const DCRTPoly &b_s , uint32_t level)
      : m_a(a), m_b(b), m_a_sign(a_s), m_b_sign(b_s), m_type(INT64), m_level(level) {}

	// Double
  explicit FPRLWECiphertextImpl(
		const vector<DCRTPoly> &a, const DCRTPoly &b,
		const vector<DCRTPoly> &a_s, const DCRTPoly &b_s, 
		const vector<DCRTPoly> &a_e, const DCRTPoly &b_e)
      : m_a(a), m_b(b), m_a_sign(a_s), m_b_sign(b_s) 
			, m_a_expo(a_e), m_b_expo(b_e), m_type(DOUBLE), m_level(2) {} 

	explicit FPRLWECiphertextImpl(
		const vector<DCRTPoly> &a, const DCRTPoly &b,
		const vector<DCRTPoly> &a_s, const DCRTPoly &b_s, 
		const vector<DCRTPoly> &a_e, const DCRTPoly &b_e, uint32_t level)
      : m_a(a), m_b(b), m_a_sign(a_s), m_b_sign(b_s) 
			, m_a_expo(a_e), m_b_expo(b_e), m_type(DOUBLE), m_level(level) {} 



	// Double & Single & Half
  explicit FPRLWECiphertextImpl(
		const vector<DCRTPoly> &a, const DCRTPoly &b,
		const vector<DCRTPoly> &a_s, const DCRTPoly &b_s, 
		const vector<DCRTPoly> &a_e, const DCRTPoly &b_e, const CTType types)
      : m_a(a), m_b(b), m_a_sign(a_s), m_b_sign(b_s) 
			, m_a_expo(a_e), m_b_expo(b_e), m_type(types), m_level(2) {} 

	explicit FPRLWECiphertextImpl(
		const vector<DCRTPoly> &a, const DCRTPoly &b,
		const vector<DCRTPoly> &a_s, const DCRTPoly &b_s, 
		const vector<DCRTPoly> &a_e, const DCRTPoly &b_e, 
		const CTType types , uint32_t level)
      : m_a(a), m_b(b), m_a_sign(a_s), m_b_sign(b_s) 
			, m_a_expo(a_e), m_b_expo(b_e), m_type(types), m_level(level) {} 


  explicit FPRLWECiphertextImpl(const FPRLWECiphertextImpl &rhs) {
    this->m_a = rhs.m_a;
    this->m_b = rhs.m_b;
		this->m_level = rhs.m_level;
		this->m_type = rhs.m_type;

		if (this->m_type == INT64 || this->m_type == DOUBLE || this->m_type == FLOAT || this->m_type == HALF) {
			this->m_a_sign = rhs.m_a_sign;
			this->m_b_sign = rhs.m_b_sign;
		}
		if (this->m_type == DOUBLE || this->m_type == FLOAT  || this->m_type == HALF) {
			this->m_a_expo = rhs.m_a_expo;
			this->m_b_expo = rhs.m_b_expo;
		}
  }

  explicit FPRLWECiphertextImpl(const FPRLWECiphertextImpl &&rhs) {
    this->m_a = std::move(rhs.m_a);
    this->m_b = std::move(rhs.m_b);
		this->m_level = std::move(rhs.m_level);
		this->m_type = std::move(rhs.m_type);
		if (this->m_type == INT64 || this->m_type == DOUBLE || this->m_type == FLOAT  || this->m_type == HALF) {
			this->m_a_sign = std::move(rhs.m_a_sign);
			this->m_b_sign = std::move(rhs.m_b_sign);
		}
		if (this->m_type == DOUBLE || this->m_type == FLOAT  || this->m_type == HALF ) {
			this->m_a_expo = std::move(rhs.m_a_expo);
			this->m_b_expo = std::move(rhs.m_b_expo);
		}
  }


	// Not Implmented
  const FPRLWECiphertextImpl &operator=(const FPRLWECiphertextImpl &rhs) {
    this->m_a = rhs.m_a;
    this->m_b = rhs.m_b;
		this->m_level = rhs.m_level;
		this->m_type = rhs.m_type;

		if (this->m_type == INT64 || this->m_type == DOUBLE || this->m_type == FLOAT  || this->m_type == HALF) {
			this->m_a_sign = rhs.m_a_sign;
			this->m_b_sign = rhs.m_b_sign;
		}
		if (this->m_type == DOUBLE || this->m_type == FLOAT  || this->m_type == HALF) {
			this->m_a_expo = rhs.m_a_expo;
			this->m_b_expo = rhs.m_b_expo;
		}
    return *this;
  }

  const FPRLWECiphertextImpl &operator=(const FPRLWECiphertextImpl &&rhs) {
    this->m_a = std::move(rhs.m_a);
    this->m_b = std::move(rhs.m_b);
		this->m_level = std::move(rhs.m_level);
		this->m_type = std::move(rhs.m_type);
		if (this->m_type == INT64 || this->m_type == DOUBLE || this->m_type == FLOAT  || this->m_type == HALF) {
			this->m_a_sign = std::move(rhs.m_a_sign);
			this->m_b_sign = std::move(rhs.m_b_sign);
		}
		if (this->m_type == DOUBLE || this->m_type == FLOAT  || this->m_type == HALF ) {
			this->m_a_expo = std::move(rhs.m_a_expo);
			this->m_b_expo = std::move(rhs.m_b_expo);
		}
    return *this;
  }
	

  const vector<DCRTPoly> &GetA()			const { return m_a; }
  const DCRTPoly &GetB()							const { return m_b; }
  const vector<DCRTPoly> &GetASign()	const { return m_a_sign; }
  const DCRTPoly &GetBSign()					const { return m_b_sign; }
  const vector<DCRTPoly> &GetAExpo()	const { return m_a_expo; }
  const DCRTPoly &GetBExpo()					const { return m_b_expo; }
  CTType GetType()	const { return m_type; }
	
	// This bug function
	void SetType(CTType types) { m_type = types; }

  void SetA(const vector<DCRTPoly> &a) { m_a = a; }
  void SetB(const DCRTPoly &b) { m_b = b; }
	uint32_t GetLevel() const {return m_level;}
	void UpLevel() {m_level +=1;}
	void DownLevel() {m_level -=1;}


	
  bool operator==(const FPRLWECiphertextImpl &other) const {
  
		bool res = ( this->m_a == other.m_a && this->m_b == other.m_b);
		res &= (this->m_level == other.m_level);
		res &= (this->m_type == other.m_type);
		if (res && this->m_type) {
			if(this->m_type == INT64 || this->m_type == DOUBLE || this->m_type == FLOAT  || this->m_type == HALF) {
				res &= (this->m_a_sign == other.m_a_sign);
				res &= (this->m_b_sign == other.m_b_sign);
			}
			if (this->m_type == DOUBLE || this->m_type == FLOAT  || this->m_type == HALF) {
				res &= (this->m_a_expo == other.m_a_expo);
				res &= (this->m_b_expo == other.m_b_expo);
			}
		}
		return res;
	}

  bool operator!=(const FPRLWECiphertextImpl &other) const {
    return !(*this == other);
  }

  template <class Archive>
  void save(Archive &ar, std::uint32_t const version) const {
    ar(::cereal::make_nvp("poly_a", m_a));
    ar(::cereal::make_nvp("poly_b", m_b));
  }
	
  template <class Archive>
  void load(Archive &ar, std::uint32_t const version) {
    if (version > SerializedVersion()) {
      PALISADE_THROW(deserialize_error,
                     "serialized object version " + std::to_string(version) +
                         " is from a later version of the library");
    }

    ar(::cereal::make_nvp("poly_a", m_a));
    ar(::cereal::make_nvp("poly_b", m_b));
  }

  std::string SerializedObjectName() const { return "FPRLWECiphertext"; }
  static uint32_t SerializedVersion() { return 1; }

 private:
  
	// Unsigned vals
	vector<DCRTPoly> m_a;
	DCRTPoly m_b;
	// sign vals
	vector<DCRTPoly> m_a_sign;
	DCRTPoly m_b_sign;
	// expo
	vector<DCRTPoly> m_a_expo;
	DCRTPoly m_b_expo;
	
	CTType m_type;
	uint32_t m_level;
};


	/**
	* @brief Class that stores the refreshing key (used in bootstrapping)
	* A three-dimensional vector of RingGSW ciphertexts
	*/
	class FPRingPackingKey : public Serializable {
		public:
			FPRingPackingKey() {}

	// dim 1 binarize p, dim2 nums n, dim3 packing num
			explicit FPRingPackingKey(uint32_t dim1, uint32_t dim2, uint32_t dim3) {
				m_key.resize(dim1);
				for (uint32_t i = 0; i < dim1; i++) {
					m_key[i].resize(dim2);
					for (uint32_t j = 0; j < dim2; j++) {m_key[i][j].resize(dim3);}
				}
			}

			explicit FPRingPackingKey(
				const std::vector<std::vector<std::vector<FPRLWECiphertextImpl>>>& key)
				: m_key(key) {}

			explicit FPRingPackingKey(const FPRingPackingKey& rhs) { this->m_key = rhs.m_key; }

			explicit FPRingPackingKey(const FPRingPackingKey&& rhs) {
				this->m_key = std::move(rhs.m_key);
			}

			const FPRingPackingKey& operator=(const FPRingPackingKey& rhs) {
				this->m_key = rhs.m_key;
				return *this;
			}

			const FPRingPackingKey& operator=(const FPRingPackingKey&& rhs) {
				this->m_key = std::move(rhs.m_key);
				return *this;
			}

			const std::vector<std::vector<std::vector<FPRLWECiphertextImpl>>>& GetElements() const {
					return m_key;
			}

			void SetElements(const std::vector<std::vector<std::vector<FPRLWECiphertextImpl>>>& key) {
				m_key = key;
			}

			std::vector<std::vector<FPRLWECiphertextImpl>>& operator[](uint32_t i) {
				return m_key[i];
			}

			const std::vector<std::vector<FPRLWECiphertextImpl>>& operator[](usint i) const {
				return m_key[i];
			}

			bool operator==(const FPRingPackingKey& other) const {
				return m_key == other.m_key;
			}

			bool operator!=(const FPRingPackingKey& other) const { return !(*this == other); }

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

			std::string SerializedObjectName() const { return "FPRingPackingKey"; }
			static uint32_t SerializedVersion() { return 1; }

		private:
			std::vector<std::vector<std::vector<FPRLWECiphertextImpl>>> m_key;
	};

	/**
 * @brief Additive LWE scheme
 */
class FPRLWEEncryptionScheme {
 public:
  FPRLWEEncryptionScheme() {}
  /**
   * Generates a secret key of dimension n using modulus q
   *
   * @param params a shared pointer to LWE scheme parameters
   * @return a shared pointer to the secret key
   */
	std::shared_ptr<FPRLWEPrivateKeyImpl> KeyGen(
      const std::shared_ptr<FPRLWECryptoParams> params,
			uint32_t h_sparse
			) const;

  /**
   * Generates a secret key of dimension N using modulus Q
   *
   * @param params a shared pointer to LWE scheme parameters
   * @return a shared pointer to the secret key
   */
  std::shared_ptr<FPRLWEPrivateKeyImpl> KeyGenN(
      const std::shared_ptr<FPRLWECryptoParams> params) const;

  /**
   * Generates a secret key of dimension packing_N using modulus Q
   *
   * @param params a shared pointer to LWE scheme parameters
   * @return a shared pointer to the secret key
   */
  std::shared_ptr<FPRLWEPrivateKeyImpl> KeyGenPackingN(
      const std::shared_ptr<FPLWECryptoParams> params) const;

	std::shared_ptr<FPRingEvaluationKey> EvalKeyGen(
			const std::shared_ptr<FPRLWECryptoParams> params,
			uint32_t EV_Bit, uint32_t EV_len_save, uint32_t EV_len_rm,
			DiscreteGaussianGeneratorImpl<NativeVector> EVDgg,
			const std::shared_ptr<const FPRLWEPrivateKeyImpl> sk
			) const;



	/**
	 * 
	 *
	 *
	 * 
	**/

	const shared_ptr<const FPRLWECiphertextImpl> DivideReduction(
		const std::shared_ptr<FPRLWECryptoParams> params,
		const std::shared_ptr<const FPRLWECiphertextImpl> ct) const;

	//	uint32_t shifter;

	//	vector<NativePoly> m_a_poly = m_a
	//}

  /**
   * Encrypts a bit using a secret key (symmetric key encryption)
   *
   * @param params a shared pointer to LWE scheme parameters
   * @param sk - the secret key
   * @param &m - the plaintext
   * @return a shared pointer to the ciphertext
   */
  const std::shared_ptr<const FPRLWECiphertextImpl> Encrypt(
      const std::shared_ptr<FPRLWECryptoParams> params,
      const std::shared_ptr<const FPRLWEPrivateKeyImpl> sk,
      DCRTPoly &m) const;

  /**
   * Decrypts the ciphertext using secret key sk
   *
   * @param params a shared pointer to LWE scheme parameters
   * @param sk the secret key
   * @param ct the ciphertext
   * @param *result plaintext result
   */
  void Decrypt(const std::shared_ptr<FPRLWECryptoParams> params,
               const std::shared_ptr<const FPRLWEPrivateKeyImpl> sk,
               const std::shared_ptr<const FPRLWECiphertextImpl> ct,
               DCRTPoly *result) const;

	const std::shared_ptr<const FPRLWECiphertextImpl> Product(	
		const std::shared_ptr<const FPRLWECiphertextImpl> ct1,
		const std::shared_ptr<const FPRLWECiphertextImpl> ct2,
		const std::shared_ptr<FPRLWECryptoParams> r_params,
    uint32_t EV_Bit,
	  uint32_t EV_len_rm,
       

		std::shared_ptr<FPRingEvaluationKey> EvalKey) const;


	std::shared_ptr<const FPLWECiphertextImpl> SampleExtract(
		const std::shared_ptr<FPRLWECryptoParams> params,
		const std::shared_ptr<const FPRLWECiphertextImpl> ct,
		uint32_t index_num) const;

	void SignedDigitDecompose(	const std::shared_ptr<FPRLWECryptoParams> params,
															const NativePoly &input,
															std::vector<NativePoly> *output) const;
  /**
   * Changes an LWE ciphertext modulo Q into an LWE ciphertext modulo q
   *
   * @param params a shared pointer to LWE scheme parameters
   * @param ctQ the input ciphertext
   * @return resulting ciphertext
   */
  std::shared_ptr<FPRLWECiphertextImpl> ModSwitch(
      const std::shared_ptr<FPRLWECryptoParams> params,
      const std::shared_ptr<const FPRLWECiphertextImpl> ctQ) const;

	};



/**
 * @brief Class that stores the LWE scheme switching key
 */
class FPRLWESwitchingKey : public Serializable {
 public:
  FPRLWESwitchingKey() {}

  explicit FPRLWESwitchingKey(
      const std::vector<std::vector<std::vector<FPRLWECiphertextImpl>>> &key)
      : m_key(key) {}

  explicit FPRLWESwitchingKey(const FPRLWESwitchingKey &rhs) {
    this->m_key = rhs.m_key;
  }


  explicit FPRLWESwitchingKey(const FPRLWESwitchingKey &&rhs) {
    this->m_key = std::move(rhs.m_key);
  }

  const FPRLWESwitchingKey &operator=(const FPRLWESwitchingKey &rhs) {
    this->m_key = rhs.m_key;
    return *this;
  }

  const FPRLWESwitchingKey &operator=(const FPRLWESwitchingKey &&rhs) {
    this->m_key = std::move(rhs.m_key);
    return *this;
  }

  const std::vector<std::vector<std::vector<FPRLWECiphertextImpl>>> &GetElements()
      const {
    return m_key;
  }

  void SetElements(
      const std::vector<std::vector<std::vector<FPRLWECiphertextImpl>>> &key) {
    m_key = key;
  }

  bool operator==(const FPRLWESwitchingKey &other) const {
    return m_key == other.m_key;
  }

  bool operator!=(const FPRLWESwitchingKey &other) const {
    return !(*this == other);
  }

  template <class Archive>
  void save(Archive &ar, std::uint32_t const version) const {
    ar(::cereal::make_nvp("k", m_key));
  }

  template <class Archive>
  void load(Archive &ar, std::uint32_t const version) {
    if (version > SerializedVersion()) {
      PALISADE_THROW(deserialize_error,
                     "serialized object version " + std::to_string(version) +
                         " is from a later version of the library");
    }

    ar(::cereal::make_nvp("k", m_key));
  }

  std::string SerializedObjectName() const { return "FPRLWEKeySwtichingKey"; }
  static uint32_t SerializedVersion() { return 1; }

 private:
	std::vector<std::vector<std::vector<FPRLWECiphertextImpl>>> m_key;
};

class FPOFUF : public Serializable {
 public:
  FPOFUF() {}

  explicit FPOFUF(
		const std::vector<std::vector<DCRTPoly>> &OFUF)
      : m_OFUF(OFUF) {}

  explicit FPOFUF(const FPOFUF &rhs) {
    this->m_OFUF = rhs.m_OFUF;
  }


  explicit FPOFUF(const FPOFUF &&rhs) {
    this->m_OFUF = std::move(rhs.m_OFUF);
  }

  const FPOFUF &operator=(const FPOFUF &rhs) {
    this->m_OFUF = rhs.m_OFUF;
    return *this;
  }

  const FPOFUF &operator=(const FPOFUF &&rhs) {
    this->m_OFUF = std::move(rhs.m_OFUF);
    return *this;
  }

  const std::vector<std::vector<DCRTPoly>> &GetElements()
      const {
    return m_OFUF;
  }

  void Adding(
			const std::vector<std::vector<DCRTPoly>> &OFUF) {
		//if (m_OFUF == nullptr) {
		if (m_OFUF.size() == 0) {
			m_OFUF.resize(OFUF.size());
			for (uint32_t i = 0; i < OFUF.size(); i++) {
				m_OFUF[i].resize(OFUF[i].size());
				for (uint32_t j = 0; j < OFUF[i].size(); j++) {
					m_OFUF[i][j] = OFUF[i][j]; 
				}
			}

		} else {
			for (uint32_t i = 0; i < OFUF.size(); i++) {
				for (uint32_t j = 0; j < OFUF[i].size(); j++) {
					m_OFUF[i][j] += OFUF[i][j]; 
				}
			}
		}
  }

  bool operator==(const FPOFUF &other) const {
    return m_OFUF == other.m_OFUF;
  }

  bool operator!=(const FPOFUF &other) const {
    return !(*this == other);
  }

  template <class Archive>
  void save(Archive &ar, std::uint32_t const version) const {
    ar(::cereal::make_nvp("OFUF", m_OFUF));
  }

  template <class Archive>
  void load(Archive &ar, std::uint32_t const version) {
    if (version > SerializedVersion()) {
      PALISADE_THROW(deserialize_error,
                     "serialized object version " + std::to_string(version) +
                         " is from a later version of the library");
    }

    ar(::cereal::make_nvp("OFUF", m_OFUF));
  }

  std::string SerializedObjectName() const { return "FPOFUF"; }
  static uint32_t SerializedVersion() { return 1; }

 private:
	std::vector<std::vector<DCRTPoly>> m_OFUF;
};






}  // namespace lbcrypto




#endif
