/*************
 * This file is modified from lwecore.h which is writen by Leo Ducas and Daniele Micciancio.
 * Full paper is listed in : eprint.iacr.org/2022/186.
 * This source is obeyed and applied following copyright law. 
 *
 * If you have any question, please contact us by the email kr3951@hanyang.ac.kr
 * Seunghwan Lee, 2022.03.04
 * *************/


// @file lwecore.h - Main Classes for Boolean circuit FHE.
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

#ifndef FPFHE_LWECORE_H
#define FPFHE_LWECORE_H

#include <string>
#include <utility>
#include <vector>

#include "math/backend.h"
#include "math/discretegaussiangenerator.h"
#include "utils/serializable.h"

namespace lbcrypto {

typedef int64_t FPLWEPlaintext;
typedef double FPLWEPlaintext_double;

/**
 * @brief Class that stores all parameters for the LWE scheme
 */
class FPLWECryptoParams : public Serializable {
 public:
  FPLWECryptoParams() : m_n(0), m_q_bit(0), m_PKS_base_bit(0){}

  /**
   * Main constructor for LWECryptoParams
   *
   * @param n lattice parameter for additive LWE scheme
   * @param N ring dimension for RingGSW/RLWE used in bootstrapping
   * @param &q modulus for additive LWE
   * @param &Q modulus for RingGSW/RLWE used in bootstrapping
   * @param std standard deviation
   * @param baseKS the base used for key switching
   */
  explicit FPLWECryptoParams(
			uint32_t n,
			uint32_t q_bit,
      double std, 
			uint32_t baseKS_bit, 
			uint32_t rm_len,
			uint32_t msg_bit)
      : m_n(n),  m_q_bit(q_bit), m_PKS_base_bit(baseKS_bit), m_PKS_l_rm(rm_len), m_msg_bit(msg_bit) {
    m_dgg.SetStd(std);

		/*
    if (Q.GetMSB() > MAX_MODULUS_SIZE) {
      std::string errMsg =
          "ERROR: Maximum size of Q supported for FHEW is 60 bits.";
      PALISADE_THROW(config_error, errMsg);
    }
		*/
    PreCompute();
  }

  /**
   * Performs precomputations based on the supplied parameters
   */
  void PreCompute() {
		m_q =1;
		m_q <<= m_q_bit;
		m_scaling = 1;
		m_scaling <<= (m_q_bit - m_msg_bit);
		// Number of digits in representing numbers mod Q
     m_PKS_l = (uint32_t)std::ceil((static_cast<double> (m_q_bit)) /
				 (static_cast<double> (m_PKS_base_bit)) );
			
		 m_PKS_l_save = m_PKS_l - m_PKS_l_rm;
 
	}

  explicit FPLWECryptoParams(const FPLWECryptoParams &rhs) {
    this->m_n = rhs.m_n;
		this->m_q_bit = rhs.m_q_bit;
    this->m_PKS_base_bit = rhs.m_PKS_base_bit;
    //this->m_digitsKS = rhs.m_digitsKS;
    this->m_dgg.SetStd(rhs.m_dgg.GetStd());
		this->m_PKS_l_rm = rhs.m_PKS_l_rm;
  }

  explicit FPLWECryptoParams(const FPLWECryptoParams &&rhs) {
    this->m_n = std::move(rhs.m_n);
		this->m_q_bit = std::move(rhs.m_q_bit);
    this->m_PKS_base_bit = std::move(rhs.m_PKS_base_bit);
		this->m_dgg.SetStd(rhs.m_dgg.GetStd());
		this->m_PKS_l_rm = std::move(rhs.m_PKS_l_rm);
  }

  const FPLWECryptoParams &operator=(const FPLWECryptoParams &rhs) {
    this->m_n = rhs.m_n;
		this->m_q_bit = rhs.m_q_bit;
    this->m_PKS_base_bit = rhs.m_PKS_base_bit;
    //this->m_digitsKS = rhs.m_digitsKS;
    this->m_dgg.SetStd(rhs.m_dgg.GetStd());
		//this->m_scaling = rhs.m_scaling;
    return *this;
  }

  const FPLWECryptoParams &operator=(const FPLWECryptoParams &&rhs) {
    this->m_n = std::move(rhs.m_n);
		this->m_q_bit = std::move(rhs.m_q_bit);
    this->m_PKS_base_bit = std::move(rhs.m_PKS_base_bit);
    //this->m_digitsKS = std::move(rhs.m_digitsKS);
    this->m_dgg.SetStd(rhs.m_dgg.GetStd());

		this->m_PKS_l_rm = std::move(rhs.m_PKS_l_rm);
    return *this;
  }

  uint32_t Getn() const { return m_n; }
  
	const uint32_t Getq_bit() const { return m_q_bit; }
	const uint64_t Getq() const { return m_q; }


  uint32_t GetBasePKSBit() const { return m_PKS_base_bit; }
  uint32_t GetPKSLen() const { return m_PKS_l; }
	uint32_t GetPKSLenSave() const { return m_PKS_l_save; }
	uint32_t GetPKSLenRm() const { return m_PKS_l_rm; }


	const DiscreteGaussianGeneratorImpl<NativeVector> &GetDgg() const {
    return m_dgg;
  }

  bool operator==(const FPLWECryptoParams &other) const {
    return m_n == other.m_n &&
					 m_q_bit == other.m_q_bit && m_dgg.GetStd() == other.m_dgg.GetStd() &&
           m_PKS_base_bit == other.m_PKS_base_bit;
				
  }

  bool operator!=(const FPLWECryptoParams &other) const {
    return !(*this == other);
  }

  template <class Archive>
  void save(Archive &ar, std::uint32_t const version) const {
    ar(::cereal::make_nvp("n", m_n));
    ar(::cereal::make_nvp("q_bit", m_q_bit));
    ar(::cereal::make_nvp("sigma", m_dgg.GetStd()));
    ar(::cereal::make_nvp("bKS_bit", m_PKS_base_bit));
		//ar(::cereal::make_nvp("scaling", m_scaling));
  }

  template <class Archive>
  void load(Archive &ar, std::uint32_t const version) {
    if (version > SerializedVersion()) {
      PALISADE_THROW(deserialize_error,
                     "serialized object version " + std::to_string(version) +
                         " is from a later version of the library");
    }

    ar(::cereal::make_nvp("n", m_n));
		 ar(::cereal::make_nvp("q_bit", m_q_bit));
		//ar(::cereal::make_nvp("scaling", m_scaling));
    double sigma;
    ar(::cereal::make_nvp("sigma", sigma));
    this->m_dgg.SetStd(sigma);
    ar(::cereal::make_nvp("bKS", m_PKS_base_bit));

    this->PreCompute();
  }

  std::string SerializedObjectName() const { return "FPLWECryptoParams"; }
  static uint32_t SerializedVersion() { return 1; }
	
	
	const uint32_t GetScalingFactor() {
		return this->m_scaling;
	}


 private:
  // lattice parameter for the additive LWE scheme
  uint32_t m_n;
	// modulus for the RingGSW/RingLWE scheme
  uint32_t m_q_bit;
  uint64_t m_q;


  // Error distribution generator
  DiscreteGaussianGeneratorImpl<NativeVector> m_dgg;
  // Base used in key switching
	// Powers of m_baseKS
	uint32_t m_PKS_base_bit;
	uint32_t m_PKS_l;
	uint32_t m_PKS_l_rm;
	uint32_t m_PKS_l_save;


  // Powers of m_baseKS
	
	// msg_bit
	uint32_t m_msg_bit;
	uint32_t m_scaling;
};

/**
 * @brief Class that stores a LWE scheme ciphertext; composed of a vector "a"
 * and integer "b"
 */


class FPLWECiphertextImpl : public Serializable {
 public:
  FPLWECiphertextImpl() {}

  explicit FPLWECiphertextImpl(const NativeVector &&a, const NativeInteger &b)
      : m_a(std::move(a)), m_b(b) {}

  explicit FPLWECiphertextImpl(const NativeVector &a, const NativeInteger &b)
      : m_a(a), m_b(b) {}

  FPLWECiphertextImpl(NativeVector &&a, NativeInteger b)
      : m_a(std::move(a)), m_b(b) {}

  explicit FPLWECiphertextImpl(const FPLWECiphertextImpl &rhs) {
    this->m_a = rhs.m_a;
    this->m_b = rhs.m_b;
  }

  explicit FPLWECiphertextImpl(const FPLWECiphertextImpl &&rhs) {
    this->m_a = std::move(rhs.m_a);
    this->m_b = std::move(rhs.m_b);
  }

  const FPLWECiphertextImpl &operator=(const FPLWECiphertextImpl &rhs) {
    this->m_a = rhs.m_a;
    this->m_b = rhs.m_b;
    return *this;
  }

  const FPLWECiphertextImpl &operator=(const FPLWECiphertextImpl &&rhs) {
    this->m_a = std::move(rhs.m_a);
    this->m_b = std::move(rhs.m_b);
    return *this;
  }

  const NativeVector &GetA() const { return m_a; }

  const NativeInteger &GetA(std::size_t i) const { return m_a[i]; }

  const NativeInteger &GetB() const { return m_b; }

  void SetA(const NativeVector &a) { m_a = a; }

  void SetB(const NativeInteger &b) { m_b = b; }

  bool operator==(const FPLWECiphertextImpl &other) const {
    return m_a == other.m_a && m_b == other.m_b;
  }

  bool operator!=(const FPLWECiphertextImpl &other) const {
    return !(*this == other);
  }

  template <class Archive>
  void save(Archive &ar, std::uint32_t const version) const {
    ar(::cereal::make_nvp("a", m_a));
    ar(::cereal::make_nvp("b", m_b));
  }

  template <class Archive>
  void load(Archive &ar, std::uint32_t const version) {
    if (version > SerializedVersion()) {
      PALISADE_THROW(deserialize_error,
                     "serialized object version " + std::to_string(version) +
                         " is from a later version of the library");
    }

    ar(::cereal::make_nvp("a", m_a));
    ar(::cereal::make_nvp("b", m_b));
  }

  std::string SerializedObjectName() const { return "FPLWECiphertext"; }
  static uint32_t SerializedVersion() { return 1; }

 private:
  NativeVector m_a;
  NativeInteger m_b;
};



/**
 * @brief Class that stores the LWE scheme secret key; contains a vector
 */


class FPLWEPrivateKeyImpl : public Serializable {
 public:
  FPLWEPrivateKeyImpl() {}


  explicit FPLWEPrivateKeyImpl(const NativeVector &s) : m_s(s) {}
  explicit FPLWEPrivateKeyImpl(const FPLWEPrivateKeyImpl &rhs) {
    this->m_s = rhs.m_s;
  }
  explicit FPLWEPrivateKeyImpl(const FPLWEPrivateKeyImpl &&rhs) {
    this->m_s = std::move(rhs.m_s);
  }
	

  const FPLWEPrivateKeyImpl &operator=(const FPLWEPrivateKeyImpl &rhs) {
    this->m_s = rhs.m_s;
    return *this;
  }

  const FPLWEPrivateKeyImpl &operator=(const FPLWEPrivateKeyImpl &&rhs) {
    this->m_s = std::move(rhs.m_s);
    return *this;
  }
  const NativeVector &GetElement() const { return m_s; }

  void SetElement(const NativeVector &s) { m_s = s; }

  bool operator==(const FPLWEPrivateKeyImpl &other) const {
    return m_s == other.m_s;
  }

  bool operator!=(const FPLWEPrivateKeyImpl &other) const {
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
  NativeVector m_s;
 };



/**
 * @brief Class that stores the LWE scheme switching key
 */
class FPLWESwitchingKey : public Serializable {
 public:
  FPLWESwitchingKey() {}

	explicit FPLWESwitchingKey(uint32_t dim0, uint32_t dim1, uint32_t dim2) {
		m_key.resize(dim0);
		for (uint32_t i = 0; i < dim0; i++) {
			m_key[i].resize(dim1);
			for (uint32_t j = 0; j < dim1; j++) {
				m_key[i][j].resize(dim2);
			}
		}
	}

  explicit FPLWESwitchingKey(
      const std::vector<std::vector<std::vector<FPLWECiphertextImpl>>> &key)
      : m_key(key) {}

  explicit FPLWESwitchingKey(const FPLWESwitchingKey &rhs) {
    this->m_key = rhs.m_key;
  }

  explicit FPLWESwitchingKey(const FPLWESwitchingKey &&rhs) {
    this->m_key = std::move(rhs.m_key);
	}

  const FPLWESwitchingKey &operator=(const FPLWESwitchingKey &rhs) {
    this->m_key = rhs.m_key;
    return *this;
  }

  const FPLWESwitchingKey &operator=(const FPLWESwitchingKey &&rhs) {
    this->m_key = std::move(rhs.m_key);
    return *this;
  }

  const std::vector<std::vector<std::vector<FPLWECiphertextImpl>>> &GetElements() const {
    return m_key;
  }

  void SetElements(
    const std::vector<std::vector<std::vector<FPLWECiphertextImpl>>> &key) {
    m_key = key;
  }

	std::vector<std::vector<FPLWECiphertextImpl>>& operator[](uint32_t i) {
		return m_key[i];
	}
 
	const std::vector<std::vector<FPLWECiphertextImpl>>& operator[](usint i) const {
		return m_key[i];
	}


  bool operator==(const FPLWESwitchingKey &other) const {
    return m_key == other.m_key;
  }

  bool operator!=(const FPLWESwitchingKey &other) const {
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

  std::string SerializedObjectName() const { return "FPLWEPrivateKey"; }
  static uint32_t SerializedVersion() { return 1; }

 private:
	std::vector<std::vector<std::vector<FPLWECiphertextImpl>>> m_key;
};

}  // namespace lbcrypto

#endif
