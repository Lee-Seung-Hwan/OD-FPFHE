/*************
 * This file is modified from lwe.cpp which is writen by Leo Ducas and Daniele Micciancio.
 * Full paper is listed in : eprint.iacr.org/2022/186.
 * This source is obeyed and applied following copyright law. 
 *
 * If you have any question, please contact us by the email kr3951@hanyang.ac.kr
 * Seunghwan Lee, 2022.03.04
 * *************/


// @file lwe.cpp - LWE Encryption Scheme implementation as described in
// https://eprint.iacr.org/2014/816 Full reference:
// @misc{cryptoeprint:2014:816,
//   author = {Leo Ducas and Daniele Micciancio},
//   title = {FHEW: Bootstrapping Homomorphic Encryption in less than a second},
//   howpublished = {Cryptology ePrint Archive, Report 2014/816},
//   year = {2014},
//   note = {\url{https://eprint.iacr.org/2014/816}},
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

#include "fplwe.h"
#include "fprlwe.h"


#include "math/binaryuniformgenerator.h"
#include "math/discreteuniformgenerator.h"
#include "math/ternaryuniformgenerator.h"

namespace lbcrypto {

std::shared_ptr<FPLWEPrivateKeyImpl> FPLWEEncryptionScheme::KeyGen(
	const std::shared_ptr<FPLWECryptoParams> params
	, uint32_t h_vals = 0
	) const {
  TernaryUniformGeneratorImpl<NativeVector> tug;
	if (h_vals == 0) {
		return std::make_shared<FPLWEPrivateKeyImpl>(  
				FPLWEPrivateKeyImpl(tug.GenerateVector(params->Getn(), params->Getq())));
	} else {
		return std::make_shared<FPLWEPrivateKeyImpl>(  
				FPLWEPrivateKeyImpl(tug.GenerateVector(params->Getn(), params->Getq(), h_vals)));
	
	}
}

// classical LWE encryption
// a is a randomly uniform vector of dimension n; with integers mod q
// b = a*s + e + m floor(q/4) is an integer mod q

std::shared_ptr<FPLWECiphertextImpl> FPLWEEncryptionScheme::Encrypt(
    const std::shared_ptr<FPLWECryptoParams> params,
    const std::shared_ptr<const FPLWEPrivateKeyImpl> sk,
    const FPLWEPlaintext &m) const {


  NativeInteger q = sk->GetElement().GetModulus();
  uint32_t n = sk->GetElement().GetLength();

  NativeInteger b =  m * NativeInteger(params->GetScalingFactor()) + params->GetDgg().GenerateInteger(q);
 	
	DiscreteUniformGeneratorImpl<NativeVector> dug;
  dug.SetModulus(q);
  NativeVector a = dug.GenerateVector(n);

  NativeInteger mu = q.ComputeMu();

	const NativeVector &s = sk->GetElement();
  for (uint32_t i = 0; i < n; ++i) {
		b += a[i].ModMulFast(s[i], q, mu);
  }
  b.ModEq(q);
  
	return std::make_shared<FPLWECiphertextImpl>(FPLWECiphertextImpl(a, b));
}



	std::shared_ptr<FPLWECiphertextImpl> FPLWEEncryptionScheme::Encrypt_WO_scaling(
    const std::shared_ptr<FPLWECryptoParams> params,
    const std::shared_ptr<const FPLWEPrivateKeyImpl> sk,
    int64_t m) const {


  NativeInteger q = sk->GetElement().GetModulus();
  uint32_t n = sk->GetElement().GetLength();

	// Roughly ...
	if (m < 0 ) { m = q.ConvertToInt()*2 + m;}   
	m = m % q.ConvertToInt();

  NativeInteger b = NativeInteger(m) + params->GetDgg().GenerateInteger(q);
 	
	DiscreteUniformGeneratorImpl<NativeVector> dug;
  dug.SetModulus(q);
  NativeVector a = dug.GenerateVector(n);

  NativeInteger mu = q.ComputeMu();

	const NativeVector &s = sk->GetElement();
  for (uint32_t i = 0; i < n; ++i) {
		b += a[i].ModMulFast(s[i], q, mu);
  }
  b.ModEq(q);
	//std::shared_ptr<FPLWECiphertextImpl> CT = std::make_shared<FPLWECiphertextImpl>(FPLWECiphertextImpl(a,b)); 
	return std::make_shared<FPLWECiphertextImpl>(FPLWECiphertextImpl(a, b));
	//return CT;


}



// classical LWE decryption
// m_result = Round(4/q * (b - a*s))

void FPLWEEncryptionScheme::Decrypt(
    const std::shared_ptr<FPLWECryptoParams> params,
    const std::shared_ptr<const FPLWEPrivateKeyImpl> sk,
    const std::shared_ptr<const FPLWECiphertextImpl> ct,
    FPLWEPlaintext_double *result) const {
  
	// TODO in the future we should add a check to make sure sk parameters match
  // the ct parameters

  // Create local variables to speed up the computations
  NativeVector a = ct->GetA();
  uint32_t n = sk->GetElement().GetLength();
  NativeVector s = sk->GetElement();
  NativeInteger q = sk->GetElement().GetModulus();

  NativeInteger mu = q.ComputeMu();

  NativeInteger inner(0);
  for (uint32_t i = 0; i < n; ++i) {
    inner += a[i].ModMulFast(s[i], q, mu);
		//std::cout << "a[" << i << "] is " << a[i] << ", and s[" << i << "] is " << s[i] << std::endl;
  }
  inner.ModEq(q);

	//std::cout << "inner result is " << inner << std::endl;
  NativeInteger r = ct->GetB();

  r.ModSubFastEq(inner, q);

  std::cout << "q is " << q <<  ", LWE Decrypting.... b - as result is " << r << std::endl;
	
  //*result = (FPLWEPlaintext_double) ( r / NativeInteger(params->GetScalingFactor()) ).ConvertToDouble();
  *result = (FPLWEPlaintext_double) ( r.ConvertToDouble() /  ((double) params ->GetScalingFactor()));

#if defined(FPFHE_DEBUG)
  double error = (4.0 * (r.ConvertToDouble() - q.ConvertToInt() / 8)) /
                     q.ConvertToDouble() -
                 static_cast<double>(*result);
  std::cerr << "error:\t" << error << std::endl;
#endif

  return;
}

// the main rounding operation used in ModSwitch (as described in Section 3 of
// https://eprint.iacr.org/2014/816) The idea is that Round(x) = 0.5 + Floor(x)
NativeInteger RoundqQ(const NativeInteger &v, const NativeInteger &q,
                      const NativeInteger &Q) {
  return NativeInteger((uint64_t)std::floor(0.5 + v.ConvertToDouble() *
                                                      q.ConvertToDouble() /
                                                      Q.ConvertToDouble())).Mod(q);
}

// Modulus switching - directly applies the scale-and-round operation RoundQ
/*
std::shared_ptr<FPLWECiphertextImpl> FPLWEEncryptionScheme::ModSwitch(
    const std::shared_ptr<FPLWECryptoParams> params,
    const std::shared_ptr<const FPLWECiphertextImpl> ctQ) const {
  NativeVector a(params->Getn(), params->Getq());

  uint32_t n = params->Getn();
  NativeInteger q = params->Getq();
  NativeInteger Q = params->GetQ();

  for (uint32_t i = 0; i < n; ++i) a[i] = RoundqQ(ctQ->GetA()[i], q, Q);

  NativeInteger b = RoundqQ(ctQ->GetB(), q, Q);

  return std::make_shared<FPLWECiphertextImpl>(FPLWECiphertextImpl(a, b));
}
*/
// Switching key as described in Section 3 of https://eprint.iacr.org/2014/816

/*
std::shared_ptr<FPLWESwitchingKey> FPLWEEncryptionScheme::KeySwitchGen(
    const std::shared_ptr<FPLWECryptoParams> params,
    const std::shared_ptr<const FPLWEPrivateKeyImpl> sk,
    const std::shared_ptr<const FPLWEPrivateKeyImpl> skN) const {
  // Create local copies of main variables
  uint32_t n = params->Getn();
  uint32_t N = params->GetN();
  NativeInteger Q = params->GetQ();
  uint32_t baseKS = params->GetBaseKS();
  std::vector<NativeInteger> digitsKS = params->GetDigitsKS();
  uint32_t expKS = digitsKS.size();

  // newSK stores negative values using modulus q
  // we need to switch to modulus Q
  NativeVector newSK = sk->GetElement();
  newSK.SwitchModulus(Q);

  NativeVector oldSK = skN->GetElement();

  DiscreteUniformGeneratorImpl<NativeVector> dug;
  dug.SetModulus(Q);

  NativeInteger mu = Q.ComputeMu();

  std::vector<std::vector<std::vector<FPLWECiphertextImpl>>> resultVec(N);

#pragma omp parallel for
  for (uint32_t i = 0; i < N; ++i) {
    std::vector<std::vector<FPLWECiphertextImpl>> vector1(baseKS);
    for (uint32_t j = 0; j < baseKS; ++j) {
      std::vector<FPLWECiphertextImpl> vector2(expKS);
      for (uint32_t k = 0; k < expKS; ++k) {
        NativeInteger b = (params->GetDgg().GenerateInteger(Q))
                              .ModAdd(oldSK[i].ModMul(j * digitsKS[k], Q), Q);

        NativeVector a = dug.GenerateVector(n);

#if NATIVEINT == 32
        for (uint32_t i = 0; i < n; ++i) {
          b.ModAddFastEq(a[i].ModMulFast(newSK[i], Q, mu), Q);
        }
#else
        for (uint32_t i = 0; i < n; ++i) {
          b += a[i].ModMulFast(newSK[i], Q, mu);
        }
        b.ModEq(Q);
#endif

        vector2[k] = FPLWECiphertextImpl(a, b);
      }
      vector1[j] = std::move(vector2);
    }
    resultVec[i] = std::move(vector1);
  }

  return std::make_shared<FPLWESwitchingKey>(FPLWESwitchingKey(resultVec));
	}
*/


// the key switching operation as described in Section 3 of
// https://eprint.iacr.org/2014/816

/*
std::shared_ptr<FPLWECiphertextImpl> FPLWEEncryptionScheme::KeySwitch(
    const std::shared_ptr<FPLWECryptoParams> params,
    const std::shared_ptr<FPLWESwitchingKey> K,
    const std::shared_ptr<const FPLWECiphertextImpl> ctQN) const {
  uint32_t n = params->Getn();
  uint32_t N = params->GetN();
  NativeInteger Q = params->GetQ();
  uint32_t baseKS = params->GetBaseKS();
  std::vector<NativeInteger> digitsKS = params->GetDigitsKS();
  uint32_t expKS = digitsKS.size();

  // creates an empty vector
  NativeVector a(n, Q);
  NativeInteger b = ctQN->GetB();
  NativeVector aOld = ctQN->GetA();

  for (uint32_t i = 0; i < N; ++i) {
    NativeInteger atmp = aOld[i];
    for (uint32_t j = 0; j < expKS; ++j, atmp /= baseKS) {
      uint64_t a0 = (atmp % baseKS).ConvertToInt();
      for (uint32_t k = 0; k < n; ++k)
        a[k].ModSubFastEq((K->GetElements()[i][a0][j]).GetA()[k], Q);
      b.ModSubFastEq((K->GetElements()[i][a0][j]).GetB(), Q);
    }
  }

  return std::make_shared<FPLWECiphertextImpl>(FPLWECiphertextImpl(a, b));
}
*/




// Switching key as described in Section 3 of https://eprint.iacr.org/2014/816
//std::shared_ptr<FPRingPackingKey> FPLWEEncryptionScheme::PackingGen(

/*
void FPLWEEncryptionScheme::PackingGen(
    const std::shared_ptr<FPRingGSWCryptoParams> params,
    const std::shared_ptr<const FPLWEPrivateKeyImpl> sk,
    const std::shared_ptr<const FPLWEPrivateKeyImpl> skPackingN) const {




	
  // Create local copies of main variables
	


	// LWE params 
	std::shared_ptr<FPLWECryptoParams> LWEparams = params->GetLWEParams();
	

	// I don't know when I use it
	//uint32_t n = LWEparams->Getn();
  //uint32_t N = LWEparams->GetpackingN();
	//std::cout<< "n is " << n << std::endl << "N is" << N << std::endl;
	
	NativeInteger Q = LWEparams->GetQ();
  //uint32_t baseKS = params->GetBaseKS();
 
	//uint32_t basePK_digit = 2;
	//std::vector<NativeInteger> digitsKS = params->GetDigitsKS();
  //uint32_t expKS = digitsKS.size();
  // newSK stores negative values using modulus q
  // we need to switch to modulus Q
  NativeVector oldLWESK = sk->GetElement();
  oldLWESK.SwitchModulus(Q);
  
	// sk to Poly
	NativePoly skNTT = NativePoly(params->GetPolyParams());
	skNTT.SetValues(skPackingN->GetElement(), Format::COEFFICIENT);
	skNTT.SetFormat(Format::EVALUATION);

	// I don't know
	//NativeVector newRingSK = skPackingN->GetElement();
  
	DiscreteUniformGeneratorImpl<NativeVector> dug;
  dug.SetModulus(Q);

	// Mu ?
  NativeInteger mu = Q.ComputeMu();
  
	std::vector<std::vector<std::vector<FPRingCiphertext>>> resultVec(N);

	
	// tempA is introduced to minimize the number of NTTs
  std::vector<NativePoly> tempA(digitsG2);

	auto result = std::make_shared<FPRingCiphertext>(digitsG2, 2); 

  for (uint32_t i = 0; i < digitsG2; ++i) {
			// 연산자 오버로딩
     (*result)[i][0] = NativePoly(dug, polyParams, Format::COEFFICIENT);
     tempA[i] = (*result)[i][0];
     (*result)[i][1] = NativePoly(params->GetLWEParams()->GetDgg(), polyParams,
                                  Format::COEFFICIENT);
   }
	// Gadget G 넣는 부분! 그니까 필요 없음
  // for (uint32_t i = 0; i < digitsG; ++i) {
  //   if (m > 0) {       // Add G Multiple
  //     (*result)[2 * i][0][0].ModAddEq(params->GetGPower()[i], Q);
  //     // [a,as+e] + G
  //     (*result)[2 * i + 1][1][0].ModAddEq(params->GetGPower()[i], Q);
  //   }
  // }
 

   // 3*digitsG2 NTTs are called
   result->SetFormat(Format::EVALUATION);
   for (uint32_t i = 0; i < digitsG2; ++i) {
     tempA[i].SetFormat(Format::EVALUATION);
     (*result)[i][1] += tempA[i];  // 곱해야됨 * skNTT;
   }
 



	// First is N, Second is decomp 2, Final is 
#pragma omp parallel for
  for (uint32_t i = 0; i < N; ++i) {
    std::vector<std::vector<FPLWECiphertextImpl>> vector1(baseKS);
    for (uint32_t j = 0; j < baseKS; ++j) {
      std::vector<FPLWECiphertextImpl> vector2(expKS);
      for (uint32_t k = 0; k < expKS; ++k) {
				// 여기서 n으로 뭐 하나 해야됨:
        NativeInteger b = (params->GetDgg().GenerateInteger(Q))
                              .ModAdd(oldSK[i].ModMul(j * digitsKS[k], Q), Q);

        NativeVector a = dug.GenerateVector(n);

#if NATIVEINT == 32
        for (uint32_t i = 0; i < n; ++i) {
          b.ModAddFastEq(a[i].ModMulFast(newSK[i], Q, mu), Q);
        }
#else
        for (uint32_t i = 0; i < n; ++i) {
          b += a[i].ModMulFast(newSK[i], Q, mu);
        }
        b.ModEq(Q);
#endif

        vector2[k] = FPLWECiphertextImpl(a, b);
      }
      vector1[j] = std::move(vector2);
    }
    resultVec[i] = std::move(vector1);
  }

  return std::make_shared<FPRingCiphertext>(FPRingPackingKey(resultVec));
	
}
*/

/*
// the key switching operation as described in Section 3 of
// https://eprint.iacr.org/2014/816
std::shared_ptr<FPRingCiphertextImpl> FPLWEEncryptionScheme::Packing(
    const std::shared_ptr<FPLWECryptoParams> params,
    const std::shared_ptr<FPLWESwitchingKey> K,
    const std::shared_ptr<const FPLWECiphertextImpl> ctQN) const {
  uint32_t n = params->Getn();
  uint32_t N = params->GetN();
  NativeInteger Q = params->GetQ();
  uint32_t baseKS = params->GetBaseKS();
  std::vector<NativeInteger> digitsKS = params->GetDigitsKS();
  uint32_t expKS = digitsKS.size();

  // creates an empty vector
  NativeVector a(n, Q);
  NativeInteger b = ctQN->GetB();
  NativeVector aOld = ctQN->GetA();

  for (uint32_t i = 0; i < N; ++i) {
    NativeInteger atmp = aOld[i];
    for (uint32_t j = 0; j < expKS; ++j, atmp /= baseKS) {
      uint64_t a0 = (atmp % baseKS).ConvertToInt();
      for (uint32_t k = 0; k < n; ++k)
        a[k].ModSubFastEq((K->GetElements()[i][a0][j]).GetA()[k], Q);
      b.ModSubFastEq((K->GetElements()[i][a0][j]).GetB(), Q);
    }
  }

  return std::make_shared<FPLWECiphertextImpl>(FPLWECiphertextImpl(a, b));
}
*/


};  // namespace lbcrypto
