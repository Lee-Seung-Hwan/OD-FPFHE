/*************
 * This file is modified from binfhecontext-ser.h which is writen by Leo Ducas and Daniele Micciancio.
 * Full paper is listed in : eprint.iacr.org/2022/186.
 * This source is obeyed and applied following copyright law. 
 *
 * If you have any question, please contact us by the email kr3951@hanyang.ac.kr
 * Seunghwan Lee, 2022.03.04
 * *************/



// @file bindhecontext-ser.h - Header file adding serialization support to
// Boolean circuit FHE
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

#ifndef FPFHE_FHECONTEXT_SER_H
#define FPFHE_FHECONTEXT_SER_H

#include "fpfhecontext.h"
#include "utils/serial.h"

// Registers types needed for serialization
CEREAL_REGISTER_TYPE(lbcrypto::FPLWECryptoParams);
CEREAL_REGISTER_TYPE(lbcrypto::FPLWECiphertextImpl);
CEREAL_REGISTER_TYPE(lbcrypto::FPLWEPrivateKeyImpl);
CEREAL_REGISTER_TYPE(lbcrypto::FPLWESwitchingKey);
CEREAL_REGISTER_TYPE(lbcrypto::FPRingGSWCryptoParams);
CEREAL_REGISTER_TYPE(lbcrypto::FPRingGSWCiphertext);
CEREAL_REGISTER_TYPE(lbcrypto::FPRingGSWBTKey);
CEREAL_REGISTER_TYPE(lbcrypto::FPBinFHEContext);

#endif
