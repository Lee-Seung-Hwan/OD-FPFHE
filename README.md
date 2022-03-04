# OD-FPFHE
Implementation of "Overflow-detectable Floating-point Fully Homomorphic Encryption".


We open all source code used in archive paper:
    Overflow-detectable Floating-point Fully Homomorphic Encryption, https://eprint.iacr.org/2022/186.

Note that our sources are modified from binfhe in PALISADE, which is written by Leo Ducas and Daniele Micciancio.

(We greatly thank them and those involved in PALISADE.)

## Installation
1. These source require PALISADE library. So, please clone PALISADE git: https://gitlab.com/palisade/palisade-release
2. We should access and modify structure NativeInteger modulo Q_0 and Q_1 in DCRTPoly (We can't find any method to modify it). So put following codes anywhere in PALISADE src/core/include/lattice/dcrtpoly.h 

    <pre><code>
    PolyType &GetElementW(uint32_t idx) { return m_vectors[idx];}; // (0)
    void  SetFormatAtIndexW(Format format, uint32_t idx) { m_vectors[idx].SetFormat(format);}; // (1)
    void  SetOutFormatW(Format format) { m_format = format;}; // (2)
    uint64_t *GetElementWForHAXL(uint32_t idx) { return reinterpret_cast<uint64_t  *> (&m_vectors[idx][0]); }; // (3)
    // The (0) must be put in dcrtpoly.h, but I can't remember (1), (2), and (3) were required. It would be deleted next version.
    </code></pre>
   
3. Copy our fpfhe folder in OD-FPFHE-src/fpfhe to PALISADE src/fpfhe
4. Copy our CMakeList.txt to PALISADE/CMakeList.txt 
    (If you already installed particular PALISADE, you can copy all script in our CMakeList.txt by searching keyword 'fpfhe'.)

5. Compile PALISADE with documentation in https://palisade-crypto.org/documentation/.

That's it :)

## How to run example codes
1. If you install palisade, then go to example and make folder 'build'
2. compile test codes by using command 'cmake ..' and 'make'.
3. Run TestCorrect, TestAmortize, TestParallel

## Explanation of test codes
If you check inside of codes, you can adjust various parameter STD128D,..., and STD192S which is explained in our paper.

We set extra parmeters TOY and it is useful to develop specific application (Note that this parameter totally insecure but fast).

The code "testCorrect" checks correctness about homomorphic addition, substraction, and multiplification.

In addition, this checks whether overflow occurs or not.

Other codes "testAmortize" and "testParallel" are used to check time comsumption while a lot of operations are run.


## Announcement
Note that some of issues are occured while developping OD-FPFHE with PALISADE.

### Issue 1. Using AVX512 and Turning the flag WITH_INTEL_HEXL=ON
For now (22.03.04), we notice that some of errors are occur while installing PALISADE with the flag WITH_INTEL_HEXL=ON as follows:

    CMake Error at cpu-features-download/cpu_features-prefix/tmp/cpu_features-gitclone.cmake:40 (message):
        Failed to checkout tag: 'master'

We think this happens due to modified git-tag from master to main google/cpu_features repository.
So, please run following lines when above errors are occurs.

1. cd <your current compile folder e.g. build>/ext_intel_hexl/src/ext_intel_hexl-build/cmake/thid-party/cpu-features/cpu-featurs-download/
2. sudo vi CMakeList.txt and replace git-tage 'master' to 'main'
3. sudo rm CMakeCache.txt
4. sudo make
5. Go back to your compile folder and run sudo make again.

### Issue 2. When WITH-INTEL-HEXL=ON is not used ...

We notice that if ring degree of both GSW ciphertext and MLWE ciphertext is not equal, the result of NTT transfrom of thoes is not correct. 

In more detail, we run following lines like implemented CKKS in PALISADE as follows:

    /* GSW precomputing for NTT transform */
    rootsQ1[i] = RootOfUnity<NativeInteger>("Cyclotomic order of GSW ciphertext", modulus Q_i);
    ChineseRemainderTransformFTT<Nativevector>::PreCompute(rootsQ1, "Cycltomic order of GSW ciphertext", "array of modulus Q_i");

    /* MLWE precomputing for NTT transform */
    rootsQ2[i] = RootOfUnity<NativeInteger>("Cyclotomic order of MLWE ciphertext", modulus Q_i);
    ChineseRemainderTransformFTT<Nativevector>::PreCompute(rootsQ2, "Cycltomic order of MLWE ciphertext", "array of modulus Q_i");

Note that CKKS requires precomputed CRT matrix with different modulo Q_i, but not the ring dimension, so this is not happen.

Howver, we use that function twice with different ring dimension, but then some results of NTT transform are not correct.

So, we have to make ring diemension of both GSW and MLWE the same. (Likewise, K of MLWE should be equal to that of GSW to achive target security level.)

If you turn off the WITH-INTEL-HEXL, parmeters of MLWE becomes that of GSW automatically.

However If WITH-INTEL-HEXL is ON, every NTT transform performs with extra library and correct as well.

Therefore, the parameter of MLWE ciphertext becomes different from that of GSW as listed in our paper when the HEXL flag is ON.

### Issue 3. How to use WITH-TCM=ON 

Somebody help me to use WITH-TCM=ON .......

### Issue 4. Using it on Mac

If memory of laptop computer is 16GB, then it is enough to run OD-FPFHE with parameters STD128D or STD128S.

We tested in both Ubuntu and Mac equipped with intel i7 CPU.

However if Mac with M1 processer is given ... we can't find any solution for running it.

Twe errors are occurs as follows:
1. The parameter NATIVE_SIZE does not work.
2. The library OPENMP cannot be found although it is installed.

## Future work

We agree that this version of OD-FPFHE have a super-slow time while homomorphic operating.

However, some breakthroughs for speeding up might be happening (or not... we definitely do not sure).

Please be waiting for more enhanced OD-FPFHE2! :D 

## Postscrpit

We will glad if these code are helpful for your research.

We keep fixing typos in our paper. (Sorry for making misleading... We expect major typos will be fixed within March, 2022.)

Seunghwan Lee. 22.03.04
