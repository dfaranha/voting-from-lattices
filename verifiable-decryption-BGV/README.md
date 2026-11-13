# verifiable-decryption-BGV

The code is forked from [tjesi/verifiable-decryption-BGV](https://github.com/tjesi/verifiable-decryption-BGV) repository. Parameters are adjusted in `parameters.h` header file to support BGV encryption with ~62 bit modulus. Make sure [NTL](https://libntl.org/doc/tour-unix.html) is installed in your system.

First, uncomment correct parameter sets in `parameters.h`, and then run `make run` to see results. For your own parameters, use `params.py` to calculate required standard deviations and bound values.

The challenges of both proofs are derived by Fiat-Shamir: SHAKE256 of the statement and the prover's first messages (`fiatshamir.cpp`, with `fips202.c` from `../cumulative_voting/lib`), recomputed by the verifier. The original drew them at random in the prover and the verifier trusted them, which let a prover pick a challenge and solve for a proof of any statement.
