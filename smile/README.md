# Improved Lattice-based Mix-Net implementation

Here, we implement SMILE protocol with `m=1` by adapting the code from [*ilmx*](https://github.com/Valeh2012/ilmx) repository. To sample polynomial coefficients with large standard deviation [*COSAC sampler*](https://github.com/raykzhao/gaussian_ac) is used.


## Run
There are, four test files from *ilmx* repository `test_rlwe, test_shortness, test_shuffle, test_mx` and one constant time gaussian sampler test `test_gaussian` from *gaussian_ac* repository. We add a new test for zero-knowledge proof of set membership (SMILE) `test_smile`. To compile and run SMILE protocol on Linux, run

```sh
make test_smile
./text_smile
```


## References

1. [*Improved Lattice-Based Mix-Nets for Electronic Voting*](https://link.springer.com/chapter/10.1007/978-3-031-08896-4_6)
2. [*SMILE: Set Membership from Ideal Lattices with Applications to Ring Signatures and Confidential Transactions*](https://link.springer.com/chapter/10.1007/978-3-030-84245-1_21)
3. [*COSAC: COmpact and Scalable Arbitrary-Centered Discrete Gaussian Sampling over Integers*](https://link.springer.com/chapter/10.1007/978-3-030-44223-1_16)


## Ballot proofs are bound to the election and the voter

`smile_prove` and `smile_proof_verify` take a context, bytes the caller chooses to name the election and the voter, and the digest of the candidate list from `smile_candidates_digest`. The first Fiat-Shamir hash absorbs a digest of the context, the BGV public key, the candidate list and the ciphertext, so a proof holds only for that voter and election. The candidate list is 32 MiB and fixed for an election, so it is hashed once, by whoever proves or verifies against it; the digest passed must be that of the list passed. `test_smile` checks that a proof is rejected for another voter and in another election, and exits non-zero otherwise.
