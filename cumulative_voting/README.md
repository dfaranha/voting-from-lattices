# Instructions

The ring arithmetic, in R_q = Z_q[X]/(X^4096 + 1) with q = 4611686018326724609,
is `include/zq.hpp`; the code no longer uses NFLlib. GMP is the only library it
needs.

## Ballot proofs are bound to the election and the voter
`vericrypt` and `verify` take a context, bytes the caller chooses (`election`
uses an election identifier and the voter's index). The first Fiat-Shamir hash
absorbs a digest of that context, the BGV public key, the commitment key and
the ciphertext, so a proof holds only for that voter, election, key and
ciphertext.

## Compile & Run
To test the verifiable encryption of cumulative ballots:
```sh
make vericrypt
./vericrypt
```

The election size is set in `include/params.hpp` (`NV` voters, `NC`
candidates, at most `LC` votes to one candidate, `LS` votes per voter) and can
be overridden at build time, e.g. `make CONFIG="-DNV=10 -DNC=8 -DLC=4 -DLS=8"`.
The proofs are large at the default `LC = 32`, so `ulimit -s unlimited` may be
needed.

## End-to-end election
`election` runs a whole election with homomorphic tallying: every voter
encrypts a random valid cumulative ballot and proves it valid, the board checks
each proof, the accepted ciphertexts are added, and the decryptor decrypts the
sum and proves the decryption correct with the proof of
`../verifiable-decryption-BGV`, which the program then checks. It also checks
that a ballot over the vote budget, a ballot with more than `LC` votes to one
candidate and a ballot changed after its proof are rejected, and that the
decryption proof does not verify for another tally or another ciphertext, that
a ballot proof does not verify for another voter or another election and the
board takes one ballot per voter, and
that a decryption proof whose challenges the prover chose itself is rejected
even though its verification equations hold. It exits non-zero if any check
fails.
```sh
make CONFIG="-DNV=4 -DNC=8 -DLC=4 -DLS=8" election
./election
```
It needs NTL, for the decryption proof. That code and this one cannot share a
translation unit, so `src/verdec_bridge.cpp` passes coefficient arrays between
them (see `include/verdec_bridge.h`); the two use the same q, the same
plaintext modulus p = 10000009 and the same form of BGV.

The decryption proof's challenges are derived by Fiat-Shamir
(`../verifiable-decryption-BGV/fiatshamir.cpp`) and recomputed by the
verifier. Its `Commitment` objects still carry their openings next to `com`;
the verifier here reads only `com`, and rebuilds the statement of the proof of
shortness from the commitment key and those `com` rather than taking the
prover's.
