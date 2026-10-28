# Instructions

The ring arithmetic, in R_q = Z_q[X]/(X^4096 + 1) with q = 4611686018326724609,
is `include/zq.hpp`; the code no longer uses NFLlib. GMP is the only library it
needs.

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
