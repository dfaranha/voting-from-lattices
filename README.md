# voting-from-lattices
Code accompannying the paper "Lattice-Based Zero-Knowledge Proofs in Action: Applications toElectronic Voting" paper, accepted for publication at "Topical Collection on Modern Zero-Knowledge Protocols", Journal of Cryptology (JoC).

## Generic Ballot Correctness Proof
Check `smile` directory for instructions.

## Ballot Correctness Proofs for Homomorphic Tallying
Check `cumulative_voting` directory for instructions.

## Verifiable Decryption
Check `verifiable-decryption-BGV` directory for more information.

## End-to-end Election with Homomorphic Tallying
`cumulative_voting/election` runs the ballot correctness proofs and the
verifiable decryption together, from the voters to a verified tally. Check the
`cumulative_voting` directory for instructions.

WARNING: This is an academic proof of concept, and in particular has not received code review. This implementation is NOT ready for any type of production use.

