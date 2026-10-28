/*
 * verdec_bridge.h
 *
 * The proof of correct decryption of ../verifiable-decryption-BGV, for the
 * tally. That code is built on NTL and its parameters.h defines N, q, p and B
 * at global scope, while params.hpp here #defines k, NV, NC and others, so the
 * two cannot share a translation unit. This interface passes plain arrays:
 * DEGREE coefficients in [0, q), in the coefficient domain (not NTT form).
 *
 * The two codebases share the ring Z_q[X]/(X^4096 + 1), q = 4611686018326724609,
 * the plaintext modulus p = 10000009 and the BGV form v - s u = m + p d, under
 * the names (u, v) = (c0, c1) and pk = (a, b) with b = a s + p e.
 */
#ifndef VERDEC_BRIDGE_H
#define VERDEC_BRIDGE_H

#include <cstddef>
#include <cstdint>

namespace verdec {

/* The ring degree the decryption proof is built for. */
size_t degree();

/* The decryptor's keys: commits to s and keeps it for proving. Seeds NTL's
 * generator, and rand(), from the operating system. */
void setup(const uint64_t a[], const uint64_t b[], const uint64_t s[]);

struct proof;

/* Proves that (u, v) decrypts to m. */
proof *prove(const uint64_t u[], const uint64_t v[], const uint64_t m[]);

/* Checks the proof of shortness of the noise and the linear relation between
 * the committed key, the committed noise, (u, v) and m. */
bool verify(const proof *pi, const uint64_t u[], const uint64_t v[],
		const uint64_t m[]);

void release(proof *pi);

/* What a prover who picks its own challenges could send: the proof of
 * shortness with C = 0 and W = A Z, or the linear proof with c = 0 and t1, t2
 * and u solved for. Either satisfies its verification equations, which
 * *equations_hold reports; only the Fiat-Shamir check can reject it. */
enum part { SHORTNESS, LINEAR };
bool verify_forged(part which, const proof *pi, const uint64_t u[],
		const uint64_t v[], const uint64_t m[], bool *equations_hold);

}

#endif
