/*
 * vericrypt.hpp
 *
 *  Created on: Dec 13, 2022
 *      Author: valeh
 */

#ifndef INCLUDE_VERICRYPT_HPP_
#define INCLUDE_VERICRYPT_HPP_

#include <iostream>
#include <string>
#include <memory>
#include <random>
#include <chrono>
#include <vector>
#include <numeric>

//using namespace std;

#include "params.hpp"
#include <BGV.hpp>
#include <BDLOP.hpp>

typedef BGV::ciphertext_t ctx_t;
typedef BGV::pk_t pk_t;
typedef BGV::sk_t sk_t;

using P = parameters::poly_t;

template<size_t KAPPA>
struct Proof{
	BDLOP::comm_t<KAPPA> t;
	array<BDLOP::comm_rnd_t<KAPPA>, k > z;
	uint8_t chash[SYMBYTES];
	ctx_t ct;
	P h;
};


#define kappa (LC+8)
using proof_t = Proof<kappa>;
typedef BDLOP::comm_pk_t<kappa> comm_key_t;

#define nmod(a, m) ((a) < 0 ? (m - (-1*(a))%m) : (a)%m)

/* ctx names the election and the voter, as bytes the caller chooses; the
 * proof holds only for that ctx, the public key, the commitment key and the
 * ciphertext it carries (see statement_digest in vericrypt.cpp). */
void vericrypt(proof_t &p, const pk_t &pubkey,const comm_key_t& commkey, const P &message,
		const uint8_t *ctx, size_t ctxlen);

bool verify(proof_t &p, const comm_key_t &commkey, const pk_t &pubkey,
		const uint8_t *ctx, size_t ctxlen);

/* A ballot proof as published bytes, always proof_bytes() long: t and ct,
 * uniform residues in NTT form, at 62 bits a coefficient; z, coefficients
 * with |z| < delta_1 - beta_1 (what verify accepts), at 33; h, whose first k
 * coefficients are zero, at 62 for the rest; then chash. proof_encode returns
 * false, writing a buffer nobody should read, for a z outside that range,
 * which no proof the verifier accepts has. proof_decode returns true only
 * for exactly such an encoding, and leaves p in the representation the
 * prover had, which is what the Fiat-Shamir hashes absorb. */
size_t proof_bytes();
bool proof_encode(uint8_t *out, const proof_t &p);
bool proof_decode(proof_t &p, const uint8_t *in, size_t len);

#endif /* INCLUDE_VERICRYPT_HPP_ */
