#ifndef SMILE_H
#define SMILE_H

#include <stdint.h>
#include "params.h"
#include "bgv.h"
#include "comm.h"

typedef struct {
  uint8_t hash[PQMX_SYMBYTES];
  uint8_t hash2[PQMX_SYMBYTES];
  poly z1[PQMX_ELL];
  commrnd z2;
  comm t;
  poly h;
} smile_proof;

#define SMILE_DIGEST_BYTES 32

/* The digest of the candidate list the proofs are bound to. The list is fixed
 * for an election and 32 MiB, so it is hashed once, by whoever proves or
 * verifies against it, and passed as cand_digest with the list itself. */
void smile_candidates_digest(uint8_t out[SMILE_DIGEST_BYTES], const poly candidates[PQMX_L]);

/* ctx names the election and the voter, as bytes the caller chooses. The
 * proof's challenges are bound to it and to the whole statement: the BGV
 * public key, the candidate list (through cand_digest, which must be the
 * digest of candidates) and the ciphertext, besides rho. */
void smile_prove(smile_proof *p, const uint8_t rho[2*PQMX_SYMBYTES], const poly bgvpk[2], poly candidates[PQMX_L],
                 const uint8_t cand_digest[SMILE_DIGEST_BYTES], bgvrnd *rnd, const bgvcp *cipher, int iota,
                 const uint8_t *ctx, size_t ctxlen);
int smile_proof_verify(smile_proof *p, const uint8_t rho[2*PQMX_SYMBYTES], const poly bgvpk[2], poly candidates[PQMX_L],
                       const uint8_t cand_digest[SMILE_DIGEST_BYTES], const bgvcp *cipher,
                       const uint8_t *ctx, size_t ctxlen);


#endif