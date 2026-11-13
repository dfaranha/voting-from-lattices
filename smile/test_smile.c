#include <string.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <time.h>
#include "randombytes.h"
#include "poly.h"
#include "comm.h"
#include "util.h"
#include "bgv.h"
#include "smile.h"

/* The context a ballot proof is bound to: an election identifier and a voter
 * index. */
static size_t context(uint8_t out[64], const uint8_t eid[16], uint32_t voter)
{
  static const char label[] = "smile test ballot";
  size_t n = sizeof(label);
  memcpy(out, label, n);
  memcpy(out + n, eid, 16);
  n += 16;
  for (int i = 0; i < 4; i++)
    out[n++] = (uint8_t)(voter >> (8 * i));
  return n;
}

/* smile_proof_verify works on the proof in place, so each check gets a deep
 * copy. */
static void *dup(const void *src, size_t len)
{
  void *d = aligned_alloc(32, len);
  memcpy(d, src, len);
  return d;
}

static smile_proof clone_proof(const smile_proof *p)
{
  smile_proof c = *p;
  c.t.t0 = (poly *) dup(p->t.t0, PQMX_MU * sizeof(poly));
  c.t.tm = (poly *) dup(p->t.tm, PQMX_M * sizeof(poly));
  c.z2.s = (poly *) dup(p->z2.s, PQMX_LAMBDA * sizeof(poly));
  c.z2.e = (poly *) dup(p->z2.e, PQMX_MU * sizeof(poly));
  c.z2.em = (poly *) dup(p->z2.em, PQMX_M * sizeof(poly));
  return c;
}

static void free_proof(smile_proof *p)
{
  free(p->t.t0);
  free(p->t.tm);
  free(p->z2.e);
  free(p->z2.em);
  free(p->z2.s);
}

/* 0 if the proof verifies under (eid, voter). */
static int verify_as(const smile_proof *p, const uint8_t seed[], const poly bgvpk[2],
                     poly *candidates, const uint8_t cand_digest[], const bgvcp *cp,
                     const uint8_t eid[16], uint32_t voter)
{
  uint8_t ctx[64];
  size_t len = context(ctx, eid, voter);
  smile_proof c = clone_proof(p);
  int res = smile_proof_verify(&c, seed, bgvpk, candidates, cand_digest, cp, ctx, len);
  free_proof(&c);
  return res;
}


int main(void)
{

  uint8_t seed[2*PQMX_SYMBYTES];

  randombytes(seed, 2*PQMX_SYMBYTES);

  poly bgvpk[2], bgvsk;
  bgv_genkey(bgvpk, &bgvsk);


  bgvcp cp;   //ballot box
  bgvrnd bgv_rnd; //bgv randomnesses

  // Randomly generate a ballot and encrpyt
  poly *candidates = (poly *) aligned_alloc(32, PQMX_L*sizeof(poly));
  poly m;
  uint8_t paddedballot[PQMX_INDCPA_MSGBYTES];
  char* choice = (char *) malloc(32);
  for(int i=0; i< PQMX_L; i++){
    memset(choice, 0, 32);
    sprintf(choice,"cand_%d", i);
    padding(paddedballot, (uint8_t*) choice, strlen(choice));
    poly_frommsg(&candidates[i], paddedballot);
  }
  
  int i = secure_random(PQMX_L);
  bgv_enc(&cp, &candidates[i], bgvpk, &bgv_rnd);
  printf("Encrypted ballot: cand_%d\n", i);
  
//   polyvec_invntt_tomont(bgvpk,2);
//   polyvec_reduce_mont(bgvpk,2);

  // generate zero-knowledge proof, bound to this election and voter 0
  uint8_t eid[16], ctx[64];
  randombytes(eid, sizeof(eid));
  size_t ctxlen = context(ctx, eid, 0);
  uint8_t cand_digest[SMILE_DIGEST_BYTES];
  clock_t tic = clock();
  smile_candidates_digest(cand_digest, candidates);
  clock_t toc = clock();
  printf("candidate list hashed in %f seconds, once per election\n", (double)(toc - tic) / CLOCKS_PER_SEC);
  smile_proof p;
  tic = clock();
  smile_prove(&p, seed, bgvpk, candidates, cand_digest, &bgv_rnd, &cp, i, ctx, ctxlen);
  toc = clock();
  printf("proof generated in %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);

  // verify the proof
  tic = clock();
  int res = verify_as(&p, seed, bgvpk, candidates, cand_digest, &cp, eid, 0);
  toc = clock();
  printf("proof %s in %f seconds\n", (res ? "rejected" : "verified"), (double)(toc - tic) / CLOCKS_PER_SEC);
  int failures = (res != 0);

  // the same proof under another voter, and in another election
  int other_voter = verify_as(&p, seed, bgvpk, candidates, cand_digest, &cp, eid, 1);
  uint8_t eid2[16];
  memcpy(eid2, eid, 16);
  eid2[0] ^= 1;
  int other_election = verify_as(&p, seed, bgvpk, candidates, cand_digest, &cp, eid2, 0);
  printf("proof for another voter %s\n", other_voter ? "rejected" : "ACCEPTED");
  printf("proof in another election %s\n", other_election ? "rejected" : "ACCEPTED");
  failures += (other_voter == 0) + (other_election == 0);

  // Start tallying
  uint8_t ballot[32], ballotOld[32];
  bgv_dec(&m, &bgvsk, &cp);
  poly_tomsg(paddedballot, &m);
  remove_padding(ballot, paddedballot);
  poly_tomsg(paddedballot, &candidates[i]);
  remove_padding(ballotOld, paddedballot);
  printf("Decrypted ballot: ");
  puts((const char*)ballot);


  printf("\n");
  if(strcmp((const char*)ballot, (const char*)ballotOld)){
    printf("Decrypted ciphertext does not match actual plaintext!\n");
    failures++;
  }

  printf("proof size: %ld\n", sizeof(smile_proof));


  free_proof(&p);

  return failures != 0;
}

