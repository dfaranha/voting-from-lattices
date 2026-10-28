#include <stdio.h>
#include <string.h>
#include <math.h>
#include "smile.h"
#include "symmetric.h"
#include "randombytes.h"
#include "util.h"
#include "rej.h"

/***
 * Name: smile_prove
 *
 * Description: Generate Zero-Knowledge proof of correct encryption using SMILE protocol (m=1)
 *
 * Arguments:
 *
 */
/* The statement the challenges are bound to, absorbed first by the prover's
 * and the verifier's first hash: the caller's context (election and voter),
 * the BGV public key, the candidate list and the ciphertext. Without it a
 * proof carries over to another voter or election, and the challenges do not
 * depend on the statement they are about. */
#define STATEMENT_BYTES 32
void smile_candidates_digest(uint8_t out[SMILE_DIGEST_BYTES], const poly candidates[PQMX_L])
{
    static const char label[] = "smile candidate list v1";
    keccak_state st;
    shake128_init(&st);
    shake128_absorb(&st, (const uint8_t *)label, sizeof(label));
    shake128_absorb(&st, (const uint8_t *)candidates, PQMX_L * sizeof(poly));
    shake128_finalize(&st);
    shake128_squeeze(out, SMILE_DIGEST_BYTES, &st);
}

static void statement_digest(uint8_t out[STATEMENT_BYTES], const uint8_t *ctx, size_t ctxlen,
                             const poly bgvpk[2], const uint8_t cand_digest[SMILE_DIGEST_BYTES],
                             const bgvcp *cipher)
{
    static const char label[] = "smile ballot proof v1";
    uint8_t len[8];
    keccak_state st;
    for (int i = 0; i < 8; i++)
        len[i] = (uint8_t)((uint64_t)ctxlen >> (8 * i));
    shake128_init(&st);
    shake128_absorb(&st, (const uint8_t *)label, sizeof(label));
    shake128_absorb(&st, len, sizeof(len));
    shake128_absorb(&st, ctx, ctxlen);
    shake128_absorb(&st, (const uint8_t *)bgvpk, 2 * sizeof(poly));
    shake128_absorb(&st, cand_digest, SMILE_DIGEST_BYTES);
    shake128_absorb(&st, (const uint8_t *)&cipher->u, sizeof(poly));
    shake128_absorb(&st, (const uint8_t *)&cipher->v, sizeof(poly));
    shake128_finalize(&st);
    shake128_squeeze(out, STATEMENT_BYTES, &st);
}

void smile_prove(smile_proof *p, const uint8_t rho[2 * PQMX_SYMBYTES], const poly bgvpk[2], poly candidates[PQMX_L],
                 const uint8_t cand_digest[SMILE_DIGEST_BYTES], bgvrnd *rnd, const bgvcp *cipher, int iota,
                 const uint8_t *ctx, size_t ctxlen)
{
    unsigned int i, j, nonce = 0;
    int ctr = 0;
    uint8_t symbuf[4 * PQMX_SYMBYTES];
    uint8_t *seed = symbuf;
    uint8_t *thash = symbuf + PQMX_SYMBYTES;
    uint8_t *chash = symbuf + 2 * PQMX_SYMBYTES;
    uint8_t *chash2 = symbuf + 3 * PQMX_SYMBYTES;

    poly two, one, poly_P = {0};
    memset(&two, 0, sizeof(poly));
    memset(&one, 0, sizeof(poly));
    two.coeffs[0] = 2;
    one.coeffs[0] = 1;
    poly_P.coeffs[0] = PQMX_P;
    poly_ntt(&two);
    poly_ntt(&one);
    poly_ntt(&poly_P);

    xof_state state;
    randombytes(seed, PQMX_SYMBYTES);

    poly *v = (poly *)aligned_alloc(32, PQMX_M * sizeof(poly));
    memset(v, 0, PQMX_M * sizeof(poly));
    v[0].coeffs[PQMX_N / PQMX_L * iota] = 1;

    commkey ck;
    commrnd r;
    comm t; 

    t.t0 = (poly *) aligned_alloc(32, PQMX_MU*sizeof(poly));
    t.tm = (poly *) aligned_alloc(32, PQMX_M*sizeof(poly));

    expand_commkey(&ck, rho);
    commit(&t, &r, v, &ck);
    // Sample g
    poly g;
    poly_uniform(&g, seed, nonce++);
    for (i = 0; i < PQMX_N / PQMX_L; i++)
    {
        g.coeffs[i] = 0;
    }
    poly_ntt(&g);
    // Commit to g
    poly_add(&t.tm[PQMX_M - 2], &t.tm[PQMX_M - 2], &g);
    poly_reduce(&t.tm[PQMX_M - 2]);

    uint8_t sdigest[STATEMENT_BYTES];
    statement_digest(sdigest, ctx, ctxlen, bgvpk, cand_digest, cipher);
    shake128_init(&state);
    shake128_absorb(&state, sdigest, STATEMENT_BYTES);
    shake128_absorb(&state, rho, PQMX_SYMBYTES);
    shake128_absorb(&state, (uint8_t *)t.t0, PQMX_MU * sizeof(poly));
    shake128_absorb(&state, (uint8_t *)t.tm, sizeof(poly));
    shake128_absorb(&state, (uint8_t *)&t.tm[PQMX_M - 2], sizeof(poly));
    shake128_finalize(&state);
    shake128_squeeze(thash, PQMX_SYMBYTES, &state);

    commrnd y, z2;
    y.s = (poly *)aligned_alloc(32, PQMX_LAMBDA * sizeof(poly));
    y.e = (poly *)aligned_alloc(32, PQMX_MU * sizeof(poly));
    y.em = (poly *)aligned_alloc(32, PQMX_M * sizeof(poly));

    z2.s = (poly *)aligned_alloc(32, PQMX_LAMBDA * sizeof(poly));
    z2.e = (poly *)aligned_alloc(32, PQMX_MU * sizeof(poly));
    z2.em = (poly *)aligned_alloc(32, PQMX_M * sizeof(poly));

    poly *P = (poly *)aligned_alloc(32, PQMX_L * sizeof(poly));
    memset(P, 0, PQMX_L * sizeof(poly));

    poly gamma0[PQMX_K], gamma1;

    poly w[PQMX_MU];
    poly by[PQMX_M];
    poly y_prime[PQMX_ELL];
    poly z_prime[PQMX_ELL];
    poly w_prime[PQMX_K];
    poly tmp;

    poly_ntt(&rnd->r1);
    poly_ntt(&rnd->r2);
    poly_ntt(&rnd->r3);

    poly omega_bin, h;
    int rej1, rej2;

    do
    {
        // Sample y using Discrete Gaussian. TODO: Integrate VCL library here
        for (i = 0; i < PQMX_LAMBDA; i++)
        {
            poly_gaussian(&y.s[i], SIGMA1);
        }
        for (i = 0; i < PQMX_MU; i++)
        {
            poly_gaussian(&y.e[i], SIGMA1);
        }
        for (i = 0; i < PQMX_M; i++)
        {
            poly_gaussian(&y.em[i], SIGMA1);
        }

        polyvec_ntt(y.s, PQMX_LAMBDA);
        polyvec_ntt(y.e, PQMX_MU);
        polyvec_ntt(y.em, PQMX_M);

        for (j = 0; j < PQMX_MU; j++)
        {
            polyvec_basemul_acc_montgomery(&w[j], &ck.b0[j * PQMX_LAMBDA], y.s, PQMX_LAMBDA);
            polyvec_basemul_acc_montgomery(&tmp, &ck.bt[j * PQMX_M], y.em, PQMX_M);
            poly_add(&w[j], &w[j], &tmp);
            poly_tomont(&w[j]);
        }
        polyvec_add(w, w, y.e, PQMX_MU);
        polyvec_reduce(w, PQMX_MU);

        for (j = 0; j < PQMX_M; j++)
        {
            polyvec_basemul_acc_montgomery(&by[j], &ck.bm[j * PQMX_LAMBDA], y.s, PQMX_LAMBDA);
            poly_tomont(&by[j]);
        }
        polyvec_add(by, by, y.em, PQMX_M);
        polyvec_reduce(by, PQMX_M);

        poly c_prime;
        do
        {
            for (i = 0; i < PQMX_ELL; i++)
            {
                poly_gaussian(&y_prime[i], SIGMA0);
            }
            polyvec_ntt(y_prime, PQMX_ELL);

            // w_prime = A * y_prime
            memset(w_prime, 0, PQMX_K * sizeof(poly));
            poly_basemul_montgomery(&tmp, &cipher->u, &y_prime[0]);
            poly_tomont(&tmp);
            poly_add(&w_prime[0], &w_prime[0], &tmp);
            poly_reduce(&w_prime[0]);
            poly_basemul_montgomery(&tmp, &bgvpk[0], &y_prime[1]);
            poly_tomont(&tmp);
            poly_add(&w_prime[0], &w_prime[0], &tmp);
            poly_reduce(&w_prime[0]);
            poly_basemul_montgomery(&tmp, &poly_P, &y_prime[2]);
            poly_tomont(&tmp);
            poly_add(&w_prime[0], &w_prime[0], &tmp);
            poly_reduce(&w_prime[0]);

            poly_basemul_montgomery(&tmp, &cipher->v, &y_prime[0]);
            poly_tomont(&tmp);
            poly_add(&w_prime[1], &w_prime[1], &tmp);
            poly_reduce(&w_prime[1]);
            poly_basemul_montgomery(&tmp, &bgvpk[1], &y_prime[1]);
            poly_tomont(&tmp);
            poly_add(&w_prime[1], &w_prime[1], &tmp);
            poly_reduce(&w_prime[1]);
            poly_basemul_montgomery(&tmp, &poly_P, &y_prime[3]);
            poly_tomont(&tmp);
            poly_add(&w_prime[1], &w_prime[1], &tmp);
            poly_reduce(&w_prime[1]);

            // commit to w_prime
            polyvec_add(&t.tm[1], &t.tm[1], w_prime, PQMX_K);
            polyvec_reduce(&t.tm[1], PQMX_K);

            // get challenge c'
            shake128_init(&state);
            shake128_absorb(&state, thash, PQMX_SYMBYTES);
            shake128_absorb(&state, (uint8_t *)w, PQMX_MU * sizeof(poly));
            shake128_absorb(&state, (uint8_t *)&t.tm[1], (PQMX_K) * sizeof(poly));
            shake128_finalize(&state);
            shake128_squeeze(chash, PQMX_SYMBYTES, &state);

            // rejection sampling 1 on z'
            // sample c
            poly_nonuniform(&c_prime, 0x80, chash, 0);
            poly_ntt(&c_prime);

            // z = y + c*r
            poly_add(&z_prime[0], &c_prime, &y_prime[0]);
            poly_reduce(&z_prime[0]);

            poly_basemul_montgomery(&tmp, &c_prime, &rnd->r1);
            poly_tomont(&tmp);
            poly_sub(&z_prime[1], &y_prime[1], &tmp);
            poly_reduce(&z_prime[1]);

            poly_basemul_montgomery(&tmp, &c_prime, &rnd->r2);
            poly_tomont(&tmp);
            poly_sub(&z_prime[2], &y_prime[2], &tmp);
            poly_reduce(&z_prime[2]);

            poly_basemul_montgomery(&tmp, &c_prime, &rnd->r3);
            poly_tomont(&tmp);
            poly_sub(&z_prime[3], &y_prime[3], &tmp);
            poly_reduce(&z_prime[3]);

            polyvec_invntt_tomont(z_prime, PQMX_ELL);
            polyvec_reduce_mont(z_prime, PQMX_ELL);
            polyvec_reduce(z_prime, PQMX_ELL);

            polyvec_invntt_tomont(y_prime, PQMX_ELL);
            polyvec_reduce_mont(y_prime, PQMX_ELL);
            polyvec_reduce(y_prime, PQMX_ELL);

            polyvec_sub(y_prime, z_prime, y_prime, PQMX_ELL);
            polyvec_reduce(y_prime, PQMX_ELL);

            rej1 = Rej0(z_prime, y_prime);
            if (rej1)
            {
                polyvec_sub(&t.tm[1], &t.tm[1], w_prime, PQMX_K);
                polyvec_reduce(&t.tm[1], PQMX_K);
                printf("Rej0 aborting...\n");
                ctr++;
            }
        } while (rej1);

        polyvec_ntt(z_prime, PQMX_ELL);

        shake128_init(&state);
        shake128_absorb(&state, chash, PQMX_SYMBYTES);
        shake128_absorb(&state, (uint8_t *)z_prime, PQMX_ELL * sizeof(poly));
        shake128_finalize(&state);
        shake128_squeeze(chash2, PQMX_SYMBYTES, &state);

        // compute P
        poly_neg(&c_prime);
        memcpy(P, candidates, PQMX_L * sizeof(poly));
        for (i = 0; i < PQMX_L; i++)
        {
            poly_ntt(&P[i]);
            poly_basemul_montgomery(&P[i], &P[i], &c_prime);
            poly_tomont(&P[i]);
        }


        // compute x1 = NTT(w' - Az') = [P[iota], 0]
        poly x1[PQMX_K];
        memset(x1, 0, PQMX_K * sizeof(poly));
        // memcpy(&x1[1], &P[iota], sizeof(poly));
        poly_basemul_montgomery(&tmp, &cipher->u, &z_prime[0]);
        poly_tomont(&tmp);
        poly_add(&x1[0], &x1[0], &tmp);
        poly_reduce(&x1[0]);
        poly_basemul_montgomery(&tmp, &bgvpk[0], &z_prime[1]);
        poly_tomont(&tmp);
        poly_add(&x1[0], &x1[0], &tmp);
        poly_reduce(&x1[0]);
        poly_basemul_montgomery(&tmp, &poly_P, &z_prime[2]);
        poly_tomont(&tmp);
        poly_add(&x1[0], &x1[0], &tmp);
        poly_reduce(&x1[0]);

        poly_basemul_montgomery(&tmp, &cipher->v, &z_prime[0]);
        poly_tomont(&tmp);
        poly_add(&x1[1], &x1[1], &tmp);
        poly_reduce(&x1[1]);
        poly_basemul_montgomery(&tmp, &bgvpk[1], &z_prime[1]);
        poly_tomont(&tmp);
        poly_add(&x1[1], &x1[1], &tmp);
        poly_reduce(&x1[1]);
        poly_basemul_montgomery(&tmp, &poly_P, &z_prime[3]);
        poly_tomont(&tmp);
        poly_add(&x1[1], &x1[1], &tmp);
        poly_reduce(&x1[1]);

        polyvec_sub(x1, w_prime, x1, PQMX_K);
        polyvec_reduce(x1, PQMX_K);

        // get gamma
        for (i = 0; i < PQMX_K; i++)
        {
            poly_uniform(&gamma0[i], chash2, i);
        }
        poly_uniform(&gamma1, chash2, PQMX_K);

        // compute x2, y1 and h
        poly x2;
        memset(&x2, 0, sizeof(poly));
        for (i = 0; i < PQMX_L; i++)
        {
            poly_basemul_montgomery(&tmp, &P[i], &gamma0[1]);
            poly_tomont(&tmp);
            poly_acc(&x2.coeffs[4 * i], &tmp);
        }

        for (i = 0; i < PQMX_L; i++)
        {
            x2.coeffs[i * 4 + 0] += gamma1.coeffs[0];
            x2.coeffs[i * 4 + 1] += gamma1.coeffs[1];
            x2.coeffs[i * 4 + 2] += gamma1.coeffs[2];
            x2.coeffs[i * 4 + 3] += gamma1.coeffs[3];
        }
        poly_reduce(&x2);

        poly y1;
        memset(&y1, 0, sizeof(poly));
        poly_basemul_montgomery(&y1, &v[0], &x2);
        poly_tomont(&y1);
        polyvec_basemul_acc_montgomery(&tmp, x1, gamma0, PQMX_K);
        poly_tomont(&tmp);
        poly_sub(&y1, &y1, &tmp);
        poly_reduce(&y1);

        y1.coeffs[0] -= gamma1.coeffs[0];
        y1.coeffs[1] -= gamma1.coeffs[1];
        y1.coeffs[2] -= gamma1.coeffs[2];
        y1.coeffs[3] -= gamma1.coeffs[3];
        poly_reduce(&y1);

        poly_add(&h, &g, &y1);
        poly_reduce(&h);

        // get alphas
        shake128_init(&state);
        shake128_absorb(&state, chash2, PQMX_SYMBYTES);
        shake128_absorb(&state, (uint8_t *)&h, sizeof(poly));
        shake128_finalize(&state);
        shake128_squeeze(chash2, PQMX_SYMBYTES, &state);

        poly alpha[2];
        poly_uniform(&alpha[0], chash2, 0);
        poly_uniform(&alpha[1], chash2, 1);

        // compute psi_sm, omega_sm, omega_bin,
        poly psi[2];
        memset(psi, 0, 2 * sizeof(poly));
        poly_basemul_montgomery(&psi[0], &x2, &by[0]);
        poly_tomont(&psi[0]);
        polyvec_basemul_acc_montgomery(&tmp, &by[1], gamma0, PQMX_K);
        poly_tomont(&tmp);
        poly_sub(&psi[0], &tmp, &psi[0]);
        poly_reduce(&psi[0]);
        poly_sub(&psi[0], &psi[0], &by[PQMX_K + 1]);
        poly_reduce(&psi[0]);

        poly_basemul_montgomery(&tmp, &two, &v[0]);
        poly_tomont(&tmp);
        poly_sub(&psi[1], &one, &tmp);
        poly_reduce(&psi[1]);
        poly_basemul_montgomery(&psi[1], &psi[1], &by[0]);
        poly_tomont(&psi[1]);

        poly remember;
        polyvec_basemul_acc_montgomery(&remember, psi, alpha, 2);
        poly_tomont(&remember);

        poly_add(&t.tm[PQMX_M - 1], &t.tm[PQMX_M - 1], &remember);
        poly_reduce(&t.tm[PQMX_M - 1]);


        poly_basemul_montgomery(&tmp, &by[0], &by[0]);
        poly_tomont(&tmp);
        poly_basemul_montgomery(&omega_bin, &tmp, &alpha[1]);
        poly_tomont(&omega_bin);
        poly_add(&omega_bin, &omega_bin, &by[PQMX_M - 1]);
        poly_reduce(&omega_bin);

        // get challenge c
        shake128_init(&state);
        shake128_absorb(&state, chash2, PQMX_SYMBYTES);
        shake128_absorb(&state, (uint8_t *)&t.tm[PQMX_M - 1], sizeof(poly));
        shake128_absorb(&state, (uint8_t *)&omega_bin, sizeof(poly));
        shake128_finalize(&state);
        shake128_squeeze(chash2, PQMX_SYMBYTES, &state);

        poly c;
        poly_nonuniform(&c, 0x80, chash2, 0);
        poly_ntt(&c);

        // z = y + c*r
        for (i = 0; i < PQMX_MU; i++)
        {
            poly_basemul_montgomery(&tmp, &c, &r.e[i]);
            poly_tomont(&tmp);
            poly_add(&z2.e[i], &tmp, &y.e[i]);
            poly_reduce(&z2.e[i]);
        }
        for (i = 0; i < PQMX_LAMBDA; i++)
        {
            poly_basemul_montgomery(&tmp, &c, &r.s[i]);
            poly_tomont(&tmp);
            poly_add(&z2.s[i], &tmp, &y.s[i]);
            poly_reduce(&z2.s[i]);
        }
        for (i = 0; i < PQMX_M; i++)
        {
            poly_basemul_montgomery(&tmp, &c, &r.em[i]);
            poly_tomont(&tmp);
            poly_add(&z2.em[i], &tmp, &y.em[i]);
            poly_reduce(&z2.em[i]);
        }

        polyvec_invntt_tomont(y.s, PQMX_LAMBDA);
        polyvec_reduce_mont(y.s, PQMX_LAMBDA);
        polyvec_reduce(y.s, PQMX_LAMBDA);

        polyvec_invntt_tomont(y.e, PQMX_MU);
        polyvec_reduce_mont(y.e, PQMX_MU);
        polyvec_reduce(y.e, PQMX_MU);

        polyvec_invntt_tomont(y.em, PQMX_M);
        polyvec_reduce_mont(y.em, PQMX_M);
        polyvec_reduce(y.em, PQMX_M);

        polyvec_invntt_tomont(z2.s, PQMX_LAMBDA);
        polyvec_reduce_mont(z2.s, PQMX_LAMBDA);
        polyvec_reduce(z2.s, PQMX_LAMBDA);

        polyvec_invntt_tomont(z2.e, PQMX_MU);
        polyvec_reduce_mont(z2.e, PQMX_MU);
        polyvec_reduce(z2.e, PQMX_MU);

        polyvec_invntt_tomont(z2.em, PQMX_M);
        polyvec_reduce_mont(z2.em, PQMX_M);
        polyvec_reduce(z2.em, PQMX_M);

        rej2 = Rej1(&z2, &y);
        if (rej2)
        {
            polyvec_sub(&t.tm[1], &t.tm[1], w_prime, PQMX_K);
            polyvec_reduce(&t.tm[1], PQMX_K);
            poly_sub(&t.tm[PQMX_M - 1], &t.tm[PQMX_M - 1], &remember);
            poly_reduce(&t.tm[PQMX_M - 1]);
            printf("Rej1 aborting...\n");
            ctr++;
        }
        // rejection sampling on z
    } while (rej2);
    printf("Rejected %d times\n", ctr);

    poly_invntt_tomont(&h);
    poly_reduce_mont(&h);
    poly_reduce(&h);

    memcpy(&p->h, &h, sizeof(poly));
    memcpy(p->hash, chash, PQMX_SYMBYTES);
    memcpy(p->hash2, chash2, PQMX_SYMBYTES);
    p->t.t0 = t.t0;
    p->t.tm = t.tm;
    p->z2.e = z2.e;
    p->z2.em = z2.em;
    p->z2.s = z2.s;
    memcpy(p->z1, z_prime, PQMX_ELL * sizeof(poly));

    free(v);
    free(P);
    free(r.e);
    free(r.em);
    free(r.s);
    free(y.e);
    free(y.em);
    free(y.s);
    free(ck.b0);
    free(ck.bm);
    free(ck.bt);
}

/***
 * Name: smile_proof_verify
 *
 * Description: Verify  SMILE proof (m=1)
 *
 * Arguments:
 *
 */
int smile_proof_verify(smile_proof *p, const uint8_t rho[2 * PQMX_SYMBYTES], const poly bgvpk[2], poly candidates[PQMX_L],
                       const uint8_t cand_digest[SMILE_DIGEST_BYTES], const bgvcp *cipher,
                       const uint8_t *ctx, size_t ctxlen)
{
    unsigned int i, j, ret = 0;
    uint8_t symbuf[4 * PQMX_SYMBYTES];
    uint8_t *seed = symbuf;
    uint8_t *thash = symbuf + PQMX_SYMBYTES;
    uint8_t *chash = symbuf + 2 * PQMX_SYMBYTES;
    uint8_t *chash2 = symbuf + 3 * PQMX_SYMBYTES;

    poly two, one, poly_P = {0};
    memset(&two, 0, sizeof(poly));
    memset(&one, 0, sizeof(poly));
    two.coeffs[0] = 2;
    one.coeffs[0] = 1;
    poly_P.coeffs[0] = PQMX_P;
    poly_ntt(&two);
    poly_ntt(&one);
    poly_ntt(&poly_P);

    xof_state state;
    randombytes(seed, PQMX_SYMBYTES);

    polyvec_invntt_tomont(p->z1, PQMX_ELL);
    polyvec_reduce_mont(p->z1, PQMX_ELL);
    polyvec_reduce(p->z1, PQMX_ELL);
    double z1_norm = polyvec_norm(p->z1, 2, PQMX_ELL);

    if (z1_norm > (SIGMA0 * sqrt(2.0 * PQMX_ELL * PQMX_N)))
    {
        ret = 1;
    }
    polyvec_ntt(p->z1, PQMX_ELL);

    double z2_norm = 0.0;
    z2_norm += pow(polyvec_norm(p->z2.s, 2, PQMX_LAMBDA), 2);
    z2_norm += pow(polyvec_norm(p->z2.e, 2, PQMX_MU), 2);
    z2_norm += pow(polyvec_norm(p->z2.em, 2, PQMX_M), 2);
    z2_norm = sqrt(z2_norm);

    if (z2_norm > (SIGMA1 * sqrt(2.0 * PQMX_KAPPA * PQMX_N)))
    {
        ret = 1;
    }

    polyvec_ntt(p->z2.s, PQMX_LAMBDA);
    polyvec_ntt(p->z2.e, PQMX_MU);
    polyvec_ntt(p->z2.em, PQMX_M);

    if (ret)
    {
        return 1;
    }

    commkey ck;
    expand_commkey(&ck, rho);

    uint8_t sdigest[STATEMENT_BYTES];
    statement_digest(sdigest, ctx, ctxlen, bgvpk, cand_digest, cipher);
    shake128_init(&state);
    shake128_absorb(&state, sdigest, STATEMENT_BYTES);
    shake128_absorb(&state, rho, PQMX_SYMBYTES);
    shake128_absorb(&state, (uint8_t *)p->t.t0, PQMX_MU * sizeof(poly));
    shake128_absorb(&state, (uint8_t *)p->t.tm, sizeof(poly));
    shake128_absorb(&state, (uint8_t *)&p->t.tm[PQMX_M - 2], sizeof(poly));
    shake128_finalize(&state);
    shake128_squeeze(thash, PQMX_SYMBYTES, &state);

    poly c, c_prime;
    poly_nonuniform(&c_prime, 0x80, p->hash, 0);
    poly_ntt(&c_prime);
    poly_nonuniform(&c, 0x80, p->hash2, 0);
    poly_ntt(&c);

    poly *P = (poly *)aligned_alloc(32, PQMX_L * sizeof(poly));
    memset(P, 0, PQMX_L * sizeof(poly));

    poly gamma0[PQMX_K], gamma1;

    poly w[PQMX_MU];
    poly tmp;

    for (j = 0; j < PQMX_MU; j++)
    {
        polyvec_basemul_acc_montgomery(&w[j], &ck.b0[j * PQMX_LAMBDA], p->z2.s, PQMX_LAMBDA);
        polyvec_basemul_acc_montgomery(&tmp, &ck.bt[j * PQMX_M], p->z2.em, PQMX_M);
        poly_add(&w[j], &w[j], &tmp);
        poly_tomont(&w[j]);
    }
    polyvec_add(w, w, p->z2.e, PQMX_MU);
    polyvec_reduce(w, PQMX_MU);

    for (i = 0; i < PQMX_MU; i++)
    {
        poly_basemul_montgomery(&tmp, &c, &p->t.t0[i]);
        poly_tomont(&tmp);
        poly_sub(&w[i], &w[i], &tmp);
        poly_reduce(&w[i]);
    }

    // w_prime = A * z_prime
    poly w_prime[PQMX_K];
    poly_basemul_montgomery(&w_prime[0], &cipher->u, &p->z1[0]);
    poly_tomont(&w_prime[0]);
    poly_basemul_montgomery(&tmp, &bgvpk[0], &p->z1[1]);
    poly_tomont(&tmp);
    poly_add(&w_prime[0], &w_prime[0], &tmp);
    poly_reduce(&w_prime[0]);
    poly_basemul_montgomery(&tmp, &poly_P, &p->z1[2]);
    poly_tomont(&tmp);
    poly_add(&w_prime[0], &w_prime[0], &tmp);
    poly_reduce(&w_prime[0]);

    poly_basemul_montgomery(&w_prime[1], &cipher->v, &p->z1[0]);
    poly_tomont(&w_prime[1]);
    poly_basemul_montgomery(&tmp, &bgvpk[1], &p->z1[1]);
    poly_tomont(&tmp);
    poly_add(&w_prime[1], &w_prime[1], &tmp);
    poly_reduce(&w_prime[1]);
    poly_basemul_montgomery(&tmp, &poly_P, &p->z1[3]);
    poly_tomont(&tmp);
    poly_add(&w_prime[1], &w_prime[1], &tmp);
    poly_reduce(&w_prime[1]);

    poly *f = (poly *)aligned_alloc(32, PQMX_M * sizeof(poly));
    poly *tm = (poly *)aligned_alloc(32, PQMX_M * sizeof(poly));
    memcpy(tm, p->t.tm, PQMX_M * sizeof(poly));

    for (i = 0; i < PQMX_K; i++)
    {
        poly_sub(&tm[1 + i], &p->t.tm[1 + i], &w_prime[i]);
        poly_reduce(&tm[1 + i]);
    }

    for (j = 0; j < PQMX_M; j++)
    {
        polyvec_basemul_acc_montgomery(&f[j], &ck.bm[j * PQMX_LAMBDA], p->z2.s, PQMX_LAMBDA);
        poly_tomont(&f[j]);
        poly_basemul_montgomery(&tmp, &c, &tm[j]);
        poly_tomont(&tmp);
        poly_sub(&f[j], &f[j], &tmp);
        poly_reduce(&f[j]);
    }
    polyvec_add(f, f, p->z2.em, PQMX_M);
    polyvec_reduce(f, PQMX_M);

    shake128_init(&state);
    shake128_absorb(&state, thash, PQMX_SYMBYTES);
    shake128_absorb(&state, (uint8_t *)w, PQMX_MU * sizeof(poly));
    shake128_absorb(&state, (uint8_t *)&p->t.tm[1], (PQMX_K) * sizeof(poly));
    shake128_finalize(&state);
    shake128_squeeze(chash, PQMX_SYMBYTES, &state);

    shake128_init(&state);
    shake128_absorb(&state, chash, PQMX_SYMBYTES);
    shake128_absorb(&state, (uint8_t *)p->z1, PQMX_ELL * sizeof(poly));
    shake128_finalize(&state);
    shake128_squeeze(chash2, PQMX_SYMBYTES, &state);

    // compute P
    poly_neg(&c_prime);
    memcpy(P, candidates, PQMX_L * sizeof(poly));
    for (i = 0; i < PQMX_L; i++)
    {
        poly_ntt(&P[i]);
        poly_basemul_montgomery(&P[i], &P[i], &c_prime);
        poly_tomont(&P[i]);
    }

    // get gamma
    for (i = 0; i < PQMX_K; i++)
    {
        poly_uniform(&gamma0[i], chash2, i);
    }
    poly_uniform(&gamma1, chash2, PQMX_K);

    // compute x2, y1 and h
    poly x2;
    memset(&x2, 0, sizeof(poly));
    for (i = 0; i < PQMX_L; i++)
    {
        poly_basemul_montgomery(&tmp, &P[i], &gamma0[1]);
        poly_tomont(&tmp);
        poly_acc(&x2.coeffs[i * 4], &tmp);
    }

    for (i = 0; i < PQMX_L; i++)
    {
        x2.coeffs[i * 4 + 0] += gamma1.coeffs[0];
        x2.coeffs[i * 4 + 1] += gamma1.coeffs[1];
        x2.coeffs[i * 4 + 2] += gamma1.coeffs[2];
        x2.coeffs[i * 4 + 3] += gamma1.coeffs[3];
    }
    poly_reduce(&x2);

    uint8_t check1 = 0;
    for (i = 0; i < PQMX_N / PQMX_L; i++)
    {
        if (p->h.coeffs[i] != 0)
        {
            check1 = 1;
        }
    }
    poly_ntt(&p->h);

    // get alphas
    shake128_init(&state);
    shake128_absorb(&state, chash2, PQMX_SYMBYTES);
    shake128_absorb(&state, (uint8_t *)&p->h, sizeof(poly));
    shake128_finalize(&state);
    shake128_squeeze(chash2, PQMX_SYMBYTES, &state);

    poly F1;
    memset(&F1, 0, sizeof(poly));
    poly_basemul_montgomery(&F1, &f[0], &x2);
    poly_tomont(&F1);
    polyvec_basemul_acc_montgomery(&tmp, &f[1], gamma0, PQMX_K);
    poly_tomont(&tmp);
    poly_sub(&F1, &tmp, &F1);
    poly_reduce(&F1);

    poly gamma010 = {0};
    gamma010.coeffs[0] = gamma1.coeffs[0];
    gamma010.coeffs[1] = gamma1.coeffs[1];
    gamma010.coeffs[2] = gamma1.coeffs[2];
    gamma010.coeffs[3] = gamma1.coeffs[3];
    poly_basemul_montgomery(&tmp, &gamma010, &c);
    poly_tomont(&tmp);
    poly_sub(&F1, &F1, &tmp);
    poly_reduce(&F1);
    poly_basemul_montgomery(&F1, &F1, &c);
    poly_tomont(&F1);

    poly alpha[2];
    poly_uniform(&alpha[0], chash2, 0);
    poly_uniform(&alpha[1], chash2, 1);

    poly F[2];
    memset(F, 0, 2 * sizeof(poly));
    poly_basemul_montgomery(&tmp, &c, &f[PQMX_K + 1]);
    poly_tomont(&tmp);
    poly_sub(&F[0], &F1, &tmp);
    poly_reduce(&F[0]);
    poly_basemul_montgomery(&tmp, &c, &p->h);
    poly_tomont(&tmp);
    poly_basemul_montgomery(&tmp, &tmp, &c);
    poly_tomont(&tmp);
    poly_sub(&F[0], &F[0], &tmp);
    poly_reduce(&F[0]);

    poly_basemul_montgomery(&F[1], &f[0], &f[0]);
    poly_tomont(&F[1]);
    poly_basemul_montgomery(&tmp, &f[0], &c);
    poly_tomont(&tmp);
    poly_add(&F[1], &F[1], &tmp);
    poly_reduce(&F[1]);

    poly omega;
    polyvec_basemul_acc_montgomery(&omega, F, alpha, 2);
    poly_tomont(&omega);
    poly_add(&omega, &omega, &f[PQMX_M - 1]);
    poly_reduce(&omega);

    // get challenge c
    shake128_init(&state);
    shake128_absorb(&state, chash2, PQMX_SYMBYTES);
    shake128_absorb(&state, (uint8_t *)&p->t.tm[PQMX_M - 1], sizeof(poly));
    shake128_absorb(&state, (uint8_t *)&omega, sizeof(poly));
    shake128_finalize(&state);
    shake128_squeeze(chash2, PQMX_SYMBYTES, &state);


    uint8_t r1 = 0, r2 = 0;
    for (j = 0; j < PQMX_SYMBYTES; j++)
    {
        r1 |= chash[j] ^ p->hash[j];
    }

    for (j = 0; j < PQMX_SYMBYTES; j++)
    {
        r2 |= chash2[j] ^ p->hash2[j];
    }

    free(tm);
    free(P);
    free(f);
    free(ck.b0);
    free(ck.bm);
    free(ck.bt);
    return check1 | (r1 | r2);
}
