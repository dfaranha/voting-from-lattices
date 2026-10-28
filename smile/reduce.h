#ifndef REDUCE_H
#define REDUCE_H

#include <stdint.h>
#include "params.h"

#define MONT 90108 // 2^64 mod q
#define MONT2 8119451664 // 2^128 mod q
#define QINV 5842042167526184961 // q^-1 mod 2^64

/*************************************************
* Name:        montgomery_reduce
*
* Description: Montgomery reduction; given a 128-bit integer a, computes
*              64-bit integer congruent to a * R^-1 mod q, where R=2^64
*
* Arguments:   - __int128 a: input integer to be reduced;
*                           has to be in {-q2^63,...,q2^63-1}
*
* Returns:     integer in {-q+1,...,q-1} congruent to a * R^-1 modulo q.
**************************************************/
#define montgomery_reduce PQMX_NAMESPACE(montgomery_reduce)
static inline int64_t montgomery_reduce(__int128 a)
{
  int64_t t;

  t = (int64_t)((uint64_t)a*QINV); // product mod 2^64; unsigned avoids signed-overflow UB
  t = (a - (__int128)t*PQMX_Q) >> 64;
  return t;
}

/*************************************************
* Name:        barrett_reduce
*
* Description: Barrett reduction; given a 64-bit integer a, computes
*              centered representative congruent to a mod q in {-(q-1)/2,...,(q-1)/2}
*
* Arguments:   - int64_t a: input integer to be reduced
*
* Returns:     integer in {-(q-1)/2,...,(q-1)/2} congruent to a modulo q.
**************************************************/
#define barrett_reduce PQMX_NAMESPACE(barrett_reduce)
static inline int64_t barrett_reduce(int64_t a) {
  int64_t t;
  const int64_t v = ( ((__int128)1 <<124 ) + PQMX_Q/2)/PQMX_Q;
  t  = ((__int128)v*a + ((__int128)1 << 123) ) >> 124;
  t *= PQMX_Q;
  return a - t;
}

/*************************************************
* Name:        fqmul
*
* Description: Multiplication followed by Montgomery reduction
*
* Arguments:   - int64_t a: first factor
*              - int64_t b: second factor
*
* Returns integer congruent to a*b*R^{-1} mod q
**************************************************/
static inline int64_t fqmul(int64_t a, int64_t b) {
  return montgomery_reduce((__int128)a*b);
}

#endif
