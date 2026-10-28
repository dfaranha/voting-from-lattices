/*
 * params.hpp
 *
 *  Created on: Sep 5, 2022
 *      Author: valeh
 */

#ifndef PARAMS_HPP_
#define PARAMS_HPP_

#include <gmpxx.h>
#include "zq.hpp"

namespace parameters{
	using poly_t = zq::poly<4096>;
	const size_t MU=1, LAMBDA=1;
}

using poly_t = zq::poly<4096>;

#define plaintextModulus 10000009  // set p next_prime after NV*LC

/* The elections benchmarked against each other share the knobs MSGS (voters),
 * NCAND (candidates), MIXERS (mix servers) and PARTIES (trustees), as
 * make CONFIG="-D...". MSGS and NCAND are NV and NC here, which still work.
 * The tally is homomorphic, with no mix servers and one decryptor, so MIXERS
 * can only be 0 and PARTIES only 1. */
#if defined(MSGS) && !defined(NV)
#define NV MSGS
#endif
#if defined(NCAND) && !defined(NC)
#define NC NCAND
#endif
#if defined(MIXERS) && MIXERS != 0
#error "cumulative_voting: the tally is homomorphic, with no mix servers, so MIXERS can only be 0"
#endif
#if defined(PARTIES) && PARTIES != 1
#error "cumulative_voting: one decryptor decrypts the tally, so PARTIES can only be 1"
#endif
#ifndef NV
#define NV 100                    // number of voters
#endif
#ifndef NC
#define NC 100                     // number of candidates
#endif
#ifndef LC
#define LC 32                      // maximum vote that can be given to a single candidate
#endif
#ifndef LS
#define LS 100                     // total amount of votes a single voter can give - vote budget
#endif
#define delta_1 (1L<<32)           // rejection sampling parameter
#define beta_1 1024                // rejection sampling and MSIS parameter
#define k 2                        // soundness boost and automorphism order parameter

#define SYMBYTES 32                // output of hash functions in bytes

#endif /* PARAMS_HPP_ */
