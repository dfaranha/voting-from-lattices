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
