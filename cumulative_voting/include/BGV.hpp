/**
 *
 * This file is part of BGV-NFLlib
 *
 * Copyright (C) 2015, 2016  CryptoExperts
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
 *
 */

#pragma once

#include <cassert>
#include <cstddef>

#include <algorithm>
#include <chrono>
#include <iostream>
#include <memory>


namespace BGV {
class sk_t;
class pk_t;
class ciphertext_t;

}  // namespace BGV

/**
 * Headers of utilitary functions implemented at the end of the document
 */
namespace BGV {
using V = typename parameters::poly_t::value_type;
using S = typename parameters::poly_t::signed_value_type;
namespace util {
inline void center(S &rop, V const &op1, V const &op2, V const &op2Div2);
inline void div_and_round(mpz_t &rop, mpz_t const &op1, mpz_t const &op2,
		mpz_t const &op2Div2);

inline void reduce(parameters::poly_t &coefficients, V multiplier,
		V const &mod_init, V const &mod_initDiv2);
//inline void lift(std::array<mpz_t, parameters::polyZ_p::degree> &coefficients,parameters::polyZ_p const &c);

}  // namespace util
}  // namespace BGV

/**
 * Class to store the secret key
 */
namespace BGV {
class sk_t {
	using P = parameters::poly_t;

public:
	/// The secret key is a polynomial
	P value; // change to binary distribution

	/// Constructor
	sk_t() {
		value = nfl::non_uniform(2);
		value.ntt_pow_phi();  // store in NTT form
	}
};
}  // namespace BGV


/**
 * Class to store the public key
 */
namespace BGV {
class pk_t {
	using P = parameters::poly_t;

public:
	/// Public key elements
	P b, a;
	P a_shoup, b_shoup;


	/// Link to evaluation key
	P::value_type qDivBy2;

	/// Interesting value
	long noise_max;

	/// Constructor
	pk_t(sk_t const &sk) {
		// Link to the evaluation key

		qDivBy2 = P::get_modulus(0) >> 1;

		// random a (already in NTT form)
		a = nfl::uniform();
		a_shoup = nfl::compute_shoup(a);

		P t = {plaintextModulus};
		t.ntt_pow_phi();


		// b = a*sk + t*e
		b = nfl::non_uniform(2);
		b.ntt_pow_phi();  // transform via NTT
		P e = a*sk.value;
		b = b*t;
		b = b + e;
		b_shoup = nfl::compute_shoup(b);

		// Set the plaintext modulus
		noise_max = mpz_sizeinbase(P::moduli_product(), 2) - 1;
	}
};
}  // namespace BGV

/**
 * Class to store a ciphertext
 */
namespace BGV {
class ciphertext_t {
	using P = parameters::poly_t;
//	using PZ = parameters::polyZ_p;

public:
	typedef P type;

	/// ciphertext_t contains two polynomials (c0, c1)
	P c0, c1;

	/// Link to public key
	pk_t *pk = nullptr;

	/// Boolean if the ciphertext is 0
	bool isnull;

	/// Constructors
	ciphertext_t() : c0(0), c1(0), pk(nullptr), isnull(true) {}
	ciphertext_t(ciphertext_t const &ct) : c0(ct.c0), c1(ct.c1) {
		if (ct.pk != nullptr) {
			pk = ct.pk;
		}
		isnull = ct.isnull;
	}
	template <typename T>
	ciphertext_t(T const &value) : c0(0), c1(0), pk(nullptr), isnull(true) {
		assert(value == 0);
	}
	template <typename T>
	ciphertext_t(pk_t &pk_in, T const &m) : c1(0), pk(&pk_in) {
		if (m == 0) {
			c0 = 0;
			isnull = true;
		} else {
			c0 = m;
			c0.ntt_pow_phi();
			//c0 = nfl::shoup(c0 * pk->delta, pk->delta_shoup);
			isnull = false;
		}
	}

	/// Assignment
	inline ciphertext_t &operator=(ciphertext_t const &ct) {
		c0 = ct.c0;
		c1 = ct.c1;
		if (ct.pk != nullptr) pk = ct.pk;
		isnull = ct.isnull;
		return *this;
	}
	template <typename Tp>
	inline ciphertext_t &operator=(Tp const &value) {
		if (value != 0) {
			assert(pk != nullptr);
			P v{value};
			v.ntt_pow_phi();
			isnull = false;
			c1 = 0;
			c0 = v; //nfl::shoup(v * pk->delta, pk->delta_shoup);
		} else {
			c0 = 0;
			c1 = 0;
			isnull = true;
		}
		return *this;
	}

	/// Additions/Substractions
	inline ciphertext_t &operator+=(ciphertext_t const &ct) {
		if (ct.isnull == false) {
			c0 = c0 + ct.c0;
			c1 = c1 + ct.c1;
			isnull = false;
		}
		return *this;
	}
	inline ciphertext_t &operator-=(ciphertext_t const &ct) {
		if (ct.isnull == false) {
			c0 = c0 - ct.c0;
			c1 = c1 - ct.c1;
			isnull = false;
		}
		return *this;
	}
	friend ciphertext_t operator+(ciphertext_t const &lhs,
			ciphertext_t const &rhs) {
		ciphertext_t ct(lhs);
		return ct += rhs;
	}
	friend ciphertext_t operator-(ciphertext_t const &lhs,
			ciphertext_t const &rhs) {
		ciphertext_t ct(lhs);
		return ct -= rhs;
	}

	inline ciphertext_t &operator+=(mpz_class const &value) {
		if (value == 0) {
			return *this;
		}
		P v{value};
		v.ntt_pow_phi();
		return *this += v;
	}
	inline ciphertext_t &operator-=(mpz_class const &value) {
		if (value == 0) {
			return *this;
		}
		P v{value};
		v.ntt_pow_phi();
		return *this -= v;
	}
	friend ciphertext_t operator+(ciphertext_t const &lhs, mpz_class const &rhs) {
		ciphertext_t ct(lhs);
		return ct += rhs;
	}
	friend ciphertext_t operator-(ciphertext_t const &lhs, mpz_class const &rhs) {
		ciphertext_t ct(lhs);
		return ct -= rhs;
	}

	inline ciphertext_t &operator+=(P::value_type const &value) {
		if (value == 0) {
			return *this;
		}
		P v{value};
		v.ntt_pow_phi();
		return *this += v;
	}
	inline ciphertext_t &operator-=(P::value_type const &value) {
		if (value == 0) {
			return *this;
		}
		P v{value};
		v.ntt_pow_phi();
		return *this -= v;
	}
	friend ciphertext_t operator+(ciphertext_t const &lhs,
			P::value_type const &rhs) {
		ciphertext_t ct(lhs);
		return ct += rhs;
	}
	friend ciphertext_t operator-(ciphertext_t const &lhs,
			P::value_type const &rhs) {
		ciphertext_t ct(lhs);
		return ct -= rhs;
	}


	/// Addition/Substraction of a polynomial
	inline ciphertext_t &operator+=(P const &p) {
		assert(pk != nullptr);
		c0 = c0 + p; // nfl::shoup(p * pk->delta, pk->delta_shoup);
		isnull = false;
		return *this;
	}
	inline ciphertext_t &operator-=(P const &p) {
		assert(pk != nullptr);
		c0 = c0 - p;// nfl::shoup(p * pk->delta, pk->delta_shoup);
		isnull = false;
		return *this;
	}
	friend ciphertext_t operator+(ciphertext_t const &lhs, P const &rhs) {
		ciphertext_t ct(lhs);
		return ct += rhs;
	}
	friend ciphertext_t operator-(ciphertext_t const &lhs, P const &rhs) {
		ciphertext_t ct(lhs);
		return ct -= rhs;
	}
};
}  // namespace BGV

/**
 * Encrypt a polynomial poly_m
 * @param ct     ciphertext (passed by reference)
 * @param pk     public key
 * @param poly_m polynomial to encrypt
 */
namespace BGV {
template <class PK, class C>
void encrypt_poly(C &ct, const PK &pk, parameters::poly_t &poly_m) {
	using P = parameters::poly_t;

	// Assume poly_m is in NTT domain
	//poly_m.ntt_pow_phi();

	P t = {plaintextModulus};
	t.ntt_pow_phi();

	// Generate a small u
	P u{nfl::non_uniform(2)};
	u.ntt_pow_phi();

	// Set the ciphertext pk
	ct.pk = (PK *)&pk;

	// Generate ct = (c0, c1)
	// where c0 = b*u + p*small error
	ct.c0 = nfl::non_uniform(2);
	ct.c0.ntt_pow_phi();
	ct.c0 = ct.c0 * t + nfl::shoup(u * pk.a, pk.a_shoup);

	// where c1 = a*u + p*small error + m
	ct.c1 = nfl::non_uniform(2);
	ct.c1.ntt_pow_phi();
	ct.c1 = ct.c1 * t + nfl::shoup(u * pk.b, pk.b_shoup) + poly_m;

	ct.isnull = false;
}
}  // namespace BGV

/**
 * Decryption of a ciphertext and recover the whole polynomial encrypted
 * @param poly_mpz pointer to the polynomial (already initialized)
 * @param sk       secret key
 * @param pk       public key
 * @param ct       ciphertext
 */
namespace BGV {
template <class SK, class PK, class C>
void decrypt_poly(parameters::poly_t &poly_m,
		const SK &sk, const PK &pk, const C &ct) {

	// Get the polynomial
	poly_m = ct.c1 - ct.c0 * sk.value;
	poly_m.invntt_pow_invphi();

	// Reduce the coefficients
	util::reduce(poly_m, plaintextModulus,
			parameters::poly_t::get_modulus(0), pk.qDivBy2);

}
}  // namespace BGV


/**
 * Return the log in base 2 of the noise in the ciphertext
 * @param  message message encrypted in the ciphertext
 * @param  sk      secret key
 * @param  pk      public key
 * @param  ct      ciphertext
 * @return         log_2(noise in ct)
 */
namespace BGV {

template <class SK, class PK, class C>
double noise(SK const &sk,PK const &pk, C const &ct) {
	using P = parameters::poly_t;

	P numerator = ct.c1 - ct.c0 * sk.value;
	numerator.invntt_pow_invphi();

	P::signed_value_type tmp, tmpMax=0, tmpMax2=0;
	double logMax = 0.0;
	double logMax2 = 0.0;

	for (auto &v : numerator) {
		util::center(tmp, v, P::get_modulus(0), pk.qDivBy2);
		logMax = std::max(logMax, std::log2(std::abs(tmp)));
		tmpMax = std::max(tmpMax, std::abs(tmp));
	}

	// std::cout << "inf noise: " << tmpMax << std::endl;

	return logMax;
}
}  // namespace BGV


namespace BGV {
using V = parameters::poly_t::value_type;
using S = typename parameters::poly_t::signed_value_type;
namespace util {

/**
 * Center op1 modulo op2
 * @param rop     result
 * @param op1     number op1 already reduced modulo op2, i.e. such that 0 <= op1
 * < op2-1
 * @param op2     modulus
 * @param op2Div2 floor(modulus/2)
 */
inline void center(S &rop, V const &op1, V const &op2,
		V const &op2Div2) {
	rop = (S)op1;
	if (op1 > op2Div2) {
		rop = rop - (S) op2;
	}
}

/**
 * Compute the quotient of op1 divided by op2 for a centered noise
 * @param rop     result
 * @param op1     numerator
 * @param op2     denominator
 * @param op2Div2 floor(denominator/2)
 */
inline void div_and_round(mpz_t &rop, mpz_t const &op1, mpz_t const &op2,
		mpz_t const &op2Div2) {
	mpz_t r;
	mpz_init2(r, mpz_size(op2) + 1);

	// Compute op1 = rop * op2 + r
	// where r has the same sign as op1
	mpz_tdiv_qr(rop, r, op1, op2);

	// Correct "rop" so that r is centered around 0
	long sgn = mpz_sgn(r);
	mpz_abs(r, r);
	if (mpz_cmp(r, op2Div2) >= 0) {
		if (sgn > 0) {
			mpz_add_ui(rop, rop, 1);
		} else {
			mpz_sub_ui(rop, rop, 1);
		}
	}

	// Clean
	mpz_clear(r);
}

/**
 * center the coefficients, multiply them by "multiplier" and then divide by the
 * divisor
 * @param coefficients pointer to the initialized coefficients
 * @param degree       number of coefficients to compute
 * @param multiplier   multiplier for the internal multiplication
 * @param divisor      denominator
 * @param divisorDiv2  floor(denominator/2)
 * @param mod_init     modulus of the initial coefficients
 * @param mod_initDiv2 floor(modulus/2)
 */

inline void reduce(parameters::poly_t &coefficients, V multiplier,
		V const &mod_init, V const &mod_initDiv2) {

	S tmp;
	for (auto &v : coefficients) {
		// Center with mod_init
		center(tmp, v, mod_init, mod_initDiv2);
		tmp = tmp % (S)multiplier;
		v =  (V)( tmp >= 0 ? tmp : (tmp+multiplier));
	}
}

/**
 * Lift the polynomial into an array of integer coefficients
 * @param coefficients pointer to the array of coefficients
 * @param c            polynomial P
 */
//inline void lift(std::array<mpz_t, parameters::polyZ_p::degree> &coefficients,
//		parameters::polyZ_p const &c) {
//	// Compute the inverse NTT
//	parameters::polyZ_p other{c};
//	other.invntt_pow_invphi();
//
//	// transform the poly into coefficients
//	other.poly2mpz(coefficients);
//}

}  // namespace util
}  // namespace BGV
