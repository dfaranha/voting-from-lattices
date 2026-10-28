/*
 * zq.hpp
 *
 * Arithmetic in R_q = Z_q[X]/(X^N + 1) for q = 4611686018326724609, the prime
 * this code took from NFLlib, which it no longer uses.
 *
 * zq::poly implements the part of NFLlib's poly that the code here calls, with
 * the same meaning, so that BGV.hpp, BDLOP.hpp and vericrypt.cpp did not have
 * to change: a polynomial is N coefficients in [0, q), kept in one 32-byte
 * aligned array, so sizeof(poly) and the bytes the hashes absorb are as before.
 * As in NFLlib, + - * act coefficient-wise, so * is the ring product only on
 * polynomials in NTT form (ntt_pow_phi).
 *
 * Two departures, both deliberate:
 *   - == and != compare every coefficient and return a bool. NFLlib's built a
 *     polynomial of per-coefficient results and converted it to bool as "any
 *     coefficient non-zero", which made != right but == true as soon as one
 *     coefficient agreed.
 *   - uniform() and non_uniform(b) sample exactly, by rejection. NFLlib masked
 *     and subtracted once, which is slightly biased.
 * The NTT is a different, equally valid transform from NFLlib's: the slots come
 * out in another order. Nothing here depends on that order: the automorphisms
 * in util.hpp act on coefficients.
 *
 * Randomness is SHAKE256 keyed with 32 bytes from getrandom(2), drawn as a
 * stream. Single-threaded.
 */
#ifndef ZQ_HPP
#define ZQ_HPP

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <algorithm>
#include <initializer_list>
#include <iterator>
#include <stdexcept>
#include <gmpxx.h>
#include <sys/random.h>

extern "C" {
#include "fips202.h"
}

/* NFLlib's sampler tags, which the code passes to construct a polynomial. */
namespace nfl {
struct uniform {};
struct non_uniform {
	uint64_t upper_bound;
	explicit non_uniform(uint64_t ub) : upper_bound(ub) {}
};
}

namespace zq {

typedef unsigned __int128 u128;

static constexpr uint64_t Q = 4611686018326724609ULL;	/* < 2^62, 2^25 | Q - 1 */

/* a * b mod Q for a, b < Q, by Barrett reduction: with k = 62 and
 * mu = floor(2^124 / Q), the estimate (x >> 61) * mu >> 63 is at most two
 * below floor(x / Q), so two conditional subtractions finish it. */
static inline uint64_t mulmod(uint64_t a, uint64_t b) {
	static const uint64_t MU = (uint64_t) ((((u128) 1) << 124) / Q);
	u128 x = (u128) a * b;
	uint64_t q1 = (uint64_t) (x >> 61);
	uint64_t q3 = (uint64_t) (((u128) q1 * MU) >> 63);
	uint64_t r = (uint64_t) x - q3 * Q;
	r -= (r >= Q) ? Q : 0;
	r -= (r >= Q) ? Q : 0;
	return r;
}

static inline uint64_t addmod(uint64_t a, uint64_t b) {
	uint64_t r = a + b;
	return r - ((r >= Q) ? Q : 0);
}

static inline uint64_t submod(uint64_t a, uint64_t b) {
	return (a >= b) ? a - b : a + Q - b;
}

static inline uint64_t powmod(uint64_t a, uint64_t e) {
	uint64_t r = 1;
	while (e) {
		if (e & 1) r = mulmod(r, a);
		a = mulmod(a, a);
		e >>= 1;
	}
	return r;
}

/* w * b mod Q for a fixed w < Q, with wp = floor(w * 2^64 / Q): Shoup. */
static inline uint64_t mulshoup(uint64_t w, uint64_t wp, uint64_t b) {
	uint64_t qh = (uint64_t) (((u128) wp * b) >> 64);
	uint64_t r = w * b - qh * Q;
	return r - ((r >= Q) ? Q : 0);
}

static inline uint64_t shoup_const(uint64_t w) {
	return (uint64_t) ((((u128) w) << 64) / Q);
}

/* The stream every sampler draws from: SHAKE256 keyed from the system. */
class prng {
	shake256incctx st;
	bool keyed = false;
	void key() {
		uint8_t seed[32];
		size_t got = 0;
		while (got < sizeof(seed)) {
			ssize_t n = getrandom(seed + got, sizeof(seed) - got, 0);
			if (n > 0) got += (size_t) n;
		}
		shake256_inc_init(&st);
		shake256_inc_absorb(&st, seed, sizeof(seed));
		shake256_inc_finalize(&st);
		keyed = true;
	}
public:
	void bytes(uint8_t *out, size_t len) {
		if (!keyed) key();
		shake256_inc_squeeze(out, len, &st);
	}
	static prng &get() {
		static prng p;
		return p;
	}
};

/* Twiddle factors of the negacyclic NTT of degree N: psi a primitive 2N-th
 * root of unity, found from the smallest quadratic non-residue, and its
 * powers in bit-reversed order with their Shoup constants. */
template <size_t N> struct ntt_tables {
	uint64_t w[N], wp[N];		/* psi^brv(i), forward */
	uint64_t iw[N], iwp[N];		/* -psi^brv(i), inverse */
	uint64_t ninv, ninvp;
	ntt_tables() {
		static_assert((N & (N - 1)) == 0, "N a power of two");
		uint64_t g = 2;
		while (powmod(g, (Q - 1) / 2) != Q - 1) g++;
		uint64_t psi = powmod(g, (Q - 1) / (2 * N));
		if (powmod(psi, N) != Q - 1) {
			throw std::runtime_error("zq: no primitive 2N-th root of unity");
		}
		unsigned logn = 0;
		while ((1UL << logn) < N) logn++;
		for (size_t i = 0; i < N; i++) {
			size_t r = 0;
			for (unsigned b = 0; b < logn; b++) r |= ((i >> b) & 1) << (logn - 1 - b);
			w[i] = powmod(psi, r);
			wp[i] = shoup_const(w[i]);
			iw[i] = (Q - w[i]) % Q;
			iwp[i] = shoup_const(iw[i]);
		}
		ninv = powmod(N, Q - 2);
		ninvp = shoup_const(ninv);
	}
	static ntt_tables const &get() {
		static ntt_tables t;
		return t;
	}
};

template <size_t Degree> class poly {
public:
	typedef uint64_t value_type;
	typedef int64_t signed_value_type;
	typedef u128 greater_value_type;
	typedef value_type *iterator;
	typedef value_type const *const_iterator;
	static constexpr size_t degree = Degree;
	static constexpr size_t nmoduli = 1;

private:
	alignas(32) value_type _data[Degree];

public:
	poly() { std::fill(begin(), end(), 0); }
	poly(nfl::uniform const &) { set_uniform(); }
	poly(nfl::non_uniform const &m) { set_bounded(m.upper_bound); }
	poly(value_type v, bool reduce_coeffs = true) { set(v, reduce_coeffs); }
	poly(std::initializer_list<value_type> values, bool reduce_coeffs = true) {
		set(values, reduce_coeffs);
	}
	/* The first coefficients from integers of any size, as NFLlib's
	 * poly(std::initializer_list<mpz_class>): reduced into [0, q). */
	poly(std::initializer_list<mpz_class> values) {
		if (values.size() > Degree) {
			throw std::runtime_error("zq: initializer longer than the degree");
		}
		size_t i = 0;
		for (auto const &v : values) {
			_data[i++] = mpz_fdiv_ui(v.get_mpz_t(), Q);
		}
		for (; i < Degree; i++) _data[i] = 0;
	}
	template <class It, class = typename std::iterator_traits<It>::value_type>
	poly(It first, It last, bool reduce_coeffs = true) {
		set(first, last, reduce_coeffs);
	}

	/* Zero, or the constant v, as NFLlib's set(value_type). */
	void set(value_type v, bool reduce_coeffs = true) {
		if (v == 0) {
			std::fill(begin(), end(), 0);
		} else {
			set({v}, reduce_coeffs);
		}
	}
	void set(std::initializer_list<value_type> values, bool reduce_coeffs = true) {
		set(values.begin(), values.end(), reduce_coeffs);
	}
	/* The first values, the rest zero; NFLlib's set(first, last). */
	template <class It> void set(It first, It last, bool reduce_coeffs = true) {
		size_t size = std::distance(first, last);
		if (size > Degree) {
			throw std::runtime_error("zq: initializer longer than the degree");
		}
		size_t i = 0;
		for (It it = first; it != last; ++it, ++i) {
			value_type v = *it;
			_data[i] = reduce_coeffs ? v % Q : v;
		}
		for (; i < Degree; i++) _data[i] = 0;
	}

	/* Uniform in [0, q). q > 2^61, so a 62-bit draw is kept with probability
	 * above 1/2; in fact q is within 2^37 of 2^62, so almost always. */
	void set_uniform() {
		value_type buf[Degree];
		size_t i = 0;
		while (i < Degree) {
			prng::get().bytes((uint8_t *) buf, sizeof(buf));
			for (size_t j = 0; j < Degree && i < Degree; j++) {
				value_type v = buf[j] & ((1ULL << 62) - 1);
				if (v < Q) _data[i++] = v;
			}
		}
	}
	/* Uniform in [-(ub - 1), ub - 1], as NFLlib's non_uniform(ub). */
	void set_bounded(uint64_t ub) {
		if (ub == 0 || ub >= Q) {
			throw std::runtime_error("zq: bound out of range");
		}
		const uint64_t range = 2 * ub - 1;
		uint64_t mask = 1;
		while (mask < range) mask = (mask << 1) | 1;
		value_type buf[Degree];
		size_t i = 0;
		while (i < Degree) {
			prng::get().bytes((uint8_t *) buf, sizeof(buf));
			for (size_t j = 0; j < Degree && i < Degree; j++) {
				uint64_t v = buf[j] & mask;
				if (v >= range) continue;
				/* v in [0, 2ub - 2]; centre it on 0. */
				_data[i++] = (v < ub) ? v : Q + v - range;
			}
		}
	}

	iterator begin() { return _data; }
	iterator end() { return _data + Degree; }
	const_iterator begin() const { return _data; }
	const_iterator end() const { return _data + Degree; }
	const_iterator cbegin() const { return _data; }
	const_iterator cend() const { return _data + Degree; }
	value_type *data() { return _data; }
	value_type const *data() const { return _data; }
	value_type &operator()(size_t, size_t i) { return _data[i]; }
	value_type const &operator()(size_t, size_t i) const { return _data[i]; }

	static constexpr value_type get_modulus(size_t) { return Q; }
	static mpz_srcptr moduli_product() {
		static struct m {
			mpz_t v;
			m() { mpz_init_set_ui(v, Q); }
		} mod;
		return mod.v;
	}

	template <class A> void poly2mpz(A &out) const {
		for (size_t i = 0; i < Degree; i++) mpz_set_ui(out[i], _data[i]);
	}
	template <class A> void mpz2poly(A const &in) {
		mpz_t t;
		mpz_init(t);
		for (size_t i = 0; i < Degree; i++) {
			mpz_mod_ui(t, in[i], Q);
			_data[i] = mpz_get_ui(t);
		}
		mpz_clear(t);
	}

	/* Forward negacyclic NTT, in place: Cooley-Tukey, psi folded in. */
	void ntt_pow_phi() {
		auto const &T = ntt_tables<Degree>::get();
		size_t k = 0;
		for (size_t len = Degree / 2; len >= 1; len >>= 1) {
			for (size_t start = 0; start < Degree; start += 2 * len) {
				k++;
				const uint64_t w = T.w[k], wp = T.wp[k];
				for (size_t j = start; j < start + len; j++) {
					uint64_t t = mulshoup(w, wp, _data[j + len]);
					_data[j + len] = submod(_data[j], t);
					_data[j] = addmod(_data[j], t);
				}
			}
		}
	}
	/* Its inverse: Gentleman-Sande, then the factor 1/N. */
	void invntt_pow_invphi() {
		auto const &T = ntt_tables<Degree>::get();
		size_t k = Degree;
		for (size_t len = 1; len < Degree; len <<= 1) {
			for (size_t start = 0; start < Degree; start += 2 * len) {
				k--;
				const uint64_t w = T.iw[k], wp = T.iwp[k];
				for (size_t j = start; j < start + len; j++) {
					uint64_t t = _data[j];
					_data[j] = addmod(t, _data[j + len]);
					_data[j + len] = mulshoup(w, wp, submod(t, _data[j + len]));
				}
			}
		}
		for (size_t j = 0; j < Degree; j++) {
			_data[j] = mulshoup(T.ninv, T.ninvp, _data[j]);
		}
	}

	poly &operator+=(poly const &o) {
		for (size_t i = 0; i < Degree; i++) _data[i] = addmod(_data[i], o._data[i]);
		return *this;
	}
	poly &operator-=(poly const &o) {
		for (size_t i = 0; i < Degree; i++) _data[i] = submod(_data[i], o._data[i]);
		return *this;
	}
	poly &operator*=(poly const &o) {
		for (size_t i = 0; i < Degree; i++) _data[i] = mulmod(_data[i], o._data[i]);
		return *this;
	}
	friend poly operator+(poly a, poly const &b) { return a += b; }
	friend poly operator-(poly a, poly const &b) { return a -= b; }
	friend poly operator*(poly a, poly const &b) { return a *= b; }
	friend poly operator-(poly const &a) { return poly() - a; }
	friend bool operator==(poly const &a, poly const &b) {
		return std::memcmp(a._data, b._data, sizeof(a._data)) == 0;
	}
	friend bool operator!=(poly const &a, poly const &b) { return !(a == b); }
};

}	/* namespace zq */

/* The rest of the NFLlib names the code calls. */
namespace nfl {

/* NFLlib's shoup(a * b, b_shoup) is a * b computed with b's precomputed
 * quotients; compute_shoup produced them. Plain products serve here. */
template <class P> inline P const &compute_shoup(P const &p) { return p; }
template <class P> inline P shoup(P const &product, P const &) { return product; }

inline void fastrandombytes(unsigned char *r, unsigned long long rlen) {
	zq::prng::get().bytes(r, (size_t) rlen);
}

}	/* namespace nfl */

#endif
