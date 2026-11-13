/*
 * verdec_bridge.cpp
 *
 * The NTL side of verdec_bridge.h: converts coefficient arrays into
 * verifiable-decryption-BGV's ZZ_pE and runs its prover and verifier, as its
 * main.cpp does, for one ciphertext (tau = 1, the homomorphic-tallying
 * parameters).
 */
#include <sys/random.h>
#include <cstdlib>
#include <stdexcept>
#include <NTL/ZZ.h>
#include "parameters.h"
#include "encryption.h"
#include "commitment.h"
#include "linearityproof.h"
#include "shortnessproof.h"
#include "verdec_bridge.h"

namespace verdec {

struct proof {
	Vec<Commitment> noiseCom;
	Vec<LinProof> linProof;
	ShortProof shortProof;
	ShortStatement shortstmt;
	proof() {
		noiseCom.SetLength(tau);
		linProof.SetLength(tau);
	}
};

static EncKey encKey;
static ComKey comKey;
static Commitment keyCom;

size_t degree() {
	return N;
}

static ZZ_pE to_ring(const uint64_t c[]) {
	ZZ_pX f;
	ZZ z;
	for (long i = 0; i < N; i++) {
		conv(z, (unsigned long) c[i]);
		SetCoeff(f, i, conv<ZZ_p>(z));
	}
	return to_ZZ_pE(f);
}

static Statement statement(const uint64_t u[], const uint64_t v[],
		const uint64_t m[]) {
	Statement stmt;
	stmt.ciphertexts[0].CTXu = to_ring(u);
	stmt.ciphertexts[0].CTXv = to_ring(v);
	stmt.messages[0] = to_ring(m);
	return stmt;
}

void setup(const uint64_t a[], const uint64_t b[], const uint64_t s[]) {
	static_assert(tau == 1, "the tally is one ciphertext");
	unsigned char seed[32];
	size_t got = 0;
	while (got < sizeof(seed)) {
		ssize_t n = getrandom(seed + got, sizeof(seed) - got, 0);
		if (n > 0) got += (size_t) n;
	}
	/* NTL's generator and rand() both start from fixed seeds otherwise, and
	 * the commitments' randomness and the proofs' masks come from them. */
	SetSeed(seed, sizeof(seed));
	unsigned int r;
	if (getrandom(&r, sizeof(r), 0) != sizeof(r)) {
		throw std::runtime_error("getrandom");
	}
	srand(r);

	ZZ_p::init(ZZ(q));
	ZZ_pE::init(ZZ_pX(INIT_MONO, N) + 1);
	encKey.PKa = to_ring(a);
	encKey.PKb = to_ring(b);
	encKey.SKs = to_ring(s);
	KeyGenCom(comKey);
	Commit(keyCom, comKey, encKey.SKs);
}

proof *prove(const uint64_t u[], const uint64_t v[], const uint64_t m[]) {
	Statement stmt = statement(u, v, m);
	proof *pi = new proof;
	ExtractNoise(pi->noiseCom, stmt, encKey, comKey);
	CreateProofStatement(pi->shortstmt, comKey, pi->noiseCom);
	ProveShortness(pi->shortProof, pi->shortstmt);
	ProveLinearRelations(pi->linProof, pi->noiseCom, comKey, keyCom,
			stmt.ciphertexts, stmt.messages);
	return pi;
}

bool verify(const proof *pi, const uint64_t u[], const uint64_t v[],
		const uint64_t m[]) {
	Statement stmt = statement(u, v, m);
	/* The statement of the proof of shortness, rebuilt from public data only:
	 * A from the commitment key, T from the noise commitments. The prover's
	 * copy also holds the openings (S), and trusting its A and T would let a
	 * prover show some other matrix short. CreateProofStatement reads the
	 * openings, so it is not used here; the commitments' com is all this
	 * function reads of them, here and in VerifyLinearRelations. */
	ShortStatement pub;
	pub.A[0][0] = comKey.B1[0]; pub.A[0][1] = comKey.B1[1];
	pub.A[0][2] = comKey.B1[2]; pub.A[0][3] = ZZ_pE(0);
	pub.A[1][0] = comKey.B2[0]; pub.A[1][1] = comKey.B2[1];
	pub.A[1][2] = comKey.B2[2]; pub.A[1][3] = ZZ_pE(1);
	for (long i = 0; i < tau; i++) {
		pub.T[0][i] = pi->noiseCom[i].com[0];
		pub.T[1][i] = pi->noiseCom[i].com[1];
	}
	bool bShort, bLin;
	VerifyShortness(bShort, pi->shortProof, pub);
	VerifyLinearRelations(bLin, pi->linProof, comKey, keyCom, pi->noiseCom,
			stmt.ciphertexts, stmt.messages);
	return bShort && bLin;
}

void release(proof *pi) {
	delete pi;
}

bool verify_forged(part which, const proof *pi, const uint64_t u[],
		const uint64_t v[], const uint64_t m[], bool *equations_hold) {
	Statement stmt = statement(u, v, m);
	proof forged = *pi;
	if (which == SHORTNESS) {
		/* A Z = T C + W with C = 0 and W = A Z: true for any T. */
		Mat<ZZ_pE> A;
		A.SetDims(2, 4);
		A[0][0] = comKey.B1[0]; A[0][1] = comKey.B1[1];
		A[0][2] = comKey.B1[2]; A[0][3] = ZZ_pE(0);
		A[1][0] = comKey.B2[0]; A[1][1] = comKey.B2[1];
		A[1][2] = comKey.B2[2]; A[1][3] = ZZ_pE(1);
		for (long i = 0; i < forged.shortProof.C.NumRows(); i++)
			for (long j = 0; j < forged.shortProof.C.NumCols(); j++)
				forged.shortProof.C[i][j] = ZZ_pE(0);
		MatrixMult(forged.shortProof.W, A, forged.shortProof.Z);
		Mat<ZZ_pE> AZ, T, TC, TCW;
		T.SetDims(2, tau);
		for (long i = 0; i < tau; i++) {
			T[0][i] = forged.noiseCom[i].com[0];
			T[1][i] = forged.noiseCom[i].com[1];
		}
		AZ.SetDims(2, lambda + 2); TC.SetDims(2, lambda + 2); TCW.SetDims(2, lambda + 2);
		MatrixMult(AZ, A, forged.shortProof.Z);
		MatrixMult(TC, T, forged.shortProof.C);
		MatrixAdd(TCW, TC, forged.shortProof.W);
		*equations_hold = (AZ == TCW);
	} else if (which == NOISE_ROW) {
		/* ProveShortness without rejection sampling, Y[3] made long. */
		ShortProof &sp = forged.shortProof;
		Mat<ZZ_pE> Y, SC, AZ, TC, TCW, C;
		Y.SetDims(4, lambda + 2); SC.SetDims(4, lambda + 2);
		SampleGaussian(Y[0], Sigma_A); SampleGaussian(Y[1], Sigma_A);
		SampleGaussian(Y[2], Sigma_A); SampleGaussian(Y[3], Sigma_A);
		ZZ_pX big;
		SetCoeff(big, 0, to_ZZ_p(power2_ZZ(37)));
		Y[3][0] += to_ZZ_pE(big);
		MatrixMult(sp.W, pi->shortstmt.A, Y);
		ShortnessChallengeFor(sp.C, pi->shortstmt, sp.W);
		MatrixMult(SC, pi->shortstmt.S, sp.C);
		MatrixAdd(sp.Z, SC, Y);
		AZ.SetDims(2, lambda + 2); TC.SetDims(2, lambda + 2); TCW.SetDims(2, lambda + 2);
		C.SetDims(tau, lambda + 2);
		MatrixMult(AZ, pi->shortstmt.A, sp.Z);
		MatrixMult(TC, pi->shortstmt.T, sp.C);
		MatrixAdd(TCW, TC, sp.W);
		ShortnessChallengeFor(C, pi->shortstmt, sp.W);
		*equations_hold = (AZ == TCW) && (C == sp.C);
	} else {
		/* With c = 0 the three equations are t1 = B1 z1, t2 = B1 z2 and
		 * u = alpha (z1[1] + B2[2] z1[2]) - (z2[1] + B2[2] z2[2]). */
		LinProof &lp = forged.linProof[0];
		const ZZ_pE alpha = -stmt.ciphertexts[0].CTXu / ZZ_p(p);
		lp.c = ZZ_pE(0);
		lp.t1 = lp.z1[0] + lp.z1[1] * comKey.B1[1] + lp.z1[2] * comKey.B1[2];
		lp.t2 = lp.z2[0] + lp.z2[1] * comKey.B1[1] + lp.z2[2] * comKey.B1[2];
		lp.u = alpha * (lp.z1[1] + lp.z1[2] * comKey.B2[2]) -
			(lp.z2[1] + lp.z2[2] * comKey.B2[2]);
		*equations_hold = true;	/* by construction, given c = 0 */
	}
	return verify(&forged, u, v, m);
}

/* ---------------------------------------------------------------------------
 * The proof as bytes. */

namespace {

const int Q_BITS = 62;			/* q < 2^62 */
const int RICE_CAP = 24;		/* unary runs past this escape to 64 bits */

struct bits_t {
	std::vector<uint8_t> *out;	/* writing */
	const uint8_t *in;			/* reading */
	size_t len, at;
	bool bad;
	void put(uint64_t v, int n) {
		for (int i = 0; i < n; i++, at++) {
			if (at / 8 >= out->size()) out->push_back(0);
			if ((v >> i) & 1) (*out)[at / 8] |= (uint8_t) (1u << (at % 8));
		}
	}
	uint64_t get(int n) {
		uint64_t v = 0;
		if ((at + n + 7) / 8 > len) { bad = true; return 0; }
		for (int i = 0; i < n; i++, at++)
			v |= (uint64_t) ((in[at / 8] >> (at % 8)) & 1) << i;
		return v;
	}
};

uint64_t residue(const ZZ_pE &a, long i) {
	return conv<unsigned long>(rep(coeff(rep(a), i)));
}

void put_uniform(bits_t &w, const ZZ_pE &a) {
	for (long i = 0; i < N; i++) w.put(residue(a, i), Q_BITS);
}

void get_uniform(bits_t &r, ZZ_pE &a) {
	ZZ_pX f;
	for (long i = 0; i < N; i++) {
		uint64_t v = r.get(Q_BITS);
		if (v >= (uint64_t) q) { r.bad = true; v = 0; }
		SetCoeff(f, i, conv<ZZ_p>(conv<ZZ>((unsigned long) v)));
	}
	a = to_ZZ_pE(f);
}

/* Rice parameter for a discrete Gaussian of width sigma: the mean of the
 * zigzagged value is about 1.6 sigma. */
int rice_k(double sigma) {
	int k = 0;
	while ((1.5957691 * sigma) / 2 >= (double) (1ULL << k)) k++;
	return k;
}

void put_gauss(bits_t &w, const ZZ_pE &a, int k) {
	for (long i = 0; i < N; i++) {
		uint64_t v = residue(a, i);
		int64_t x = (v > (uint64_t) q / 2) ? -(int64_t) ((uint64_t) q - v) : (int64_t) v;
		uint64_t u = ((uint64_t) x << 1) ^ (uint64_t) (x >> 63);
		uint64_t hi = u >> k;
		if (hi >= (uint64_t) RICE_CAP) {
			w.put((1ULL << RICE_CAP) - 1, RICE_CAP);
			w.put(u, 64);
			continue;
		}
		w.put((1ULL << hi) - 1, (int) hi);
		w.put(0, 1);
		w.put(u & ((1ULL << k) - 1), k);
	}
}

void get_gauss(bits_t &r, ZZ_pE &a, int k) {
	ZZ_pX f;
	for (long i = 0; i < N; i++) {
		int hi = 0;
		uint64_t u;
		while (hi < RICE_CAP && r.get(1) == 1) hi++;
		if (hi == RICE_CAP) {
			u = r.get(64);
			if ((u >> k) < (uint64_t) RICE_CAP) r.bad = true;	/* needless escape */
		} else {
			u = ((uint64_t) hi << k) | r.get(k);
		}
		int64_t x = (int64_t) (u >> 1) ^ -(int64_t) (u & 1);
		if (x > (int64_t) (q / 2) || x < -(int64_t) (q / 2)) { r.bad = true; x = 0; }
		uint64_t v = (x < 0) ? (uint64_t) q - (uint64_t) (-x) : (uint64_t) x;
		SetCoeff(f, i, conv<ZZ_p>(conv<ZZ>((unsigned long) v)));
	}
	a = to_ZZ_pE(f);
}

bool put_chall(bits_t &w, const ZZ_pE &c) {
	int set = 0;
	for (long i = 0; i < N; i++) {
		uint64_t v = residue(c, i);
		if (v == 0) continue;
		if (v != 1 && v != (uint64_t) q - 1) return false;
		w.put((uint64_t) i, 12);
		w.put(v == 1 ? 1 : 0, 1);
		set++;
	}
	return set == 36;
}

void get_chall(bits_t &r, ZZ_pE &c) {
	ZZ_pX f;
	long last = -1;
	for (int j = 0; j < 36; j++) {
		long pos = (long) r.get(12);
		long sign = r.get(1) ? 1 : -1;
		if (pos <= last) { r.bad = true; continue; }
		last = pos;
		SetCoeff(f, pos, sign);
	}
	c = to_ZZ_pE(f);
}

}

std::vector<uint8_t> encode(const proof *pi) {
	std::vector<uint8_t> out;
	bits_t w{&out, nullptr, 0, 0, false};
	const int kc = rice_k(Sigma_C), ka = rice_k(Sigma_A);
	bool ok = true;

	static_assert(N == 4096, "12-bit positions");
	for (long i = 0; i < tau; i++) {
		put_uniform(w, pi->noiseCom[i].com[0]);
		put_uniform(w, pi->noiseCom[i].com[1]);
	}
	for (long i = 0; i < tau; i++) {
		const LinProof &lp = pi->linProof[i];
		put_uniform(w, lp.t1);
		put_uniform(w, lp.t2);
		put_uniform(w, lp.u);
		ok &= put_chall(w, lp.c);
		for (long j = 0; j < 3; j++) put_gauss(w, lp.z1[j], kc);
		for (long j = 0; j < 3; j++) put_gauss(w, lp.z2[j], kc);
	}
	const ShortProof &sp = pi->shortProof;
	for (long i = 0; i < sp.W.NumRows(); i++)
		for (long j = 0; j < sp.W.NumCols(); j++) put_uniform(w, sp.W[i][j]);
	for (long i = 0; i < sp.C.NumRows(); i++)
		for (long j = 0; j < sp.C.NumCols(); j++) {
			uint64_t b = residue(sp.C[i][j], 0);
			ok &= deg(rep(sp.C[i][j])) <= 0 && b <= 1;
			w.put(b, 1);
		}
	for (long i = 0; i < sp.Z.NumRows(); i++)
		for (long j = 0; j < sp.Z.NumCols(); j++) put_gauss(w, sp.Z[i][j], ka);
	if (!ok) out.clear();
	return out;
}

proof *decode(const uint8_t *buf, size_t len) {
	bits_t r{nullptr, buf, len, 0, false};
	const int kc = rice_k(Sigma_C), ka = rice_k(Sigma_A);
	proof *pi = new proof;

	for (long i = 0; i < tau; i++) {
		get_uniform(r, pi->noiseCom[i].com[0]);
		get_uniform(r, pi->noiseCom[i].com[1]);
	}
	for (long i = 0; i < tau; i++) {
		LinProof &lp = pi->linProof[i];
		get_uniform(r, lp.t1);
		get_uniform(r, lp.t2);
		get_uniform(r, lp.u);
		get_chall(r, lp.c);
		for (long j = 0; j < 3; j++) get_gauss(r, lp.z1[j], kc);
		for (long j = 0; j < 3; j++) get_gauss(r, lp.z2[j], kc);
	}
	ShortProof &sp = pi->shortProof;
	for (long i = 0; i < sp.W.NumRows(); i++)
		for (long j = 0; j < sp.W.NumCols(); j++) get_uniform(r, sp.W[i][j]);
	for (long i = 0; i < sp.C.NumRows(); i++)
		for (long j = 0; j < sp.C.NumCols(); j++) sp.C[i][j] = ZZ_pE((long) r.get(1));
	for (long i = 0; i < sp.Z.NumRows(); i++)
		for (long j = 0; j < sp.Z.NumCols(); j++) get_gauss(r, sp.Z[i][j], ka);
	if (!r.bad && ((r.at + 7) / 8 != len ||
			(r.at % 8 != 0 && (buf[r.at / 8] >> (r.at % 8)) != 0))) r.bad = true;
	if (r.bad) {
		delete pi;
		return nullptr;
	}
	return pi;
}

}
