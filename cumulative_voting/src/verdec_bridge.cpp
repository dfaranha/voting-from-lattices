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

}
