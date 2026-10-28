/*
 * election.cpp
 *
 * A cumulative-voting election with homomorphic tallying, end to end: NV
 * voters each spread LS votes over NC candidates, at most LC to one; each
 * encrypts the ballot with BGV and proves it valid (vericrypt); the board
 * checks every proof; the accepted ciphertexts are added; the decryptor
 * decrypts the sum and proves the decryption correct with
 * ../verifiable-decryption-BGV; and anyone checks that proof and reads the
 * tally. It also checks that invalid ballots, a ballot changed after its
 * proof, and a decryption proof moved to another tally or ciphertext are all
 * rejected. Exits non-zero if any check fails.
 *
 * Build and run with make election && ./election; the sizes are params.hpp's
 * NV, NC, LC and LS, which make CONFIG="-DNV=... " overrides.
 */
#include <array>
#include <chrono>
#include <cstring>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <vector>

/* vericrypt.hpp and the headers it includes name std:: types unqualified,
 * as main.cpp does. */
using namespace std;

#include "vericrypt.hpp"
#include "verdec_bridge.h"

static_assert(LS < NC * LC, "a valid ballot, and one vote over the budget, must exist");
static_assert(NC <= P::degree, "one coefficient per candidate");

/* ---------------------------------------------------------------------------
 * Reporting
 */
static int failures = 0, checks = 0;

static void check(const string &what, bool ok) {
	string line = "Testing if " + what + "...";
	if (line.size() < 78) line.append(78 - line.size(), ' ');
	cout << line << (ok ? "[PASS]" : "[FAIL]") << endl;
	failures += !ok;
	checks++;
}

static double seconds_since(chrono::steady_clock::time_point t0) {
	return chrono::duration<double>(chrono::steady_clock::now() - t0).count();
}

/* ---------------------------------------------------------------------------
 * Ballots: one coefficient per candidate, in the coefficient domain, which is
 * how vericrypt takes the message.
 */
typedef array<uint64_t, NC> votes_t;

static mt19937_64 &rng() {
	static mt19937_64 g{random_device{}()};
	return g;
}

/* LS votes, one at a time, each to a candidate that has fewer than LC. */
static votes_t valid_votes() {
	votes_t v{};
	uniform_int_distribution<size_t> pick(0, NC - 1);
	for (int cast = 0; cast < LS; cast++) {
		size_t j;
		do {
			j = pick(rng());
		} while (v[j] >= LC);
		v[j]++;
	}
	return v;
}

static P ballot(const votes_t &v) {
	P m{0};
	copy(v.begin(), v.end(), m.begin());
	return m;
}

/* Proves a ballot and checks a copy of the proof, since verify() moves the
 * proof's responses into NTT form in place. */
static bool cast_and_check(proof_t &out, const votes_t &v, const pk_t &pk,
		const comm_key_t &ck) {
	vericrypt(out, pk, ck, ballot(v));
	auto copy = make_unique<proof_t>(out);
	return verify(*copy, ck, pk);
}

/* ---------------------------------------------------------------------------
 * Coefficient arrays for the decryption proof, in [0, q).
 */
typedef vector<uint64_t> coeffs_t;

static coeffs_t coeffs(P p, bool ntt_form) {
	if (ntt_form) p.invntt_pow_invphi();
	return coeffs_t(p.begin(), p.end());
}

int main() {
	cout << "\n** End-to-end test of cumulative voting with homomorphic tallying:\n\n";
	cout << "    " << NV << " voters, " << NC << " candidates, at most " << LC
		<< " votes to one, " << LS << " votes each\n";
	cout << "    ballot proof: " << sizeof(proof_t) / 1024.0 << " KiB in memory\n\n";

	if (verdec::degree() != P::degree) {
		cerr << "the decryption proof is for another degree" << endl;
		return 1;
	}

	/* Keys. */
	auto t0 = chrono::steady_clock::now();
	sk_t sk;
	pk_t pk(sk);
	auto ck = make_unique<comm_key_t>();
	verdec::setup(coeffs(pk.a, true).data(), coeffs(pk.b, true).data(),
			coeffs(sk.value, true).data());
	double t_keys = seconds_since(t0);

	/* Voting and the board's check. */
	double t_cast = 0, t_check = 0;
	vector<uint64_t> expected(NC, 0);
	BGV::ciphertext_t tally;
	bool all_accepted = true;
	auto pf = make_unique<proof_t>();
	for (int i = 0; i < NV; i++) {
		votes_t v = valid_votes();
		t0 = chrono::steady_clock::now();
		vericrypt(*pf, pk, *ck, ballot(v));
		t_cast += seconds_since(t0);

		t0 = chrono::steady_clock::now();
		auto copy = make_unique<proof_t>(*pf);
		bool ok = verify(*copy, *ck, pk);
		t_check += seconds_since(t0);

		all_accepted &= ok;
		if (ok) {
			tally += pf->ct;
			for (size_t j = 0; j < NC; j++) expected[j] += v[j];
		}
	}
	check("the board accepts every voter's ballot proof", all_accepted);

	/* Invalid ballots: the prover runs, but the board must refuse the proof. */
	{
		votes_t v = valid_votes();
		/* One more vote than the budget, to a candidate with room for it. */
		size_t j = 0;
		while (v[j] >= LC) j++;
		v[j]++;
		check("a ballot over the vote budget is rejected",
				!cast_and_check(*pf, v, pk, *ck));
	}
	{
		/* LC + 1 votes to one candidate, the budget kept where it can be. */
		votes_t v{};
		v[0] = LC + 1;
		int left = LS - (LC + 1);
		for (size_t j = 1; j < NC && left > 0; j++) {
			int give = min(left, LC);
			v[j] = give;
			left -= give;
		}
		check("a ballot with more than LC votes to one candidate is rejected",
				!cast_and_check(*pf, v, pk, *ck));
	}
	{
		/* A valid proof, then its ciphertext changed: one vote moved. */
		votes_t v = valid_votes();
		bool honest = cast_and_check(*pf, v, pk, *ck);
		P delta{0};
		delta.data()[0] = 1;
		delta.ntt_pow_phi();
		pf->ct.c1 = pf->ct.c1 + delta;
		auto copy = make_unique<proof_t>(*pf);
		check("a ballot changed after its proof is rejected",
				honest && !verify(*copy, *ck, pk));
	}

	/* Decryption. */
	t0 = chrono::steady_clock::now();
	P result;
	BGV::decrypt_poly(result, sk, pk, tally);
	double t_decrypt = seconds_since(t0);
	bool tally_ok = true;
	for (size_t j = 0; j < P::degree; j++) {
		tally_ok &= result.data()[j] == (j < NC ? expected[j] : 0);
	}
	check("the decrypted tally is the votes cast", tally_ok);

	coeffs_t u = coeffs(tally.c0, true), w = coeffs(tally.c1, true);
	coeffs_t m = coeffs(result, false);

	t0 = chrono::steady_clock::now();
	verdec::proof *dp = verdec::prove(u.data(), w.data(), m.data());
	double t_dprove = seconds_since(t0);

	t0 = chrono::steady_clock::now();
	bool dec_ok = verdec::verify(dp, u.data(), w.data(), m.data());
	double t_dverify = seconds_since(t0);
	check("the proof of correct decryption verifies", dec_ok);

	{
		coeffs_t other = m;
		other[0] = (other[0] + 1) % plaintextModulus;
		check("the decryption proof does not verify for another tally",
				!verdec::verify(dp, u.data(), w.data(), other.data()));
	}
	{
		/* The same plaintext under another ciphertext: the tally plus a fresh
		 * encryption of zero. */
		BGV::ciphertext_t zero;
		P none{0};
		BGV::encrypt_poly(zero, pk, none);
		BGV::ciphertext_t moved = tally;
		moved += zero;
		coeffs_t u2 = coeffs(moved.c0, true), w2 = coeffs(moved.c1, true);
		check("the decryption proof does not verify for another ciphertext",
				!verdec::verify(dp, u2.data(), w2.data(), m.data()));
	}
	{
		bool eq;
		bool accepted = verdec::verify_forged(verdec::SHORTNESS, dp, u.data(),
				w.data(), m.data(), &eq);
		check("a shortness proof with a challenge of its own is rejected",
				eq && !accepted);
		accepted = verdec::verify_forged(verdec::LINEAR, dp, u.data(), w.data(),
				m.data(), &eq);
		check("a linear proof with a challenge of its own is rejected",
				eq && !accepted);
	}
	verdec::release(dp);

	cout << "\n    tally:";
	for (size_t j = 0; j < NC; j++) cout << " " << result.data()[j];
	cout << "\n    log2 of the tally's noise: " << BGV::noise(sk, pk, tally) << "\n";
	cout << "    timing, wall clock: keys " << t_keys << " s; voting "
		<< t_cast << " s, " << 1000 * t_cast / NV << " ms per ballot; board check "
		<< t_check << " s, " << 1000 * t_check / NV << " ms per ballot\n";
	cout << "    timing: decryption " << 1000 * t_decrypt << " ms; its proof "
		<< t_dprove << " s to prove, " << t_dverify << " s to verify\n\n";
	cout << (failures ? "FAILED: " : "All checks passed: ") << checks - failures
		<< " of " << checks << "\n";
	return failures != 0;
}
