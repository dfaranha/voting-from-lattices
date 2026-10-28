/*
 * election.cpp
 *
 * A cumulative-voting election with homomorphic tallying, end to end: NV
 * voters each spread LS votes over NC candidates, at most LC to one; each
 * encrypts the ballot with BGV and proves it valid (vericrypt); the board
 * checks every proof; the accepted ciphertexts are added; the decryptor
 * decrypts the sum and proves the decryption correct with
 * ../verifiable-decryption-BGV; and anyone checks that proof and reads the
 * tally. Each ballot proof is bound to the election and the voter. It also
 * checks that invalid ballots, a ballot changed after its proof, a ballot
 * proof moved to another voter or election, a second ballot from one voter,
 * a decryption proof moved to another tally or ciphertext, and decryption
 * proofs whose challenges the prover chose itself are all rejected. Exits non-zero if any check fails.
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

/* ---------------------------------------------------------------------------
 * The context a ballot proof is bound to: the election's identifier and the
 * voter's index on the roll.
 */
typedef array<uint8_t, 16> election_id_t;

static election_id_t new_election_id() {
	election_id_t id;
	random_device rd;
	for (auto &b : id) b = (uint8_t) rd();
	return id;
}

static vector<uint8_t> context(const election_id_t &eid, uint32_t voter) {
	static const char label[] = "cumulative voting ballot";
	vector<uint8_t> c(label, label + sizeof(label));
	c.insert(c.end(), eid.begin(), eid.end());
	for (int i = 0; i < 4; i++) c.push_back((uint8_t) (voter >> (8 * i)));
	return c;
}

/* Checks a copy of the proof, since verify() moves the proof's responses into
 * NTT form in place. */
static bool check_proof(const proof_t &pf, const pk_t &pk, const comm_key_t &ck,
		const election_id_t &eid, uint32_t voter) {
	auto copy = make_unique<proof_t>(pf);
	vector<uint8_t> c = context(eid, voter);
	return verify(*copy, ck, pk, c.data(), c.size());
}

static void prove(proof_t &out, const votes_t &v, const pk_t &pk,
		const comm_key_t &ck, const election_id_t &eid, uint32_t voter) {
	vector<uint8_t> c = context(eid, voter);
	vericrypt(out, pk, ck, ballot(v), c.data(), c.size());
}

/* The bulletin board: one ballot per voter on the roll, each checked against
 * that voter's context, added to the tally if accepted. */
struct board_t {
	election_id_t eid;
	vector<bool> voted = vector<bool>(NV, false);
	BGV::ciphertext_t tally;

	bool cast(uint32_t voter, const proof_t &pf, const pk_t &pk,
			const comm_key_t &ck) {
		if (voter >= NV || voted[voter]) return false;
		if (!check_proof(pf, pk, ck, eid, voter)) return false;
		voted[voter] = true;
		tally += pf.ct;
		return true;
	}
};

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
	board_t board;
	board.eid = new_election_id();
	bool all_accepted = true;
	auto pf = make_unique<proof_t>();
	for (uint32_t i = 0; i < NV; i++) {
		votes_t v = valid_votes();
		t0 = chrono::steady_clock::now();
		prove(*pf, v, pk, *ck, board.eid, i);
		t_cast += seconds_since(t0);

		t0 = chrono::steady_clock::now();
		bool ok = board.cast(i, *pf, pk, *ck);
		t_check += seconds_since(t0);

		all_accepted &= ok;
		if (ok) {
			for (size_t j = 0; j < NC; j++) expected[j] += v[j];
		}
	}
	check("the board accepts every voter's ballot proof", all_accepted);
	BGV::ciphertext_t &tally = board.tally;

	/* The last voter's proof, elsewhere. */
	{
		check("a ballot proof moved to another voter is rejected",
				!check_proof(*pf, pk, *ck, board.eid, NV > 1 ? 0 : 1));
		election_id_t other = board.eid;
		other[0] ^= 1;
		check("a ballot proof moved to another election is rejected",
				!check_proof(*pf, pk, *ck, other, NV - 1));
		/* A second, valid ballot from a voter who has voted. */
		auto again = make_unique<proof_t>();
		prove(*again, valid_votes(), pk, *ck, board.eid, 0);
		check("the board rejects a second ballot from the same voter",
				check_proof(*again, pk, *ck, board.eid, 0) &&
				!board.cast(0, *again, pk, *ck));
	}

	/* Invalid ballots: the prover runs, but the board must refuse the proof. */
	{
		votes_t v = valid_votes();
		/* One more vote than the budget, to a candidate with room for it. */
		size_t j = 0;
		while (v[j] >= LC) j++;
		v[j]++;
		prove(*pf, v, pk, *ck, board.eid, NV);
		check("a ballot over the vote budget is rejected",
				!check_proof(*pf, pk, *ck, board.eid, NV));
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
		prove(*pf, v, pk, *ck, board.eid, NV);
		check("a ballot with more than LC votes to one candidate is rejected",
				!check_proof(*pf, pk, *ck, board.eid, NV));
	}
	{
		/* A valid proof, then its ciphertext changed: one vote moved. */
		votes_t v = valid_votes();
		prove(*pf, v, pk, *ck, board.eid, NV);
		bool honest = check_proof(*pf, pk, *ck, board.eid, NV);
		P delta{0};
		delta.data()[0] = 1;
		delta.ntt_pow_phi();
		pf->ct.c1 = pf->ct.c1 + delta;
		check("a ballot changed after its proof is rejected",
				honest && !check_proof(*pf, pk, *ck, board.eid, NV));
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
