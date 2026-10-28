#include <array>
#include <iostream>
#include <string>
#include <memory>
#include <chrono>
#include <vector>

using namespace std;

#include "vericrypt.hpp"

using P = parameters::poly_t;
using V = typename P::value_type;

array<V, NC> ballot_generator(){
	bool flag = false;
	int lim;
	array<V, NC> tmp;

	while(!flag){
		fill(tmp.begin(), tmp.end(), 0);
		lim = LS;
		for(int i=0; i<NC;i++){
			if(lim >0){
				auto b = rand() % (LC+1);
				tmp[i] = (V) b;
				lim -= b;
			}
		}
		auto acc = accumulate(tmp.begin(), tmp.end(), 0);
		flag = (acc == (V) LS);
	}
	return tmp;
}

void random_ballot(P& out, array<V, NC> ballot){
	shuffle(ballot.begin(), ballot.end(), default_random_engine(chrono::system_clock::now().time_since_epoch().count()));
	memcpy(out.begin(), ballot.begin(), NC*sizeof(V));
}

int main(){
	cout << "Number of voters:        " << NV << endl;
	cout << "Number of candiates:     " << NC << endl;
	cout << "Max vote per candidate:  " << LC << endl;
	cout << "Total vote budget(l_s):  " << LS << endl;
	cout << "sizeof proof:            " << sizeof(proof_t) << "B\n";

	srand(time(NULL));

	sk_t secret_key;
	pk_t public_key(secret_key);

	comm_key_t commkey = comm_key_t();

	auto total_duration_voter=0, total_duration_verifier=0, total_duration_tally=0;
	int rejected_ballots = 0;

	array<V, NC> valid_ballot = ballot_generator();
	BGV::ciphertext_t tally;
	P real_tally = {0};
	int i = 0;
	for(i=0;i<NV;i++){
		P m = {0};
		random_ballot(m, valid_ballot);
//		cout << "ballot: ";
//		for(int j=0; j< NC; j++){
//			cout << m.data()[j] << ", ";
//		}
//		cout << endl;

		proof_t proof;
		auto start = chrono::steady_clock::now();
		string ctx = "vericrypt demo, voter " + to_string(i);
		vericrypt(proof, public_key, commkey, m, (const uint8_t *) ctx.data(), ctx.size());
		auto end = chrono::steady_clock::now();
		auto duration_voter = chrono::duration_cast<chrono::milliseconds>(end-start).count();
		// cout << "verifiable encryption time: " << duration_voter << "(ms)" << endl;

		start = chrono::steady_clock::now();
		bool res = verify(proof, commkey, public_key, (const uint8_t *) ctx.data(), ctx.size());
		end = chrono::steady_clock::now();
		auto duration_verifier = chrono::duration_cast<chrono::milliseconds>(end-start).count();

		if(!res) {
			cout << "rejected ballot #" << (rejected_ballots++) << endl;
			continue;
		}

		total_duration_voter += duration_voter;
		total_duration_verifier += duration_verifier;
		start = chrono::steady_clock::now();
		tally += proof.ct;
		end = chrono::steady_clock::now();
		auto duration_tally = chrono::duration_cast<chrono::microseconds>(end-start).count();
		total_duration_tally += duration_tally;
		real_tally = real_tally + m;
		// cout << "single encryption noise: " << BGV::noise(secret_key, public_key, proof.ct) << endl;

	}

	cout << "total encryption time:   " << total_duration_voter << "(ms)" << endl;
	cout << "total verifier time:     " << total_duration_verifier << "(ms)" << endl;
	cout << "total tallying time:     " << total_duration_tally << "(μs)" << endl;
	cout << "rejected ballots:        " << rejected_ballots << endl;

	P result;
	auto start = chrono::steady_clock::now();
	BGV::decrypt_poly(result, secret_key, public_key, tally);
	auto end = chrono::steady_clock::now();
	auto duration = chrono::duration_cast<chrono::microseconds>(end-start).count();
	cout << "Decryption time:         " << duration << "(μs)" << endl;
	cout << "Decrypted tally is " << (result == real_tally ? "" : "not") << "equal to actual tally" << endl;
	//cout << result << endl;

	auto noise = BGV::noise(secret_key, public_key, tally);
	cout << "log2 noise:              " << noise << endl;
	cout << "noise max:               " << public_key.noise_max << endl;

	return 0;
}
