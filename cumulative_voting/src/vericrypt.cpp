#include <iostream>
#include <string>
#include <memory>
#include <random>
#include <chrono>
#include <vector>
#include <numeric>

using namespace std;

template <class T, size_t Align, class... Args>
T* alloc_aligned(size_t n, Args&& ... args)
{
	T* ret;
	if (posix_memalign((void**) &ret, Align, sizeof(T)*n) != 0) {
		throw std::bad_alloc();
	}
	for (size_t i = 0; i < n; i++) {
		new (&ret[i]) T(std::forward<Args>(args)...);
	}
	return ret;
}


//#include "params.hpp"
#include "vericrypt.hpp"
#include "util.hpp"


using V = typename P::value_type;

/* The statement the Fiat-Shamir challenges are bound to: the caller's context
 * (election and voter), the BGV public key, the commitment key and the
 * ciphertext. The first hash of the prover and of the verifier absorbs it
 * before anything else, so every later challenge depends on it: a proof does
 * not carry over to another voter, election, key or ciphertext, and the
 * statement cannot be chosen after the challenges. */
static void statement_digest(uint8_t out[SYMBYTES], const uint8_t *ctx, size_t ctxlen,
		const pk_t &pubkey, const comm_key_t &commkey, const ctx_t &ct){
	static const char label[] = "cumulative_voting vericrypt v1";
	uint8_t len[8];
	for(int i=0;i<8;i++) len[i] = (uint8_t)((uint64_t)ctxlen >> (8*i));
	shake128incctx st;
	shake128_inc_init(&st);
	shake128_inc_absorb(&st, (const uint8_t *)label, sizeof(label));
	shake128_inc_absorb(&st, len, sizeof(len));
	shake128_inc_absorb(&st, ctx, ctxlen);
	shake128_inc_absorb(&st, (const uint8_t *)pubkey.a.data(), sizeof(P));
	shake128_inc_absorb(&st, (const uint8_t *)pubkey.b.data(), sizeof(P));
	for(auto &row : commkey.B0l) for(auto &e : row) shake128_inc_absorb(&st, (const uint8_t *)e.data(), sizeof(P));
	for(auto &row : commkey.B0k) for(auto &e : row) shake128_inc_absorb(&st, (const uint8_t *)e.data(), sizeof(P));
	for(auto &row : commkey.b) for(auto &e : row) shake128_inc_absorb(&st, (const uint8_t *)e.data(), sizeof(P));
	shake128_inc_absorb(&st, (const uint8_t *)ct.c0.data(), sizeof(P));
	shake128_inc_absorb(&st, (const uint8_t *)ct.c1.data(), sizeof(P));
	shake128_inc_finalize(&st);
	shake128_inc_squeeze(out, SYMBYTES, &st);
	shake128_inc_ctx_release(&st);
}

void vericrypt(proof_t &p, const pk_t &pubkey,const comm_key_t& commkey, const P &message,
		const uint8_t *ctx, size_t ctxlen){

	uint8_t symbuf[4*SYMBYTES];
	memset(symbuf, 0, 4*SYMBYTES);
	uint8_t *thash = symbuf; // to store hash of public variables before while loop
	uint8_t *whash = symbuf+2*SYMBYTES; // to store hash of variables withing while loop i.e depends on w
	uint8_t *chash = symbuf+3*SYMBYTES;

	shake128incctx state;

	P poly_m = message;
	poly_m.ntt_pow_phi();

	P t = {plaintextModulus};
	t.ntt_pow_phi();

	P r_enc = nfl::non_uniform(2);
	P e_u = nfl::non_uniform(2);
	P e_v = nfl::non_uniform(2);

	r_enc.ntt_pow_phi();

	p.ct.c0 = e_u;
	p.ct.c0.ntt_pow_phi();
	p.ct.c0 = p.ct.c0 * t + r_enc * pubkey.a;

	p.ct.c1 = e_v;
	p.ct.c1.ntt_pow_phi();
	p.ct.c1 = p.ct.c1 * t + r_enc * pubkey.b + poly_m;

	p.ct.isnull = false;


	vector<P> shat;
	shat.resize(0);
	r_enc.invntt_pow_invphi();
	shat.push_back(r_enc);
	shat.push_back(e_u);
	shat.push_back(e_v);
	shat.push_back(message);

	P tmp = message;
	P tmp2{0}, tmp3{0};
	for(int i=1; i<= LC;i++){
		fill(tmp2.begin(), tmp2.end(), (V) i);
		tmp2 = message - tmp2;
		tmp = tmp * tmp2;
		shat.push_back(tmp);
	}

	BDLOP::comm_rnd_t<kappa> r_comm;
	p.t = BDLOP::comm_t<kappa>(commkey, r_comm);

	BDLOP::commit_polyvec(p.t, shat, 0, true);

	P g = nfl::uniform();
	fill(g.begin(), g.begin() + k, 0);
	g.ntt_pow_phi();

	BDLOP::commit_poly(p.t, g, shat.size(), true);

	uint8_t sdigest[SYMBYTES];
	statement_digest(sdigest, ctx, ctxlen, pubkey, commkey, p.ct);
	shake128_inc_init(&state);
	shake128_inc_absorb(&state, sdigest, SYMBYTES);
	shake128_inc_absorb(&state, (uint8_t *)p.t.t0.begin(), parameters::MU*sizeof(P));
	shake128_inc_absorb(&state, (uint8_t *)p.t.t.begin(),  (shat.size() + 1)*sizeof(P));
	shake128_inc_finalize(&state);
	shake128_inc_squeeze(thash,SYMBYTES, &state);

	array<BDLOP::comm_rnd_t<kappa>, k> y;

	// start rejection sampling
	P* alphas = alloc_aligned<P, 32>(4*k, P{0});
	P* alphas_prime = alloc_aligned<P, 32>(LC*k, P{0});
	P* gammas = alloc_aligned<P, 32>(3*k, P{0});

	bool rej = false;
	auto counter = 0;

	p.ct.c0.invntt_pow_invphi();
	p.ct.c1.invntt_pow_invphi();

	P t5,t6,t7;

	while(!rej){
		counter++;
//		cout << "trial " << counter++ << endl;
		for(auto &el : y){
			for(auto &yi : el.re){
				yi = nfl::non_uniform(delta_1);
				yi.ntt_pow_phi();
			}
			for(auto &yi : el.rs){
				yi = nfl::non_uniform(delta_1);
				yi.ntt_pow_phi();
			}
			for(auto &yi : el.rem){
				yi = nfl::non_uniform(delta_1);
				yi.ntt_pow_phi();
			}
		}


		array<array<P, parameters::MU>, k> w;
		for(int i=0; i<k;i++){
			for(size_t j =0; j<parameters::MU;j++){
				w[i][j] = inner_product(commkey.B0l[j].begin(), commkey.B0l[j].end(), y[i].rem.begin(), P{0});
				tmp = inner_product(commkey.B0k[j].begin(), commkey.B0k[j].end(), y[i].re.begin(), P{0});
				w[i][j] = w[i][j] + tmp + y[i].rs[j];
			}
		}

		shake128_inc_init(&state);
		shake128_inc_absorb(&state, thash,  SYMBYTES);
		shake128_inc_absorb(&state, (uint8_t *)w.begin(), k*parameters::MU*sizeof(P));
		shake128_inc_finalize(&state);
		shake128_inc_squeeze(whash,SYMBYTES, &state);

		P one{1};fill(one.begin(), one.end(), (P::value_type) 1);

		utils::sample_uniform(alphas, whash, 4*k);

		tmp = inner_product(commkey.b[LC+6].begin(), commkey.b[LC+6].end(), y[0].rem.begin(), P{0});
		tmp = tmp + y[0].re[LC+6];

		P v = inner_product(commkey.b[LC+5].begin(), commkey.b[LC+5].end(), y[0].rem.begin(), P{0});
		v = v + y[0].re[LC+5];

		P three, bjyi;
		fill(three.begin(), three.end(), (P::value_type) 3);
		tmp3 = {0};
		for(int i=0;i<k;i++){
			for(size_t j=0; j<3; j++){
				// <b_j, y_i>
				bjyi = inner_product(commkey.b[j].begin(), commkey.b[j].end(), y[i].rem.begin(), P{0});
				bjyi = bjyi + y[i].re[j];

				// m_nd+1 = <b_nd+2, y_0> - sum_ij alphas_ij * sigma^(-i)((3s_j)<b_j, y_i>^2))
				tmp2 = bjyi * bjyi;
				tmp2 = three * shat[j] * tmp2;
				for(int ii=0; ii<i; ii++) utils::inv_galois_transform(tmp2, tmp2, 2*P::degree/k + 1);
				tmp2 = alphas[3*i+j] * tmp2;
				tmp  = tmp - tmp2;

				// m_nd+2 = sum_ij alphas_ij * sigma^(-i)((3s_j^2 - 1)<b_j, y_i>))
				tmp2 = three * shat[j] * shat[j];
				tmp2 = tmp2 - one;
				tmp2 = tmp2 * bjyi;
				for(int ii=0; ii<i; ii++) utils::inv_galois_transform(tmp2, tmp2, 2*P::degree/k + 1);
				tmp2 = alphas[3*i+j] * tmp2;
				tmp3  = tmp3 +tmp2;

				// v = <b_nd+2, y_0> + sum_ij alphas_ij * sigma^(-i)(<b_j, y_i>^3))
				tmp2 = bjyi * bjyi * bjyi;
				for(int ii=0; ii<i; ii++) utils::inv_galois_transform(tmp2, tmp2, 2*P::degree/k + 1);
				tmp2 = alphas[3*i+j] * tmp2;
				v  = v + tmp2;
			}
		}

		t5 = p.t.t[LC+5] + tmp; //BDLOP::commit_poly(p.t, tmp, LC+5, true);
		t6 = p.t.t[LC+6] + tmp3; //BDLOP::commit_poly(p.t, tmp3, LC+6, true);

		P vprime;
		tmp = {0};
		for(size_t i=0;i<k;i++){
			tmp = inner_product(commkey.b[LC+3].begin(), commkey.b[LC+3].end(), y[i].rem.begin(), P{0});
			tmp = tmp + y[i].re[LC+3];

			for(size_t ii=0;ii<i;ii++) utils::inv_galois_transform(tmp, tmp, 2*P::degree/k + 1);
			tmp = alphas[3*k+i]*tmp;
			vprime = vprime + tmp;
		}



		shake128_inc_init(&state);
		shake128_inc_absorb(&state, whash,  SYMBYTES);
		shake128_inc_absorb(&state, (uint8_t *)t5.data(), sizeof(P));
		shake128_inc_absorb(&state, (uint8_t *)t6.data(), sizeof(P));
		shake128_inc_absorb(&state, (uint8_t *)v.data(), sizeof(P));
		shake128_inc_absorb(&state, (uint8_t *)vprime.data(), sizeof(P));
		shake128_inc_finalize(&state);
		shake128_inc_squeeze(whash, SYMBYTES, &state);

		utils::sample_uniform(alphas_prime, whash, LC*k);

		P v2 = inner_product(commkey.b[LC+7].begin(), commkey.b[LC+7].end(), y[0].rem.begin(), P{0});
		v2 = v2 + y[0].re[LC+7];

		tmp = {0};
		for(int i=0;i<k;i++){
			P b3yi = inner_product(commkey.b[3].begin(), commkey.b[3].end(), y[i].rem.begin(), P{0});
			b3yi = b3yi + y[i].re[3];
			for(size_t j=0;j<LC;j++){
				tmp2 = inner_product(commkey.b[3+j].begin(), commkey.b[3+j].end(), y[i].rem.begin(), P{0});
				tmp2 = tmp2 + y[i].re[3+j];

				tmp3 = b3yi * tmp2;
				for(int ii=0; ii<i; ii++) utils::inv_galois_transform(tmp3, tmp3, 2*P::degree/k + 1);
				tmp3 = tmp3 * alphas_prime[i*LC+j];
				v2 = v2 + tmp3;

				fill(tmp3.begin(), tmp3.end(), (P::value_type) (j+1));
				tmp3 = shat[3] - tmp3;
				tmp2 = tmp2 * tmp3;
				tmp2 = tmp2 + b3yi * shat[3+j];
				tmp3 = inner_product(commkey.b[3+j+1].begin(), commkey.b[3+j+1].end(), y[i].rem.begin(), P{0});
				tmp3 = tmp3 + y[i].re[3+j+1];
				tmp2 = tmp3 - tmp2;
				for(int ii=0; ii<i; ii++) utils::inv_galois_transform(tmp2, tmp2, 2*P::degree/k + 1);
				tmp2 = alphas_prime[i*LC+j] * tmp2;
				tmp  = tmp + tmp2;
			}
		}
		t7 = p.t.t[LC+7] + tmp; //BDLOP::commit_poly(p.t, tmp, LC+7, true);

		shake128_inc_init(&state);
		shake128_inc_absorb(&state, whash,  SYMBYTES);
		shake128_inc_absorb(&state, (uint8_t *) t7.data(), sizeof(P));
		shake128_inc_absorb(&state, (uint8_t *) v2.data(), sizeof(P));
		shake128_inc_finalize(&state);
		shake128_inc_squeeze(whash,2*SYMBYTES, &state);

		utils::sample_uniform(gammas,whash, 3*k);

		array<P, 4> pis;
		P pka = pubkey.a;
		pka.invntt_pow_invphi();
		utils::poly_transpose(pka);
		pka.ntt_pow_phi();

		P pkb = pubkey.b;
		pkb.invntt_pow_invphi();
		utils::poly_transpose(pkb);
		pkb.ntt_pow_phi();

		V kinv = (V) utils::modInverse((P::signed_value_type) k, (P::signed_value_type) P::get_modulus(0)); // k * kinv = 1 mod q

		array<P, k> vs;
		for(size_t i=0;i<k;i++){
			vs[i] = inner_product(commkey.b[LC+4].begin(), commkey.b[LC+4].end(), y[i].rem.begin(), P{0});
			vs[i] = vs[i] + y[i].re[LC+4];
		}

		P const_d, ls;
		fill(const_d.begin(), const_d.end(), (V) P::degree);
		fill(ls.begin(), ls.end(), (V)LS);
		mpz_t ls_mpz, gamma3_mpz;
		mpz_init(ls_mpz);
		mpz_init(gamma3_mpz);
		mpz_set_ui(ls_mpz, (V)LS);

		p.h = g;
		for(int _mu = 0; _mu<k;_mu++){
			P Xmu;
			Xmu.data()[_mu] = kinv;
			Xmu.ntt_pow_phi(); // find faster way of getting X^mu

			pis[0] = pka * gammas[3*_mu] + pkb * gammas[3*_mu+1];
			pis[1] = t * gammas[3*_mu];
			pis[2] = t * gammas[3*_mu+1];
			pis[3] = gammas[3*_mu+1];// + P{gammas[3*_mu + 2].data()[0]};


			for(auto &pi : pis){
				pi.invntt_pow_invphi();
			}

			mpz_t ug;
			mpz_init(ug);

			gammas[3*_mu].invntt_pow_invphi();
			gammas[3*_mu+1].invntt_pow_invphi();
			gammas[3*_mu+2].invntt_pow_invphi();
			for(int i=0; i < P::degree; i++){
				pis[3].data()[i] += gammas[3*_mu + 2].data()[0];
			}
			mpz_set_ui(gamma3_mpz, gammas[3*_mu + 2].data()[0]);

			utils::poly_inner_product_Zq(ug, p.ct.c0, gammas[3*_mu]);
			utils::poly_inner_product_Zq(ug, p.ct.c1, gammas[3*_mu+1]);
			mpz_addmul(ug, gamma3_mpz, ls_mpz);
			mpz_mod_ui(ug, ug, P::get_modulus(0));


			tmp = inner_product(pis.begin(), pis.end(), shat.begin(), P{0});
			tmp = const_d * tmp;
			tmp.invntt_pow_invphi();

			tmp = tmp - P{(V)mpz_get_ui(ug)};
			tmp2 = tmp;
			tmp3 = {0};
			for(int i=0;i<k;i++){
				tmp = tmp2;
				for(int j=0;j<i;j++) utils::galois_transform(tmp, tmp, 2*P::degree/k + 1);
				tmp3 = tmp3 + tmp;
			}

			tmp3.ntt_pow_phi();
			tmp3 = tmp3 * Xmu;
			p.h = p.h + tmp3;


			for(auto i=0;i<k;i++){
				tmp3 = {0};
				for(auto _nu=0;_nu<k;_nu++){
					tmp2 = {0};
					for(auto j=0; j<4;j++){
						tmp = inner_product( commkey.b[j].begin(), commkey.b[j].end(), y[nmod(i-_nu, k)].rem.begin(), P{0});
						tmp = tmp + y[nmod(i-_nu, k)].re[j];
						tmp = tmp * pis[j];
						tmp = tmp * const_d;
						tmp2 = tmp2 + tmp;
					}

					tmp2.invntt_pow_invphi();
					tmp = tmp2;
					for(int j=0;j<_nu;j++) utils::galois_transform(tmp, tmp, 2*P::degree/k + 1);
					tmp3 = tmp3 + tmp;
				}
				tmp3.ntt_pow_phi();
				tmp3 = tmp3 * Xmu;
				vs[i] = vs[i] + tmp3;
			}
		}
		p.h.invntt_pow_invphi();

		shake128_inc_init(&state);
		shake128_inc_absorb(&state, whash,  2*SYMBYTES);
		shake128_inc_absorb(&state, (uint8_t *) p.h.data(), sizeof(P));
		shake128_inc_absorb(&state, (uint8_t *) vs.data(), k*sizeof(P));
		shake128_inc_finalize(&state);
		shake128_inc_squeeze(chash, SYMBYTES, &state);

		P challenge = utils::sample_ternary(0x80, chash);

		for(int i=0;i<k; i++){
			tmp = challenge;
			for(int j=0;j<i;j++) utils::galois_transform(tmp, tmp, 2*P::degree/k + 1);
			tmp.ntt_pow_phi();
			for(size_t j=0;j<parameters::LAMBDA;j++) p.z[i].rem[j] = y[i].rem[j] + tmp * r_comm.rem[j];
			for(size_t j=0;j<parameters::MU;j++) p.z[i].rs[j] = y[i].rs[j] + tmp * r_comm.rs[j];
			for(size_t j=0;j<kappa;j++) p.z[i].re[j] = y[i].re[j] + tmp * r_comm.re[j];
		}

		memcpy(p.chash, chash, SYMBYTES);

		for(int i=0; i<k; i++){
			for(size_t j=0;j<parameters::LAMBDA; j++) p.z[i].rem[j].invntt_pow_invphi();
			for(size_t j=0;j<parameters::MU; j++) p.z[i].rs[j].invntt_pow_invphi();
			for(size_t j=0;j<kappa; j++) p.z[i].re[j].invntt_pow_invphi();
		}

		rej = true;
		for(int i=0; i<k; i++){
			auto nrm = utils::inf_norm(p.z[i].rem);
			if( nrm >=  delta_1 - beta_1 ){
				rej &= false;
			}
			nrm = utils::inf_norm(p.z[i].rs);
			if( nrm >= delta_1 - beta_1){
				rej &= false;
			}
			nrm = utils::inf_norm(p.z[i].re);
			if( nrm >= delta_1 - beta_1){
				rej &= false;
			}
		}
	}

	p.t.t[LC+5] = t5;
	p.t.t[LC+6] = t6;
	p.t.t[LC+7] = t7;

	p.ct.c0.ntt_pow_phi();
	p.ct.c1.ntt_pow_phi();

//	for(uint8_t *c = chash; c < chash+SYMBYTES; c++){
//		printf("%02x", *c);
//	}printf("\n");

	free(alphas);
	free(alphas_prime);
	free(gammas);

}

bool verify(proof_t &p, const comm_key_t &commkey, const pk_t &pubkey,
		const uint8_t *ctx, size_t ctxlen){


	bool  flag = true;
	for(int i=0; i<k; i++){
		auto nrm = utils::inf_norm(p.z[i].rem);

		if( nrm >=  delta_1 - beta_1 ){
			flag &= false;
		}
		nrm = utils::inf_norm(p.z[i].rs);

		if( nrm >= delta_1 - beta_1){
			flag &= false;
		}
		nrm = utils::inf_norm(p.z[i].re);

		if( nrm >= delta_1 - beta_1){
			flag &= false;
		}
	}

	for(int i=0; i<k; i++){
		for(size_t j=0;j<parameters::LAMBDA; j++) p.z[i].rem[j].ntt_pow_phi();
		for(size_t j=0;j<parameters::MU; j++) p.z[i].rs[j].ntt_pow_phi();
		for(size_t j=0;j<kappa; j++) p.z[i].re[j].ntt_pow_phi();
	}

	if(!flag) return false;
//	cout << "check  1 ✅\n";

	uint8_t symbuf[4*SYMBYTES];
	memset(symbuf, 0, 4*SYMBYTES);
	uint8_t *thash = symbuf; // to store hash of public variables before while loop
	uint8_t *whash = symbuf+2*SYMBYTES; // to store hash of variables withing while loop i.e depends on w
	uint8_t *chash = symbuf+3*SYMBYTES;

	shake128incctx state;

	uint8_t sdigest[SYMBYTES];
	statement_digest(sdigest, ctx, ctxlen, pubkey, commkey, p.ct);
	shake128_inc_init(&state);
	shake128_inc_absorb(&state, sdigest, SYMBYTES);
	shake128_inc_absorb(&state, (uint8_t *)p.t.t0.data(), parameters::MU*sizeof(P));
	shake128_inc_absorb(&state, (uint8_t *)p.t.t.data(),  (4+LC+1)*sizeof(P));
	shake128_inc_finalize(&state);
	shake128_inc_squeeze(thash,SYMBYTES, &state);

	array<P, k> sigma_c;
	P tmp;

	P c = utils::sample_ternary(0x80, p.chash);

//	p.c.invntt_pow_invphi();
	for(int i=0;i<k;i++){
		tmp = c;
		for(int j=0;j<i;j++) utils::galois_transform(tmp, tmp, 2*P::degree/k + 1);
		tmp.ntt_pow_phi();
		sigma_c[i] = tmp;
	}
	c.ntt_pow_phi();

	array<array<P, parameters::MU>, k> w;
	for(int i=0; i<k;i++){
		for(size_t j=0; j<parameters::MU;j++){
			w[i][j] = inner_product(commkey.B0l[j].begin(), commkey.B0l[j].end(), p.z[i].rem.begin(), P{0});
			tmp = inner_product(commkey.B0k[j].begin(), commkey.B0k[j].end(), p.z[i].re.begin(), P{0});
			w[i][j] = w[i][j] + tmp + p.z[i].rs[j];
			w[i][j] = w[i][j] - sigma_c[i]* p.t.t0[j];
		}
	}

	shake128_inc_init(&state);
	shake128_inc_absorb(&state, thash,  SYMBYTES);
	shake128_inc_absorb(&state, (uint8_t *)w.begin(), k*parameters::MU*sizeof(P));
	shake128_inc_finalize(&state);
	shake128_inc_squeeze(whash, SYMBYTES, &state);


	P* alphas = alloc_aligned<P, 32>(4*k, P{0});
	P* alphas_prime = alloc_aligned<P, 32>(LC*k, P{0});
	P* gammas = alloc_aligned<P, 32>(3*k, P{0});

	utils::sample_uniform(alphas, whash, 4*k);

	P f2, f3, f4, fji;
	f2 = inner_product(commkey.b[LC+5].begin(), commkey.b[LC+5].end(), p.z[0].rem.begin(), P{0});
	f2 = f2 + p.z[0].re[LC+5];
	f2 = f2 - c * p.t.t[LC+5];

	f3 = inner_product(commkey.b[LC+6].begin(), commkey.b[LC+6].end(), p.z[0].rem.begin(), P{0});
	f3 = f3 + p.z[0].re[LC+6];
	f3 = f3 - c * p.t.t[LC+6];

	P tmp2, tmp3;
	tmp = P{0};

	for(int i=0; i<k;i++){
		for(size_t j=0;j<3;j++){
			fji = inner_product(commkey.b[j].begin(), commkey.b[j].end(), p.z[i].rem.begin(), P{0});
			fji = fji + p.z[i].re[j];
			fji = fji - sigma_c[i] * p.t.t[j];
			tmp2 = fji;
			tmp3 = fji - sigma_c[i];
			tmp2 = tmp2 * tmp3;
			tmp3 = fji + sigma_c[i];
			tmp2 = tmp2 * tmp3;
			for(int ii=0; ii<i; ii++) utils::inv_galois_transform(tmp2, tmp2, 2*P::degree/k + 1);
			tmp2 = tmp2 * alphas[3*i+j];
			tmp  = tmp + tmp2;
		}
	}

	P v = tmp + f2;
	v = v + c*f3;

	P vprime;
	for(size_t i=0; i<k;i++){
		tmp = inner_product(commkey.b[LC+3].begin(), commkey.b[LC+3].end(), p.z[i].rem.begin(), P{0});
		tmp = tmp + p.z[i].re[LC+3];
		tmp = tmp - sigma_c[i] * p.t.t[LC+3];
		tmp2 = tmp;
		for(size_t j=0;j<i;j++) utils::inv_galois_transform(tmp2, tmp2, 2*P::degree/k + 1);
		tmp2 = alphas[3*k+i]*tmp2;
		vprime = vprime + tmp2;
	}

	shake128_inc_init(&state);
	shake128_inc_absorb(&state, whash,  SYMBYTES);
	shake128_inc_absorb(&state, (uint8_t *)p.t.t[LC+5].data(), sizeof(P));
	shake128_inc_absorb(&state, (uint8_t *)p.t.t[LC+6].data(), sizeof(P));
	shake128_inc_absorb(&state, (uint8_t *)v.data(), sizeof(P));
	shake128_inc_absorb(&state, (uint8_t *)vprime.data(), sizeof(P));
	shake128_inc_finalize(&state);
	shake128_inc_squeeze(whash,SYMBYTES, &state);

	utils::sample_uniform(alphas_prime, whash, LC*k);


	f4 = inner_product(commkey.b[LC+7].begin(), commkey.b[LC+7].end(), p.z[0].rem.begin(), P{0});
	f4 = f4 + p.z[0].re[LC+7];
	f4 = f4 - c * p.t.t[LC+7];

	tmp = {0};
	P f3i;
	for(int i=0; i<k;i++){
		f3i = inner_product(commkey.b[3].begin(), commkey.b[3].end(), p.z[i].rem.begin(), P{0});
		f3i = f3i + p.z[i].re[3];
		f3i = f3i - sigma_c[i] * p.t.t[3];
		for(size_t j=0;j<LC;j++){
			fji = inner_product(commkey.b[3+j].begin(), commkey.b[3+j].end(), p.z[i].rem.begin(), P{0});
			fji = fji + p.z[i].re[3+j];
			fji = fji - sigma_c[i] * p.t.t[3+j];
			fill(tmp2.begin(), tmp2.end(), (P::value_type) (j+1));
			tmp2 = sigma_c[i] * tmp2;
			tmp2 = f3i + tmp2;
			tmp2 = fji * tmp2;
			fji = inner_product(commkey.b[3+j+1].begin(), commkey.b[3+j+1].end(), p.z[i].rem.begin(), P{0});
			fji = fji + p.z[i].re[3+j+1];
			fji = fji - sigma_c[i] * p.t.t[3+j+1];
			tmp2 = tmp2 + sigma_c[i] * fji;
			for(int ii=0; ii<i; ii++) utils::inv_galois_transform(tmp2, tmp2, 2*P::degree/k + 1);
			tmp2 = tmp2 * alphas_prime[i*LC + j];
			tmp  = tmp + tmp2;
		}
	}
	P v2 = tmp + f4;

	shake128_inc_init(&state);
	shake128_inc_absorb(&state, whash,  SYMBYTES);
	shake128_inc_absorb(&state, (uint8_t *)p.t.t[LC+7].data(), sizeof(P));
	shake128_inc_absorb(&state, (uint8_t *)v2.data(), sizeof(P));
	shake128_inc_finalize(&state);
	shake128_inc_squeeze(whash,2*SYMBYTES, &state);

	utils::sample_uniform(gammas, whash, 3*k);

	array<P, k> vs;
	for(size_t i=0; i<k;i++){
		vs[i] = inner_product(commkey.b[LC+4].begin(), commkey.b[LC+4].end(), p.z[i].rem.begin(), P{0});
		vs[i] = vs[i] + p.z[i].re[LC+4];
	}


	for(int i=0;i<k;i++){
		if(p.h.data()[i] != 0) return false;
	}
//	cout << "check  2 ✅\n";

	array<P, 4> pis;
	P pka = pubkey.a;
	pka.invntt_pow_invphi();
	utils::poly_transpose(pka);
	pka.ntt_pow_phi();

	P pkb = pubkey.b;
	pkb.invntt_pow_invphi();
	utils::poly_transpose(pkb);
	pkb.ntt_pow_phi();

	//P Xmu, X{0,1};
	//X.ntt_pow_phi();
	//Xmu = one;
	V kinv = (V) utils::modInverse((P::signed_value_type) k, (P::signed_value_type) P::get_modulus(0)); // k * kinv = 1 mod q
	p.ct.c0.invntt_pow_invphi();
	p.ct.c1.invntt_pow_invphi();

	P const_d, ls;
	fill(const_d.begin(), const_d.end(), (V) P::degree);
	fill(ls.begin(), ls.end(), (V) LS);

	mpz_t ls_mpz, gamma3_mpz;
	mpz_init(ls_mpz);
	mpz_init(gamma3_mpz);
	mpz_set_ui(ls_mpz, (V)LS);

	P t = {plaintextModulus};
	t.ntt_pow_phi();
	P tau;
	for(int _mu = 0; _mu<k;_mu++){
		P Xmu;
		Xmu.data()[_mu] = kinv;
		Xmu.ntt_pow_phi(); // find faster way of getting X^mu

		pis[0] = pka * gammas[3*_mu];
		pis[0] = pis[0] + pkb * gammas[3*_mu+1];
		pis[1] = t * gammas[3*_mu];
		pis[2] = t * gammas[3*_mu+1];
		pis[3] = gammas[3*_mu+1]; // + P{gammas[3*_mu + 2].data()[0]};

		for(auto &pi : pis){
			pi.invntt_pow_invphi();
		}

		mpz_t ug;
		mpz_init(ug);

		gammas[3*_mu].invntt_pow_invphi();
		gammas[3*_mu+1].invntt_pow_invphi();
		gammas[3*_mu+2].invntt_pow_invphi();
		for(int i=0; i < P::degree; i++){
			pis[3].data()[i] += gammas[3*_mu + 2].data()[0];
		}
		mpz_set_ui(gamma3_mpz, gammas[3*_mu + 2].data()[0]);

		utils::poly_inner_product_Zq(ug, p.ct.c0, gammas[3*_mu]);
		utils::poly_inner_product_Zq(ug, p.ct.c1, gammas[3*_mu+1]);
		mpz_addmul(ug, gamma3_mpz, ls_mpz);
		mpz_mod_ui(ug, ug, (unsigned long int) P::get_modulus(0));

		tmp = inner_product(pis.begin(), pis.end(), p.t.t.begin(), P{0});
		tmp = const_d * tmp;
//		tmp2 = const_d*p.t.t[3] - ls;
//		tmp2 = gammas_prime[_mu]* tmp2;
//		tmp = tmp + tmp2;
		tmp.invntt_pow_invphi();

		tmp = tmp - P{(V)mpz_get_ui(ug)};
		tmp2 = tmp;
		tmp3 = {0};
		for(int i=0;i<k;i++){
			tmp = tmp2;
			for(int j=0;j<i;j++) utils::galois_transform(tmp, tmp, 2*P::degree/k + 1);
			tmp3 = tmp3 + tmp;
		}
		tmp3.ntt_pow_phi();
		tmp3 = tmp3 * Xmu;
		tau = tau + tmp3;

		for(auto i=0;i<k;i++){
			tmp3 = {0};

			for(auto _nu=0;_nu<k;_nu++){
				tmp2 = {0};
				for(auto j=0; j<4;j++){
					tmp = inner_product( commkey.b[j].begin(), commkey.b[j].end(), p.z[nmod(i-_nu, k)].rem.begin(), P{0});
					tmp = tmp + p.z[nmod(i-_nu, k)].re[j];
					tmp = tmp * pis[j];
					tmp = tmp * const_d;
					tmp2 = tmp2 + tmp;
				}

//				tmp = inner_product( commkey.b[3].begin(), commkey.b[3].end(), p.z[nmod(i-_nu, k)].rem.begin(), P{0});
//				tmp = tmp + p.z[nmod(i-_nu, k)].re[3];
//				tmp = const_d*tmp;
//				tmp = tmp *gammas_prime[_mu];
//				tmp2 = tmp2 + tmp;

				tmp2.invntt_pow_invphi();
				tmp = tmp2;
				for(int j=0;j<_nu;j++) utils::galois_transform(tmp, tmp, 2*P::degree/k + 1);
				tmp3 = tmp3 + tmp;
			}
			tmp3.ntt_pow_phi();
			tmp3 = tmp3 * Xmu;
			vs[i] = vs[i] + tmp3;
		}
	}

	p.h.ntt_pow_phi();
	tmp2 = tau + p.t.t[LC+4];
	tmp2 = tmp2 - p.h;
	for(auto i=0; i<k;i++){
		vs[i] = vs[i] - sigma_c[i] * tmp2;
	}
	p.h.invntt_pow_invphi();

	p.ct.c0.ntt_pow_phi();
	p.ct.c1.ntt_pow_phi();

	shake128_inc_init(&state);
	shake128_inc_absorb(&state, whash,  2*SYMBYTES);
	shake128_inc_absorb(&state, (uint8_t *) p.h.begin(), sizeof(P));
	shake128_inc_absorb(&state, (uint8_t *) vs.begin(), k*sizeof(P));
	shake128_inc_finalize(&state);
	shake128_inc_squeeze(chash, SYMBYTES, &state);

//	for(uint8_t *c = chash; c < chash+SYMBYTES; c++){
//		printf("%02x", *c);
//	}printf("\n");


	P challenge = utils::sample_ternary(0x80, chash);
	challenge.ntt_pow_phi();
	if(challenge != c) return false;
//	cout << "check  3 ✅\n";

	free(alphas);
	free(alphas_prime);
	free(gammas);

	return true;
}


