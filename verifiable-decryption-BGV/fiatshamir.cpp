#include <cstring>
#include <stdexcept>
#include <vector>
#include "fiatshamir.h"

Transcript::Transcript(const char *label) : squeezing(false) {
  shake256_inc_init(&state);
  shake256_inc_absorb(&state, (const uint8_t *) label, strlen(label) + 1);}

Transcript::~Transcript(){
  shake256_inc_ctx_release(&state);}

void Transcript::Absorb(const ZZ_pE &x){
  if (squeezing) throw std::logic_error("absorb after squeeze");
  const ZZ_pX &f = rep(x);
  unsigned char buf[8*N];
  for (long i = 0; i < N; i++){
    unsigned long c = to_ulong(rep(coeff(f,i)));
    for (int b = 0; b < 8; b++) buf[8*i+b] = (unsigned char) (c >> (8*b));}
  shake256_inc_absorb(&state, buf, sizeof(buf));}

void Transcript::Absorb(const Vec<ZZ_pE> &v){
  for (long i = 0; i < v.length(); i++) Absorb(v[i]);}

void Transcript::Absorb(const Mat<ZZ_pE> &m){
  for (long i = 0; i < m.NumRows(); i++) Absorb(m[i]);}

void Transcript::Squeeze(unsigned char *out, size_t len){
  if (!squeezing){shake256_inc_finalize(&state); squeezing = true;}
  shake256_inc_squeeze(out, len, &state);}

void LinearChallenge(ZZ_pE &c, Transcript &t){
  static_assert((N & (N-1)) == 0, "N a power of two");
  ZZ_pX f; std::vector<bool> used(N, false);
  int set = 0; unsigned char buf[3];
  while (set < 36){
    t.Squeeze(buf, sizeof(buf));
    long pos = (buf[0] | ((long) buf[1] << 8)) & (N-1);
    if (used[pos]) continue;
    used[pos] = true; set++;
    SetCoeff(f, pos, (buf[2] & 1) ? 1 : -1);}
  c = to_ZZ_pE(f);}

void ShortnessChallenge(Mat<ZZ_pE> &C, Transcript &t){
  unsigned char byte = 0; int left = 0;
  for (long i = 0; i < C.NumRows(); i++){
    for (long j = 0; j < C.NumCols(); j++){
      if (left == 0){t.Squeeze(&byte, 1); left = 8;}
      C[i][j] = ZZ_pE(byte & 1); byte >>= 1; left--;}}}
