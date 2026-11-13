
#ifndef HEADER_FS
#define HEADER_FS

#include <NTL/ZZ_pE.h>
#include <NTL/vec_ZZ_pE.h>
#include <NTL/mat_ZZ_pE.h>
#include "parameters.h"

extern "C" {
#include "fips202.h"
}

using namespace NTL;

/* The challenges of both proofs, by Fiat-Shamir: SHAKE256 of a label, the
 * statement and the prover's first messages. A ring element is absorbed as
 * its N coefficients in [0, q), 8 bytes each, little-endian. */
class Transcript {
  shake256incctx state;
  bool squeezing;
public:
  explicit Transcript(const char *label);
  ~Transcript();
  Transcript(const Transcript &) = delete;
  Transcript &operator=(const Transcript &) = delete;
  void Absorb(const ZZ_pE &);
  void Absorb(const Vec<ZZ_pE> &);
  void Absorb(const Mat<ZZ_pE> &);
  void Squeeze(unsigned char *, size_t);
};

/* c with exactly 36 coefficients +-1 at distinct positions, the rest 0. */
void LinearChallenge(ZZ_pE &c, Transcript &);

/* C with entries in {0, 1}, of the dimensions it already has. */
void ShortnessChallenge(Mat<ZZ_pE> &C, Transcript &);

#endif
