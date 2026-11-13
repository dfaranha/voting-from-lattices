
#include "linearityproof.h"
#include "fiatshamir.h"

/* The challenge of one linear relation, from its statement (the commitment
 * key, the two commitments, alpha and beta) and the prover's first messages
 * t1, t2 and u. */
static void Challenge(ZZ_pE &c, const ComKey &comKey, const Vec<ZZ_pE> &X1,
  const Vec<ZZ_pE> &X2, const ZZ_pE &alpha, const ZZ_pE &beta,
  const ZZ_pE &t1, const ZZ_pE &t2, const ZZ_pE &u){
  Transcript t("verifiable-decryption-BGV linear relation");
  t.Absorb(comKey.B1); t.Absorb(comKey.B2); t.Absorb(X1); t.Absorb(X2);
  t.Absorb(alpha); t.Absorb(beta); t.Absorb(t1); t.Absorb(t2); t.Absorb(u);
  LinearChallenge(c, t);}

void ProveLinearRelation(LinProof &proof, const ComKey &comKey, const
  Commitment &X1, const Commitment &X2, const ZZ_pE &alpha, const ZZ_pE &beta){
    Vec<ZZ_pE> y1, y2, cr1, cr2; y1.SetLength(3),
    y2.SetLength(3), cr1.SetLength(3), cr2.SetLength(3);
    bool success = false, success1 = false, success2 = false;
    while (!success){
      SampleGaussian(y1,Sigma_C), SampleGaussian(y2,Sigma_C);
      proof.t1 = y1[0] + (y1[1]*comKey.B1[1]) + (y1[2]*comKey.B1[2]);
      proof.t2 = y2[0] + (y2[1]*comKey.B1[1]) + (y2[2]*comKey.B1[2]);
      proof.u = alpha * (y1[1] + y1[2]*comKey.B2[2]) -
        (y2[1] + y2[2]*comKey.B2[2]);
      Challenge(proof.c, comKey, X1.com, X2.com, alpha, beta, proof.t1,
        proof.t2, proof.u);
    for(int j = 0; j < 3; j++){
      cr1[j] = proof.c * X1.r[j]; cr2[j] = proof.c * X2.r[j];
      proof.z1[j] = y1[j] + cr1[j]; proof.z2[j] = y2[j] + cr2[j];}
    RejectionSamplingLinearity(success1, proof.z1, cr1);
    RejectionSamplingLinearity(success2, proof.z2, cr2);
    success = success1 and success2;}
}

void ProveLinearRelations(Vec<LinProof> &linProof, const Vec<Commitment>
  &noiseCom, const ComKey comKey, const Commitment keyCom, const
  Vec<Ciphertext> ciphertexts, const Vec<ZZ_pE> messages){
  for(int j = 0; j < tau; j++){
    ProveLinearRelation(linProof[j], comKey, keyCom, noiseCom[j],
    -ciphertexts[j].CTXu/ZZ_p(p),(ciphertexts[j].CTXv-messages[j])/ZZ_p(p));
}}

void VerifyLinearRelation(bool &b, const LinProof &proof, const ComKey &comKey,
  const Vec<ZZ_pE> &X1, const Vec<ZZ_pE> &X2, const ZZ_pE &alpha, const ZZ_pE
  &beta){b = true;
  ZZ_pE c;
  Challenge(c, comKey, X1, X2, alpha, beta, proof.t1, proof.t2, proof.u);
  if (c != proof.c){b = false;}
  if (!CheckNorm(proof.z1,Bound_C) || !CheckNorm(proof.z2,Bound_C)){b = false;}
  if ((proof.z1[0] + proof.z1[1]*comKey.B1[1] + proof.z1[2]*comKey.B1[2])
  != (proof.t1 + proof.c*X1[0])){b = false;}
  if ((proof.z2[0] + proof.z2[1]*comKey.B1[1] + proof.z2[2]*comKey.B1[2])
  != (proof.t2 + proof.c*X2[0])){b = false;}
  if ((alpha*(proof.z1[1] + proof.z1[2]*comKey.B2[2]) - (proof.z2[1] +
    proof.z2[2]*comKey.B2[2])) != ((alpha*X1[1]+beta-X2[1]))*proof.c + proof.u)
    {b = false;}
}

void VerifyLinearRelations(bool &b, const Vec<LinProof> &linProof, const ComKey
  &comKey, const Commitment &keyCom, const Vec<Commitment> &noiseCom,
  const Vec<Ciphertext> &ciphertexts, const Vec<ZZ_pE> &messages){
  bool bTemp; b = true;
  for(int j = 0; j < tau; j++){
    VerifyLinearRelation(bTemp,linProof[j],comKey,keyCom.com, noiseCom[j].com,
    - ciphertexts[j].CTXu/ZZ_p(p), (ciphertexts[j].CTXv -messages[j])/ZZ_p(p));
    b = b && bTemp;
}}
