/*
 * Check 6 of stage A (design section 12) for the C code: the Fock matrices of AtomHF.cpp
 * against finite differences of the energy E_av (HFEnergy). For unnormalized coefficients
 * E_av is a polynomial of degree 4, and dE_av/dc_a = 2 qbar_a F_a c_a exactly; the central
 * difference along a random direction, Richardson-extrapolated from h = 1e-2 and 5e-3, is
 * exact for a quartic, so what is left is rounding (|E| eps / h).
 * Build from the repository root (the harness includes Quanty's single compilation unit):
 *   g++ -I. Doc/Claude/QuantyDFT_StageA/QuantyDFT_StageA_proto/check6_fock_gradient.cpp \
 *       -D LUA_USE_MACOSX -DBLASunderscore -O2 -o check6 -framework Accelerate -lreadline
 * Result (2026-11-11): relative agreement 2.5e-13 .. 2.8e-10 for Ni 3d8 4s2 (Dirac, ten
 * subshells), O 2p4 and Fe 3d6 4s2 (Dirac and Schroedinger); absolute <= 3.3e-10. The
 * last run (Fe, Schroedinger, jjAngular) prints E_av with the j-subshells (check 5).
 */
#define main QuantyMain
#include "LuaInterface/main.cpp"
#undef main

static double Eof(DiracHFType *hf)
{
  double E, a, b, c, d, e;
  HFOrbitalsOnPoints(hf);
  HFEnergy(hf, &E, &a, &b, &c, &d);
  return E;
}

static void run(double Z, const char *conf, int rel, int jj)
{
  RadialBasisType basis;
  DiracAtomType atom;
  DiracHFType hf;
  uint32_t nShell, q, a, i, Nf, gr;
  AtomShellType shell[40];
  double cc = 137.035999084;
  nShell = 0;
  if (AtomParseConfiguration(conf, 40, &nShell, shell)) { printf("parse failed\n"); return; }
  RadialBasisCreate(8, 120, 1e-6/Z, 40., 13, &basis);
  double *V = (double *) malloc(basis.nPoint*sizeof(double));
  for (q = 0; q < basis.nPoint; q++) V[q] = -Z/basis.rq[q];
  if (AtomSetup(&atom, &basis, V, Z, rel, cc, XCLDAVWN, 1, 1e-30, nShell, shell)) { printf("setup failed\n"); return; }
  AtomSCF(&atom, 200, 8, 0.5, 1e-10, 1e-8, 1, 0);
  if (HFSetup(&hf, &atom, jj)) { printf("HFSetup failed\n"); return; }
  double E, Ek, En, Ed, Ex;
  HFOrbitalsOnPoints(&hf);
  HFEnergy(&hf, &E, &Ek, &En, &Ed, &Ex);
  printf("%s Z=%g %s %s: nSub %u kMax %u  E_HF(LDA orbitals) = %.10f (kin %.6f nuc %.6f dir %.6f exch %.6f), LDA Etot %.10f\n",
         conf, Z, rel ? "Dirac" : "Schroedinger", jj ? "jj" : "", hf.nSub, hf.kMax, E, Ek, En, Ed, Ex, atom.Etot);
  if (jj) { HFFree(&hf); AtomFree(&atom); RadialBasisFree(&basis); free(V); return; }
  Nf = atom.Nf;
  double **F = (double **) malloc(hf.nSub*sizeof(double *));
  for (a = 0; a < hf.nSub; a++) F[a] = (double *) malloc(Nf*Nf*sizeof(double));
  HFFockMatrices(&hf, F);
  // directional derivatives, a random direction per subshell on its active functions
  unsigned long s = 12345;
  for (a = 0; a < hf.nSub; a++)
  {
    double *c = hf.sub[a].coef, *dc = (double *) calloc(Nf, sizeof(double)), *c0 = (double *) malloc(Nf*sizeof(double));
    gr = atom.orb[hf.sub[a].orbital].group;
    for (i = 0; i < atom.group[gr].N; i++) { s = s*6364136223846793005UL + 1442695040888963407UL; dc[atom.group[gr].active[i]] = ((s >> 11)*(1.0/9007199254740992.0) - 0.5)*fabs(c[atom.group[gr].active[i]]+1e-3); }
    for (i = 0; i < Nf; i++) c0[i] = c[i];
    double an = 0.;
    for (i = 0; i < Nf; i++) { uint32_t j; for (j = 0; j < Nf; j++) an += dc[i]*F[a][i + j*Nf]*c0[j]; }
    an *= 2.*hf.sub[a].qbar;
    double D[2];
    for (int t = 0; t < 2; t++)
    {
      double h = (t == 0) ? 1e-2 : 5e-3, Ep, Em;
      for (i = 0; i < Nf; i++) c[i] = c0[i] + h*dc[i];
      Ep = Eof(&hf);
      for (i = 0; i < Nf; i++) c[i] = c0[i] - h*dc[i];
      Em = Eof(&hf);
      D[t] = (Ep - Em)/(2*h);
    }
    for (i = 0; i < Nf; i++) c[i] = c0[i];
    double fd = (4*D[1] - D[0])/3;
    printf("  sub %u n%d l%d kappa %d qbar %.4f: analytic %.12e  FD %.12e  rel %.1e\n", a, hf.sub[a].n, hf.sub[a].l, hf.sub[a].kappa, hf.sub[a].qbar, an, fd, fabs(an-fd)/fabs(an));
    free(dc); free(c0);
  }
  for (a = 0; a < hf.nSub; a++) free(F[a]);
  free(F);
  HFFree(&hf); AtomFree(&atom); RadialBasisFree(&basis); free(V);
}

int main()
{
  setvbuf(stdout, NULL, _IONBF, 0);
  run(28, "1s2 2s2 2p6 3s2 3p6 3d8 4s2", 1, 0);
  run(8, "1s2 2s2 2p4", 0, 0);
  run(8, "1s2 2s2 2p4", 1, 0);
  run(26, "1s2 2s2 2p6 3s2 3p6 3d6 4s2", 1, 0);
  run(26, "1s2 2s2 2p6 3s2 3p6 3d6 4s2", 0, 0);
  run(26, "1s2 2s2 2p6 3s2 3p6 3d6 4s2", 0, 1);
  return 0;
}
