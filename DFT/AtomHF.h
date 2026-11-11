/*
 * AtomHF.h
 *
 * Configuration-average (Dirac-)Hartree-Fock atom (DFT stage A4), see AtomHF.cpp.
 */

#ifndef DFT_ATOMHF_H_
#define DFT_ATOMHF_H_

// one subshell of the energy functional: its angular labels, shell and radial function
typedef struct {
  int n, l, kappa;           // kappa = 0: nonrelativistic shell (Cowan's angular factors)
  uint32_t shell;            // index of its nl shell
  double g, qbar;            // degeneracy and mean occupation
  double *coef;              // radial coefficients (not owned)
  uint32_t orbital;          // the DiracAtomType orbital it belongs to
  int closed;                // its shell is full
  int empty;                 // its shell is empty (q = 0): not in the SCF, Fock operator of V^N
  double energy;             // eps_a = c_a^T F_a c_a (HFSCF)
} HFSubshellType;

typedef struct {
  DiracAtomType *atom;       // orbitals, groups, basis (not owned)
  uint32_t nSub, nShell, kMax, nPois;
  int ks;                    // 1: Kohn-Sham operators and energy (atom->method) instead of Hartree-Fock
  HFSubshellType *sub;
  double *shellQ, *shellG;   // occupation q_s and degeneracy G_s (4l+2, or 2j+1 for an n kappa subshell)
  double *Gam;               // Gamma^k_ab at [(a*nSub + b)*(kMax+1) + k]
  double *Rself;             // 2 N_aa / qbar_a
  double *Rpair;             // N_ab / qbar_a at [a*nSub + b]
  double *Npair;             // N_ab (N_aa on the diagonal)
  double **band;             // factorized Poisson matrices, k = 0 .. kMax
  double *P, *Q;             // the radial functions on the quadrature points, [a*nPoint + q]
  // HFSCF results
  uint32_t iterations;
  int converged;
  double Etot, Ekin, Eenuc, Edirect, Eexchange, Eidentity, Eextra, dE, gradient, pairGradient, residual;
} DiracHFType;

double HFGamma(int la, int kappaA, int lb, int kappaB, int k);
int HFSetup(DiracHFType *hf, DiracAtomType *atom, int jjAngular);
void HFFree(DiracHFType *hf);
void HFOrbitalsOnPoints(DiracHFType *hf);
int RadialYkOnPoints(const RadialBasisType *b, uint32_t k, const double *band, const double *rho, double *Yr);
int HFYk(DiracHFType *hf, uint32_t k, const double *rho, double *Yr);
double HFIntegral(const RadialBasisType *b, const double *rho, const double *Yr);
int HFEnergy(DiracHFType *hf, double *Etot, double *Ekin, double *Eenuc, double *Edirect, double *Eexchange);
double HFCoupling(const DiracHFType *hf, uint32_t a, uint32_t b);
int HFShared(const DiracHFType *hf, uint32_t a);
int HFRedundantPair(const DiracHFType *hf, uint32_t a, uint32_t b);
int HFFockMatrices(DiracHFType *hf, double **F);
int HFExchangeMatrix(DiracHFType *hf, uint32_t group, uint32_t b, uint32_t k, double *K);
int HFGradient(DiracHFType *hf, double **F, double *gradient, double *pairGradient);
int HFSCF(DiracHFType *hf, uint32_t maxIterations, uint32_t history, double mixing, double energyTolerance,
          double residualTolerance, int rayleigh);

#endif /* DFT_ATOMHF_H_ */
