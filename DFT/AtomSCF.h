/*
 * AtomSCF.h
 *
 * The self-consistent Dirac-Kohn-Sham atom (DFT stage A3), see AtomSCF.cpp.
 */

#ifndef DFT_ATOMSCF_H_
#define DFT_ATOMSCF_H_

// a shell of the configuration: nl (kappa = 0) or one relativistic subshell n kappa
typedef struct {
  int n, l, kappa;
  double q;                  // its occupation
} AtomShellType;

// one subshell: Dirac (n, kappa), or Schroedinger (n, l) with kappa = 0
typedef struct {
  int n, l, kappa;
  double q;                  // mean occupation qbar
  uint32_t shell;            // its shell of the configuration
  uint32_t group, root;      // symmetry group, root index n-l-1 in its electronic branch
  double energy;             // eigenvalue
  double *coef;              // Nf coefficients (DKB, or B-spline)
  double *extra;             // NULL, or an Nf x Nf one-electron term of this orbital only (A5)
} DiracOrbitalType;

// the orbitals of one kappa (or l)
typedef struct {
  int kappa, l;
  uint32_t N, *active;       // the functions in use (DiracActiveFunctions)
  double *H0, *S, *Hfix;     // kinetic part, overlap, H0 + M[V_N]; Nf x Nf
} DiracGroupType;

typedef struct {
  RadialBasisType *basis;    // not owned
  double Z, c, nElectron;
  int relativistic, method;
  DFTReal cutoff;            // density below which e_xc = 0
  uint32_t Nf;               // coefficients per orbital: 2n (Dirac), n (Schroedinger)
  uint32_t nOrb, nGroup, nShell;
  DiracOrbitalType *orb;
  AtomShellType *shell;      // the shells of the configuration
  DiracGroupType *group;
  DFTReal *VN;               // nuclear potential on the quadrature points
  DFTReal *f, *df;           // 4 pi r^2 n and its derivative
  DFTReal *v, *g;            // the last input potential (local part and GGA coupling)
  double *poisson;           // factorized Poisson matrix, k = 0
  // results
  uint32_t iterations;
  int converged;
  double residual, dE;      // residual of the mixed quantity, change of Etot
  double potentialResidual;  // |v_out - v_in| weighted with the density, in both modes
  double Etot, Ekin, Eenuc, Ecoul, Exc, Eeig;
} DiracAtomType;

int AtomParseConfiguration(const char *configuration, uint32_t maxShell, uint32_t *nShell, AtomShellType *shell);
double AtomShellDegeneracy(const AtomShellType *shell);
void AtomOrbitalOnPoints(const RadialBasisType *basis, int relativistic, int kappa, double c, const double *coef,
                         double *P, double *dP, double *Q, double *dQ);
uint32_t AtomLocalFunctions(const RadialBasisType *basis, int relativistic, int kappa, double c, uint32_t q,
                            uint32_t *idx, double *P, double *dP, double *Q, double *dQ);
void AtomPotentialMatrix(const RadialBasisType *basis, int relativistic, int kappa, double c, const double *v, const double *g, double *M);
int AtomKineticMatrices(const RadialBasisType *basis, int relativistic, int kappa, int l, double c, double *H0, double *S);
double AtomBilinear(const double *x, const double *A, const double *y, uint32_t Nf);
DFTReal AtomThomasFermiScreening(DFTReal Z, DFTReal r);
void AtomDensity(DiracAtomType *atom);
int AtomPotential(DiracAtomType *atom, const double *f, const double *df, double *v, double *g, double *Ecoul, double *Exc);
int AtomDensityPotential(DiracAtomType *atom, double *v, double *g, double *Ecoul, double *Exc);
int AtomSetup(DiracAtomType *atom, RadialBasisType *basis, const DFTReal *VN, double Z, int relativistic, double c, int method,
              int originBoundary, DFTReal cutoff, uint32_t nShell, const AtomShellType *shell);
void AtomFree(DiracAtomType *atom);
void AtomSignConvention(DiracAtomType *atom, uint32_t o);
int AtomSolveOrbitals(DiracAtomType *atom, const double *v, const double *g, int rayleigh);
int AtomSCF(DiracAtomType *atom, uint32_t maxIterations, uint32_t history, double mixing, double energyTolerance,
            double residualTolerance, int rayleigh, int mixDensity);

#endif /* DFT_ATOMSCF_H_ */
