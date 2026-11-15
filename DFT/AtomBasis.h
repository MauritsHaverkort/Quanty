/*
 * AtomBasis.h
 *
 * Atomic basis functions (DFT stage A5): generating calculations per orbital, see
 * AtomBasis.cpp.
 */

#ifndef DFT_ATOMBASIS_H_
#define DFT_ATOMBASIS_H_

// the settings of the self-consistent fields
typedef struct {
  uint32_t maxIterations, history;
  double mixing, energyTolerance;
  double residualTolerance;      // Kohn-Sham ground state (AtomSCF)
  double orbitalTolerance;       // the orbital step of HFSCF (Hartree-Fock, and every generating SCF)
  int rayleigh, mixDensity;
} AtomSCFOptionsType;

// the results per orbital and of the ground state
typedef struct {
  double *coef;                  // nOrb x Nf: the basis functions (before AtomOrthonormalizeBasis: as generated)
  double *energy;                // eigenvalue of its generating equation (with Z+Q and V_beta)
  double *energyGround;          // <a|F_ground|a> of the final (orthonormalized) function
  double *beta;                  // <beta V_beta> of its confinement, int V_beta (P^2 - Q^2) dr
  double *Egenerating;           // total energy of its generating calculation (ground state if none)
  int *generated;                // 1 if it has its own generating calculation
  int *converged;                // its generating SCF converged (1 for the others)
  double **Fground;              // the Fock (or Kohn-Sham) matrices of the ground state, per orbital
  // ground state
  double Etot, Ekin, Eenuc, Edirect, Eexchange, dE, residual, gradient, pairGradient;
  double Eeig;                   // Etot through the eigenvalues (Kohn-Sham: AtomSCF; Hartree-Fock: the identity of HFSCF)
  double potentialResidual;      // Kohn-Sham: the density-weighted potential residual of AtomSCF
  double EavJJ;                  // Schroedinger Hartree-Fock: E_av with the angular factors of the j subshells
  uint32_t iterations;
  int groundConverged;
} AtomBasisType;

int AtomTargetMatrices(DiracAtomType *atom, uint32_t o, double Q, const DFTReal *Vbeta, double *extra, double *beta);
int AtomRootSolve(DiracAtomType *atom, uint32_t group, const double *F, uint32_t root, int rayleigh, double *coef, double *energy);
int AtomBasisFunctions(DiracAtomType *atom, int hartreeFock, const double *Q, const double *rc, double Nconf,
                       const AtomSCFOptionsType *opt, AtomBasisType *res);
int AtomOrthonormalizeBasis(DiracAtomType *atom, double *coef, int method, double redundancy, double *minEigenvalue, double *orthoError);
void AtomBasisFree(AtomBasisType *res, uint32_t nOrb);

#endif /* DFT_ATOMBASIS_H_ */
