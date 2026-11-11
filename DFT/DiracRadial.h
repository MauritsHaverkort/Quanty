/*
 * DiracRadial.h
 *
 * Radial Dirac equation on a dual-kinetic-balance (DKB) B-spline basis (DFT stage A,
 * design sections 2-5): radial basis and quadrature, nucleus models, DKB matrices, the
 * generalized eigenproblem and root selection. See DiracRadial.cpp.
 */

#ifndef DFT_DIRACRADIAL_H_
#define DFT_DIRACRADIAL_H_

// Gauss points per knot interval at most: the size of the 34-digit table in
// Core/BasicMath/Interpolation/GaussLegendre; above it the weights are good to 4e-10 only
#define DiracMaxQuadraturePoints 127

// on-site orthonormalization (design section 6.1)
#define DiracGramSchmidt 0
#define DiracLowdin      1

// nucleus models (design section 3)
#define NucleusPoint    0
#define NucleusGaussian 1
#define NucleusFermi    2
#define NucleusFunction 3

// knots, quadrature points and weights of one radial B-spline basis
typedef struct {
  uint32_t degree;          // polynomial degree d (order d+1)
  uint32_t nInterval;       // knot intervals, [0, r1] included
  uint32_t nKnot, nSpline;  // nSpline = nKnot - degree - 1; the DKB set uses B_1 .. B_{nSpline-2}
  uint32_t nQ, nPoint;      // Gauss points per interval, all quadrature points
  double *knot, *rq, *wq;
} RadialBasisType;

typedef struct {
  int model;                // NucleusPoint, NucleusGaussian, NucleusFermi, NucleusFunction
  DFTReal Z;                // charge
  DFTReal A;                // mass number (0 if not given)
  DFTReal Rrms;             // root-mean-square radius (bohr); Gaussian and Fermi
  DFTReal cFermi, aFermi;   // Fermi half-density radius and diffuseness (bohr)
} NucleusType;

int RadialBasisCreate(uint32_t degree, uint32_t nInterval, double r1, double rMax, uint32_t nQ, RadialBasisType *basis);
int RadialBasisFromKnots(const double *knot, uint32_t nKnot, uint32_t degree, uint32_t nQ, RadialBasisType *basis);
void RadialBasisFree(RadialBasisType *basis);

int NucleusFermiParameters(DFTReal Rrms, DFTReal t, DFTReal *cFermi, DFTReal *aFermi);
int NucleusDensity(const NucleusType *nucleus, const RadialBasisType *basis, DFTReal *rho);
DFTReal NucleusCharge(const RadialBasisType *basis, const DFTReal *rho);
int NucleusPotential(const NucleusType *nucleus, const RadialBasisType *basis, const DFTReal *rho, DFTReal *V);

int DiracKappaL(int kappa);
int DiracDKBMatrices(const RadialBasisType *basis, int kappa, double c, const DFTReal *V, const DFTReal *Vbeta, double *H, double *S);
void DiracActiveFunctions(const RadialBasisType *basis, int kappa, int originBoundary, uint32_t *active, uint32_t *nActive, uint32_t *nLarge);
int DiracRadialEigen(const RadialBasisType *basis, int kappa, double c, const DFTReal *V, const DFTReal *Vbeta, int rayleigh,
                     int originBoundary, double *E, double *C, uint32_t *nE, double *residual, double *Eeigensolver);
int DiracGeneralizedEigen(const double *Hf, const double *Sf, uint32_t Nf, const uint32_t *active, uint32_t N, int rayleigh,
                          double *E, double *C, double *residual, double *Eeigensolver);
uint32_t DiracElectronicFirst(const double *E, uint32_t nE, double c);
int DiracRootIndex(uint32_t n, int kappa, uint32_t *index);
int DiracComponents(const RadialBasisType *basis, int kappa, double c, const double *coef, const double *r, uint32_t nr,
                    double *P, double *Q);
int DiracSignConvention(const RadialBasisType *basis, int kappa, double c, double *coef);

int ConfinementPotential(const RadialBasisType *basis, DFTReal rc, DFTReal N, DFTReal *Vbeta);
DFTReal ConfinementRMax(DFTReal rc, DFTReal N, DFTReal Vmax);
DFTReal DiracExpectation(const RadialBasisType *basis, int kappa, double c, const double *coef, const DFTReal *f, int sign);
int DiracDKBOverlap(const RadialBasisType *basis, int kappa, double c, double *S);
int DiracOrthonormalize(const double *S, uint32_t N, double *C, uint32_t m, int method, double redundancy, double *minEigenvalue);
int DiracOrthonormalizePass(const double *S, uint32_t N, double *C, uint32_t m, int method, double redundancy, double *minEigenvalue);
double DiracOrthonormalityError(const double *S, uint32_t N, const double *C, uint32_t m);

#endif /* DFT_DIRACRADIAL_H_ */
