/*
 * NucleusMultipole.cpp
 *
 * Non-spherical nuclear charge densities, "Function3D" (DFT stage A5c, design section 3.1).
 * The density rho_N(x, y, z) (charge per volume) is expanded in real spherical harmonics of
 * even L <= Lmax,
 *     rho_N(r) = sum_LM rho_LM(r) S_LM(r^),   rho_LM(r) = int dOmega rho_N S_LM,
 * by a product quadrature on every radial quadrature point: Gauss-Legendre in cos(theta)
 * and uniform in phi, exact for the S_LM up to the order of the grids. The potential energy
 * of an electron in the field of multipole LM (design section 3.1, from section 7):
 *     V_LM(r) = -4 pi/(2L+1) Y^L[r^2 rho_LM](r)/r,   V_N(r) = sum_LM V_LM(r) S_LM(r^).
 * L = 0 is the spherical nucleus of the atomic SCF (rho_00 S_00 = rho_00/sqrt(4 pi)); L >= 2
 * is the one-particle operator NewOperator("NuclearMultipole", ...).
 *
 * The real spherical harmonics (no Condon-Shortley phase in them), with the complex Y_LM
 * of Condon and Shortley:
 *     S_L0 = Y_L0,
 *     S_LM  = sqrt2 N_LM P_L^M(cos theta) cos(M phi) = (Y_L,-M + (-1)^M Y_LM)/sqrt2,      M > 0,
 *     S_L,-M = sqrt2 N_LM P_L^M(cos theta) sin(M phi) = i (Y_L,-M - (-1)^M Y_LM)/sqrt2,   M > 0,
 * N_LM = sqrt((2L+1)/(4 pi) (L-M)!/(L+M)!), P_L^M the associated Legendre function without the
 * (-1)^M. The same S_LM define rho_LM and the operator, so the result does not depend on how
 * other parts of Quanty label tesseral harmonics.
 *
 * All routines return 0 on success and 1 on error and do not print.
 */

#include "NucleusMultipole.h"

// S_LM at cos(theta) = x and phi
double RealSphericalHarmonic(int L, int M, double x, double phi)
{
  int m = abs(M), l;
  double pmm = 1., pm1, pl = 0., s = sqrt(fmax(0., (1. - x)*(1. + x))), norm, f;
  if (m > L) return 0.;
  for (l = 1; l <= m; l++) pmm *= (2.*l - 1.)*s;          // P_m^m = (2m-1)!! (1-x^2)^(m/2)
  if (L == m) pl = pmm;
  else
  {
    pm1 = x*(2.*m + 1.)*pmm;                               // P_{m+1}^m
    if (L == m + 1) pl = pm1;
    else
      for (l = m + 2; l <= L; l++)
      {
        pl = ((2.*l - 1.)*x*pm1 - (l + m - 1.)*pmm)/(l - m);
        pmm = pm1;
        pm1 = pl;
      }
  }
  f = 1.;
  for (l = L - m + 1; l <= L + m; l++) f *= l;             // (L+M)!/(L-M)!
  norm = sqrt((2.*L + 1.)/(4.*M_PI)/f);
  if (M == 0) return norm*pl;
  return sqrt(2.)*norm*pl*((M > 0) ? cos(m*phi) : sin(m*phi));
}

// the number of (L, M) with even L <= Lmax; their order: L = 0, 2, 4 .., M = -L .. L
uint32_t NucleusMultipoleCount(uint32_t Lmax)
{
  uint32_t L, n = 0;
  for (L = 0; L <= Lmax; L += 2) n += 2*L + 1;
  return n;
}

void NucleusMultipoleLM(uint32_t i, int *L, int *M)
{
  int l = 0;
  while (i >= (uint32_t) (2*l + 1)) { i -= (uint32_t) (2*l + 1); l += 2; }
  *L = l;
  *M = (int) i - l;
}

/* rho_LM on the radial points: rho[(q*nTheta + i)*nPhi + j] the density at r_q, cos(theta) =
 * xTheta[i] (Gauss-Legendre weights wTheta) and phi = 2 pi j/nPhi. rhoLM[lm*nP + q]. */
int NucleusMultipoleProject(const double *rho, uint32_t nP, uint32_t nTheta, uint32_t nPhi, const double *xTheta,
                            const double *wTheta, uint32_t Lmax, double *rhoLM)
{
  uint32_t nLM = NucleusMultipoleCount(Lmax), lm, q, i, j;
  int L, M;
  double *S = (double *) malloc(nTheta*nPhi*sizeof(double)), w = 2.*M_PI/nPhi;
  for (lm = 0; lm < nLM; lm++)
  {
    NucleusMultipoleLM(lm, &L, &M);
    for (i = 0; i < nTheta; i++)
      for (j = 0; j < nPhi; j++) S[i*nPhi + j] = wTheta[i]*w*RealSphericalHarmonic(L, M, xTheta[i], w*j);
    for (q = 0; q < nP; q++)
    {
      DFTReal sum = 0., comp = 0.;
      for (i = 0; i < nTheta*nPhi; i++) NeumaierAdd(&sum, &comp, S[i]*rho[q*nTheta*nPhi + i]);
      rhoLM[lm*nP + q] = sum + comp;
    }
  }
  free(S);
  return 0;
}

/* V_LM = -4 pi/(2L+1) Y^L[r^2 rho_LM]/r on the quadrature points of b, every (L, M) of
 * NucleusMultipoleCount(Lmax) (VLM[lm*nPoint + q]; L = 0 included). Returns 1 if the Poisson
 * equation fails. */
int NucleusMultipolePotentials(const RadialBasisType *b, uint32_t Lmax, const double *rhoLM, double *VLM)
{
  uint32_t nP = b->nPoint, nLM = NucleusMultipoleCount(Lmax), lm, q, ld = b->degree + 1;
  int L, M, Lband = -1, error = 0;
  double *band = (double *) malloc(ld*(b->nSpline - 1)*sizeof(double)), *f = (double *) malloc(nP*sizeof(double));
  double *Yr = (double *) malloc(nP*sizeof(double));
  for (lm = 0; lm < nLM && !error; lm++)
  {
    NucleusMultipoleLM(lm, &L, &M);
    if (L != Lband)
    {
      error = BSplinePoissonSetup(b->knot, b->nKnot, b->degree, b->rq, b->wq, nP, (uint32_t) L, band, ld);
      Lband = L;
    }
    for (q = 0; q < nP; q++) f[q] = b->rq[q]*b->rq[q]*rhoLM[lm*nP + q];
    if (!error) error = RadialYkOnPoints(b, (uint32_t) L, band, f, Yr);
    for (q = 0; q < nP; q++) VLM[lm*nP + q] = -4.*M_PI/(2.*L + 1.)*Yr[q];
  }
  free(band); free(f); free(Yr);
  return error;
}

// the multipole moment q_LM = int r^L S_LM rho_N d^3r = int r^(L+2) rho_LM dr
double NucleusMultipoleMoment(const RadialBasisType *b, int L, const double *rhoLM)
{
  DFTReal sum = 0., comp = 0.;
  uint32_t q;
  for (q = 0; q < b->nPoint; q++) NeumaierAdd(&sum, &comp, b->wq[q]*pow(b->rq[q], L + 2)*rhoLM[q]);
  return sum + comp;
}

/* <kappaA muA| S_LM |kappaB muB> (muA2 = 2 muA, muB2 = 2 muB) of the spin-angular functions
 * Omega_kappa,mu, identical for the large (Omega_kappa) and the small (Omega_-kappa)
 * components: with S_LM = sum_q u_q Y_Lq (file header) and the Wigner-Eckart theorem in
 * Grant's convention,
 *     <kappa mu|C^L_q|kappa' mu'> = (-1)^(j-mu) (j L j'; -mu q mu') <kappa||C^L||kappa'>,
 *     <kappa||C^L||kappa'> = (-1)^(j+1/2) sqrt((2j+1)(2j'+1)) (j L j'; 1/2 0 -1/2) pi(l+l'+L),
 * which depends on l only through the parity, the same for l and lbar; Y_Lq = sqrt((2L+1)/
 * (4 pi)) C^L_q. */
void RealHarmonicMatrixElement(int kappaA, int muA2, int L, int M, int kappaB, int muB2, double *re, double *im)
{
  int jA2 = 2*abs(kappaA) - 1, jB2 = 2*abs(kappaB) - 1, lA = DiracKappaL(kappaA), lB = DiracKappaL(kappaB), q2 = muA2 - muB2, m = abs(M);
  double red, c, ur = 0., ui = 0., s2 = 1./sqrt(2.);
  *re = *im = 0.;
  if ((lA + lB + L) % 2 != 0 || q2 % 2 != 0) return;
  int q = q2/2;
  if (abs(q) != m) return;
  red = (((jA2 + 1)/2) % 2 == 0 ? 1. : -1.)*sqrt((jA2 + 1.)*(jB2 + 1.))*ThreeJSymbolHalf(jA2, 1, 2*L, 0, jB2, -1);
  c = ((((jA2 - muA2)/2) % 2 == 0) ? 1. : -1.)*ThreeJSymbolHalf(jA2, -muA2, 2*L, q2, jB2, muB2)*red*sqrt((2.*L + 1.)/(4.*M_PI));
  // u_q of S_LM
  if (M == 0) ur = 1.;
  else if (M > 0) ur = (q == -m) ? s2 : ((m % 2 == 0) ? s2 : -s2);
  else ui = (q == -m) ? s2 : ((m % 2 == 0) ? -s2 : s2);
  *re = ur*c;
  *im = ui*c;
}
