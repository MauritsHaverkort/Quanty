/*
 * DiracRadial.cpp
 *
 * The radial Dirac equation for one electron in a spherical potential V (vector) and
 * V_beta (Lorentz scalar), on a dual-kinetic-balance (DKB) B-spline basis. DFT/HF add-on,
 * stage A; design Doc/Claude/QuantyDFT_StageA/QuantyDFT_StageA_DiracAtom.md, sections
 * 2-5. Hartree atomic units, energies measured from mc^2, c an argument.
 *
 * Radial spinor psi = (1/r) (P Omega_kappa,mu ; i Q Omega_-kappa,mu), int (P^2+Q^2) dr = 1,
 * and
 *     ( V + V_beta              c (-d/dr + kappa/r) ) (P)     (P)
 *     ( c (d/dr + kappa/r)      V - V_beta - 2c^2   ) (Q) = E (Q).
 *
 * All routines return 0 on success and 1 on error and do not print.
 */

#include "DiracRadial.h"

// compensated (Neumaier) summation: sum + comp is the sum of everything added
void NeumaierAdd(DFTReal *sum, DFTReal *comp, DFTReal x)
{
  DFTReal t = *sum + x;
  if (fabs(*sum) >= fabs(x)) *comp += (*sum - t) + x;
  else *comp += (x - t) + *sum;
  *sum = t;
}

/* The radial basis: exponential knots on [0, rMax] (BSplineKnotsExponential) with
 * nInterval intervals, the first [0, r1], and nQ Gauss points per interval. The caller
 * frees it with RadialBasisFree. Returns 1 for nQ = 0 or nQ > DiracMaxQuadraturePoints,
 * for invalid knots, or for fewer than 3 B-splines. */
int RadialBasisCreate(uint32_t degree, uint32_t nInterval, double r1, double rMax, uint32_t nQ, RadialBasisType *basis)
{
  double *knot;
  uint32_t nKnot;
  int error;
  basis->knot = NULL;
  basis->rq = NULL;
  basis->wq = NULL;
  if (degree == 0 || degree > BSplineMaxDegree) return 1;
  if (BSplineKnotsExponential(degree, nInterval, r1, rMax, &knot, &nKnot)) return 1;
  error = RadialBasisFromKnots(knot, nKnot, degree, nQ, basis);
  free(knot);
  return error;
}

/* The radial basis for a given knot vector (copied), which must start with degree+1
 * knots at 0, with nQ Gauss points per non-empty interval. The caller frees it with
 * RadialBasisFree. Returns 1 for nQ = 0 or nQ > DiracMaxQuadraturePoints, knots that do
 * not start at 0, or fewer than 3 B-splines. */
int RadialBasisFromKnots(const double *knot, uint32_t nKnot, uint32_t degree, uint32_t nQ, RadialBasisType *basis)
{
  uint32_t i;
  basis->knot = NULL;
  basis->rq = NULL;
  basis->wq = NULL;
  if (nQ == 0 || nQ > DiracMaxQuadraturePoints || degree == 0 || degree > BSplineMaxDegree || nKnot < degree + 4) return 1;
  if (knot[0] != 0. || knot[degree] != 0.) return 1;
  basis->knot = (double *) malloc(nKnot*sizeof(double));
  for (i = 0; i < nKnot; i++) basis->knot[i] = knot[i];
  basis->nKnot = nKnot;
  basis->degree = degree;
  basis->nInterval = 0;
  for (i = 0; i+1 < nKnot; i++) if (knot[i] < knot[i+1]) basis->nInterval++;
  basis->nSpline = nKnot - degree - 1;
  basis->nQ = nQ;
  if (BSplineQuadratureSetup(basis->knot, basis->nKnot, nQ, &(basis->rq), &(basis->wq), &(basis->nPoint)))
  {
    RadialBasisFree(basis);
    return 1;
  }
  return 0;
}

void RadialBasisFree(RadialBasisType *basis)
{
  free(basis->knot);
  free(basis->rq);
  free(basis->wq);
  basis->knot = NULL;
  basis->rq = NULL;
  basis->wq = NULL;
}

/* Fermi nucleus rho_0 / (1 + exp((r - cFermi)/aFermi)) with skin thickness t (10% to 90%
 * of the density): aFermi = t/(4 ln 3), and cFermi from the root-mean-square radius,
 *     Rrms^2 = 3/5 cFermi^2 + 7/5 pi^2 aFermi^2
 * (Visscher and Dyall 1997; the terms of order exp(-cFermi/aFermi) are dropped). Returns 1
 * if Rrms^2 <= 7/5 pi^2 aFermi^2, i.e. for nuclei too small for the skin thickness. */
int NucleusFermiParameters(DFTReal Rrms, DFTReal t, DFTReal *cFermi, DFTReal *aFermi)
{
  DFTReal a = t/(4.*log(3.)), c2;
  c2 = 5./3.*(Rrms*Rrms - 7./5.*M_PI*M_PI*a*a);
  if (!(t > 0.) || !(c2 > 0.)) return 1;
  *cFermi = sqrt(c2);
  *aFermi = a;
  return 0;
}

// 4 pi int r^2 rho dr on the quadrature points, compensated summation
DFTReal NucleusCharge(const RadialBasisType *basis, const DFTReal *rho)
{
  DFTReal sum = 0., comp = 0.;
  for (uint32_t q = 0; q < basis->nPoint; q++) NeumaierAdd(&sum, &comp, basis->wq[q]*basis->rq[q]*basis->rq[q]*rho[q]);
  return 4.*M_PI*(sum + comp);
}

/* The nuclear charge density (per volume) on the quadrature points, normalized to Z:
 *     Gaussian  Z (zeta/pi)^(3/2) exp(-zeta r^2),   zeta = 3/(2 Rrms^2)
 *     Fermi     rho_0 / (1 + exp((r - cFermi)/aFermi)), rho_0 from 4 pi int r^2 rho = Z
 * Returns 1 for the point nucleus and for "Function" (whose density the caller has). */
int NucleusDensity(const NucleusType *nucleus, const RadialBasisType *basis, DFTReal *rho)
{
  uint32_t q;
  DFTReal zeta, norm;
  if (nucleus->model == NucleusGaussian)
  {
    zeta = 1.5/(nucleus->Rrms*nucleus->Rrms);
    for (q = 0; q < basis->nPoint; q++) rho[q] = nucleus->Z*pow(zeta/M_PI, 1.5)*exp(-zeta*basis->rq[q]*basis->rq[q]);
    return 0;
  }
  if (nucleus->model == NucleusUniform)
  {
    DFTReal R = sqrt(5./3.)*nucleus->Rrms;
    for (q = 0; q < basis->nPoint; q++) rho[q] = (basis->rq[q] < R) ? 3.*nucleus->Z/(4.*M_PI*R*R*R) : 0.;
    return 0;
  }
  if (nucleus->model == NucleusFermi)
  {
    for (q = 0; q < basis->nPoint; q++) rho[q] = 1./(1. + exp((basis->rq[q] - nucleus->cFermi)/nucleus->aFermi));
    norm = NucleusCharge(basis, rho);
    if (!(norm > 0.)) return 1;
    for (q = 0; q < basis->nPoint; q++) rho[q] *= nucleus->Z/norm;
    return 0;
  }
  return 1;
}

/* The potential energy V_N of an electron in the field of the nucleus, on the quadrature
 * points (design section 3):
 *     Point     -Z/r
 *     Gaussian  -Z erf(sqrt(zeta) r)/r,   zeta = 3/(2 Rrms^2)
 *     Uniform   a homogeneously charged sphere, R = sqrt(5/3) Rrms (FPLO's "constant density"):
 *               -Z (3 - r^2/R^2)/(2R) inside, -Z/r outside (A5e)
 *     Fermi     -Y^0[4 pi r^2 rho](r)/r   (radial Poisson equation, design section 7)
 *     Function  the same, with the density rho (per volume) given on the quadrature
 *               points; it should be normalized to Z, the caller checks that.
 *     Function3D the same with the spherical part rho_00/sqrt(4 pi) (NucleusMultipole.cpp).
 * rho is used only for "Function" and "Function3D" (NULL otherwise). Returns 1 for an unknown model or if
 * the Poisson solver fails. */
int NucleusPotential(const NucleusType *nucleus, const RadialBasisType *basis, const DFTReal *rho, DFTReal *V)
{
  uint32_t q, ld = basis->degree + 1, n = basis->nSpline - 1;
  DFTReal zeta, *f, *band, *Y;
  const DFTReal *rhoN = rho;
  int error = 0;
  if (nucleus->model == NucleusPoint)
  {
    for (q = 0; q < basis->nPoint; q++) V[q] = -nucleus->Z/basis->rq[q];
    return 0;
  }
  if (nucleus->model == NucleusGaussian)
  {
    zeta = 1.5/(nucleus->Rrms*nucleus->Rrms);
    for (q = 0; q < basis->nPoint; q++) V[q] = -nucleus->Z*erf(sqrt(zeta)*basis->rq[q])/basis->rq[q];
    return 0;
  }
  if (nucleus->model == NucleusUniform)
  {
    DFTReal R = sqrt(5./3.)*nucleus->Rrms, r;
    for (q = 0; q < basis->nPoint; q++)
    {
      r = basis->rq[q];
      V[q] = (r < R) ? -nucleus->Z*(3. - r*r/(R*R))/(2.*R) : -nucleus->Z/r;
    }
    return 0;
  }
  if (nucleus->model != NucleusFermi && nucleus->model != NucleusFunction && nucleus->model != NucleusFunction3D) return 1;
  if (nucleus->model != NucleusFermi && rho == NULL) return 1;
  f = (DFTReal *) malloc(basis->nPoint*sizeof(DFTReal));
  if (nucleus->model == NucleusFermi)
  {
    if (NucleusDensity(nucleus, basis, f)) { free(f); return 1; }
    rhoN = f;
  }
  // the radial charge density 4 pi r^2 rho (in dr), overwriting f
  for (q = 0; q < basis->nPoint; q++) f[q] = 4.*M_PI*basis->rq[q]*basis->rq[q]*rhoN[q];
  band = (DFTReal *) malloc(ld*n*sizeof(DFTReal));
  Y = (DFTReal *) malloc(basis->nSpline*sizeof(DFTReal));
  error = BSplinePoissonSetup(basis->knot, basis->nKnot, basis->degree, basis->rq, basis->wq, basis->nPoint, 0, band, ld);
  if (!error) error = BSplinePoissonSolve(basis->knot, basis->nKnot, basis->degree, basis->rq, basis->wq, basis->nPoint, 0, band, ld, f, Y);
  if (!error)
    for (q = 0; q < basis->nPoint; q++) V[q] = -BSplineDerivative(basis->rq[q], basis->knot, basis->nKnot, basis->degree, Y, 0)/basis->rq[q];
  free(f);
  free(band);
  free(Y);
  return error;
}

// l of the large component: kappa for kappa > 0, -kappa-1 for kappa < 0
int DiracKappaL(int kappa)
{
  return (kappa > 0) ? kappa : -kappa - 1;
}

/* The dense matrix M[i + j*n] = int B_{1+i}^(a) f B_{1+j}^(b) dr, i, j = 0 .. n-1, over the
 * working B-splines B_1 .. B_{nSpline-2}, with f on the quadrature points (NULL for 1).
 * BSplineMatrixBanded stores i >= j only; for a != b the upper triangle comes from the
 * call with a and b swapped, M^(a,b)_ij = M^(b,a)_ji. */
int DiracDenseMatrix(const RadialBasisType *basis, const double *fq, uint32_t a, uint32_t b, double *M)
{
  uint32_t n = basis->nSpline - 2, ld = basis->degree + 1, i, j;
  double *band = (double *) malloc(ld*n*sizeof(double));
  int error = BSplineMatrixBanded(basis->knot, basis->nKnot, basis->degree, basis->rq, basis->wq, basis->nPoint, fq, a, b, 1, n, band, ld);
  for (i = 0; i < n*n; i++) M[i] = 0.;
  for (j = 0; j < n && !error; j++)
    for (i = j; i < n && i <= j + basis->degree; i++)
    {
      M[i + j*n] = band[(i-j) + j*ld];
      if (a == b) M[j + i*n] = band[(i-j) + j*ld];
    }
  if (a != b && !error)
  {
    error = BSplineMatrixBanded(basis->knot, basis->nKnot, basis->degree, basis->rq, basis->wq, basis->nPoint, fq, b, a, 1, n, band, ld);
    for (j = 0; j < n && !error; j++)
      for (i = j+1; i < n && i <= j + basis->degree; i++) M[j + i*n] = band[(i-j) + j*ld];
  }
  free(band);
  return error;
}

// block (row0, col0) of the N x N matrix H += coef M (or coef M^T); M is n x n
void DiracAddBlock(double *H, uint32_t N, uint32_t row0, uint32_t col0, uint32_t n, double coef, const double *M, int transpose)
{
  uint32_t i, j;
  for (j = 0; j < n; j++)
    for (i = 0; i < n; i++)
      H[(row0 + i) + (col0 + j)*N] += coef*(transpose ? M[j + i*n] : M[i + j*n]);
}

/* The DKB matrices H and S of size N = 2n, n = nSpline-2, dense and column major (both
 * triangles), in the basis u_i^(1) (indices 0 .. n-1) and u_i^(2) (n .. 2n-1) of design
 * section 5:
 *     u_i^(1) = ( B_i , h (B_i' + kappa B_i/r) ),   u_i^(2) = ( h (B_i' - kappa B_i/r) , B_i ),
 * h = 1/(2c), with
 *     H_ab = int [ P_a (V+Vb) P_b + Q_a (V-Vb-2c^2) Q_b + c (P_a' Q_b + Q_a P_b')
 *                  + c kappa/r (P_a Q_b + Q_a P_b) ] dr,      S_ab = int (P_a P_b + Q_a Q_b) dr.
 * V and Vbeta are given on the quadrature points (Vbeta NULL for 0).
 *
 * Expanding P, Q, P', Q' in B, B', B'' and powers of 1/r point by point (no integration by
 * parts, so the quadrature sums equal those of the direct evaluation) gives, with
 * M^(a,b)_m[g] = int B_i^(a) g r^-m B_j^(b) dr, Wp = V+Vb, Wm = V-Vb,
 *     K+-[g] = M^(1,1)_0 +- kappa (M^(1,0)_1 + M^(0,1)_1) + kappa^2 M^(0,0)_2     [g],
 *     H11 = M^(0,0)_0[Wp] + h^2 K+[Wm] + 1/2 K+[1]
 *     H22 = h^2 K-[Wp] + M^(0,0)_0[Wm] - 2c^2 M^(0,0)_0[1] + 1/2 (M^(2,0)_0 + M^(0,2)_0)[1]
 *           - kappa(kappa-1) M^(0,0)_2[1]
 *     H12 = h (M^(0,1)_0 - kappa M^(0,0)_1)[Wp] + h (M^(1,0)_0 + kappa M^(0,0)_1)[Wm]
 *           + h/2 (M^(1,2)_0 + kappa M^(0,2)_1 - kappa(kappa-1) M^(1,0)_2
 *                  - kappa^2(kappa-1) M^(0,0)_3)[1]
 *     S11 = M^(0,0)_0[1] + h^2 K+[1],   S22 = M^(0,0)_0[1] + h^2 K-[1],
 *     S12 = h (M^(0,1)_0 + M^(1,0)_0)[1],   H21 = H12^T,   S21 = S12^T.
 * In H12 the terms -c (B_i' + kappa B_i/r) B_j of the kinetic and of the -2c^2 part
 * cancel exactly and are left out. This list was checked against the direct evaluation
 * of the prototype (dkb.py) to 1e-16 relative for kappa = -4 .. 3. The integrands with
 * r^-3 and r^-2 B' are not integrable at r = 0 for the first working spline (B_1 ~ r) when
 * |kappa| >= 2; the Gauss sum is finite, and the eigenvalues are insensitive to it.
 * DiracRadialEigen leaves these functions out by default (DiracActiveFunctions). */
int DiracDKBMatrices(const RadialBasisType *basis, int kappa, double c, const DFTReal *V, const DFTReal *Vbeta, double *H, double *S)
{
  uint32_t n = basis->nSpline - 2, N = 2*n, q, i, j, m;
  double h = 1./(2.*c), k = (double) kappa;
  double *M, *f[3][4], *g;
  int error = 0, w;
  if (basis->nSpline < 3 || !(c > 0.)) return 1;
  M = (double *) malloc(n*n*sizeof(double));
  // f[w][m] = g_w / r^m with g_0 = 1, g_1 = Wp = V+Vb, g_2 = Wm = V-Vb
  for (w = 0; w < 3; w++)
    for (m = 0; m < 4; m++)
    {
      f[w][m] = (double *) malloc(basis->nPoint*sizeof(double));
      for (q = 0; q < basis->nPoint; q++)
      {
        g = f[w][m];
        if (w == 0) g[q] = 1.;
        else if (w == 1) g[q] = V[q] + ((Vbeta != NULL) ? Vbeta[q] : 0.);
        else g[q] = V[q] - ((Vbeta != NULL) ? Vbeta[q] : 0.);
        g[q] /= pow(basis->rq[q], (double) m);
      }
    }
  for (i = 0; i < N*N; i++) { H[i] = 0.; S[i] = 0.; }

  // g = 1
  error |= DiracDenseMatrix(basis, f[0][0], 0, 0, M);
  DiracAddBlock(S, N, 0, 0, n, 1., M, 0);
  DiracAddBlock(S, N, n, n, n, 1., M, 0);
  DiracAddBlock(H, N, n, n, n, -2.*c*c, M, 0);
  error |= DiracDenseMatrix(basis, f[0][0], 1, 1, M);                  // K+-: M^(1,1)_0
  DiracAddBlock(S, N, 0, 0, n, h*h, M, 0);
  DiracAddBlock(S, N, n, n, n, h*h, M, 0);
  DiracAddBlock(H, N, 0, 0, n, 0.5, M, 0);
  error |= DiracDenseMatrix(basis, f[0][1], 1, 0, M);                  // K+-: M^(1,0)_1 + M^(0,1)_1
  DiracAddBlock(S, N, 0, 0, n, h*h*k, M, 0);
  DiracAddBlock(S, N, 0, 0, n, h*h*k, M, 1);
  DiracAddBlock(S, N, n, n, n, -h*h*k, M, 0);
  DiracAddBlock(S, N, n, n, n, -h*h*k, M, 1);
  DiracAddBlock(H, N, 0, 0, n, 0.5*k, M, 0);
  DiracAddBlock(H, N, 0, 0, n, 0.5*k, M, 1);
  error |= DiracDenseMatrix(basis, f[0][2], 0, 0, M);                  // K+-: M^(0,0)_2
  DiracAddBlock(S, N, 0, 0, n, h*h*k*k, M, 0);
  DiracAddBlock(S, N, n, n, n, h*h*k*k, M, 0);
  DiracAddBlock(H, N, 0, 0, n, 0.5*k*k, M, 0);
  DiracAddBlock(H, N, n, n, n, -k*(k-1.), M, 0);
  error |= DiracDenseMatrix(basis, f[0][0], 2, 0, M);                  // M^(2,0)_0 + M^(0,2)_0
  DiracAddBlock(H, N, n, n, n, 0.5, M, 0);
  DiracAddBlock(H, N, n, n, n, 0.5, M, 1);
  error |= DiracDenseMatrix(basis, f[0][0], 0, 1, M);                  // S12
  DiracAddBlock(S, N, 0, n, n, h, M, 0);
  DiracAddBlock(S, N, 0, n, n, h, M, 1);
  error |= DiracDenseMatrix(basis, f[0][0], 1, 2, M);                  // H12 kinetic
  DiracAddBlock(H, N, 0, n, n, 0.5*h, M, 0);
  error |= DiracDenseMatrix(basis, f[0][1], 0, 2, M);
  DiracAddBlock(H, N, 0, n, n, 0.5*h*k, M, 0);
  error |= DiracDenseMatrix(basis, f[0][2], 1, 0, M);
  DiracAddBlock(H, N, 0, n, n, -0.5*h*k*(k-1.), M, 0);
  error |= DiracDenseMatrix(basis, f[0][3], 0, 0, M);
  DiracAddBlock(H, N, 0, n, n, -0.5*h*k*k*(k-1.), M, 0);

  // g = Wp = V + Vb
  error |= DiracDenseMatrix(basis, f[1][0], 0, 0, M);
  DiracAddBlock(H, N, 0, 0, n, 1., M, 0);
  error |= DiracDenseMatrix(basis, f[1][0], 1, 1, M);                  // h^2 K-[Wp] in H22
  DiracAddBlock(H, N, n, n, n, h*h, M, 0);
  error |= DiracDenseMatrix(basis, f[1][1], 1, 0, M);
  DiracAddBlock(H, N, n, n, n, -h*h*k, M, 0);
  DiracAddBlock(H, N, n, n, n, -h*h*k, M, 1);
  error |= DiracDenseMatrix(basis, f[1][2], 0, 0, M);
  DiracAddBlock(H, N, n, n, n, h*h*k*k, M, 0);
  error |= DiracDenseMatrix(basis, f[1][0], 0, 1, M);                  // H12
  DiracAddBlock(H, N, 0, n, n, h, M, 0);
  error |= DiracDenseMatrix(basis, f[1][1], 0, 0, M);
  DiracAddBlock(H, N, 0, n, n, -h*k, M, 0);

  // g = Wm = V - Vb
  error |= DiracDenseMatrix(basis, f[2][0], 0, 0, M);
  DiracAddBlock(H, N, n, n, n, 1., M, 0);
  error |= DiracDenseMatrix(basis, f[2][0], 1, 1, M);                  // h^2 K+[Wm] in H11
  DiracAddBlock(H, N, 0, 0, n, h*h, M, 0);
  error |= DiracDenseMatrix(basis, f[2][1], 1, 0, M);
  DiracAddBlock(H, N, 0, 0, n, h*h*k, M, 0);
  DiracAddBlock(H, N, 0, 0, n, h*h*k, M, 1);
  error |= DiracDenseMatrix(basis, f[2][2], 0, 0, M);
  DiracAddBlock(H, N, 0, 0, n, h*h*k*k, M, 0);
  error |= DiracDenseMatrix(basis, f[2][0], 1, 0, M);                  // H12
  DiracAddBlock(H, N, 0, n, n, h, M, 0);
  error |= DiracDenseMatrix(basis, f[2][1], 0, 0, M);
  DiracAddBlock(H, N, 0, n, n, h*k, M, 0);

  // the lower-left blocks are the transposes of the upper-right ones
  for (j = 0; j < n; j++)
    for (i = 0; i < n; i++)
    {
      H[(n + j) + i*N] = H[i + (n + j)*N];
      S[(n + j) + i*N] = S[i + (n + j)*N];
    }
  for (w = 0; w < 3; w++)
    for (m = 0; m < 4; m++) free(f[w][m]);
  free(M);
  return error;
}

/* The DKB functions in use. On the first working spline B_1 (~ r), u_1^(1) has
 * Q(0) = h a B_1'(0)(1+kappa) and u_1^(2) has P(0) = h b B_1'(0)(1-kappa). With
 * originBoundary != 0 (default, decided by Maurits, 2026-11-11):
 *   kappa = +1     u_1^(1) is left out: its Q(0) != 0 is bound without limit by the -Z/r
 *                  of a point nucleus and makes a spurious state (-16542.7 Ha for U);
 *   |kappa| >= 2   both are left out: P(0) = Q(0) = 0, and the integrands with B_1^2/r^3
 *                  etc. that are logarithmically divergent are gone;
 *   kappa = -1     both are kept: u_1^(2), with Q = B_1 ~ r, is the only independent linear
 *                  term of Q at the origin, which a point nucleus needs where -Z/r exceeds
 *                  2c^2 and kinetic balance fails (without it Q' of hydrogen 1s is 500 times
 *                  worse at r = 1e-3); kappa = -1 has no spurious state.
 * active[0 .. *nActive-1] receives the indices (0 .. 2n-1, u^(1) first) of the functions
 * in use, and *nLarge the number of them of type u^(1). */
void DiracActiveFunctions(const RadialBasisType *basis, int kappa, int originBoundary, uint32_t *active, uint32_t *nActive, uint32_t *nLarge)
{
  uint32_t n = basis->nSpline - 2, i, m = 0;
  *nLarge = 0;
  for (i = 0; i < 2*n; i++)
  {
    if (originBoundary && i == 0 && kappa != -1) continue;
    if (originBoundary && i == n && kappa != 1 && kappa != -1) continue;
    active[m++] = i;
    if (i < n) (*nLarge)++;
  }
  *nActive = m;
}

/* The eigenvalues E (ascending) and eigenvectors of H c = E S c (DiracGeneralizedEigen) in the DKB
 * functions in use (DiracActiveFunctions); *nE receives their number, at most 2n. The
 * eigenvectors are returned with all 2n DKB coefficients (0 for the functions left out),
 * column j at C[i + j*2n], normalized as C^T S C = 1; E and C must have room for 2n and
 * (2n)^2 values. Electronic states are those with E > -c^2 (DiracElectronicFirst); below
 * them lie the states of the negative continuum, one per function u^(2) in use.
 *
 * rayleigh != 0 replaces every eigenvalue by the Rayleigh quotient c^T H c / c^T S c of
 * its eigenvector, with the H and S as assembled. The eigenvalues of a dense solver carry
 * a rounding error of about DBL_EPSILON times the largest |E|, which is ~1e12 Ha for a
 * point nucleus (r1 = 1e-6/Z), i.e. 1e-8 Ha; the Rayleigh quotient is second order in the
 * eigenvector error and brings them to the discretization error (Z = 80: 1e-10, Z = 1:
 * 1e-13). If Eeigensolver is not NULL it receives the eigenvalues of the eigensolver itself.
 * If residual is not NULL it receives max_i |(H c_j - E_j S c_j)_i| for every j, with the
 * returned E_j. */
int DiracRadialEigen(const RadialBasisType *basis, int kappa, double c, const DFTReal *V, const DFTReal *Vbeta, int rayleigh,
                     int originBoundary, double *E, double *C, uint32_t *nE, double *residual, double *Eeigensolver)
{
  uint32_t Nf = 2*(basis->nSpline - 2), nLarge, *active;
  double *Hf, *Sf;
  int error;
  Hf = (double *) malloc(Nf*Nf*sizeof(double));
  Sf = (double *) malloc(Nf*Nf*sizeof(double));
  if (DiracDKBMatrices(basis, kappa, c, V, Vbeta, Hf, Sf)) { free(Hf); free(Sf); return 1; }
  active = (uint32_t *) malloc(Nf*sizeof(uint32_t));
  DiracActiveFunctions(basis, kappa, originBoundary, active, nE, &nLarge);
  error = DiracGeneralizedEigen(Hf, Sf, Nf, active, *nE, rayleigh, E, C, residual, Eeigensolver);
  free(Hf);
  free(Sf);
  free(active);
  return error;
}

/* H c = E S c for the N functions active[0 .. N-1] of the Nf x Nf matrices Hf, Sf (dense,
 * column major, the lower triangle is used). E (ascending, N values) and the eigenvectors
 * with all Nf coefficients (0 outside the active functions), column j at C[i + j*Nf],
 * normalized as C^T S C = 1; C must have room for Nf*N values. rayleigh, residual and
 * Eeigensolver as in DiracRadialEigen. Returns 1 if LAPACK fails.
 * dsygv (QR iteration), not dsygvd (divide and conquer): the reduced tridiagonal matrix is
 * strongly graded (largest eigenvalue 2e14 for the Schroedinger basis with r1 = 1e-6,
 * 1e12 for the DKB basis of a point nucleus), and divide and conquer gives the small
 * eigenvalues only to DBL_EPSILON max |E| absolute: hydrogen 1s = -0.4602 instead of -0.5
 * (Schroedinger; the same number with Accelerate's and with scipy's LAPACK), and 1e-8 Ha
 * noise for the DKB basis. QR iteration keeps them: -0.5 to 1e-12, and 5e-10 instead of
 * 1e-8 Ha for the DKB basis. Diagonal scaling of S does not help (the spectrum stays). */
int DiracGeneralizedEigen(const double *Hf, const double *Sf, uint32_t Nf, const uint32_t *active, uint32_t N, int rayleigh,
                          double *E, double *C, double *residual, double *Eeigensolver)
{
  uint32_t lwork, itype = 1, i, j;
  double *H, *S, *Cr, *Sc, *work, *HC, *SC, one = 1., zero = 0., chc, csc;
  int info = 0;
  char jobz = 'V', uplo = 'L', nt = 'N';
  H = (double *) malloc(N*N*sizeof(double));
  S = (double *) malloc(N*N*sizeof(double));
  for (j = 0; j < N; j++)
    for (i = 0; i < N; i++)
    {
      H[i + j*N] = Hf[active[i] + active[j]*Nf];
      S[i + j*N] = Sf[active[i] + active[j]*Nf];
    }
  Cr = (double *) malloc(N*N*sizeof(double));
  for (i = 0; i < N*N; i++) Cr[i] = H[i];
  Sc = (double *) malloc(N*N*sizeof(double));         // dsygv overwrites it with its Cholesky factor
  for (i = 0; i < N*N; i++) Sc[i] = S[i];
  lwork = 64*N;
  work = (double *) malloc(lwork*sizeof(double));
  dsygv(&itype, &jobz, &uplo, &N, Cr, &N, Sc, &N, E, work, &lwork, &info);
  free(work);
  free(Sc);
  if (info == 0 && Eeigensolver != NULL) for (j = 0; j < N; j++) Eeigensolver[j] = E[j];
  if (info == 0 && (rayleigh || residual != NULL))
  {
    // H C and S C with the matrices as assembled (both triangles)
    HC = (double *) malloc(N*N*sizeof(double));
    SC = (double *) malloc(N*N*sizeof(double));
    dgemm(&nt, &nt, &N, &N, &N, &one, H, &N, Cr, &N, &zero, HC, &N);
    dgemm(&nt, &nt, &N, &N, &N, &one, S, &N, Cr, &N, &zero, SC, &N);
    for (j = 0; j < N; j++)
    {
      if (rayleigh)
      {
        chc = 0.;
        csc = 0.;
        for (i = 0; i < N; i++)
        {
          chc += Cr[i + j*N]*HC[i + j*N];
          csc += Cr[i + j*N]*SC[i + j*N];
        }
        E[j] = chc/csc;
      }
      if (residual != NULL)
      {
        residual[j] = 0.;
        for (i = 0; i < N; i++) residual[j] = fmax(residual[j], fabs(HC[i + j*N] - E[j]*SC[i + j*N]));
      }
    }
    free(HC);
    free(SC);
  }
  // the eigenvectors with all Nf coefficients
  for (i = 0; i < Nf*N; i++) C[i] = 0.;
  for (j = 0; j < N; j++)
    for (i = 0; i < N; i++) C[active[i] + j*Nf] = Cr[i + j*N];
  free(Cr);
  free(H);
  free(S);
  return (info != 0);
}

// index of the first electronic eigenvalue, E > -c^2 (nE if there is none)
uint32_t DiracElectronicFirst(const double *E, uint32_t nE, double c)
{
  uint32_t i = 0;
  while (i < nE && !(E[i] > -c*c)) i++;
  return i;
}

/* Root selection (design section 5): subshell (n, kappa) is the electronic eigenvalue
 * with index n - l(kappa) - 1, counted from DiracElectronicFirst. Returns 1 for n <= l. */
int DiracRootIndex(uint32_t n, int kappa, uint32_t *index)
{
  int l = DiracKappaL(kappa);
  if (kappa == 0 || (int) n <= l) return 1;
  *index = n - (uint32_t) l - 1;
  return 0;
}

/* P and Q of the spinor with DKB coefficients coef (2n: a_i for u_i^(1), then b_i for
 * u_i^(2)) at the points r[0 .. nr-1] (r > 0); 0 outside the knots:
 *     P = sum_i a_i B_i + h b_i (B_i' - kappa B_i/r),   Q = sum_i h a_i (B_i' + kappa B_i/r) + b_i B_i.
 * Either of P and Q may be NULL. */
int DiracComponents(const RadialBasisType *basis, int kappa, double c, const double *coef, const double *r, uint32_t nr,
                    double *P, double *Q)
{
  uint32_t n = basis->nSpline - 2, d = basis->degree, p, mu, j;
  int64_t i;
  double val[2*(BSplineMaxDegree+1)], h = 1./(2.*c), a, b, B, dB;
  for (p = 0; p < nr; p++)
  {
    if (P != NULL) P[p] = 0.;
    if (Q != NULL) Q[p] = 0.;
    mu = BSplineFindInterval(r[p], basis->knot, basis->nKnot);
    if (mu >= basis->nKnot || !(r[p] > 0.)) continue;
    BSplineValues(r[p], mu, basis->knot, basis->nKnot, d, 1, val);
    for (j = 0; j <= d; j++)
    {
      i = (int64_t) mu - d + j - 1;        // index in the working set B_1 .. B_n
      if (i < 0 || i >= n) continue;
      a = coef[i];
      b = coef[n + i];
      B = val[j];
      dB = val[d + 1 + j];
      if (P != NULL) P[p] += a*B + h*b*(dB - kappa*B/r[p]);
      if (Q != NULL) Q[p] += h*a*(dB + kappa*B/r[p]) + b*B;
    }
  }
  return 0;
}

/* The sign of P on its innermost lobe, P given on nP points of increasing r: the first local
 * extremum of |P| with |P| > DBL_EPSILON^(1/3) max |P| (6e-6). The floor skips what is not a
 * lobe near r = 0: rounding noise of the coefficients for large l (below sqrt(DBL_EPSILON)
 * max |P|), and a tiny linear term that the DKB basis leaves in P of kappa = +1 near the
 * origin, 1.6e-7 max |P| with a uniform-sphere nucleus (with which the sign of Re 2p1/2
 * came out opposite to FPLO's; A5e). Genuine innermost lobes are 1e-3 max |P| or more.
 * Returns +1, -1, or 0 if P vanishes. Used by DiracSignConvention and AtomSignConvention. */
double InnermostLobeSign(const double *P, uint32_t nP)
{
  uint32_t q;
  double Pmax = 0., floor;
  for (q = 0; q < nP; q++) if (fabs(P[q]) > Pmax) Pmax = fabs(P[q]);
  floor = cbrt(DBL_EPSILON)*Pmax;
  for (q = 0; q + 1 < nP; q++)
    if (fabs(P[q]) > floor && fabs(P[q]) >= fabs(P[q+1])) return (P[q] > 0.) ? 1. : -1.;
  return 0.;
}

/* Sign convention: P > 0 on its innermost lobe (InnermostLobeSign), as in Grant's book and in
 * FPLO's basis files (checked against FPLO's Re functions, A5e). Flips the sign of coef if
 * needed. */
int DiracSignConvention(const RadialBasisType *basis, int kappa, double c, double *coef)
{
  uint32_t q, nP = basis->nPoint, n2 = 2*(basis->nSpline - 2);
  double *P = (double *) malloc(nP*sizeof(double)), sign;
  DiracComponents(basis, kappa, c, coef, basis->rq, nP, P, NULL);
  sign = InnermostLobeSign(P, nP);
  free(P);
  if (sign < 0.) for (q = 0; q < n2; q++) coef[q] = -coef[q];
  return (sign != 0.) ? 0 : 1;
}

/* Lorentz-scalar confinement (design section 6) on the quadrature points:
 *     V_beta(r) = (r/rc)^N,
 * rc the confinement radius (P_nl P_at rho in the design, combined by the caller). Enters
 * the Dirac equation as beta V_beta, a position-dependent mass c^2 + V_beta. Returns 1 for
 * rc <= 0 or N <= 0. */
int ConfinementPotential(const RadialBasisType *basis, DFTReal rc, DFTReal N, DFTReal *Vbeta)
{
  uint32_t q;
  if (!(rc > 0.) || !(N > 0.)) return 1;
  for (q = 0; q < basis->nPoint; q++) Vbeta[q] = pow(basis->rq[q]/rc, N);
  return 0;
}

/* The end of the knots for a confinement (r/rc)^N: where it reaches Vmax,
 *     rMax = rc Vmax^(1/N).
 * The largest eigenvalue of the DKB problem, and with it the rounding error of every
 * eigenvalue, grows with V_beta(rMax); Vmax = 1e5 Ha gives rMax = 4.2 rc for N = 8,
 * enough for the confined orbitals to vanish, and max |E| ~ 1e10 instead of 1e15 Ha for
 * rMax = 40, rc = 0.5. */
DFTReal ConfinementRMax(DFTReal rc, DFTReal N, DFTReal Vmax)
{
  return rc*pow(Vmax, 1./N);
}

/* int f (P^2 - Q^2) dr of the spinor with DKB coefficients coef, f on the quadrature
 * points; with f = V_beta the confinement energy <beta V_beta> (design section 6). With
 * sign = +1 instead the expectation value int f (P^2 + Q^2) dr. Compensated summation. */
DFTReal DiracExpectation(const RadialBasisType *basis, int kappa, double c, const double *coef, const DFTReal *f, int sign)
{
  uint32_t q;
  double *P = (double *) malloc(basis->nPoint*sizeof(double)), *Q = (double *) malloc(basis->nPoint*sizeof(double));
  DFTReal sum = 0., comp = 0.;
  DiracComponents(basis, kappa, c, coef, basis->rq, basis->nPoint, P, Q);
  for (q = 0; q < basis->nPoint; q++) NeumaierAdd(&sum, &comp, basis->wq[q]*f[q]*(P[q]*P[q] + sign*Q[q]*Q[q]));
  free(P);
  free(Q);
  return sum + comp;
}

// the DKB overlap matrix S (2n x 2n, dense, column major) of DiracDKBMatrices; it does not depend on V
int DiracDKBOverlap(const RadialBasisType *basis, int kappa, double c, double *S)
{
  uint32_t N = 2*(basis->nSpline - 2), q;
  double *H = (double *) malloc(N*N*sizeof(double)), *V = (double *) malloc(basis->nPoint*sizeof(double));
  int error;
  for (q = 0; q < basis->nPoint; q++) V[q] = 0.;
  error = DiracDKBMatrices(basis, kappa, c, V, NULL, H, S);
  free(H);
  free(V);
  return error;
}

/* On-site orthonormalization (design section 6.1) of m orbitals of one kappa on one knot
 * set: the columns of C (N x m, column major, the DKB coefficient vectors), S the DKB
 * overlap matrix (N x N). With the overlap O = C^T S C (m x m):
 *     DiracGramSchmidt  O = L L^T (Cholesky), C <- C L^-T: ordered Gram-Schmidt in the
 *                       order of the columns; the first column is only normalized, every
 *                       other one changes least with respect to those before it;
 *     DiracLowdin       C <- C O^(-1/2) (symmetric).
 * The step is exact in exact arithmetic, but in floating point it leaves max |C^T S C - 1|
 * of about DBL_EPSILON cond(O): 3e-10 for two 3d orbitals of Z = 25 and 26 (smallest
 * eigenvalue of O 2.8e-7). It is therefore done twice ("twice is enough" as long as
 * DBL_EPSILON cond(O) < 1, which the redundancy check guarantees): the second pass, on an
 * O = 1 + O(1e-10), brings the error to 1e-16 and keeps the Gram-Schmidt order (the
 * first orbital stays the normalized input) and the symmetry of Lowdin (decided by
 * Maurits, 2026-11-11; the design said one pass suffices).
 * *minEigenvalue receives the smallest eigenvalue of the O of the input. If it is below
 * redundancy, the orbitals are (nearly) linearly dependent; C is left unchanged and 1 is
 * returned. Returns 1 also if a LAPACK routine fails or the method is unknown. */
int DiracOrthonormalize(const double *S, uint32_t N, double *C, uint32_t m, int method, double redundancy, double *minEigenvalue)
{
  double second;
  if (DiracOrthonormalizePass(S, N, C, m, method, redundancy, minEigenvalue)) return 1;
  return DiracOrthonormalizePass(S, N, C, m, method, 0., &second);
}

// one pass of DiracOrthonormalize, with the same arguments
int DiracOrthonormalizePass(const double *S, uint32_t N, double *C, uint32_t m, int method, double redundancy, double *minEigenvalue)
{
  double *SC, *O, *U, *w, *work, *T, one = 1., zero = 0.;
  uint32_t i, j, l, lwork = 3*m + 1 + m*m;
  int info = 0, error = 0;
  char uplo = 'L', jobz = 'V', side = 'R', trans = 'T', diag = 'N', nt = 'N';
  if (m == 0) return 1;
  SC = (double *) malloc(N*m*sizeof(double));
  O = (double *) malloc(m*m*sizeof(double));
  U = (double *) malloc(m*m*sizeof(double));
  w = (double *) malloc(m*sizeof(double));
  work = (double *) malloc(lwork*sizeof(double));
  // O = C^T S C
  dgemm(&nt, &nt, &N, &m, &N, &one, (double *) S, &N, C, &N, &zero, SC, &N);
  dgemm(&trans, &nt, &m, &m, &N, &one, C, &N, SC, &N, &zero, O, &m);
  for (j = 0; j < m; j++)
    for (i = j+1; i < m; i++) O[i + j*m] = O[j + i*m] = 0.5*(O[i + j*m] + O[j + i*m]);
  for (i = 0; i < m*m; i++) U[i] = O[i];
  dsyev(&jobz, &uplo, &m, U, &m, w, work, &lwork, &info);
  *minEigenvalue = (info == 0) ? w[0] : NAN;
  if (info != 0 || !(w[0] >= redundancy)) error = 1;
  else if (method == DiracGramSchmidt)
  {
    dpotrf(&uplo, &m, O, &m, &info);
    if (info != 0) error = 1;
    else dtrsm(&side, &uplo, &trans, &diag, &N, &m, &one, O, &m, C, &N);
  }
  else if (method == DiracLowdin)
  {
    // O^(-1/2) = U diag(w^(-1/2)) U^T
    T = (double *) malloc(m*m*sizeof(double));
    for (j = 0; j < m; j++)
      for (i = 0; i < m; i++)
      {
        double sum = 0.;
        for (l = 0; l < m; l++) sum += U[i + l*m]*U[j + l*m]/sqrt(w[l]);
        T[i + j*m] = sum;
      }
    for (i = 0; i < N*m; i++) SC[i] = C[i];
    dgemm(&nt, &nt, &N, &m, &m, &one, SC, &N, T, &m, &zero, C, &N);
    free(T);
  }
  else error = 1;
  free(SC);
  free(O);
  free(U);
  free(w);
  free(work);
  return error;
}

// max_ij |(C^T S C - 1)_ij| for the m columns of C (N x m), S (N x N)
double DiracOrthonormalityError(const double *S, uint32_t N, const double *C, uint32_t m)
{
  double *SC = (double *) malloc(N*m*sizeof(double)), *O = (double *) malloc(m*m*sizeof(double)), one = 1., zero = 0., err = 0.;
  uint32_t i, j;
  char nt = 'N', trans = 'T';
  dgemm(&nt, &nt, &N, &m, &N, &one, (double *) S, &N, (double *) C, &N, &zero, SC, &N);
  dgemm(&trans, &nt, &m, &m, &N, &one, (double *) C, &N, SC, &N, &zero, O, &m);
  for (j = 0; j < m; j++)
    for (i = 0; i < m; i++) err = fmax(err, fabs(O[i + j*m] - ((i == j) ? 1. : 0.)));
  free(SC);
  free(O);
  return err;
}
