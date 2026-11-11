/*
 * FourierBessel.cpp
 *
 * The Fourier-Bessel transforms of the radial functions (DFT stage A5d, design section 10,
 * item 3; prepared for stage C):
 *     P~(k) = sqrt(2/pi) int r P(r) j_l(kr) dr,     Q~(k) = sqrt(2/pi) int r Q(r) j_lbar(kr) dr,
 * the radial part of the momentum-space spinor of psi = (P Omega_kappa,mu, i Q Omega_-kappa,mu)/r,
 * normalized as int k^2 (P~^2 + Q~^2) dk = int (P^2 + Q^2) dr. The integral runs over the
 * knot intervals with Gauss-Legendre points, n_q = (degree + 1) + ceil(k_max dr/pi) per
 * interval (design 10.3: the polynomial order plus the number of half periods of j_l),
 * split into pieces of at most 127 points (the 34-digit table); P and Q are evaluated
 * exactly (DKBSplineEvaluate, BSplineDerivative). j_l: SphericalBesselJ (the recurrences).
 * Returns 0 on success, 1 on error; does not print.
 */

#include "FourierBessel.h"

/* P~ and Q~ at the nk points k[] of the orbital with coefficients coef (Dirac: DKB, kappa;
 * Schroedinger: B-splines, l, Q~ = 0). */
int FourierBesselTransform(const RadialBasisType *b, int relativistic, int kappa, int l, double c, const double *coef,
                           uint32_t nk, const double *k, double *Pk, double *Qk)
{
  uint32_t i, j, s, nPiece, n, nq, ns = b->nSpline, nw = ns - 2, iq;
  int lP = relativistic ? DiracKappaL(kappa) : l, lQ = relativistic ? DiracKappaL(-kappa) : 0;
  double kMax = 0., *A = (double *) calloc(ns, sizeof(double)), *B = (double *) calloc(ns, sizeof(double));
  double x[DiracMaxQuadraturePoints], w[DiracMaxQuadraturePoints], pref = sqrt(2./M_PI);
  DFTReal *sP = (DFTReal *) calloc(nk, sizeof(DFTReal)), *cP = (DFTReal *) calloc(nk, sizeof(DFTReal));
  DFTReal *sQ = (DFTReal *) calloc(nk, sizeof(DFTReal)), *cQ = (DFTReal *) calloc(nk, sizeof(DFTReal));
  int error = 0;
  for (i = 0; i < nk; i++) if (k[i] > kMax) kMax = k[i];
  for (i = 0; i < nw; i++) { A[i+1] = coef[i]; if (relativistic) B[i+1] = coef[nw + i]; }
  for (iq = 0; iq + 1 < b->nKnot && !error; iq++)
  {
    double t0 = b->knot[iq], t1 = b->knot[iq+1];
    if (!(t1 > t0)) continue;
    nq = (b->degree + 1) + (uint32_t) ceil(kMax*(t1 - t0)/M_PI);
    nPiece = (nq + DiracMaxQuadraturePoints - 1)/DiracMaxQuadraturePoints;
    n = (nq + nPiece - 1)/nPiece;
    for (s = 0; s < nPiece && !error; s++)
    {
      double a = t0 + (t1 - t0)*s/nPiece, e = t0 + (t1 - t0)*(s + 1)/nPiece;
      error = GaussLegendreNodesWeights(n, a, e, x, w);
      for (j = 0; j < n; j++)
      {
        double P, Q = 0., r = x[j];
        if (relativistic)
        {
          P = DKBSplineEvaluate(r, b->knot, b->nKnot, b->degree, kappa, 0, c, A, B, 0);
          Q = DKBSplineEvaluate(r, b->knot, b->nKnot, b->degree, kappa, 1, c, A, B, 0);
        }
        else P = BSplineDerivative(r, b->knot, b->nKnot, b->degree, A, 0);
        for (i = 0; i < nk; i++)
        {
          NeumaierAdd(sP + i, cP + i, w[j]*r*P*SphericalBesselJ(k[i]*r, (uint32_t) lP));
          if (relativistic) NeumaierAdd(sQ + i, cQ + i, w[j]*r*Q*SphericalBesselJ(k[i]*r, (uint32_t) lQ));
        }
      }
    }
  }
  for (i = 0; i < nk; i++)
  {
    Pk[i] = pref*(sP[i] + cP[i]);
    Qk[i] = relativistic ? pref*(sQ[i] + cQ[i]) : 0.;
  }
  free(A); free(B); free(sP); free(cP); free(sQ); free(cQ);
  return error;
}
