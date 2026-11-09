/*
 * BSplineMatrix.cpp
 *
 * Integrals over B-splines, see BSplineMatrix.h. Conventions as in BSpline.h; the cell
 * iv is the knot interval [knot[iv], knot[iv+1]], ks = degree + 1.
 */

#include "BSplineMatrix.h"

/* Gauss-Legendre points on every non-empty knot interval, nQ per interval, in the order
 * of the intervals. Allocates *rq and *wq (the caller frees them); *nPoint is their
 * length. nQ points integrate a polynomial of degree 2nQ-1 on each interval exactly, so
 * nQ = degree+1 integrates the product of two B-splines exactly. Returns 0, or 1 for
 * fewer than 2 knots or nQ = 0. */
int BSplineQuadratureSetup(const double *knot, uint32_t nKnot, uint32_t nQ, double **rq, double **wq, uint32_t *nPoint)
{
  uint32_t iv, nInterval = 0, p = 0;
  double *r, *w;
  if (nKnot < 2 || nQ == 0) return 1;
  for (iv = 0; iv+1 < nKnot; iv++) if (knot[iv] < knot[iv+1]) nInterval++;
  r = (double *) malloc((nInterval*nQ + 1)*sizeof(double));
  w = (double *) malloc((nInterval*nQ + 1)*sizeof(double));
  for (iv = 0; iv+1 < nKnot; iv++)
  {
    if (!(knot[iv] < knot[iv+1])) continue;
    GaussLegendreNodesWeights(nQ, knot[iv], knot[iv+1], r+p, w+p);
    p += nQ;
  }
  *rq = r;
  *wq = w;
  *nPoint = p;
  return 0;
}

/* The banded matrix
 *     M_ij = int B_i^(a)(r) f(r) B_j^(b)(r) dr  ~  sum_q wq[q] fq[q] B_i^(a)(rq[q]) B_j^(b)(rq[q])
 * for the basis functions first .. first+n-1 (matrix indices 0 .. n-1), from the
 * quadrature points rq, wq (BSplineQuadratureSetup) and the values fq of f on them
 * (fq = NULL for f = 1). a and b are derivative orders.
 * Only i-j = 0 .. degree can be nonzero. The lower band is stored as in LAPACK ('L'):
 *     band[(i-j) + j*ldBand] = M_ij,   i >= j,   ldBand >= degree+1.
 * For a = b the matrix is symmetric; for a != b, M_ji is the same call with a and b
 * swapped. Returns 0, or 1 for an invalid degree or ldBand. */
int BSplineMatrixBanded(const double *knot, uint32_t nKnot, uint32_t degree,
                        const double *rq, const double *wq, uint32_t nPoint, const double *fq,
                        uint32_t a, uint32_t b, uint32_t first, uint32_t n, double *band, uint32_t ldBand)
{
  double val[(BSplineMaxDegree+1)*(BSplineMaxDegree+1)], wf;
  uint32_t q, ii, jj, mu, nDeriv = (a > b) ? a : b;
  int64_t i, j;
  if (degree > BSplineMaxDegree || nDeriv > BSplineMaxDegree || ldBand < degree+1) return 1;
  for (j = 0; j < n; j++)
    for (ii = 0; ii <= degree; ii++) band[ii + j*ldBand] = 0.;
  for (q = 0; q < nPoint; q++)
  {
    mu = BSplineFindInterval(rq[q], knot, nKnot);
    if (mu >= nKnot) continue;
    BSplineValues(rq[q], mu, knot, nKnot, degree, nDeriv, val);
    wf = (fq != NULL) ? wq[q]*fq[q] : wq[q];
    for (jj = 0; jj <= degree; jj++)
    {
      j = (int64_t) mu - degree + jj - first;
      if (j < 0 || j >= n) continue;
      for (ii = jj; ii <= degree; ii++)
      {
        i = (int64_t) mu - degree + ii - first;
        if (i >= n) break;
        band[(ii-jj) + j*ldBand] += wf*val[a*(degree+1)+ii]*val[b*(degree+1)+jj];
      }
    }
  }
  return 0;
}

/* The basis functions on nQ Gauss points of every cell iv = degree .. nKnot-2:
 *     BSval[iv*ks*nQ + r*nQ + q] = B_{iv-degree+r}(x_q),   r = 0 .. degree,
 * x_q the q-th Gauss point of the cell. Empty cells and functions beyond the last are
 * 0. BSval has room for (nKnot-1)*ks*nQ values. (Martin's PutBSplinesOnQuadGrid.)
 * Returns 0, or 1 for an invalid degree. */
int BSplineCellValues(const double *knot, uint32_t nKnot, uint32_t degree, uint32_t nQ, double *BSval)
{
  double val[BSplineMaxDegree+1], *x, *w;
  uint32_t ks = degree+1, nv = nKnot-1, nB = nKnot-degree-1, iv, q, r;
  if (degree > BSplineMaxDegree || nKnot < degree+2) return 1;
  for (iv = 0; iv < nv*ks*nQ; iv++) BSval[iv] = 0.;
  x = (double *) malloc(nQ*sizeof(double));
  w = (double *) malloc(nQ*sizeof(double));
  for (iv = degree; iv < nv; iv++)
  {
    if (!(knot[iv] < knot[iv+1])) continue;
    GaussLegendreNodesWeights(nQ, knot[iv], knot[iv+1], x, w);
    for (q = 0; q < nQ; q++)
    {
      BSplineValues(x[q], iv, knot, nKnot, degree, 0, val);
      for (r = 0; r <= degree && iv-degree+r < nB; r++) BSval[iv*ks*nQ + r*nQ + q] = val[r];
    }
  }
  free(x);
  free(w);
  return 0;
}

/* Cell tables for the radial Slater integrals R^k (Martin Brass' cell method), from the
 * basis functions on the Gauss points of the cells (BSplineCellValues, same nQ). With
 * i, j, i2, j2 counted from the first basis function of cell iv (iv-degree):
 *     rk1[iv*ks^2 + i*ks + j] = int_cell r^k        B_i B_j dr
 *     rk2[iv*ks^2 + i*ks + j] = int_cell r^(-k-1)   B_i B_j dr
 *     rkD[iv*ks^4 + i1*ks^3 + j1*ks^2 + i2*ks + j2]
 *         = int_cell dr1 B_i1(r1) B_i2(r1) r1^(-k-1) int_{knot[iv]}^{r1} dr2 r2^k B_j1(r2) B_j2(r2)
 * the last with an inner Gauss rule of nQ points on [knot[iv], r1]. rk1 and rk2 have room
 * for (nKnot-1)*ks^2 values, rkD for (nKnot-1)*ks^4; any of them may be NULL. None of
 * the integrals is exact (r^(-k-1) is no polynomial); the error falls with nQ.
 * Returns 0, or 1 for an invalid degree. */
int BSplineSlaterCellTables(const double *knot, uint32_t nKnot, uint32_t degree, uint32_t nQ, const double *BSval,
                            int k, double *rk1, double *rk2, double *rkD)
{
  uint32_t ks = degree+1, ks2 = ks*ks, ks3 = ks2*ks, ks4 = ks3*ks, nv = nKnot-1, nB = nKnot-degree-1, i;
  int64_t iv;
  if (degree > BSplineMaxDegree || nKnot < degree+2) return 1;
  if (rk1 != NULL) for (i = 0; i < nv*ks2; i++) rk1[i] = 0.;
  if (rk2 != NULL) for (i = 0; i < nv*ks2; i++) rk2[i] = 0.;
  if (rkD != NULL) for (i = 0; i < nv*ks4; i++) rkD[i] = 0.;

#pragma omp parallel for
  for (iv = degree; iv < (int64_t) nv; iv++)
  {
    double val[BSplineMaxDegree+1], Itot, I2, s1, s2;
    double *u1, *w1, *u2, *w2, *BSval2;
    const double *B = BSval + iv*ks*nQ;
    uint32_t m, i1, j1, i2, j2, q, idx1, idx2, r;
    if (!(knot[iv] < knot[iv+1])) continue;
    m = (uint32_t)((iv < (int64_t) nB-1 ? iv : (int64_t) nB-1) - (iv - degree) + 1);   /* splines in this cell */
    u1 = (double *) malloc(nQ*sizeof(double));
    w1 = (double *) malloc(nQ*sizeof(double));
    u2 = (double *) malloc(nQ*sizeof(double));
    w2 = (double *) malloc(nQ*sizeof(double));
    BSval2 = (double *) malloc(nQ*ks*sizeof(double));
    GaussLegendreNodesWeights(nQ, knot[iv], knot[iv+1], u1, w1);

    if (rk1 != NULL || rk2 != NULL)
    {
      for (i1 = 0; i1 < m; i1++)
        for (j1 = i1; j1 < m; j1++)
        {
          s1 = 0.;
          s2 = 0.;
          for (q = 0; q < nQ; q++)
          {
            s1 += pow(u1[q], k)    * w1[q] * B[i1*nQ + q] * B[j1*nQ + q];
            s2 += pow(u1[q], -k-1) * w1[q] * B[i1*nQ + q] * B[j1*nQ + q];
          }
          if (rk1 != NULL) rk1[iv*ks2 + i1*ks + j1] = rk1[iv*ks2 + j1*ks + i1] = s1;
          if (rk2 != NULL) rk2[iv*ks2 + i1*ks + j1] = rk2[iv*ks2 + j1*ks + i1] = s2;
        }
    }

    if (rkD != NULL)
    {
      for (idx1 = 0; idx1 < nQ; idx1++)
      {
        GaussLegendreNodesWeights(nQ, knot[iv], u1[idx1], u2, w2);
        for (idx2 = 0; idx2 < nQ; idx2++)
        {
          BSplineValues(u2[idx2], (uint32_t) iv, knot, nKnot, degree, 0, val);
          for (r = 0; r < m; r++) BSval2[r*nQ + idx2] = val[r];
        }
        for (i1 = 0; i1 < m; i1++)
          for (j1 = 0; j1 < m; j1++)
            for (i2 = i1; i2 < m; i2++)
              for (j2 = j1; j2 < m; j2++)
              {
                I2 = 0.;
                for (idx2 = 0; idx2 < nQ; idx2++)
                  I2 += w2[idx2] * pow(u2[idx2], k) * BSval2[j1*nQ + idx2] * BSval2[j2*nQ + idx2];
                Itot = I2 * w1[idx1] * pow(u1[idx1], -k-1) * B[i1*nQ + idx1] * B[i2*nQ + idx1];
                rkD[iv*ks4 + i1*ks3 + j1*ks2 + i2*ks + j2] += Itot;
                if (i1 != i2) rkD[iv*ks4 + i2*ks3 + j1*ks2 + i1*ks + j2] += Itot;
                if (j1 != j2) rkD[iv*ks4 + i1*ks3 + j2*ks2 + i2*ks + j1] += Itot;
                if (i1 != i2 && j1 != j2) rkD[iv*ks4 + i2*ks3 + j2*ks2 + i1*ks + j1] += Itot;
              }
      }
    }
    free(u1);
    free(w1);
    free(u2);
    free(w2);
    free(BSval2);
  }
  return 0;
}

/* The radial Slater integral
 *     R^k = int int dr1 dr2  f_i1(r1) f_i2(r1)  r<^k / r>^(k+1)  f_j1(r2) f_j2(r2)
 * of four functions given by their B-spline coefficients ci1, cj1, ci2, cj2
 * (nKnot-degree-1 each), from the cell tables of BSplineSlaterCellTables for this k.
 * Cells iv < jv contribute rk1(iv) rk2(jv) and the reverse; equal cells contribute rkD.
 * (Martin's AssembleSlaterInts2.) */
double BSplineSlaterIntegral(const double *ci1, const double *cj1, const double *ci2, const double *cj2,
                             const double *rk1, const double *rk2, const double *rkD, uint32_t nKnot, uint32_t degree)
{
  double val = 0., I11, I21, I12, I22;
  uint32_t n = degree, ks = degree+1, ks2 = ks*ks, ks3 = ks2*ks, ks4 = ks3*ks;
  uint32_t nv = nKnot-1, Nspline = nv - n;
  uint32_t iv, jv, lowerIdx, lowerIdxj, iMax, jMax, i1, j1, i2, j2;

  for (iv = 0; iv < nv; iv++)
  {
    lowerIdx = iv > n ? iv-n : 0;
    iMax = lowerIdx + ks <= Nspline ? ks : Nspline - lowerIdx;
    for (i1 = 0; i1 < iMax; i1++)
      for (j1 = 0; j1 < iMax; j1++)
        for (i2 = 0; i2 < iMax; i2++)
          for (j2 = 0; j2 < iMax; j2++)
            val += (rkD[iv*ks4 + i1*ks3 + j1*ks2 + i2*ks + j2] + rkD[iv*ks4 + j1*ks3 + i1*ks2 + j2*ks + i2])
                   * ci1[i1+lowerIdx] * ci2[i2+lowerIdx] * cj1[j1+lowerIdx] * cj2[j2+lowerIdx];
  }

  for (iv = 0; iv < nv; iv++)
  {
    lowerIdx = iv > n ? iv-n : 0;
    iMax = lowerIdx + ks <= Nspline ? ks : Nspline - lowerIdx;
    I11 = 0.;
    I21 = 0.;
    for (i1 = 0; i1 < iMax; i1++)
      for (i2 = 0; i2 < iMax; i2++)
      {
        I11 += rk1[iv*ks2 + i1*ks + i2] * ci1[i1+lowerIdx] * ci2[i2+lowerIdx];
        I21 += rk1[iv*ks2 + i1*ks + i2] * cj1[i1+lowerIdx] * cj2[i2+lowerIdx];
      }
    I12 = 0.;
    I22 = 0.;
    for (jv = iv+1; jv < nv; jv++)
    {
      lowerIdxj = jv > n ? jv-n : 0;
      jMax = lowerIdxj + ks <= Nspline ? ks : Nspline - lowerIdxj;
      for (j1 = 0; j1 < jMax; j1++)
        for (j2 = 0; j2 < jMax; j2++)
        {
          I12 += rk2[jv*ks2 + j1*ks + j2] * cj1[j1+lowerIdxj] * cj2[j2+lowerIdxj];
          I22 += rk2[jv*ks2 + j1*ks + j2] * ci1[j1+lowerIdxj] * ci2[j2+lowerIdxj];
        }
    }
    val += I11*I12 + I21*I22;
  }
  return val;
}
