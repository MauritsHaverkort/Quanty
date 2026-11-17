/*
 * DKBSpline.cpp
 *
 * Radial Dirac spinor components in the dual-kinetic-balance B-spline basis, see
 * DKBSpline.h (DFT stage A, design sections 5 and 10).
 */

#include "DKBSpline.h"

/* The k-th derivative at t of one component of a DKB spinor (see DKBSpline.h),
 *     f = u + h (v' + s v/r),   h = 1/(2c),
 * with (u, v, s) = (S_a, S_b, -kappa) for the large component P (component 0) and
 * (S_b, S_a, kappa) for the small one Q (component 1); a and b are the coefficients of
 * S_a and S_b over all nKnot-degree-1 B-splines. Then
 *     f^(k) = u^(k) + h ( v^(k+1) + s (v/r)^(k) ).
 * The derivatives of the splines come from BSplineDerivative (coefficient differencing).
 * v(0) = 0, because the first B-spline, the only one nonzero at 0, has coefficient 0, so
 *     (v/r)^(k)(t) = int_0^1 s^k v^(k+1)(t s) ds = t^(-k-1) int_0^t x^k v^(k+1)(x) dx.
 * For k >= 1 this integral is used: its integrand is a polynomial of degree degree-1 on
 * every knot interval, which ceil(degree/2) Gauss-Legendre points integrate exactly, and
 * it has no cancellation. (The Leibniz sum sum_j k!/j! v^(j) (-1)^(k-j) / t^(k-j+1) has
 * terms of order v/t^(k+1) and loses up to 4 digits for k = 2 at t = 1e-6.) For k = 0,
 * v/t is used directly; at t = 0 the limit v'(0).
 * 0 outside [knot[0], knot[nKnot-1]], value and derivatives. The knots must start with
 * degree+1 knots at 0. */
double DKBSplineEvaluate(double t, const double *knot, uint32_t nKnot, uint32_t degree, int kappa, uint32_t component, double c,
                         const double *a, const double *b, uint32_t k)
{
  const double *u = (component == 0) ? a : b, *v = (component == 0) ? b : a;
  double s = (component == 0) ? -kappa : kappa, h = 1./(2.*c), vr = 0., lo, hi;
  double xg[BSplineMaxDegree/2+2], wg[BSplineMaxDegree/2+2];
  uint32_t nG = degree/2 + 1, iv, q;
  if (!(t >= knot[0]) || !(t <= knot[nKnot-1])) return 0.;
  if (k == 0)
  {
    vr = (t > 0.) ? BSplineDerivative(t, knot, nKnot, degree, v, 0)/t : BSplineDerivative(0., knot, nKnot, degree, v, 1);
  }
  else if (t > 0.)
  {
    for (iv = 0; iv+1 < nKnot && knot[iv] < t; iv++)
    {
      lo = knot[iv];
      hi = (knot[iv+1] < t) ? knot[iv+1] : t;
      if (!(lo < hi)) continue;
      GaussLegendreNodesWeights(nG, lo, hi, xg, wg);
      for (q = 0; q < nG; q++) vr += wg[q]*pow(xg[q], (double) k)*BSplineDerivative(xg[q], knot, nKnot, degree, v, k+1);
    }
    vr /= pow(t, (double) (k+1));
  }
  else
  {
    vr = BSplineDerivative(0., knot, nKnot, degree, v, k+1)/(k+1);
  }
  return BSplineDerivative(t, knot, nKnot, degree, u, k) + h*(BSplineDerivative(t, knot, nKnot, degree, v, k+1) + s*vr);
}

namespace Interpolation {

DKBSpline::DKBSpline() {
	type = "DKBSpline";
	degree = 0;
	kappa = 0;
	component = 0;
	speedOfLight = 0.;
	nQ = 0;
}

DKBSpline::~DKBSpline() {
}

/* The component (0 = P, 1 = Q) of the spinor with DKB coefficients coef[0 .. 2n-1]
 * (a_i for u_i^(1), then b_i for u_i^(2), over the working B-splines B_1 .. B_n,
 * n = nKnot-degree-3), with nQ Gauss points per knot interval for x and Integrate.
 * Returns 0, or 1 if the knots do not start with degree+1 knots at 0. */
int DKBSpline::Create(const double *knot, uint32_t nKnot, uint32_t d, int k, uint32_t comp, double c, const double *coef, uint32_t nQuad)
{
	uint32_t nSpline = nKnot - d - 1, n = nSpline - 2, nPoint, i;
	double *rq, *wq;
	if (nKnot < d + 4 || knot[0] != 0. || knot[d] != 0.) return 1;
	knots = std::vector<double>(knot, knot + nKnot);
	degree = d;
	kappa = k;
	component = comp;
	speedOfLight = c;
	nQ = nQuad;
	a = std::vector<double>(nSpline, 0.);
	b = std::vector<double>(nSpline, 0.);
	for (i = 0; i < n; i++)
	{
		a[i+1] = coef[i];
		b[i+1] = coef[n+i];
	}
	if (BSplineQuadratureSetup(knot, nKnot, nQ, &rq, &wq, &nPoint)) return 1;
	x = std::vector<double>(rq, rq + nPoint);
	free(rq);
	free(wq);
	y = (*this)(x);
	type = "DKBSpline";
	return 0;
}

double DKBSpline::operator()(double t)
{
	return DKBSplineEvaluate(t, knots.data(), (uint32_t) knots.size(), degree, kappa, component, speedOfLight, a.data(), b.data(), 0);
}

double DKBSpline::operator()(double t, unsigned int k)
{
	return DKBSplineEvaluate(t, knots.data(), (uint32_t) knots.size(), degree, kappa, component, speedOfLight, a.data(), b.data(), k);
}

std::vector<double> DKBSpline::operator()(const std::vector<double> &t)
{
	std::vector<double> val = std::vector<double>(t.size());
	for (uint32_t j = 0; j < (uint32_t) t.size(); j++) val[j] = (*this)(t[j]);
	return val;
}

std::vector<double> DKBSpline::operator()(const std::vector<double> &t, unsigned int k)
{
	std::vector<double> val = std::vector<double>(t.size());
	for (uint32_t j = 0; j < (uint32_t) t.size(); j++) val[j] = (*this)(t[j], k);
	return val;
}

/* int_lo^hi f dr with nQ Gauss-Legendre points on every knot interval inside [lo, hi]
 * (clipped at lo and hi); 0 outside the knots. On the first interval f is a polynomial
 * and the sum is exact for nQ >= degree/2+1; elsewhere f contains B/r and the sum
 * converges exponentially with nQ. */
double DKBSplineIntegrate(DKBSpline &f, double lo, double hi)
{
	double sum = 0., sign = 1., l, r, *xg, *wg;
	uint32_t iv, q, nKnot = (uint32_t) f.knots.size();
	if (hi < lo) { l = lo; lo = hi; hi = l; sign = -1.; }
	xg = (double *) malloc(f.nQ*sizeof(double));
	wg = (double *) malloc(f.nQ*sizeof(double));
	for (iv = 0; iv+1 < nKnot; iv++)
	{
		l = (f.knots[iv] > lo) ? f.knots[iv] : lo;
		r = (f.knots[iv+1] < hi) ? f.knots[iv+1] : hi;
		if (!(l < r)) continue;
		GaussLegendreNodesWeights(f.nQ, l, r, xg, wg);
		for (q = 0; q < f.nQ; q++) sum += wg[q]*f(xg[q]);
	}
	free(xg);
	free(wg);
	return sign*sum;
}

} /* namespace Interpolation */
