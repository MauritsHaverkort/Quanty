/*
 * BSpline.cpp
 *
 * B-splines in plain C, see BSpline.h for the conventions.
 *
 * Martin Brass wrote the first version (2019) in C++ with a recursive Cox-de Boor
 * evaluation, O(2^degree) per basis function. Rewritten in 2026 with plain arrays and
 * de Boor's algorithm, O(degree^2) per point for all degree+1 nonzero basis functions
 * and their derivatives. The recursive functions are kept below as EvaluateOld and
 * EvalDerivativeOld; they are compiled only when BSPLINE_COMPILE_OLD is defined, which
 * only the comparison against the new code does.
 */

#include "BSpline.h"
#include "TBSpline.cpp"

/* knot[i], with i clamped to [0, nKnot-1]. Beyond its ends the knot vector is
 * continued with copies of its end knots; this changes none of the existing basis
 * functions and lets BSplineValues work near the ends of knots that are not clamped. */
static double BSplineKnot(const double *knot, uint32_t nKnot, int64_t i)
{
  if (i < 0) return knot[0];
  if (i >= (int64_t) nKnot) return knot[nKnot-1];
  return knot[i];
}

/* The interval mu with knot[mu] <= t < knot[mu+1]. A point on an interior knot gets
 * the interval to its right, the last knot the last non-empty interval. Returns nKnot
 * when t lies outside [knot[0], knot[nKnot-1]], is NaN, or all knots coincide.
 * Binary search, O(log nKnot). */
uint32_t BSplineFindInterval(double t, const double *knot, uint32_t nKnot)
{
  uint32_t lo, hi, mid;
  if (nKnot < 2) return nKnot;
  if (!(t >= knot[0]) || !(t <= knot[nKnot-1]) || knot[0] == knot[nKnot-1]) return nKnot;
  if (t == knot[nKnot-1])
  {
    hi = nKnot-1;
    while (knot[hi-1] == knot[nKnot-1]) hi--;
    return hi-1;
  }
  lo = 0;
  hi = nKnot-1;                    /* knot[lo] <= t < knot[hi] */
  while (hi - lo > 1)
  {
    mid = lo + (hi - lo)/2;
    if (knot[mid] <= t) lo = mid;
    else hi = mid;
  }
  return lo;
}

/* The degree+1 basis functions that can be nonzero on the knot interval
 * [knot[interval], knot[interval+1]], B_{interval-degree} .. B_{interval}, and their
 * derivatives, at t:
 *     val[k*(degree+1) + j] = d^k/dt^k B_{interval-degree+j}(t),   k = 0 .. nDeriv.
 * The interval must be non-empty (from BSplineFindInterval); t is not checked to lie in
 * it. Derivatives above the degree are 0.
 * de Boor's triangular recurrence for the values, with the derivatives from differences
 * of its coefficients (Piegl and Tiller, The NURBS Book, algorithm A2.3). Only sums of
 * positive terms enter the values, and every denominator is a knot distance that
 * contains the interval, so it is stable for repeated knots too. O(degree^2) per point.
 * Returns 0, or 1 if the degree exceeds BSplineMaxDegree or the interval is not one. */
int BSplineValues(double t, uint32_t interval, const double *knot, uint32_t nKnot, uint32_t degree,
                  uint32_t nDeriv, double *val)
{
  double ndu[BSplineMaxDegree+1][BSplineMaxDegree+1];
  double a[2][BSplineMaxDegree+1];
  double left[BSplineMaxDegree+1], right[BSplineMaxDegree+1];
  double saved, temp, d;
  int p = (int) degree, n, j, r, k, s1, s2, rk, pk, j1, j2;
  int64_t mu = interval;
  if (degree > BSplineMaxDegree || interval+1 >= nKnot) return 1;
  n = (nDeriv < degree) ? (int) nDeriv : p;

  ndu[0][0] = 1.;
  for (j = 1; j <= p; j++)
  {
    left[j]  = t - BSplineKnot(knot, nKnot, mu+1-j);
    right[j] = BSplineKnot(knot, nKnot, mu+j) - t;
    saved = 0.;
    for (r = 0; r < j; r++)
    {
      ndu[j][r] = right[r+1] + left[j-r];          /* knot distance */
      temp = ndu[r][j-1]/ndu[j][r];
      ndu[r][j] = saved + right[r+1]*temp;          /* basis function of degree j */
      saved = left[j-r]*temp;
    }
    ndu[j][j] = saved;
  }
  for (j = 0; j <= p; j++) val[j] = ndu[j][p];

  for (r = 0; r <= p; r++)
  {
    s1 = 0;
    s2 = 1;
    a[0][0] = 1.;
    for (k = 1; k <= n; k++)
    {
      d = 0.;
      rk = r - k;
      pk = p - k;
      if (r >= k)
      {
        a[s2][0] = a[s1][0]/ndu[pk+1][rk];
        d = a[s2][0]*ndu[rk][pk];
      }
      j1 = (rk >= -1) ? 1 : -rk;
      j2 = (r-1 <= pk) ? k-1 : p-r;
      for (j = j1; j <= j2; j++)
      {
        a[s2][j] = (a[s1][j] - a[s1][j-1])/ndu[pk+1][rk+j];
        d += a[s2][j]*ndu[rk+j][pk];
      }
      if (r <= pk)
      {
        a[s2][k] = -a[s1][k-1]/ndu[pk+1][r];
        d += a[s2][k]*ndu[r][pk];
      }
      val[k*(p+1)+r] = d;
      j = s1; s1 = s2; s2 = j;
    }
  }
  r = p;
  for (k = 1; k <= n; k++)
  {
    for (j = 0; j <= p; j++) val[k*(p+1)+j] *= r;
    r *= (p - k);
  }
  for (k = n+1; k <= (int) nDeriv; k++)
    for (j = 0; j <= p; j++) val[k*(p+1)+j] = 0.;
  return 0;
}

/* The k-th derivative of basis function B_i of the given degree at t; 0 outside its
 * support, outside the knots, and for k > degree. */
double BSplineEvaluateDerivative(uint32_t i, uint32_t degree, const double *knot, uint32_t nKnot, double t, uint32_t k)
{
  double val[(BSplineMaxDegree+1)*(BSplineMaxDegree+1)];
  uint32_t mu;
  if (k > degree || degree > BSplineMaxDegree || i + degree + 1 >= nKnot) return 0.;
  mu = BSplineFindInterval(t, knot, nKnot);
  if (mu >= nKnot || i > mu || i + degree < mu) return 0.;
  BSplineValues(t, mu, knot, nKnot, degree, k, val);
  return val[k*(degree+1) + i + degree - mu];
}

/* Basis function B_i of the given degree at t; 0 outside its support and outside the
 * knots. */
double BSplineEvaluate(uint32_t i, uint32_t degree, const double *knot, uint32_t nKnot, double t)
{
  return BSplineEvaluateDerivative(i, degree, knot, nKnot, t, 0);
}

/* The k-th derivative of the spline sum_i c[i] B_i(t), c[0 .. nKnot-degree-2]. 0 outside
 * the knots and for k > degree.
 * The k-th derivative is a spline of degree - k whose coefficients are differences of
 * the c (de Boor, A Practical Guide to Splines, bvalue):
 *     c'_i = p (c_i - c_{i-1}) / (t_{i+p} - t_i),   p the degree before the step.
 * Only the degree+1 coefficients that matter on the interval of t are differenced, and
 * the result is evaluated with the (positive) basis functions of degree - k. Summing
 * c_i B_i^(k)(t) instead would cancel terms of order (degree/h)^k to a result of order
 * 1 and lose up to 3 digits for the third derivative. */
double BSplineDerivative(double t, const double *knot, uint32_t nKnot, uint32_t degree, const double *c, uint32_t k)
{
  double a[BSplineMaxDegree+1], val[BSplineMaxDegree+1], sum = 0.;
  uint32_t mu, j, m, nB;
  int64_t i, p = degree, imu;
  if (k > degree || degree > BSplineMaxDegree || nKnot < degree + 2) return 0.;
  nB = nKnot - degree - 1;
  mu = BSplineFindInterval(t, knot, nKnot);
  if (mu >= nKnot) return 0.;
  imu = mu;
  for (j = 0; j <= degree; j++)
  {
    i = imu - p + j;
    a[j] = (i >= 0 && i < nB) ? c[i] : 0.;
  }
  for (m = 1; m <= k; m++)
    for (j = degree; j >= m; j--)
      a[j] = (p - m + 1)*(a[j] - a[j-1])
             /(BSplineKnot(knot, nKnot, imu + j - m + 1) - BSplineKnot(knot, nKnot, imu - p + j));
  BSplineValues(t, mu, knot, nKnot, degree - k, 0, val);
  for (j = k; j <= degree; j++) sum += a[j]*val[j-k];
  return sum;
}

/* The spline sum_i c[i] B_i(t) by de Boor's algorithm on knot interval idx, which must
 * satisfy degree <= idx and knot[idx] <= t <= knot[idx+1]. Prints an error and returns
 * NaN for idx < degree. */
double BSplineDeBoor(uint32_t idx, double t, const double *knot, uint32_t degree, const double *c)
{
  double d[BSplineMaxDegree+1], alpha;
  uint32_t i, j;
  if (idx < degree)
  {
    printf("Error deBoor: idx < n\n");
    return NAN;
  }
  if (degree > BSplineMaxDegree) return NAN;
  for (i = 0; i <= degree; i++) d[i] = c[i + idx - degree];
  for (i = 1; i <= degree; i++)
  {
    for (j = degree; j > i-1; j--)
    {
      alpha = (t - knot[j+idx-degree])/(knot[j+1+idx-i] - knot[j+idx-degree]);
      d[j] = (1. - alpha)*d[j-1] + alpha*d[j];
    }
  }
  return d[degree];
}

/* Knots for interpolation at the nt points t[0] < .. < t[nt-1] (knot averaging): the
 * end points degree+1 times, in between the averages of degree consecutive points.
 * knot has room for nt+degree+1 values. Returns 0, 1 for degree 0, 2 for
 * nt < degree+1. */
int BSplineKnotsAveraged(const double *t, uint32_t nt, uint32_t degree, double *knot)
{
  uint32_t i, j;
  if (degree == 0) return 1;
  if (nt < degree + 1) return 2;
  for (i = 0; i <= degree; i++)
  {
    knot[i] = t[0];
    knot[i + nt] = t[nt-1];
  }
  for (i = 0; i < nt - (degree+1); i++)
  {
    knot[i + degree+1] = 0.;
    for (j = i+1; j < i+1+degree; j++) knot[i + degree+1] += t[j];
    knot[i + degree+1] /= degree;
  }
  return 0;
}

/* Collocation matrix for interpolation at the nt points t: B[i + j*nt] = B_j(t[i]),
 * column major, nt x nt. Needs nKnot = nt + degree + 1. Returns 0, or 1 for a
 * mismatch. */
int BSplineCollocationMatrix(const double *knot, uint32_t nKnot, uint32_t degree, const double *t, uint32_t nt, double *B)
{
  double val[BSplineMaxDegree+1];
  uint32_t i, j, mu;
  int64_t col;
  if (nKnot != nt + degree + 1 || degree > BSplineMaxDegree) return 1;
  for (i = 0; i < nt*nt; i++) B[i] = 0.;
  for (i = 0; i < nt; i++)
  {
    mu = BSplineFindInterval(t[i], knot, nKnot);
    if (mu >= nKnot) continue;
    BSplineValues(t[i], mu, knot, nKnot, degree, 0, val);
    for (j = 0; j <= degree; j++)
    {
      col = (int64_t) mu - degree + j;
      if (col >= 0 && col < nt) B[i + col*nt] = val[j];
    }
  }
  return 0;
}

/* The collocation matrix of BSplineCollocationMatrix in LAPACK general band storage, as
 * dgbsv and dgbtrf want it: B_j(t_i) at AB[(kl+ku+i-j) + j*ldab], ldab = 2*kl+ku+1, the
 * first kl rows left free for the LU factors. kl and ku, the numbers of sub- and
 * superdiagonals, follow from the points (degree or less for knot averaging). Solving
 * with it costs O(nt degree^2) instead of O(nt^3) for the dense matrix. Allocates *AB
 * (the caller frees it). Returns 0, or 1 for a mismatch of nKnot and nt. */
int BSplineCollocationBand(const double *knot, uint32_t nKnot, uint32_t degree, const double *t, uint32_t nt,
                           uint32_t *kl, uint32_t *ku, uint32_t *ldab, double **AB)
{
  double val[BSplineMaxDegree+1], *ab;
  uint32_t i, j, mu, l = 0, u = 0, ld;
  int64_t col;
  if (nKnot != nt + degree + 1 || degree > BSplineMaxDegree) return 1;
  for (i = 0; i < nt; i++)
  {
    mu = BSplineFindInterval(t[i], knot, nKnot);
    if (mu >= nKnot) continue;
    for (j = 0; j <= degree; j++)
    {
      col = (int64_t) mu - degree + j;
      if (col < 0 || col >= nt) continue;
      if (col < i && i - col > l) l = (uint32_t)(i - col);
      if (col > i && col - i > u) u = (uint32_t)(col - i);
    }
  }
  ld = 2*l + u + 1;
  ab = (double *) calloc((size_t) ld*nt, sizeof(double));
  for (i = 0; i < nt; i++)
  {
    mu = BSplineFindInterval(t[i], knot, nKnot);
    if (mu >= nKnot) continue;
    BSplineValues(t[i], mu, knot, nKnot, degree, 0, val);
    for (j = 0; j <= degree; j++)
    {
      col = (int64_t) mu - degree + j;
      if (col >= 0 && col < nt) ab[(l + u + i - col) + col*ld] = val[j];
    }
  }
  *kl = l;
  *ku = u;
  *ldab = ld;
  *AB = ab;
  return 0;
}

/* Knots for radial functions after Froese Fischer: degree+1 knots at 0, then m+1 steps
 * of h = 2^-m, then a geometric grid r -> r(1+h) while h r < hMax, then steps of hMax up
 * to rMax, and the last knot degree+1 times; everything divided by Z. Allocates *knot
 * (the caller frees it). Returns 0. */
int BSplineKnotsFischer(int m, double hMax, double rMax, double Z, uint32_t degree, double **knot, uint32_t *nKnot)
{
  double h = pow(2., -m), last, *x;
  uint32_t n = 0, size = 64, i;
  int j;
  x = (double *) malloc(size*sizeof(double));
#define BSplineKnotsFischerAppend(value) \
  { if (n == size) { size *= 2; x = (double *) realloc(x, size*sizeof(double)); } x[n] = (value); n++; }
  for (i = 0; i <= degree; i++) BSplineKnotsFischerAppend(0.);
  for (j = 0; j <= m; j++) BSplineKnotsFischerAppend(x[n-1] + h);
  while (h*x[n-1] < hMax) BSplineKnotsFischerAppend(x[n-1]*(1. + h));
  while (x[n-1] < rMax) BSplineKnotsFischerAppend(x[n-1] + hMax);
  last = x[n-1];
  for (i = 1; i <= degree; i++) BSplineKnotsFischerAppend(last);
#undef BSplineKnotsFischerAppend
  for (i = 0; i < n; i++) x[i] /= Z;
  *knot = x;
  *nKnot = n;
  return 0;
}


/* Martin Brass' recursive Cox-de Boor evaluation (2019), replaced by BSplineEvaluate and
 * BSplineEvaluateDerivative. Not compiled into Quanty; kept for the record and for the
 * comparison of old and new code, which defines BSPLINE_COMPILE_OLD. */
#ifdef BSPLINE_COMPILE_OLD

double EvaluateOld(uint32_t i, uint32_t n, const vector<double> &x, double t)
{
	double s1 = 0.;
	double s2 = 0.;

	if (n==0)
	{
		if ( (x[i] <= t && t < x[i+1]) || (t == x[i+1] && t == x[x.size()-1] ) ) //TODO: Find a better way for this
		{
			return 1.;
		}
		else
		{
			return 0.;
		}
	}
	else
	{
		if (x[i+n] != x[i])
		{
			s1 = (t - x[i])/(x[i+n]-x[i]) * EvaluateOld(i, n-1,x,t);
		}
		if (x[i+n+1] != x[i+1])
		{
			s2 = (x[i+n+1] - t)/(x[i+n+1]-x[i+1]) * EvaluateOld(i+1,n-1,x,t);
		}
		return s1 + s2;
	}
}

double EvalDerivativeOld(uint32_t i, uint32_t n, const vector<double> &x,  double t, uint32_t k)
{
	if (k > 0)
	{
		double E1 =  EvalDerivativeOld(i,n-1,x,t,k-1);
		double E2 = EvalDerivativeOld(i+1,n-1,x,t,k-1);
		double dx1 = (x[i + n] -x[i]);
		double dx2 = (x[i + n +1] -x[i+1]);
		if (dx1 == 0)
		{
			if (dx2 == 0)
			{
				return 0;
			}
			else
			{
				return - (n / dx2 * E2);
			}
		}
		else
		{
			if (dx2 == 0)
			{
				return n / dx1 * E1;
			}
			else
			{
				return n / dx1 * E1 - (n / dx2 * E2);
			}
		}

	}
	else
	{
		return EvaluateOld(i,n,x,t);
	}
}

#endif /* BSPLINE_COMPILE_OLD */


/* Martin Brass' overlap matrix and projection onto B-splines (2019). They had no callers
 * when the B-spline code was rewritten; BSplineMatrixBanded replaces them. Kept for the
 * record, not compiled. */
#if 0

double OverlapOld(uint32_t i, uint32_t j, uint32_t ks, const vector<double> &x)
{
	uint32_t N = (uint32_t) x.size();
	uint32_t ia = std::max(ks,std::max(i,j));
	uint32_t ib = std::min(N - ks -1, std::min(i,j) + ks+1);

	double integral = 0;
	for (uint32_t iv = ia; iv < ib; iv++)
	{
		double a = x[iv  ];
		double b = x[iv+1];
		vector<double> u, w;
		Quad::GaussianQuadratureWeights(ks+1,u,w,a,b);
		for (uint32_t m = 0; m < ks+1; m++)
		{
			integral += EvaluateOld(i,ks,x,u[m])*EvaluateOld(j,ks,x,u[m])*w[m];
		}
	}
	return integral;
}

void OverlapMatrixOld(uint32_t ks, vector<double> x, vector<double> &S)
{
	uint32_t Nspline = (uint32_t)x.size() - ks - 1 -2;
	S = vector<double>(Nspline*(ks+1));
	//S = vector<double>(Nspline*Nspline);
/*	for (int k = 0; k < Nspline*Nspline; k++)
	{
		S[k] = 100;
	}*/
	for (uint32_t j = 0; j < Nspline; j++)
	{
		for ( uint32_t i = j; i <= std::min(Nspline-1,j+ks); i++)  // n = #sub-diagonals
		{
			//S[(i-j)*Nspline + j] = OverlapOld(i+1,j+1,ks,x);
			S[(i-j) + j*(ks+1)] = OverlapOld(i+1,j+1,ks,x);
		}
	}
}

void ComponentsOld(const vector<Interpolation::InterpolatingFunction*> &funs, uint32_t ks, const vector<double> &x, vector<double> &S, vector<double> &c)
{
	uint32_t N = (uint32_t) x.size();
	uint32_t Nfuns = (uint32_t) funs.size();
	uint32_t Nspline = N - ks-3;
	c = vector<double> (Nspline*Nfuns);
	for (uint32_t i = 0; i < Nspline; i++)
	{
		for (uint32_t j = 0; j < Nfuns; j++)
		{
			c[i + j*Nspline] = 0;
			for (uint32_t iv = std::max(ks,i+1); iv < std::min(N-ks-1,i+1+ks); iv++)
			{
				double a = x[iv];
				double b = x[iv+1];
				vector<double> u,w;
				Quad::GaussianQuadratureWeights(ks+1,u,w,a,b);
				for (uint32_t m = 0; m < ks+1; m++)
				{
					c[i + j*Nspline] += (*funs[j])(u[m]) * EvaluateOld(i,ks,x,u[m]) * w[m];
				}
			}
		}
	}
	printf("c = {");
	for (uint32_t i = 0; i < Nspline*Nfuns; i++)
	{
		printf("%.16f, ",c[i]);
	}
	printf("}\n");
  /*
	uint32_t eins = 1;
	double one = 1.;
	double null = 0.;
  */
	const char *l = "L";
	int info = 0;
	uint32_t LDA = ks+1;
/*	dpbsv(l,&Nspline,&ks,&Nfuns,&S[0],&Nspline,&c[0],&Nspline,&info);
	if (info)
	{
		printf("error in BSpline::Components: \n\tLAPACK dpbsv (positive banded symmetric linear equation) info = %d\n",info);
		fflush(stdout);
	}*/
	printf("???\n");
	dpbsv(l,&Nspline,&ks,&Nfuns,&S[0],&LDA,&c[0],&Nspline,&info);
	if (info)
	{
		printf("error in BSpline::Components: \n\tLAPACK dpbsv (positive banded symmetric linear equation) info = %d\n",info);
		fflush(stdout);
	}


}

#endif
