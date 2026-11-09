/*
 * GaussLegendre.cpp
 *
 * Gauss-Legendre quadrature, see GaussLegendre.h.
 */

/*  returns a pointer to an array of length 2*n contaning the n
 *  Quadrature nodes and weights
 *  GaussLegendreQuadratureCompact defined in BasicMath_StandardFunctions.h
 */
const double *GaussLegendreQuadrature(int n)
{
  return &(GaussLegendreQuadratureCompact[n*(n-1)]);
}
const double *GaussLegendreQuadratureNode(int n, int i)
{
 return &(GaussLegendreQuadratureCompact[n*(n-1)+(2*i)]);
}
const double *GaussLegendreQuadratureWeight(int n, int i)
{
 return &(GaussLegendreQuadratureCompact[n*(n-1)+(2*i+1)]);
}

/* Gauss-Legendre quadrature with N points on [a,b]: nodes x[N] in ascending order and
 * weights w[N], allocated by the caller. For N <= 127 they come from the 34-digit table
 * GaussLegendreQuadratureCompact, above that from the Golub-Welsch eigenvalue problem.
 * Returns 0, or 1 for N = 0. */
int GaussLegendreNodesWeights(uint32_t N, double a, double b, double *x, double *w)
{
  uint32_t i;
  const double *q;
  if (N == 0) return 1;
  if (N > 127) return GaussLegendreNodesWeightsGolubWelsch(N, a, b, x, w);
  q = GaussLegendreQuadrature((int)N);
  for (i = 0; i < N; i++)
  {
    x[i] = (b*(q[2*i]+1.) + a*(1.-q[2*i]))/2.;
    w[i] = q[2*i+1]*(b-a)/2.;
  }
  return 0;
}

/* Gauss-Legendre quadrature with N points on [a,b] from the eigenvalues of the Jacobi
 * matrix of the Legendre polynomials (Golub and Welsch), and the weights from
 * P_{N+1}. Used for N > 127, where there is no table. Compared with the table for
 * N <= 127 the nodes agree to 2e-15, the weights only to 4e-10 relative (the smallest
 * weights, near the ends of the interval). Newton steps on the nodes would make the
 * weights exact to rounding, if a quadrature this large is ever needed.
 * Martin Brass (2019), rewritten with plain arrays. Returns 0, 1 for N = 0, 2 if LAPACK
 * fails. */
int GaussLegendreNodesWeightsGolubWelsch(uint32_t N, double a, double b, double *x, double *w)
{
  uint32_t j, one = 1;
  char compz = 'N';
  int info = 0;
  double p, *subdiag, *dummy = NULL;
  if (N == 0) return 1;
  subdiag = (double*) malloc(N*sizeof(double));
  for (j = 0; j < N; j++) x[j] = 0.;
  for (j = 1; j < N; j++) subdiag[j-1] = j/sqrt(4.*j*j - 1.);
  dsteqr(&compz, &N, x, subdiag, dummy, &one, dummy, &info);
  free(subdiag);
  if (info != 0) return 2;
  for (j = 0; j < N; j++)
  {
    p = Quad::LegendreP(x[j], N+1);
    w[j] = 2.*(1. - x[j]*x[j])/(N+1)/(N+1)/p/p*(b-a)/2.;
    x[j] = (b*(x[j]+1.) + a*(1.-x[j]))/2.;
  }
  return 0;
}

// namespace Quad: Martin Brass, May 2019 (Martin/Math/Quad.cpp)
namespace Quad
{

// N points on [a,b], see GaussLegendreNodesWeights
void GaussianQuadratureWeights(uint32_t N, vector<double> &x, vector<double> &w, double a, double b)
{
	x = vector<double>(N);
	w = vector<double>(N);
	if (N > 0) GaussLegendreNodesWeights(N, a, b, &x[0], &w[0]);
}
void GaussianQuadratureWeights(uint32_t N, const vector<double> &x0, vector<double> &x, vector<double> &w)
{
	uint32_t n = (uint32_t) x0.size() -1;
	x = vector<double>(n*N);
	w = vector<double>(n*N);
	vector<double> x1, w1;
	for (uint32_t i = 0; i < n; i++)
	{
		GaussianQuadratureWeights(N,x1,w1,x0[i],x0[i+1]);
		for (uint32_t j = 0; j < N; j++)
		{
			x[i*N + j] = x1[j];
			w[i*N + j] = w1[j];
		}
	}
}

double LegendreP(double x, uint32_t n)			// Recursion formula for Legendre polynomials
{
	if (n == 0) return 0;
	double p0 = 1;
	double p1 = x;
	double p = 1;

	for (uint32_t i = 2; i <= n; i++)
	{
		p = ((2 * i - 1) * x * p1 - (i - 1) * p0 ) / i;
		p0 = p1;
		p1 = p;
	}
	return p;
}

template <typename T>
double Integrate(T &f, double a, double b, uint32_t N)
{
	double val = 0;
	vector<double> x,w;
	GaussianQuadratureWeights(N,x,w,a,b);
	for (uint32_t j = 0; j < N; j++)
	{
		val += f(x[j])*w[j];
	}
	return val;
}

} // namespace Quad


