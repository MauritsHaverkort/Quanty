/*
 * TBSpline.cpp
 *
 *  Created on: Sep 10, 2019
 *      Author: brass
 */

#include "TBSpline.h"

namespace BSpline {

TBSpline::TBSpline() {
	type = "BSpline";
	degree = 0;
}

TBSpline::~TBSpline() {
	// TODO Auto-generated destructor stub
}

TBSpline::TBSpline(const std::vector<double> &u,const std::vector<double> &k,const std::vector<double> &v, uint32_t n)
{
	y = v;
	degree = n;
	knots = k;
	x = u;
	points = u;
	type = "BSpline";
}

void TBSpline::Create(const std::vector<double> &u, const std::vector<double> &v,uint32_t n)
{
	uint32_t N = (uint32_t) u.size();
	y = v;
	degree = n;
	points = u;
	knots = std::vector<double>(N + n + 1);
	int error = BSplineKnotsAveraged(&u[0], N, n, knots.data());
	if (error == 1)
	{
		printf("Error BSpline: Interpolation with n = 0 not implemented\n");
		return;
	}
	if (error == 2)
	{
		printf("Error BSpline: number of points = %u < %u = n+1\n", N, n+1);
		return;
	}
	uint32_t nx = (uint32_t) knots.size()-2*n;
	x = std::vector<double>(nx);
	for (uint32_t i = 0; i < nx; i++)
	{
		x[i] = knots[i+n];
	}
	// the coefficients solve the banded collocation system
	uint32_t kl, ku, ldab;
	double *AB;
	BSplineCollocationBand(knots.data(), (uint32_t) knots.size(), n, &u[0], N, &kl, &ku, &ldab, &AB);
	uint32_t NRHS = 1;
	vector<int> IPIV = vector<int>(N);
	int info = 0;
	dgbsv(&N,&kl,&ku,&NRHS,AB,&ldab,&(IPIV[0]),&(y[0]),&N,&info);
	free(AB);
	if (info)
	{
		printf("Error TBSpline.Create: dgbsv info = %i\n",info);
	}
	type = "BSpline";
}

// The spline at t; 0 outside the knots
double TBSpline::operator()(double t)
{
	return BSplineDerivative(t, knots.data(), (uint32_t) knots.size(), degree, y.data(), 0);
}

// The k-th derivative of the spline at t; 0 outside the knots
double TBSpline::operator()(double t, unsigned int k)
{
	if (k > degree)
	{
		printf("Warning: you just took a derivative that is higher than the degree of the interpolating function.\n This clearly is zero and I will return zero, but I'm wondering if this result ismeaningfull...\n");
		return 0;
	}
	return BSplineDerivative(t, knots.data(), (uint32_t) knots.size(), degree, y.data(), k);
}

std::vector<double> TBSpline::operator()(const std::vector<double> &t)
{
	std::vector<double> val = std::vector<double>(t.size());
	for (uint32_t j = 0; j < (uint32_t) t.size(); j++)
	{
		val[j] = BSplineDerivative(t[j], knots.data(), (uint32_t) knots.size(), degree, y.data(), 0);
	}
	return val;
}

std::vector<double> TBSpline::operator()(const std::vector<double> &t, unsigned int k)
{
	std::vector<double> val = std::vector<double>(t.size());
	if (k > degree)
	{
		printf("Warning: you just took a derivative that is higher than the degree of the interpolating function.\nThis clearly is zero and I will return zero, but I'm wondering if this result is meaningfull...\n");
		return val;
	}
	for (uint32_t j = 0; j < (uint32_t) t.size(); j++)
	{
		val[j] = BSplineDerivative(t[j], knots.data(), (uint32_t) knots.size(), degree, y.data(), k);
	}
	return val;
}

/* Arithmetic on B-splines. As for the other interpolating functions the result is a new
 * function of the type of the first, through op(f(p), g(p)) on the union of the points
 * of f and g: for a B-spline the points it interpolates (its data), for the other types
 * their grid. The result has the degree of f. -f, and f + g and f - g for two B-splines
 * of the same degree and knots, act on the coefficients and are exact everywhere. */

// the points on which arithmetic samples a function
static std::vector<double> SamplePoints(Interpolation::InterpolatingFunction &f)
{
	if (f.type == "BSpline" && ((TBSpline &) f).points.size() > 0) return ((TBSpline &) f).points;
	return f.GetX();
}

static void Arithmetic(TBSpline &f, Interpolation::InterpolatingFunction &g, char op, TBSpline &r)
{
	if ((op == '+' || op == '-') && g.type == "BSpline")
	{
		TBSpline &h = (TBSpline &) g;
		if (h.degree == f.degree && h.knots == f.knots && h.y.size() == f.y.size())
		{
			r.degree = f.degree;
			r.knots = f.knots;
			r.x = f.x;
			r.points = f.points;
			r.y = std::vector<double>(f.y.size());
			for (uint32_t i = 0; i < (uint32_t) f.y.size(); i++) r.y[i] = (op == '+') ? f.y[i] + h.y[i] : f.y[i] - h.y[i];
			r.type = "BSpline";
			return;
		}
	}
	std::vector<double> p, v;
	Interpolation::MergeGrids(SamplePoints(f), SamplePoints(g), p);
	v = std::vector<double>(p.size());
	for (uint32_t i = 0; i < (uint32_t) p.size(); i++)
	{
		double a = f(p[i]), b = g(p[i]);
		switch (op)
		{
			case '+': v[i] = a + b; break;
			case '-': v[i] = a - b; break;
			case '*': v[i] = a * b; break;
			default:  v[i] = a / b; break;
		}
	}
	r.Create(p, v, f.degree);
}

void Add(TBSpline &f, Interpolation::InterpolatingFunction &g, TBSpline &r) { Arithmetic(f, g, '+', r); }
void Sub(TBSpline &f, Interpolation::InterpolatingFunction &g, TBSpline &r) { Arithmetic(f, g, '-', r); }
void Mul(TBSpline &f, Interpolation::InterpolatingFunction &g, TBSpline &r) { Arithmetic(f, g, '*', r); }
void Div(TBSpline &f, Interpolation::InterpolatingFunction &g, TBSpline &r) { Arithmetic(f, g, '/', r); }

// -f: the coefficients change sign, exact
void UnaryMinus(TBSpline &f, TBSpline &g)
{
	g.degree = f.degree;
	g.knots = f.knots;
	g.x = f.x;
	g.points = f.points;
	g.y = std::vector<double>(f.y.size());
	for (uint32_t i = 0; i < (uint32_t) f.y.size(); i++) g.y[i] = -f.y[i];
	g.type = "BSpline";
}

// f^n through f(p)^n on the points of f
void Pow(TBSpline &f, double n, TBSpline &r)
{
	std::vector<double> p = SamplePoints(f), v = std::vector<double>(p.size());
	for (uint32_t i = 0; i < (uint32_t) p.size(); i++) v[i] = pow(f(p[i]), n);
	r.Create(p, v, f.degree);
}

// f..g through f(g(p)) on the points of f
void Concat(TBSpline &f, Interpolation::InterpolatingFunction &g, TBSpline &r)
{
	std::vector<double> p = SamplePoints(f), v = std::vector<double>(p.size());
	for (uint32_t i = 0; i < (uint32_t) p.size(); i++) v[i] = f(g(p[i]));
	r.Create(p, v, f.degree);
}

// the same degree, knots and coefficients, to 1e-14 relative
bool operator==(TBSpline &f, TBSpline &g)
{
	const double eps = 1E-14;
	if (f.type != g.type || f.degree != g.degree || f.knots.size() != g.knots.size() || f.y.size() != g.y.size()) return false;
	for (uint32_t i = 0; i < (uint32_t) f.knots.size(); i++)
		if (fabs(f.knots[i] - g.knots[i]) > eps*fabs(f.knots[i])) return false;
	for (uint32_t i = 0; i < (uint32_t) f.y.size(); i++)
		if (fabs(f.y[i] - g.y[i]) > eps*fabs(f.y[i])) return false;
	return true;
}

} /* namespace BSpline */
