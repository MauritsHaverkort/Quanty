/*
 * Quad.cpp
 *
 *  Created on: May 2, 2019
 *      Author: brass
 */

#include "Quad.h"

namespace Quad
{

void GaussianQuadratureWeights(uint32_t N, vector<double> &x, vector<double> &w, double a, double b)
{
	x = vector<double>(N);
	w = vector<double>(N);
	vector<double> subdiag = vector<double>(N-1);
	for (uint32_t j = 1; j < N; j++)
	{
		subdiag[j-1] = j / sqrt(4 * j*j -1);
	}
	char compz = 'N';
	int info;
	uint32_t bar = 1;
	double *foo=NULL;
	dsteqr(&compz,&N,&x[0],&subdiag[0],foo,&bar,foo,&info);
	for (uint32_t j = 0; j < N; j++)
	{
		double p = LegendreP(x[j],N+1);
		w[j] = 2 * (1 - x[j]*x[j]) / (N+1) / (N+1) / p / p * (b-a)/2;
		x[j] = (b*(x[j]+1.)+a*(1.-x[j]))/2.;
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


