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
	type = "BSpline";
}

void TBSpline::Create(const std::vector<double> &u, const std::vector<double> &v,uint32_t n)
{
	uint32_t N = (uint32_t) u.size();
	y = v;
	degree = n;
	knots = std::vector<double>(N + n + 1);
	int error = BSplineKnotsAveraged(&u[0], N, n, &knots[0]);
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
	std::vector<double> B = std::vector<double>(N*N);
	BSplineCollocationMatrix(&knots[0], (uint32_t) knots.size(), n, &u[0], N, &B[0]);
	uint32_t NRHS = 1;
	vector<int> IPIV = vector<int>(N);
	int info = 0;
	dgesv(&N,&NRHS,&(B[0]),&N,&(IPIV[0]),&(y[0]),&N,&info);
	if (info)
	{
		printf("Error TBSpline.Create: dgesv info = %i\n",info);
	}
	type = "BSpline";
}

// The spline at t; 0 outside the knots
double TBSpline::operator()(double t)
{
	return BSplineDerivative(t, &knots[0], (uint32_t) knots.size(), degree, &y[0], 0);
}

// The k-th derivative of the spline at t; 0 outside the knots
double TBSpline::operator()(double t, unsigned int k)
{
	if (k > degree)
	{
		printf("Warning: you just took a derivative that is higher than the degree of the interpolating function.\n This clearly is zero and I will return zero, but I'm wondering if this result ismeaningfull...\n");
		return 0;
	}
	return BSplineDerivative(t, &knots[0], (uint32_t) knots.size(), degree, &y[0], k);
}

std::vector<double> TBSpline::operator()(const std::vector<double> &t)
{
	std::vector<double> val = std::vector<double>(t.size());
	for (uint32_t j = 0; j < (uint32_t) t.size(); j++)
	{
		val[j] = BSplineDerivative(t[j], &knots[0], (uint32_t) knots.size(), degree, &y[0], 0);
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
		val[j] = BSplineDerivative(t[j], &knots[0], (uint32_t) knots.size(), degree, &y[0], k);
	}
	return val;
}

} /* namespace BSpline */
