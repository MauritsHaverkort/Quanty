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
	y = v;
	degree = n;
	GenerateGrid(u,knots,n);
	uint32_t nx = (uint32_t) knots.size()-2*n;
	x = std::vector<double>(nx);
	for (uint32_t i = 0; i < nx; i++)
	{
		x[i] = knots[i+n];
	}
	std::vector<double> B;
	GetKnotMatrix(knots,n,u,B);
	uint32_t N = (uint32_t) u.size();
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

double TBSpline::operator()(double t)
{
	uint32_t idx = degree;
	uint32_t N = (uint32_t) knots.size();
	if ( t < knots[0] || t > knots[N-1])
	{
		printf("Error in TBSpline(): %E is out of range\n",t);
		return NAN;
	}
	while ( idx < N-1 && t > knots[idx+1])
	{
		idx++;
	}
	return deBoor(idx,t,knots,degree,y);
}
double TBSpline::operator()(double t, unsigned int k)
{
	if (k > degree)
	{
		printf("Warning: you just took a derivative that is higher than the degree of the interpolating function.\n This clearly is zero and I will return zero, but I'm wondering if this result ismeaningfull...\n");
		return 0;
	}
	uint32_t idx = degree;
	uint32_t N = (uint32_t) knots.size();
	if ( t < knots[0] || t > knots[N-1])
	{
		printf("Error in TBSpline(): %E is out of range\n",t);
		return NAN;
	}
	while ( idx < N-1 && t > knots[idx+1])
	{
		idx++;
	}
	return Derivative(idx,t,knots,degree,y,k);
}
std::vector<double> TBSpline::operator()(const std::vector<double> &t)
{
	std::vector<double> val = std::vector<double>(t.size());
	uint32_t i = degree;
	for (uint32_t j = 0; j < (uint32_t) t.size(); j++)
	{
		while (i < (uint32_t) knots.size()-1 && t[j] > knots[i+1])
		{
			i++;
		}

		val[j] = deBoor(i,t[j],knots,degree,y);;
		if(j < (uint32_t) t.size()-1) //TODO improve that
		{
			if (t[j+1] < t[j])
			{
				i = 1;
			}
		}
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
	uint32_t i = degree;
	for (uint32_t j = 0; j < (uint32_t) t.size(); j++)
	{
		while (i < (uint32_t) knots.size()-1 && t[j] > knots[i+1])
		{
			i++;
		}

		val[j] = Derivative(i,t[j],knots,degree,y,k);
		if(j < (uint32_t) t.size()-1) //TODO improve that
		{
			if (t[j+1] < t[j])
			{
				i = 1;
			}
		}
	}
	return val;
}

} /* namespace BSpline */
