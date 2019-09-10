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

void TBSpline::Create(const std::vector<double> &u, const std::vector<double> &v,uint32_t n)
{
	x = u;
	y = v;
	degree = n;
	GenerateGrid(x,knots,n);
	std::vector<double> B;
	GetKnotMatrix(knots,n,x,B);
	uint32_t N = (uint32_t) x.size();
	uint32_t NRHS = 1;
	vector<int> IPIV = vector<int>(N);
	int info = 0;
	dgesv(&N,&NRHS,&(B[0]),&N,&(IPIV[0]),&(y[0]),&N,&info);

}

double TBSpline::operator()(double t)
{
	uint32_t idx = 0;
	uint32_t N = (uint32_t) x.size();
	if ( t < x[0] || t > x[N-1])
	{
		printf("Error in TBSpline(): %E is out of range\n",t);
		return NAN;
	}
	while ( idx < N-1 && t >= x[idx+1])
	{
		idx++;
	}
	return deBoor(idx,t,knots,degree,y);
}
double TBSpline::operator()(double t, unsigned int k)
{
	return 0;
}
std::vector<double> TBSpline::operator()(const std::vector<double> &t)
{
	std::vector<double> val = std::vector<double>(t.size());
	return val;
}
std::vector<double> TBSpline::operator()(const std::vector<double> &t, unsigned int k)
{
	std::vector<double> val = std::vector<double>(t.size());
	return val;
}

} /* namespace BSpline */
