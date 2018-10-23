/*
 * CubicHermiteSpline.cpp
 *
 *  Created on: Oct 23, 2018
 *      Author: brass
 */

#include "CubicHermiteSpline.h"

namespace Interpolation {

CubicHermiteSpline::CubicHermiteSpline() {
	// TODO Auto-generated constructor stub

}

CubicHermiteSpline::~CubicHermiteSpline() {
	// TODO Auto-generated destructor stub
}

CubicHermiteSpline::CubicHermiteSpline(const std::vector<double> &u, const std::vector<double> &v, const std::vector<double> &w)
{
	Create(u,v,w);
}

void CubicHermiteSpline::Create(const std::vector<double> &u, const std::vector<double> &v, const std::vector<double> &w)
{
	x = u;
	y = v;
	m = w;
}

double CubicHermiteSpline::operator()(double t)
{
	for (unsigned int i = 1; i < x.size(); i++)
	{
		if (t <= x[i])
		{
			double d  = x[i]-x[i-1];
			double tx = (t-x[i-1])/d;
			return y[i-1]  + tx * (m[i-1]*d + tx * (3*(y[i]-y[i-1]) - d*(2*m[i-1]+m[i]) + tx * (2*(y[i-1]-y[i]) + d*(m[i-1]+m[i]))));
		}
	}

	std::cout << "error: value out of interpolation range" << std::endl;
	return NAN;
}

double CubicHermiteSpline::operator()(double t, unsigned int k)
{
	return NAN;
}
std::vector<double> CubicHermiteSpline::operator()(const std::vector<double> &t)
{
	std::vector<double> val;
	return val;
}

} /* namespace Interpolation */
