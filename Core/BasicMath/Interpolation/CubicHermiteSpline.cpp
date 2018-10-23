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

CubicHermiteSpline::CubicHermiteSpline(const std::vector<double> &u, const std::vector<double> &v)
{
	Create(u,v);
}
CubicHermiteSpline::CubicHermiteSpline(const std::vector<double> &u, const std::vector<double> &v, const std::vector<double> &w)
{
	Create(u,v,w);
}

void CubicHermiteSpline::Create(const std::vector<double> &u, const std::vector<double> &v)
{
	unsigned long n = u.size();
	std::vector<double> delta;
	std::vector<double> m;
	delta.reserve(n-1);
	m = std::vector<double>(n);
	for(unsigned long k = 0; k < n-1; k++)
	{
		delta[k] = (v[k+1] - v[k])/(u[k+1] - u[k]);
	}
	m[0  ] = delta[0  ];
	m[n-1] = delta[n-2];
	for(unsigned long k = 1; k < n-1; k++)
	{
		m[k] = (delta[k-1] + delta[k]) / 2.;
	}
	for(unsigned long k = 0; k < n-1; k++)
	{
		if (delta[k] == 0)
		{
			m[k  ] = 0;
			m[k+1] = 0;
		}
		else
		{
			double alpha = m[k  ]/delta[k];
			double beta  = m[k+1]/delta[k];
			double gamma = 0;
			if (k > 0) gamma = m[k]/delta[k-1];
			if (alpha < 0 || gamma < 0)
			{
				m[k] = 0;
			}
			else
			{
				double r2 = alpha*alpha + beta*beta;
				if (r2 > 9)
				{
					double tau = 3./sqrt(r2);
					m[k  ] *= tau;
					m[k+1] *= tau;
				}
			}
		}
	}
	Create(u,v,m);
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
