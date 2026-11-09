/*
 * CubicHermiteSpline.cpp
 *
 *  Created on: Oct 23, 2018
 *      Author: brass
 */

#include "CubicHermiteSpline.h"

namespace Interpolation {

CubicHermiteSpline::CubicHermiteSpline() {
	type = "CubicHermiteSpline";

}

CubicHermiteSpline::~CubicHermiteSpline() {
	// TODO Auto-generated destructor stub
}

CubicHermiteSpline::CubicHermiteSpline(const std::vector<double> &u, const std::vector<double> &v)
{
	Create(u,v);
    type = "CubicHermiteSpline";
}
CubicHermiteSpline::CubicHermiteSpline(const std::vector<double> &u, const std::vector<double> &v, const std::vector<double> &w)
{
	Create(u,v,w);
    type = "CubicHermiteSpline";
}

void CubicHermiteSpline::Create(const std::vector<double> &u, const std::vector<double> &v)
{
	unsigned long n = u.size();
	std::vector<double> delta;
	std::vector<double> m;
	m = std::vector<double>(n);
	if (n < 2)
	{
		Create(u,v,m);
		return;
	}
	delta = std::vector<double>(n-1);
	for(unsigned long k = 0; k < n-1; k++)
	{
		delta[k] = (v[k+1] - v[k])/(u[k+1] - u[k]);
	}
	m[0  ] = delta[0  ];
	m[n-1] = delta[n-2];
	for(unsigned long k = 1; k < n-1; k++)
	{
		m[k] = (delta[k-1] + delta[k]) / 2;//((u[k+1]-u[k])*delta[k-1] + (u[k]-u[k-1])*delta[k]) / (2 * (u[k+1]-u[k-1]) );
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
			//double gamma = 0;
			//if (k > 0) gamma = m[k]/delta[k-1];
			if (alpha < 0 ||beta < 0) // wikipedia says gamma but original paper says beta
			{
				m[k] = 0;
				m[k+1] = 0; // according to original paper
			}
			else
			{
				//double r2 = alpha*alpha + beta*beta;
				//if (r2 > 9)
				if (alpha+beta > 3)
				{
					//double tau = 3./sqrt(r2);
					double tau = 3./(alpha+beta);
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

// The spline at t; 0 outside the grid
double CubicHermiteSpline::operator()(double t)
{
	return (*this)(t, 0);
}

// The k-th derivative of the spline at t; 0 outside the grid and for k > 3
double CubicHermiteSpline::operator()(double t, unsigned int k)
{
	uint32_t i = GridInterval(x, t);
	if (i == 0) return 0.;
	double d  = x[i]-x[i-1];
	double tx = (t-x[i-1])/d;
	switch (k) {
		case 0:
			return y[i-1]  + tx * (m[i-1]*d + tx * (3*(y[i]-y[i-1]) - d*(2*m[i-1]+m[i]) + tx * (2*(y[i-1]-y[i]) + d*(m[i-1]+m[i]))));
		case 1:
			return m[i-1] + tx * (6*(y[i]-y[i-1]) - 2*d*(2*m[i-1]+m[i]) + 3*tx * (2*(y[i-1]-y[i]) + d*(m[i-1]+m[i])))/d;
		case 2:
			return (6*(y[i]-y[i-1]) - 2*d*(2*m[i-1]+m[i]) + 6*tx * (2*(y[i-1]-y[i]) + d*(m[i-1]+m[i])))/d/d;
		case 3:
			return 6 * (2*(y[i-1]-y[i]) + d*(m[i-1]+m[i]))/d/d/d;
		default:
			return 0.;
	}
}

std::vector<double> CubicHermiteSpline::operator()(const std::vector<double> &t)
{
	std::vector<double> val = std::vector<double>(t.size());
	for (uint32_t j = 0; j < (uint32_t) t.size(); j++) val[j] = (*this)(t[j], 0);
	return val;
}

std::vector<double> CubicHermiteSpline::operator()(const std::vector<double> &t, unsigned int k)
{
	std::vector<double> val = std::vector<double>(t.size());
	for (uint32_t j = 0; j < (uint32_t) t.size(); j++) val[j] = (*this)(t[j], k);
	return val;
}

void UnaryMinus(CubicHermiteSpline &f,CubicHermiteSpline &g)
{
    uint32_t n = (uint32_t) f.x.size();
    g.x = f.x;
    g.y = std::vector<double>(n);
    g.m = std::vector<double>(n);
    for (uint32_t i = 0; i < n; i++)
    {
    	g.y[i]  = -f.y[i];
    	g.m[i]  = -f.m[i];
    }
}

bool operator==(CubicHermiteSpline &f,CubicHermiteSpline &g)
{
    const double eps = 1E-14;
    vector<double> x1 = f.GetX();
    vector<double> x2 = g.GetX();
    vector<double> y1 = f.GetY();
    vector<double> y2 = g.GetY();
    if (x1.size() != x2.size() || f.type != g.type) return false;
    for (uint32_t i = 0; i < (uint32_t) x1.size(); i++)
    {
        if ( fabs(x1[i] - x2[i]) > eps*fabs(x1[i]) || fabs(y1[i] - y2[i]) > eps*fabs(y1[i])|| fabs(f.m[i] - g.m[i]) > eps*fabs(f.m[i])) return false;
    }
    return true;
}

} /* namespace Interpolation */
