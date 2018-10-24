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
	for (unsigned int i = 1; i < x.size(); i++)
	{
		if (t <= x[i])
		{
			double d  = x[i]-x[i-1];
			double tx = (t-x[i-1])/d;
			switch (k) {
				case 0:
					return y[i-1]  + tx * (m[i-1]*d + tx * (3*(y[i]-y[i-1]) - d*(2*m[i-1]+m[i]) + tx * (2*(y[i-1]-y[i]) + d*(m[i-1]+m[i]))));
					break;
				case 1:
					return m[i-1] + tx * (6*(y[i]-y[i-1]) - 2*d*(2*m[i-1]+m[i]) + 3*tx * (2*(y[i-1]-y[i]) + d*(m[i-1]+m[i])))/d;
					break;
				case 2:
					return (6*(y[i]-y[i-1]) - 2*d*(2*m[i-1]+m[i]) + 6*tx * (2*(y[i-1]-y[i]) + d*(m[i-1]+m[i])))/d/d;
					break;
				case 3:
					return 6 * (2*(y[i-1]-y[i]) + d*(m[i-1]+m[i]))/d/d/d;
					break;
				default:
					return 0;
					break;
			}

		}
	}
	std::cout << "error: value out of interpolation range" << std::endl;
	return NAN;
}
std::vector<double> CubicHermiteSpline::operator()(const std::vector<double> &t)
{
	std::vector<double> val = std::vector<double>(t.size());
	uint32_t i = 1;
	for (uint32_t j = 0; j < t.size(); j++)
	{
		while (t[j] > x[i] && i < x.size()-1)
		{
			i++;
		}
		double d  = x[i]-x[i-1];
		double tx = (t[j]-x[i-1])/d;
		val[j] = y[i-1]  + tx * (m[i-1]*d + tx * (3*(y[i]-y[i-1]) - d*(2*m[i-1]+m[i]) + tx * (2*(y[i-1]-y[i]) + d*(m[i-1]+m[i]))));
		if(j < t.size()-1) //TODO improve that
		{
			if (t[j+1] < t[j])
			{
				i = 1;
			}
		}
	}
	return val;
}
std::vector<double> CubicHermiteSpline::operator()(const std::vector<double> &t, unsigned int k)
{
	std::vector<double> val = std::vector<double>(t.size());
	uint32_t i = 1;
	for (uint32_t j = 0; j < t.size(); j++)
	{
		while (t[j] > x[i] && i < x.size()-1)
		{
			i++;
		}
		double d  = x[i]-x[i-1];
		double tx = (t[j]-x[i-1])/d;
		switch (k) {
			case 0:
				val[j] = y[i-1]  + tx * (m[i-1]*d + tx * (3*(y[i]-y[i-1]) - d*(2*m[i-1]+m[i]) + tx * (2*(y[i-1]-y[i]) + d*(m[i-1]+m[i]))));
				break;
			case 1:
				val[j] = m[i-1] + tx * (6*(y[i]-y[i-1]) - 2*d*(2*m[i-1]+m[i]) + 3*tx * (2*(y[i-1]-y[i]) + d*(m[i-1]+m[i])))/d;
				break;
			case 2:
				val[j] = (6*(y[i]-y[i-1]) - 2*d*(2*m[i-1]+m[i]) + 6*tx * (2*(y[i-1]-y[i]) + d*(m[i-1]+m[i])))/d/d;
				break;
			case 3:
				val[j] = 6 * (2*(y[i-1]-y[i]) + d*(m[i-1]+m[i]))/d/d/d;
				break;
			default:
				val[j] = 0;
				break;
		}
		if(j < t.size()-1) //TODO improve that
		{
			if (t[j+1] < t[j])
			{
				i = 1;
			}
		}
	}
	return val;
}

} /* namespace Interpolation */
