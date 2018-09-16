/*
 * CubicSpline.cpp
 *
 *  Created on: Sep 8, 2018
 *      Author: brass
 */

#include <iostream>

namespace Interpolation {

CubicSpline::CubicSpline() {
	// TODO Auto-generated constructor stub

}

CubicSpline::~CubicSpline() {
	// TODO Auto-generated destructor stub
}

CubicSpline::CubicSpline(std::vector<double> u, std::vector<double> v)
{
	Create(u,v);
}

void CubicSpline::Create(std::vector<double> u, std::vector<double> v)
{
	x = u;
	y = v;
	int n = (int)x.size()-1;
	y2 = std::vector<double>(n+1);
	std::vector<double> cc,dd;
	cc = std::vector<double>(n-1);
	dd = std::vector<double>(n-1);
	double a,b,c, d;
	double h0 = x[1]-x[0];
	double h1 = x[2]-x[1];
	b = 2/3.*(h0+h1);
	c = h1/3.;
	d = (y[2] - y[1])/h1 -(y[1]-y[0])/h0;
	cc[0] = c/b;
	dd[0] = d/b;
	for (int i = 1; i <= n-2; i++)
	{
		h0 = x[i+1] - x[i];
		h1 = x[i+2] - x[i+1];
		a = h0 / 3.;
		b = 2. / 3. * (h1+h0);
		c = h1/3.;
		d = (y[i+2]-y[i+1])/h1 - (y[i+1]-y[i])/h0;
		cc[i] = c/(b-a*cc[i-1]);
		dd[i] = (d - dd[i-1] * a) / (b -a * cc[i-1]);
	}
	y2[0] = 0;
	y2[n-1] = dd[n-2];
	for (int i = n-2; i > 0; i--)
	{
		y2[i] = dd[i-1] - cc[i-1] * y2[i+1];
	}
	y2[n] = 0;
}

double CubicSpline::operator()(double t)
{
	for (unsigned int i = 1; i < x.size(); i++)
	{
		if (t <= x[i])
		{
			double a,b,c,h, xx;
			xx = t - x[i-1];
			h = x[i] - x[i-1];
			b = y2[i-1];
			a = (y2[i] -y2[i-1]) / h / 3.;
			c = (y[i]-y[i-1])/h -h/3. * (y2[i] + 2*y2[i-1]);

			return y[i-1] + xx*(c +xx*(b + xx * a));
		}
	}

	std::cout << "error: value out of interpolation range" << std::endl;
	return NAN;
}

double CubicSpline::operator()(double t,unsigned int k)
{
	for (unsigned int i = 1; i < x.size(); i++)
	{
		if (t <= x[i])
		{
			double a,b,c,h, xx;
			xx = t - x[i-1];
			h = x[i] - x[i-1];
			b = y2[i-1];
			a = (y2[i] -y2[i-1]) / h / 3.;
			c = (y[i]-y[i-1])/h -h/3. * (y2[i] + 2*y2[i-1]);
			switch (k) {
				case 0:
					return y[i-1] + xx*(c +xx*(b + xx * a));
					break;
				case 1:
					return c +xx*(2*b + xx * a *3);
					break;
				case 2:
					return 2*b + xx * a *6;
					break;
				default:
					std::cout << "error: order of derivative of interpolated function at most 2" << std::endl;
					return NAN;
					break;
			}

		}
	}

	std::cout << "error: value out of interpolation range" << std::endl;
	return NAN;
}

std::vector<double> CubicSpline::operator()(std::vector<double> t)
{
	std::vector<double> val = std::vector<double>(t.size());
	uint32_t i = 1;
	for (uint32_t j = 0; j < t.size(); j++)
	{
		while (t[j] > x[i] && i < x.size()-1)
		{
			i++;
		}
		double a,b,c,h, xx;
		xx = t[j] - x[i-1];
		h = x[i] - x[i-1];
		b = y2[i-1];
		a = (y2[i] -y2[i-1]) / h / 3.;
		c = (y[i]-y[i-1])/h -h/3. * (y2[i] + 2*y2[i-1]);

		val[j] = y[i-1] + xx*(c +xx*(b + xx * a));
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
