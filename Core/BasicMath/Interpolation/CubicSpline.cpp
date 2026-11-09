/*
 * CubicSpline.cpp
 *
 *  Created on: Sep 8, 2018
 *      Author: brass
 */

#include <iostream>

namespace Interpolation {

CubicSpline::CubicSpline() {
	type = "CubicSpline";

}

CubicSpline::~CubicSpline() {
	// TODO Auto-generated destructor stub
}

CubicSpline::CubicSpline(const std::vector<double> &u, const std::vector<double> &v)
{
	Create(u,v);
    type = "CubicSpline";
}

void CubicSpline::Create(const std::vector<double> &u, const std::vector<double> &v)
{
	x = u;
	y = v;
	int n = (int)x.size()-1;
	y2 = std::vector<double>(n+1 > 0 ? n+1 : 0);
	if (n < 2) return;  // 2 points: a straight line, y'' = 0
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

void CubicSpline::Create(const std::vector<double> &u, const std::vector<double> &v, double yp0, double ypn)
{
    x = u;
	y = v;
	int n = (int)x.size()-1;
    if (n < 1)
    {
        y2 = std::vector<double>(n+1 > 0 ? n+1 : 0);
        return;
    }
    std::vector<double> a = std::vector<double>(n);
    std::vector<double> b = std::vector<double>(n+1);
    std::vector<double> c = std::vector<double>(n);
    std::vector<double> d = std::vector<double>(n+1);
    double h0 = x[1] - x[0];
    double h1 = x[n] -x[n-1];
    d[0] = (y[1]-y[0])/h0 - yp0;
    d[n] = ypn - (y[n]-y[n-1])/h1;
    c[0]   = 2*h0 / 6.;
    a[n-1] = 2*h1 / 6.;
    b[0]   = 2*h0 / 3.;
    b[n]   = 2*h1 / 3.;
    for (uint32_t i = 1; i < n; i++)
    {
        h0     = x[i] - x[i-1];
        h1     = x[i+1] - x[i];
        a[i-1] = 2*h0 / 6.;
        b[i]   = 2*(h0+h1) / 3.;
        c[i]   = 2*h1 / 6.;
        d[i]   = (y[i+1] - y[i])/h1 - (y[i] - y[i-1])/h0;
    }
    ThomasAlgorithm(a,b,c,d,y2);
}

// The spline at t; 0 outside the grid
double CubicSpline::operator()(double t)
{
	return (*this)(t, 0);
}

// The k-th derivative of the spline at t; 0 outside the grid and for k > 3
double CubicSpline::operator()(double t,unsigned int k)
{
	uint32_t i = GridInterval(x, t);
	if (i == 0) return 0.;
	double a,b,c,h, xx;
	xx = t - x[i-1];
	h = x[i] - x[i-1];
	b = y2[i-1];
	a = (y2[i] -y2[i-1]) / h / 3.;
	c = (y[i]-y[i-1])/h -h/3. * (y2[i] + 2*y2[i-1]);
	switch (k) {
		case 0:
			return y[i-1] + xx*(c +xx*(b + xx * a));
		case 1:
			return c +xx*(2*b + xx * a *3);
		case 2:
			return 2*b + xx * a *6;
		case 3:
			return a*6;
		default:
			return 0.;
	}
}

std::vector<double> CubicSpline::operator()(const std::vector<double> &t)
{
	std::vector<double> val = std::vector<double>(t.size());
	for (uint32_t j = 0; j < (uint32_t) t.size(); j++) val[j] = (*this)(t[j], 0);
	return val;
}

std::vector<double> CubicSpline::operator()(const std::vector<double> &t, unsigned int k)
{
	std::vector<double> val = std::vector<double>(t.size());
	for (uint32_t j = 0; j < (uint32_t) t.size(); j++) val[j] = (*this)(t[j], k);
	return val;
}

std::vector<double> CubicSpline::GetY2()
{
	return y2;
}
double CubicSpline::GetY2(uint32_t i)
{
	return y2[i];
}

void Pow(InterpolatingFunction &f, double n, CubicSpline &g)
{
    std::vector<double> x,y;
    x = f.GetX();
    uint32_t N = (uint32_t) x.size();
    y = f(x);
    for (uint32_t i = 0; i < N; i++)
    {
        y[i] = pow(y[i],n);
    }
    g.Create(x,y);
}

void Concat(InterpolatingFunction &f, InterpolatingFunction &g, CubicSpline &fg)
{
    std::vector<double> x,y;
    x = f.GetX();
    uint32_t n = (uint32_t) x.size();
    y = g(x);
    for (uint32_t i = 0; i < n; i++)
    {
        y[i] = f(y[i]);
    }
    fg.Create(x,y);
}

void UnaryMinus(CubicSpline &f,CubicSpline &g)
{
    uint32_t n = (uint32_t) f.x.size();
    g.x = f.x;
    g.y = std::vector<double>(n);
    g.y2 = std::vector<double>(n);
    for (uint32_t i = 0; i < n; i++)
    {
    	g.y[i]  = -f.y[i];
    	g.y2[i]  = -f.y2[i];
    }
}

void Add(InterpolatingFunction &f, InterpolatingFunction &g, CubicSpline &fplusg)
{
    MergeGrids(f.x,g.x,fplusg.x);
    uint32_t n = (uint32_t) fplusg.x.size();
    fplusg.y = std::vector<double>(n);
    fplusg.y2 = std::vector<double>(n);
    for (uint32_t i = 0; i < n; i++)
    {
    	fplusg.y[i]  = f(fplusg.x[i])   + g(fplusg.x[i]);
    	fplusg.y2[i] = (f(fplusg.x[i],2) + g(fplusg.x[i],2))/2;
    }

    /*
    double yp0 = f(x[0],1) + g(x[0],1);
    double ypn = f(x[n-1],1) + g(x[n-1],1);
    fplusg.Create(x,y,yp0,ypn);*/
}

void Mul(InterpolatingFunction &f, InterpolatingFunction &g, CubicSpline &fplusg)
{
    std::vector<double> x,y,x1,x2;
    x1 = f.GetX();
    x2 = g.GetX();
    MergeGrids(x1,x2,x);
    uint32_t n = (uint32_t) x.size();
    y = std::vector<double>(n);
    for (uint32_t i = 0; i < n; i++)
    {
        y[i] = f(x[i]) * g(x[i]);
    }
    fplusg.Create(x,y);
}

void Sub(InterpolatingFunction &f, InterpolatingFunction &g, CubicSpline &fplusg)
{
	MergeGrids(f.x,g.x,fplusg.x);
	uint32_t n = (uint32_t) fplusg.x.size();
	fplusg.y = std::vector<double>(n);
	fplusg.y2 = std::vector<double>(n);
	for (uint32_t i = 0; i < n; i++)
	{
		fplusg.y[i]  = f(fplusg.x[i])   - g(fplusg.x[i]);
		fplusg.y2[i] = (f(fplusg.x[i],2) - g(fplusg.x[i],2))/2;
	}
}

void Div(InterpolatingFunction &f, InterpolatingFunction &g, CubicSpline &fplusg)
{
    std::vector<double> x,y,x1,x2;
    x1 = f.GetX();
    x2 = g.GetX();
    MergeGrids(x1,x2,x);
    uint32_t n = (uint32_t) x.size();
    y = std::vector<double>(n);
    for (uint32_t i = 0; i < n; i++)
    {
        y[i] = f(x[i]) / g(x[i]);
    }
    fplusg.Create(x,y);
}

bool operator==(CubicSpline &f,CubicSpline &g)
{
    const double eps = 1E-14;
    vector<double> x1 = f.GetX();
    vector<double> x2 = g.GetX();
    vector<double> y1 = f.GetY();
    vector<double> y2 = g.GetY();
    if (x1.size() != x2.size() || f.type != g.type) return false;
    for (uint32_t i = 0; i < (uint32_t) x1.size(); i++)
    {
        if ( fabs(x1[i] - x2[i]) > eps*fabs(x1[i]) || fabs(y1[i] - y2[i]) > eps*fabs(y1[i])|| fabs(f.y2[i] - g.y2[i]) > eps*fabs(f.y2[i])) return false;
    }
    return true;
}


} /* namespace Interpolation */

void ThomasAlgorithm(const std::vector<double> &a,const std::vector<double> &b,const std::vector<double> &c,const std::vector<double> &d, std::vector<double> &x)
{
    uint32_t n = (uint32_t) b.size();
    std::vector<double> cc = std::vector<double>(n-1);
    std::vector<double> dd = std::vector<double>(n);
    x = std::vector<double>(n);
    cc[0] = c[0]/b[0];
    dd[0] = d[0]/b[0];
    
    for (uint32_t i = 1; i < n-1; i++)
    {
        double denom = b[i] - cc[i-1]*a[i-1];
        cc[i] = c[i]/denom;
        dd[i] = (d[i] - dd[i-1] * a[i-1]) / denom;
    }
    x[n-1] = (d[n-1] - dd[n-2] * a[n-2]) / (b[n-1] - cc[n-2]*a[n-2]);
    for (int i = n-2; i >= 0; i--)
    {
        x[i] = dd[i] - cc[i] * x[i+1];
    }
}















































