/*
 * InterpolatingFunction.cpp
 *
 *  Created on: Sep 8, 2018
 *      Author: brass
 */

namespace Interpolation {

InterpolatingFunction::InterpolatingFunction() {
	// TODO Auto-generated constructor stub

}

InterpolatingFunction::~InterpolatingFunction() {
	// TODO Auto-generated destructor stub
}

double InterpolatingFunction::operator()(double t)
{
	printf("error: call of virtual InterpolatingFunction::operator()\n");
	return NAN;
}
double InterpolatingFunction::operator()(double t, unsigned int n)
{
	printf("error: call of virtual InterpolatingFunction::operator()\n");
	return NAN;
}
std::vector<double> InterpolatingFunction::operator()(const std::vector<double> &t)
{
	printf("error: call of virtual InterpolatingFunction::operator()\n");
	return std::vector<double>();
}
std::vector<double> InterpolatingFunction::operator()(const std::vector<double> &t, unsigned int k)
{
	printf("error: call of virtual InterpolatingFunction::operator()\n");
	return std::vector<double>();
}

void InterpolatingFunction::Create(const std::vector<double> &u, const std::vector<double> &v)
{
    printf("you should not see this message\n");
}

std::vector<double> InterpolatingFunction::GetX()
{
	return x;
}

std::vector<double> InterpolatingFunction::GetY()
{
	return y;
}

double InterpolatingFunction::GetX(uint32_t i)
{
	return x[i];
}

double InterpolatingFunction::GetY(uint32_t i)
{
	return y[i];
}

void Pow(InterpolatingFunction &f, double n, InterpolatingFunction &g)
{
    std::vector<double> x,y;
    x = f.GetX();
    y = f(x);
    for (uint32_t i = 0; i < x.size(); i++)
    {
        y[i] = pow(y[i],n);
    }
    g.Create(x,y);
}

void Concat(InterpolatingFunction &f, InterpolatingFunction &g, InterpolatingFunction &fg)
{
    std::vector<double> x,y;
    x = f.GetX();
    y = g(x);
    for (uint32_t i = 0; i < x.size(); i++)
    {
        y[i] = f(y[i]);
    }
    fg.Create(x,y);
}

void Add(InterpolatingFunction &f, InterpolatingFunction &g, InterpolatingFunction &fplusg)
{
    std::vector<double> x,y,x1,x2;
    x1 = f.GetX();
    x2 = g.GetX();
    MergeGrids(x1,x2,x);
    y = std::vector<double>(x.size());
    for (uint32_t i = 0; i < x.size(); i++)
    {
        y[i] = f(x[i]) + g(x[i]);
    }
    fplusg.Create(x,y);
}

void Mul(InterpolatingFunction &f, InterpolatingFunction &g, InterpolatingFunction &fplusg)
{
    std::vector<double> x,y,x1,x2;
    x1 = f.GetX();
    x2 = g.GetX();
    MergeGrids(x1,x2,x);
    y = std::vector<double>(x.size());
    for (uint32_t i = 0; i < x.size(); i++)
    {
        y[i] = f(x[i]) * g(x[i]);
    }
    fplusg.Create(x,y);
}

void Sub(InterpolatingFunction &f, InterpolatingFunction &g, InterpolatingFunction &fplusg)
{
    std::vector<double> x,y,x1,x2;
    x1 = f.GetX();
    x2 = g.GetX();
    MergeGrids(x1,x2,x);
    y = std::vector<double>(x.size());
    for (uint32_t i = 0; i < x.size(); i++)
    {
        y[i] = f(x[i]) - g(x[i]);
    }
    fplusg.Create(x,y);
}

void Div(InterpolatingFunction &f, InterpolatingFunction &g, InterpolatingFunction &fplusg)
{
    std::vector<double> x,y,x1,x2;
    x1 = f.GetX();
    x2 = g.GetX();
    MergeGrids(x1,x2,x);
    y = std::vector<double>(x.size());
    for (uint32_t i = 0; i < x.size(); i++)
    {
        y[i] = f(x[i]) / g(x[i]);
    }
    fplusg.Create(x,y);
}

bool operator==(InterpolatingFunction &f,InterpolatingFunction &g)
{
    const double eps = 1E-14;
    vector<double> x1 = f.GetX();
    vector<double> x2 = g.GetX();
    vector<double> y1 = f.GetY();
    vector<double> y2 = g.GetY();
    if (x1.size() != x2.size() || f.type != g.type) return false;
    for (uint32_t i = 0; i < (uint32_t) x1.size(); i++)
    {
        if ( fabs(x1[i] - x2[i]) > eps*fabs(x1[i]) || fabs(y1[i] - y2[i]) > eps*fabs(y1[i])) return false;
    }
    return true;
}

/* Gauss points per interval of f.x that integrate f exactly. Two suffice for the
 * cubic splines; a B-spline of degree d needs (d+2)/2, at least 2. */
static uint32_t IntegrationPoints(InterpolatingFunction &f)
{
	if (f.type == "BSpline")
	{
		uint32_t degree = ((BSpline::TBSpline &) f).degree;
		return std::max((uint32_t) 2, (degree + 2)/2);
	}
	return 2;
}

double Integrate(InterpolatingFunction &f)
{
	uint32_t N = IntegrationPoints(f);
	double val = 0;
	vector<double> x,w,y;
	for (uint32_t i = 0; i < (uint32_t) f.x.size()-1; i++)
	{
		Quad::GaussianQuadratureWeights(N,x,w,f.x[i],f.x[i+1]);
		y = f(x);
		for (uint32_t j = 0; j < N; j++)
		{
			val += y[j]*w[j];
		}

	}
	return val;
}

double IntegrateInterval(InterpolatingFunction &f, uint32_t idx)
{
	uint32_t N = IntegrationPoints(f);
	double val = 0;
	vector<double> x,w,y;
	Quad::GaussianQuadratureWeights(N,x,w,f.x[idx],f.x[idx+1]);
	y = f(x);
	for (uint32_t j = 0; j < N; j++)
	{
		val += y[j]*w[j];
	}
	return val;
}

double Integrate(InterpolatingFunction &f, double a, double b)
{
	const uint32_t N = IntegrationPoints(f);
	double val = 0;
	double sign = 1;
	uint32_t i = 0;
	if (a > b)
	{
		sign = a;
		a = b;
		b = sign;
		sign = -1;
	}
	while(i < (uint32_t) f.x.size() && a > f.x[i] )
		 i++;
	if (b <= f.x[i])
	{
		val = Quad::Integrate(f,a,b,N);
	}
	else
	{
		val = Quad::Integrate(f,a,f.x[i],N);
		while(i < (uint32_t) f.x.size() -1 && b > f.x[i+1] )
		{
			val += Quad::Integrate(f,f.x[i],f.x[i+1],N);
			i++;
		}
		val += Quad::Integrate(f,f.x[i],b,N);
	}
	return sign*val;
}

/* The interval i = 1 .. n-1 of the grid x[0] < .. < x[n-1] with x[i-1] <= t <= x[i];
 * a grid point belongs to the interval on its left. 0 when t lies outside the grid, is
 * NaN, or the grid has fewer than 2 points. Binary search. */
uint32_t GridInterval(const std::vector<double> &x, double t)
{
	uint32_t n = (uint32_t) x.size(), lo = 1, hi, mid;
	if (n < 2 || !(t >= x[0]) || !(t <= x[n-1])) return 0;
	hi = n-1;                       // x[hi] >= t
	while (lo < hi)
	{
		mid = lo + (hi - lo)/2;
		if (x[mid] >= t) hi = mid;
		else lo = mid + 1;
	}
	return lo;
}

/* The sorted union of two grids; points that differ by less than 1e-14 relative count
 * as one, the smallest of them is kept. */
void MergeGrids(const std::vector<double> &x1,const std::vector<double> &x2, std::vector<double> &xout )
{
    const double eps = 1E-14;
    std::vector<double> all = std::vector<double>(x1.size() + x2.size());
    for (uint32_t i = 0; i < (uint32_t) x1.size(); i++) all[i] = x1[i];
    for (uint32_t i = 0; i < (uint32_t) x2.size(); i++) all[i + x1.size()] = x2[i];
    if (all.size() > 0) qsort(&(all[0]), all.size(), sizeof(double), &qsortcomparedoubles);
    xout.clear();
    for (uint32_t i = 0; i < (uint32_t) all.size(); i++)
    {
        if (xout.size() == 0 || !(fabs(all[i] - xout.back()) <= eps * (fabs(all[i]) + eps)))
        {
            xout.push_back(all[i]);
        }
    }
}

} /* namespace Interpolation */
