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

double Integrate(InterpolatingFunction &f)
{
	uint32_t N = 2;
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
	uint32_t N = 2;
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
	const uint32_t N = 2;
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

void MergeGrids(const std::vector<double> &x1,const std::vector<double> &x2, std::vector<double> &xout )
{
    uint32_t n1 = (uint32_t) x1.size();
    uint32_t n2 = (uint32_t) x2.size();
    uint32_t nout = n1+n2;
    xout = std::vector<double>(nout);
    const double eps = 1E-14;
    for (uint32_t i = 0; i < n1; i++)
    {
        xout[i] = x1[i];
    }
    for (uint32_t i = 0; i < n2; i++)
    {
        xout[i+n1] = x2[i];
    }    
    qsort(&(xout[0]), nout, sizeof(double), &qsortcomparedoubles);
    for (std::vector<double>::iterator it = xout.end()-1; it != xout.begin(); it--)
    {
        if (fabs(((*it)-(*(it-1)))) <= eps * (fabs(*it)+eps))
        {
            xout.erase(it);
        }
    }
}

} /* namespace Interpolation */
