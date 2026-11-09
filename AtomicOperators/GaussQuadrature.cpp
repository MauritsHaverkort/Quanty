/*
 * GaussQuadrature.cpp
 *
 *  Created on: Sep 11, 2018
 *      Author: brass
 */

namespace Quad {

using std::vector;

// k Gauss-Legendre nodes on [-1,1], ascending, from the 34-digit table in
// Core/BasicMath/Interpolation/GaussLegendre (k = 1 .. 127)
vector<double> GetKnots(uint32_t k)
{
	vector<double> knots;
	if (k == 0 || k > 127)
	{
		printf("error in GaussQuad::GetKnots: number of knots must be 1 .. 127; got %u\n",k);
		fflush(stdout);
		return knots;
	}
	knots = vector<double>(k);
	vector<double> weights = vector<double>(k);
	GaussLegendreNodesWeights(k, -1., 1., &knots[0], &weights[0]);
	return knots;
}

// the k Gauss-Legendre weights on [-1,1] that belong to GetKnots(k)
vector<double> GetWeights(uint32_t k)
{
	vector<double> weights;
	if (k == 0 || k > 127)
	{
		printf("error in GaussQuad::GetWeights: number of knots must be 1 .. 127; got %u\n",k);
		fflush(stdout);
		return weights;
	}
	vector<double> knots = vector<double>(k);
	weights = vector<double>(k);
	GaussLegendreNodesWeights(k, -1., 1., &knots[0], &weights[0]);
	return weights;
}

// Martin Brass' tables (2018), to 6 digits, for k = 1 .. 7; replaced by the 34-digit
// table above. Kept for the record, not compiled.
#if 0
vector<double> GetKnotsOld(uint32_t k)
{
	vector<double> knots;
	switch (k) {
		case 1:
			knots.push_back(0.);
			break;
		case 2:
			knots.push_back(-0.57735);
			knots.push_back( 0.57735);
			break;
		case 3:
			knots.push_back(-0.774597);
			knots.push_back(0.);
			knots.push_back(0.774597);
			break;
		case 4:
			knots.push_back(-0.861136);
			knots.push_back(-0.339981);
			knots.push_back(0.339981);
			knots.push_back(0.861136);
			break;
		case 5:
			knots.push_back(-0.90618);
			knots.push_back(-0.538469);
			knots.push_back(0.);
			knots.push_back(0.538469);
			knots.push_back(0.90618);
			break;
		case 6:
			knots.push_back(-0.93247);
			knots.push_back(-0.661209);
			knots.push_back(-0.238619);
			knots.push_back(0.238619);
			knots.push_back(0.661209);
			knots.push_back(0.93247);
			break;
		case 7:
			knots.push_back(-0.949108);
			knots.push_back(-0.741531);
			knots.push_back(-0.405845);
			knots.push_back(0.);
			knots.push_back(0.405845);
			knots.push_back(0.741531);
			knots.push_back(0.949108);
			break;
		default:
			printf("error in GaussQuad::GetKnots: max number of knots = 7; got %u\n",k);
			fflush(stdout);
			break;
	}
	return knots;
}

vector<double> GetWeightsOld(uint32_t k)
{
	vector<double> weights;
	switch (k) {
		case 1:
			weights.push_back(2.);
			break;
		case 2:
			weights.push_back(1.);
			weights.push_back(1.);
			break;
		case 3:
			weights.push_back(0.555556);
			weights.push_back(0.888889);
			weights.push_back(0.555556);
			break;
		case 4:
			weights.push_back(0.347855);
			weights.push_back(0.652145);
			weights.push_back(0.652145);
			weights.push_back(0.347855);
			break;
		case 5:
			weights.push_back(0.236927);
			weights.push_back(0.478629);
			weights.push_back(0.568889);
			weights.push_back(0.478629);
			weights.push_back(0.236927);
			break;
		case 6:
			weights.push_back(0.171324);
			weights.push_back(0.360762);
			weights.push_back(0.467914);
			weights.push_back(0.467914);
			weights.push_back(0.360762);
			weights.push_back(0.171324);
			break;
		case 7:
			weights.push_back(0.129485);
			weights.push_back(0.279705);
			weights.push_back(0.38183);
			weights.push_back(0.417959);
			weights.push_back(0.38183);
			weights.push_back(0.279705);
			weights.push_back(0.129485);
			break;
		default:
			printf("error in GaussQuad::GetWeights: max number of knots = 7; got %u\n",k);
			fflush(stdout);
			break;
	}
	return weights;
}

#endif

vector<double> GetWeights(uint32_t k, vector<double> x)
{
	vector<double> w = GetWeights(k);
	unsigned long  n = x.size();
	vector<double> weights = vector<double>(k*(n-1));
	for (unsigned long i = 0; i < n-1; i++)
	{
		for(uint32_t j = 0; j < k; j++)
		{
			weights[j + k*i] = (x[i+1]-x[i])/2 * w[j];
		}
	}
	return weights;
}

vector<double> TrafoXpoints(vector<double> x, double a, double b)
{
	vector<double> z = vector<double>(x.size());
	for(int i = 0; i < x.size(); i++)
	{
		z[i] = (b*(x[i]+1.)+a*(1.-x[i]))/2.;
	}
	return z;
}

vector<double> GetKnots(uint32_t k, double a, double b)
{
	return TrafoXpoints(GetKnots(k),a,b);
}

vector<double> GetKnots(uint32_t k, vector<double> grid)
{
	vector<double> result;
	result.reserve(k * (grid.size()-1));
	for (unsigned long i = 0; i < grid.size()-1; i++)
	{
		vector<double> knots = GetKnots(k,grid[i],grid[i+1]);
		result.insert(result.end(),knots.begin(),knots.end());
	}

	return result;
}

double Integrate(vector<double> f,double a, double b)
{
	if ( a == b)
		return 0.;
	uint32_t n = (uint32_t) f.size();
	if (n > 127)
	{
		cout << "error: in GaussQuad::Integrate: number of knots too big" << endl;
		return 0;
	}
	else
	{
		vector<double> weights = GetWeights(n);
		double result = 0;
		for (uint32_t i = 0; i < n; i++)
		{
			result += (b-a)/2*weights[i]*f[i];
		}
		return result;
	}
}

double Integrate(uint32_t k, vector<double> x, vector<double> f)
{
	uint32_t n = (uint32_t) x.size();
	if (k > 127)
	{
		cout << "error: in GaussQuad::Integrate: number of knots too big" << endl;
		return 0;
	}
	else
	{
		vector<double> weights = GetWeights(k);
		double result = 0;
		for (uint32_t j = 0; j < n; j++)
		{
			for (uint32_t i = 0; i < k; i++)
			{
				result += (x[j+1]-x[j])/2*weights[i]*f[i + k*j];
			}
		}
		return result;
	}
}

double GetMultipoleIntegral(Interpolation::CubicSpline f1, Interpolation::CubicSpline f2, uint32_t i, int k)
{
	double a1,a2,b1,b2,c1,c2,y1,y2, h,x0,x1;

	x0 = f1.GetX(i);
	x1 = f1.GetX(i+1);
	h = x1 - x0;

	b1 = f1.GetY2(i);
	a1 = (f1.GetY2(i+1) -f1.GetY2(i)) / h / 3.;
	c1 = (f1.GetY(i+1)-f1.GetY(i))/h -h/3. * (f1.GetY2(i+1) + 2*f1.GetY2(i));
	y1 = f1.GetY(i);

	b2 = f2.GetY2(i);
	a2 = (f2.GetY2(i+1) -f2.GetY2(i)) / h / 3.;
	c2 = (f2.GetY(i+1)-f2.GetY(i))/h -h/3. * (f2.GetY2(i+1) + 2*f2.GetY2(i));
	y2 = f2.GetY(i);

	double px02 = x0*x0;
	double px03 = px02*x0;
	double px04 = px03*x0;
	double px05 = px04*x0;
	double px06 = px05*x0;
	double px01k = pow(x0,1 + k);
	double px11k = pow(x1,1 + k);
	double px02k = px01k*x0;
	double px12k = px11k*x0;
	double px03k = px02k*x0;
	double px13k = px12k*x0;
	double px04k = px03k*x0;
	double px14k = px13k*x0;
	double px05k = px04k*x0;
	double px15k = px14k*x0;
	double px06k = px05k*x0;
	double px16k = px15k*x0;
	double px07k = px06k*x0;
	double px17k = px16k*x0;

	return (b2*c1*px03*(px01k - px11k))/(1 + k) +
			   (b1*c2*px03*(px01k - px11k))/(1 + k) +
			   (a2*b1*px05*(px01k - px11k))/(1 + k) +
			   (a1*b2*px05*(px01k - px11k))/(1 + k) +
			   (c1*c2*px02*(-px01k + px11k))/(1 + k) +
			   (b1*b2*px04*(-px01k + px11k))/(1 + k) +
			   (a2*c1*px04*(-px01k + px11k))/(1 + k) +
			   (a1*c2*px04*(-px01k + px11k))/(1 + k) +
			   (a1*a2*px06*(-px01k + px11k))/(1 + k) +
			   (2*c1*c2*x0*(px02k - px12k))/(2 + k) +
			   (4*b1*b2*px03*(px02k - px12k))/(2 + k) +
			   (4*a2*c1*px03*(px02k - px12k))/(2 + k) +
			   (4*a1*c2*px03*(px02k - px12k))/(2 + k) +
			   (6*a1*a2*px05*(px02k - px12k))/(2 + k) +
			   (3*b2*c1*px02*(-px02k + px12k))/(2 + k) +
			   (3*b1*c2*px02*(-px02k + px12k))/(2 + k) +
			   (5*a2*b1*px04*(-px02k + px12k))/(2 + k) +
			   (5*a1*b2*px04*(-px02k + px12k))/(2 + k) +
			   (3*b2*c1*x0*(px03k - px13k))/(3 + k) +
			   (3*b1*c2*x0*(px03k - px13k))/(3 + k) +
			   (10*a2*b1*px03*(px03k - px13k))/(3 + k) +
			   (10*a1*b2*px03*(px03k - px13k))/(3 + k) +
			   (c1*c2*(-px03k + px13k))/(3 + k) +
			   (6*b1*b2*px02*(-px03k + px13k))/(3 + k) +
			   (6*a2*c1*px02*(-px03k + px13k))/(3 + k) +
			   (6*a1*c2*px02*(-px03k + px13k))/(3 + k) +
			   (15*a1*a2*px04*(-px03k + px13k))/(3 + k) +
			   (4*b1*b2*x0*(px04k - px14k))/(4 + k) +
			   (4*a2*c1*x0*(px04k - px14k))/(4 + k) +
			   (4*a1*c2*x0*(px04k - px14k))/(4 + k) +
			   (20*a1*a2*px03*(px04k - px14k))/(4 + k) +
			   (b2*c1*(-px04k + px14k))/(4 + k) +
			   (b1*c2*(-px04k + px14k))/(4 + k) +
			   (10*a2*b1*px02*(-px04k + px14k))/(4 + k) +
			   (10*a1*b2*px02*(-px04k + px14k))/(4 + k) +
			   (5*a2*b1*x0*(px05k - px15k))/(5 + k) +
			   (5*a1*b2*x0*(px05k - px15k))/(5 + k) +
			   (b1*b2*(-px05k + px15k))/(5 + k) +
			   (a2*c1*(-px05k + px15k))/(5 + k) +
			   (a1*c2*(-px05k + px15k))/(5 + k) +
			   (15*a1*a2*px02*(-px05k + px15k))/(5 + k) +
			   (6*a1*a2*x0*(px06k - px16k))/(6 + k) +
			   (a2*b1*(-px06k + px16k))/(6 + k) +
			   (a1*b2*(-px06k + px16k))/(6 + k) +
			   (a1*a2*(-px07k + px17k))/(7 + k) +
			   (c2*x0*(px01k - px11k)*y1)/(1 + k) +
			   (a2*px03*(px01k - px11k)*y1)/(1 + k) +
			   (b2*px02*(-px01k + px11k)*y1)/(1 + k) +
			   (2*b2*x0*(px02k - px12k)*y1)/(2 + k) +
			   (c2*(-px02k + px12k)*y1)/(2 + k) +
			   (3*a2*px02*(-px02k + px12k)*y1)/(2 + k) +
			   (3*a2*x0*(px03k - px13k)*y1)/(3 + k) +
			   (b2*(-px03k + px13k)*y1)/(3 + k) +
			   (a2*(-px04k + px14k)*y1)/(4 + k) +
			   (c1*x0*(px01k - px11k)*y2)/(1 + k) +
			   (a1*px03*(px01k - px11k)*y2)/(1 + k) +
			   (b1*px02*(-px01k + px11k)*y2)/(1 + k) +
			   (2*b1*x0*(px02k - px12k)*y2)/(2 + k) +
			   (c1*(-px02k + px12k)*y2)/(2 + k) +
			   (3*a1*px02*(-px02k + px12k)*y2)/(2 + k) +
			   (3*a1*x0*(px03k - px13k)*y2)/(3 + k) +
			   (b1*(-px03k + px13k)*y2)/(3 + k) +
			   (a1*(-px04k + px14k)*y2)/(4 + k) +
			   ((-px01k + px11k)*y1*y2)/(1 + k);
}

double GetMultipoleIntegral(Interpolation::CubicSpline f1, Interpolation::CubicSpline f2, int k)
{
	double result = 0;
	uint32_t n = (uint32_t)f1.GetX().size();
	if (n != f2.GetX().size())
	{
		std::cout << "error in GetMultipoleIntegral: grids are not the same." << std::endl;
		return NAN;
	}
	for (uint32_t i = 0; i < n-1; i++)
	{
		result += GetMultipoleIntegral(f1,f2,i,k);
	}
	return result;
}

double GetMultipoleIntegral(Interpolation::InterpolatingFunction &f1, Interpolation::InterpolatingFunction &f2, int k, uint32_t NknotsPerInterval)
{
	double res = 0;
	vector<double> x = f1.GetX();
	vector<double> knots = Quad::GetKnots(NknotsPerInterval,x);
	vector<double> w = Quad::GetWeights(NknotsPerInterval,x);
	vector<double> f1val = f1(knots);
	vector<double> f2val = f2(knots);
	for (uint32_t i = 0; i < knots.size(); i++)
	{
		res += f1val[i] * f2val[i] * pow(knots[i],k) * w[i];
	}
	return res;
}

double GetMultipoleIntegral(Interpolation::InterpolatingFunction &f1, Interpolation::InterpolatingFunction &f2, int k, double q, uint32_t NknotsPerInterval)
{
	double res = 0;
	vector<double> x = f1.GetX();
	vector<double> knots = Quad::GetKnots(NknotsPerInterval,x);
	vector<double> w = Quad::GetWeights(NknotsPerInterval,x);
	vector<double> f1val = f1(knots);
	vector<double> f2val = f2(knots);
	for (uint32_t i = 0; i < knots.size(); i++)
	{
		res += f1val[i] * f2val[i] * SphericalBesselJ(q*knots[i],k) * w[i];
	}
	return res;
}

} /* namespace Quad */


