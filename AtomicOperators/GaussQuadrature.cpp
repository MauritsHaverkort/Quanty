/*
 * GaussQuadrature.cpp
 *
 *  Created on: Sep 11, 2018
 *      Author: brass
 */

namespace GaussQuad {

using std::vector;

vector<double> GetKnots(uint32_t k)
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

vector<double> GetWeights(uint32_t k)
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
	if (n > 7)
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
	if (k > 7)
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

} /* namespace GaussQuad */


