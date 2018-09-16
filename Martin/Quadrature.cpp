/*
 * Quadrature.cpp
 *
 *  Created on: Jun 28, 2018
 *      Author: brass
 */

#include "Quadrature.h"
#include <iostream>

using std::vector;
using std::cout;
using std::endl;

Quadrature::Quadrature() {
	weights = vector<vector<double> >(7);
	weights[0].push_back(2.);
	weights[1].push_back(1.);
	weights[1].push_back(1.);
	weights[2].push_back(0.555556);
	weights[2].push_back(0.888889);
	weights[2].push_back(0.555556);
	weights[3].push_back(0.347855);
	weights[3].push_back(0.652145);
	weights[3].push_back(0.652145);
	weights[3].push_back(0.347855);
	weights[4].push_back(0.236927);
	weights[4].push_back(0.478629);
	weights[4].push_back(0.568889);
	weights[4].push_back(0.478629);
	weights[4].push_back(0.236927);
	weights[5].push_back(0.171324);
	weights[5].push_back(0.360762);
	weights[5].push_back(0.467914);
	weights[5].push_back(0.467914);
	weights[5].push_back(0.360762);
	weights[5].push_back(0.171324);
	weights[6].push_back(0.129485);
	weights[6].push_back(0.279705);
	weights[6].push_back(0.38183);
	weights[6].push_back(0.417959);
	weights[6].push_back(0.38183);
	weights[6].push_back(0.279705);
	weights[6].push_back(0.129485);

	knots = vector<vector<double> >(7);
	knots[0].push_back(0.);
	knots[1].push_back(-0.57735);
	knots[1].push_back( 0.57735);
	knots[2].push_back(-0.774597);
	knots[2].push_back(0.);
	knots[2].push_back(0.774597);
	knots[3].push_back(-0.861136);
	knots[3].push_back(-0.339981);
	knots[3].push_back(0.339981);
	knots[3].push_back(0.861136);
	knots[4].push_back(-0.90618);
	knots[4].push_back(-0.538469);
	knots[4].push_back(0.);
	knots[4].push_back(0.538469);
	knots[4].push_back(0.90618);
	knots[5].push_back(-0.93247);
	knots[5].push_back(-0.661209);
	knots[5].push_back(-0.238619);
	knots[5].push_back(0.238619);
	knots[5].push_back(0.661209);
	knots[5].push_back(0.93247);
	knots[6].push_back(-0.949108);
	knots[6].push_back(-0.741531);
	knots[6].push_back(-0.405845);
	knots[6].push_back(0.);
	knots[6].push_back(0.405845);
	knots[6].push_back(0.741531);
	knots[6].push_back(0.949108);

}

Quadrature::~Quadrature() {
	// TODO Auto-generated destructor stub
}

vector<double> Quadrature::GetKnots(int numberOfKnots)
{
	if ((numberOfKnots <= NKnotsMax) && (numberOfKnots > 0))
	{
		return knots[numberOfKnots-1];
	}
	else
	{
		cout << "error: in Quadrature::GetKnots: numberOfKnots too big or smaller than 1" << endl;
		return {};
	}
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
vector<double> Quadrature::GetKnots(int numberOfKnots,double a, double b)
{
	if ((numberOfKnots <= NKnotsMax) && (numberOfKnots > 0))
	{
		return TrafoXpoints(knots[numberOfKnots-1],a,b);
	}
	else
	{
		cout << "error: in Quadrature::GetKnots: numberOfKnots too big or smaller than 1" << endl;
		return {};
	}
}

double Quadrature::Integrate(vector<double> f)
{
	unsigned long n = f.size();
	if (n > NKnotsMax)
	{
		cout << "error: in Quadrature::Integrate: numberOfKnots too big" << endl;
		return 0;
	}
	else
	{
		double result = 0;
		for (unsigned long i = 0; i < n; i++)
		{
			result += weights[n-1][i]*f[i];
		}
		return result;
	}
}

double Quadrature::Integrate(vector<double> f,double a, double b)
{
	if ( a == b)
		return 0.;
	unsigned long n = f.size();
	if (n > NKnotsMax)
	{
		cout << "error: in Quadrature::Integrate: numberOfKnots too big" << endl;
		return 0;
	}
	else
	{
		double result = 0;
		for (unsigned long i = 0; i < n; i++)
		{
			result += (b-a)/2*weights[n-1][i]*f[i];
		}
		return result;
	}
}
