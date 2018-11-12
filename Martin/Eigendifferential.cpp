/*
 * Eigendifferential.cpp
 *
 *  Created on: Sep 3, 2018
 *      Author: brass
 */

#include "Quadrature.h"
#include "Eigendifferential.h"

double EvaluateEigendifferential(double (*f)(double, int, double), double x, int kappa, double q1, double q2, int Nknots=4)
{
	Quadrature quad;
	vector<double> knots = quad.GetKnots(Nknots,q1,q2);
	vector<double> fval = vector<double>(knots.size());
	for (int i = 0; i < fval.size(); i++)
	{
		fval[i] = f(x,kappa,knots[i]);
	}
	return quad.Integrate(fval,q1,q2)/sqrt(q2-q1);
}

double WaveFunFreeParticleG(double x, int kappa, double q)
{
	double E = speedoflight * sqrt(q*q + speedoflight*speedoflight);
	double cf = Sign(kappa) * speedoflight * q / (E + speedoflight*speedoflight);
	double a = sqrt(2 * q*q / M_PI / (1 + cf*cf));
	int l;
	if (kappa > 0)
	{
		l = kappa;
	}
	else
	{
		l = -kappa - 1;
	}
	return a * SphericalBesselJ(q*x,l);
}

double WaveFunFreeParticleNonRelativistic(double x, int l, double q)
{
	return SphericalBesselJ(q*x,l);
}

double WaveFunFreeParticleF(double x, int kappa, double q)
{
	double E = speedoflight * sqrt(q*q + speedoflight*speedoflight);
	double cf = Sign(kappa) * speedoflight * q / (E + speedoflight*speedoflight);
	double a = sqrt(2 * q*q / M_PI / (1 + cf*cf));
	int l;
	if (kappa > 0)
	{
		l = kappa-1;
	}
	else
	{
		l = -kappa;
	}
	return a * cf * SphericalBesselJ(q*x,l);
}

vector<double> GetEigenDiffOnGrid(vector<double> x, int kappa, double q1, double q2)
{
	int NGridPoints = (int) x.size();
	vector<double> eigendiff = vector<double>(2*NGridPoints);
	for (int i = 0; i < NGridPoints; i++)
	{
		eigendiff[i            ] = x[i] * EvaluateEigendifferential(&WaveFunFreeParticleG,x[i],kappa,q1,q2);
		eigendiff[i+NGridPoints] = x[i] * EvaluateEigendifferential(&WaveFunFreeParticleF,x[i],kappa,q1,q2);
	}
	return eigendiff;
}

// output: eigendiff - must be initialized !
void GetEigenDiffOnGrid(const vector<double> &x, int kappa, double q1, double q2,double (*f)(double, int, double),vector<double> &eigendiff)
{
	if (x.size() != eigendiff.size())
	{
		printf("eigendiff not properly initialized");
	}
	else
	{
		for (int i = 0; i < eigendiff.size(); i++)
		{
			eigendiff[i] = x[i] * EvaluateEigendifferential(f,x[i],kappa,q1,q2);
		}
	}
}

vector<vector<double> > CreateEigenDiffsOnGrid(const vector<double> &x, int kappa, double Emin, double dE, int NeigenDiffs)
{
	vector<vector<double> > eigendiffs = vector<vector<double> >(NeigenDiffs);
	for (int i = 0; i < NeigenDiffs; i++)
	{
		double q1 = sqrt(2*(Emin +  i    * dE));//sqrt((Emin + i     * dE)*(Emin + i     * dE)/speedoflight/speedoflight - speedoflight*speedoflight);
		double q2 = sqrt(2*(Emin + (i+1) * dE));//sqrt((Emin + (i+1) * dE)*(Emin + (i+1) * dE)/speedoflight/speedoflight - speedoflight*speedoflight);

		eigendiffs[i] = GetEigenDiffOnGrid(x,kappa,q1,q2);
	}
	return eigendiffs;
}

// Non Relativistic Eigen Differentials
void CreateEigenDiffs(const vector<double> &x, int l, double Emin, double dE, int NeigenDiffs, vector<Interpolation::CubicSpline> &eigendifferentials)
{
	vector<double> vals = vector<double>(x.size());
	eigendifferentials = vector<Interpolation::CubicSpline>(NeigenDiffs);
	for (int i = 0; i < NeigenDiffs; i++)
	{
		double q1 = sqrt(2*(Emin +  i    * dE));//sqrt((Emin + i     * dE)*(Emin + i     * dE)/speedoflight/speedoflight - speedoflight*speedoflight);
		double q2 = sqrt(2*(Emin + (i+1) * dE));//sqrt((Emin + (i+1) * dE)*(Emin + (i+1) * dE)/speedoflight/speedoflight - speedoflight*speedoflight);

		GetEigenDiffOnGrid(x,l,q1,q2,&WaveFunFreeParticleNonRelativistic,vals);
		eigendifferentials[i].Create(x,vals);
	}
}
// Relativistic Eigen Differentials
void CreateEigenDiffs(const vector<double> &x, int kappa, double Emin, double dE, int NeigenDiffs, vector<Interpolation::CubicSpline> &eigendifferentialsG, vector<Interpolation::CubicSpline> &eigendifferentialsF)
{
	vector<double> vals = vector<double>(x.size());
	eigendifferentialsG = vector<Interpolation::CubicSpline>(NeigenDiffs);
	eigendifferentialsF = vector<Interpolation::CubicSpline>(NeigenDiffs);
	for (int i = 0; i < NeigenDiffs; i++)
	{
		double q1 = sqrt(2*(Emin +  i    * dE));//sqrt((Emin + i     * dE)*(Emin + i     * dE)/speedoflight/speedoflight - speedoflight*speedoflight);
		double q2 = sqrt(2*(Emin + (i+1) * dE));//sqrt((Emin + (i+1) * dE)*(Emin + (i+1) * dE)/speedoflight/speedoflight - speedoflight*speedoflight);

		GetEigenDiffOnGrid(x,kappa,q1,q2,&WaveFunFreeParticleG,vals);
		eigendifferentialsG[i].Create(x,vals);

		GetEigenDiffOnGrid(x,kappa,q1,q2,&WaveFunFreeParticleF,vals);
		eigendifferentialsF[i].Create(x,vals);
	}
}
