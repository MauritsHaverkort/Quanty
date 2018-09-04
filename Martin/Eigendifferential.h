/*
 * Eigendifferential.h
 *
 *  Created on: Sep 3, 2018
 *      Author: brass
 */

#ifndef MARTIN_EIGENDIFFERENTIAL_H_
#define MARTIN_EIGENDIFFERENTIAL_H_

double EvaluateEigendifferential(double (*f)(double, int, double), double x, int kappa, double q1, double q2, int Nknots);
double WaveFunFreeParticleG(double x, int kappa, double q);
double WaveFunFreeParticleF(double x, int kappa, double q);
vector<double> GetEigenDiffOnGrid(vector<double> x, int kappa, double q1, double q2);
vector<vector<double> > CreateEigenDiffsOnGrid(vector<double> x, int kappa, double Emin, double dE, int NeigenDiffs);



#endif /* MARTIN_EIGENDIFFERENTIAL_H_ */
