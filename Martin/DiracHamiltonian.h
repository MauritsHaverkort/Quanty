/*
 * DiracHamiltonian.h
 *
 *  Created on: Aug 31, 2018
 *      Author: brass
 */

#ifndef MARTIN_DIRACHAMILTONIAN_H_
#define MARTIN_DIRACHAMILTONIAN_H_

double DiracH(double *basisfun,double *dualfun,double *derivfun,const double *x, int NGridPoints, const int kappa, const int Z, const double c);

double DiracH(vector<double> basisfun,vector<double> dualfun,vector<double> derivfun,vector<double> x, int NGridPoints, const int kappa, const int Z, const double c);

double **CalculateDiracMatrixElements(double **basisfuns, double **derivativefuns, double **dualfuns, const unsigned int *indices, const unsigned int Nindices, const double *x,
		unsigned int NGridPoints, const int kappa, const int Z, const double c);

vector<vector<double> > CalculateDiracMatrixElements(vector<vector<double> > basisfuns, vector<vector<double> > derivativefuns, vector<vector<double> > dualfuns,
		const unsigned int *indices, const unsigned int Nindices, vector<double> x,	unsigned int NGridPoints, const int kappa, const int Z, const double c);

#endif /* MARTIN_DIRACHAMILTONIAN_H_ */
