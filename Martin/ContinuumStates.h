/*
 * ContinuumStates.h
 *
 *  Created on: Aug 31, 2018
 *      Author: brass
 */

#ifndef MARTIN_CONTINUUMSTATES_H_
#define MARTIN_CONTINUUMSTATES_H_

double Gamma(double const x); // Martin Brass copied this from Martin Sundermann and Martin Rotter edit by MS

double Binomial(double x, double y);

double Phi(double n, double a, double b);

double EnergyDO(int n, int k, int l, double w);

double NormDO(int n, int k, int l, double w);

double C1DO(int n, int k, int l, double w);

double C2DO(int n, int k, int l, double w);

double EKinDO(int n, int m, int k, int l, double w);

double **CalculateContinuumBoundDiracElements(double **basisfuns, double **derivativefuns, double **dualcont, const int NcontFuns, const unsigned int *indices,
		const unsigned int Nindices, const double *x, unsigned int NGridPoints, const int kappa, const int Z);

vector<vector<double> > CalculateContinuumBoundDiracElements(vector<vector<double> > basisfuns, vector<vector<double> > derivativefuns,
		vector<vector<double> > dualcont, const int NcontFuns, const unsigned int *indices,
		const unsigned int Nindices, vector<double> x, unsigned int NGridPoints, const int kappa, const int Z);

double *KineticEnergyCB(double **basisfuns, double **derivativefuns, double **dualfuns,double **dualcont, const int NcontFuns, const unsigned int *indices,
		const unsigned int Nindices, const double *x, unsigned int NGridPoints, const int kappa,double w);

vector<double> KineticEnergyCB(vector<vector<double> > basisfuns, vector<vector<double> > derivativefuns, vector<vector<double> > dualfuns,
		vector<vector<double> > dualcont, const int NcontFuns, const unsigned int *indices,
		const unsigned int Nindices, vector<double> x, unsigned int NGridPoints, const int kappa,double w);

double *MetricCB(double **basisfuns, double **dualfuns, const unsigned int *indices,
		const unsigned int Nindices, double **dualcont, const int NcontFuns, unsigned int NGridPoints);

vector<double> MetricCB(vector<vector<double> > basisfuns, vector<vector<double> > dualfuns, const unsigned int *indices,
		const unsigned int Nindices, vector<vector<double> > dualcont, const int NcontFuns, unsigned int NGridPoints);

double *SubMatrix(double *M, int dim, int n1, int n2, int m1, int m2 );

vector<double> SubMatrix(vector<double> M, int dim, int n1, int n2, int m1, int m2 );

void vAdd(const double *a,const double aleph, double* b,const int n);

void vAdd(vector<double> a,const double aleph, vector<double> &b);

void svMul(const double aleph, double *a, const int n);

void svMul(const double aleph, vector<double> &a);

void CalcContinuumWFs(double **basisfuns, const unsigned int *indices, const unsigned int Nindices, double **contfuns,
		const int Ncont,const int NGridpoints, double *U, double *R);

void CalcContinuumWFs(vector<vector<double> > basisfuns, const unsigned int *indices, const unsigned int Nindices,
		vector<vector<double> > &contfuns, const int Ncont,const int NGridpoints, vector<double> U, vector<double> R);

double *CreateBaseForContinuum(vector<vector<double> > basisfuns,vector<vector<double> > derivfuns,vector<double> x,
		SKappaGroup group, vector<vector<double> > &contfuns, vector<vector<double> > h0Cont);

vector<double> CreateEKinCB(vector<vector<double> > basisfuns, vector<vector<double> > derivativefuns, vector<vector<double> > dualfuns,
		vector<vector<double> > dualcont, const unsigned int *indices, const unsigned int Nindices, vector<double> x, const int kappa,vector<vector<double> > h0cont);

vector<vector<double> > CreateEigenDifferentialEnergies(uint32_t n, double Emin, double dE);

#endif /* MARTIN_CONTINUUMSTATES_H_ */
