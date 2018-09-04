/*
 * DiracHamiltonian.c
 *
 *  Created on: May 4, 2018
 *      Author: brass
 */

#include "DiracHamiltonian.h"

#ifndef speedoflight
#define speedoflight 137.0
#endif

//double ddot();

double DiracH(double *basisfun,double *dualfun,double *derivfun,const double *x, int NGridPoints, const int kappa, const int Z)
{
	double *H11G, *H22F,*H12F, *H21G, result, r;
	//const double speedoflight = 137;
	unsigned int i, eins;
	eins = 1;

	H11G = (double *) malloc(NGridPoints * sizeof(double));
	H22F = (double *) malloc(NGridPoints * sizeof(double));
	H12F = (double *) malloc(NGridPoints * sizeof(double));
	H21G = (double *) malloc(NGridPoints * sizeof(double));

	if ( x[0] == 0. )
	{
		H11G[0] = 0.;
		H22F[0] = 0.;
		H12F[0] = 0.;
		H21G[0] = 0.;
		i = 1;
	}
	else
	{
		i = 0;
	}
	while ( i < NGridPoints)
	{
		r = x[i];
		H11G[i] = (-Z)/r * basisfun[i];
		H22F[i] = (-2*speedoflight*speedoflight-Z/r) * basisfun[NGridPoints+i];
		H12F[i] = kappa/r * basisfun[NGridPoints+i] - derivfun[NGridPoints+i];
		H21G[i] = kappa/r * basisfun[i] + derivfun[i];
		i++;
	}
	result = ddot((unsigned int *) &NGridPoints,dualfun,&eins,H11G,&eins) + speedoflight * ddot((unsigned int *) &NGridPoints,dualfun,&eins,H12F,&eins)
			+ddot((unsigned int *) &NGridPoints,&(dualfun[NGridPoints]),&eins,H22F,&eins) + speedoflight * ddot((unsigned int *) &NGridPoints,&(dualfun[NGridPoints]),&eins,H21G,&eins);
	free(H11G);
	free(H12F);
	free(H21G);
	free(H22F);

	return result;
}

double DiracH(vector<double> basisfun,vector<double> dualfun,vector<double> derivfun,vector<double> x, int NGridPoints, const int kappa, const int Z)
{
	double *H11G, *H22F,*H12F, *H21G, result, r;
	//const double speedoflight = 137;
	unsigned int i, eins;
	eins = 1;

	H11G = (double *) malloc(NGridPoints * sizeof(double));
	H22F = (double *) malloc(NGridPoints * sizeof(double));
	H12F = (double *) malloc(NGridPoints * sizeof(double));
	H21G = (double *) malloc(NGridPoints * sizeof(double));

	if ( x[0] == 0. )
	{
		H11G[0] = 0.;
		H22F[0] = 0.;
		H12F[0] = 0.;
		H21G[0] = 0.;
		i = 1;
	}
	else
	{
		i = 0;
	}
	while ( i < NGridPoints)
	{
		r = x[i];
		H11G[i] = (-Z)/r * basisfun[i];
		H22F[i] = (-2*speedoflight*speedoflight-Z/r) * basisfun[NGridPoints+i];
		H12F[i] = kappa/r * basisfun[NGridPoints+i] - derivfun[NGridPoints+i];
		H21G[i] = kappa/r * basisfun[i] + derivfun[i];
		i++;
	}
	result = ddot((unsigned int *) &NGridPoints,&(dualfun[0]),&eins,H11G,&eins) + speedoflight * ddot((unsigned int *) &NGridPoints,&(dualfun[0]),&eins,H12F,&eins)
			+ddot((unsigned int *) &NGridPoints,&(dualfun[NGridPoints]),&eins,H22F,&eins) + speedoflight * ddot((unsigned int *) &NGridPoints,&(dualfun[NGridPoints]),&eins,H21G,&eins);
	free(H11G);
	free(H12F);
	free(H21G);
	free(H22F);
	//printf("result = %f\n",result);
	return result;
}
//void PrintMartix2();
double **CalculateDiracMatrixElements(double **basisfuns, double **derivativefuns, double **dualfuns, const unsigned int *indices, const unsigned int Nindices, const double *x,
		unsigned int NGridPoints, const int kappa, const int Z)
{
	double **result = (double **) malloc(Nindices * sizeof(double *));

	unsigned int i, j, n, m;
	for ( i = 0; i < Nindices; i++ )
	{
		result[i] = (double *) malloc(Nindices * sizeof(double));
		n = indices[i];
		for ( j = i; j < Nindices; j++)
		{
			m = indices[j];
			result[i][j] = DiracH(basisfuns[m],dualfuns[n],derivativefuns[m], x, NGridPoints, kappa, Z);
		}
	}

	return result;

}

vector<vector<double> > CalculateDiracMatrixElements(vector<vector<double> > basisfuns, vector<vector<double> > derivativefuns, vector<vector<double> > dualfuns,
		const unsigned int *indices, const unsigned int Nindices, vector<double> x,	unsigned int NGridPoints, const int kappa, const int Z)
{
	vector<vector<double> > result;
	vector<double> help;

	unsigned int i, j, n, m;
	for ( i = 0; i < Nindices; i++ )
	{
		result.insert(result.end(),help);
		for ( j = 0; j < Nindices; j++)
		{
			result[i].insert(result[i].end(),0.); //TODO make that better
		}
	}
	for ( i = 0; i < result.size(); i++ )
	{
		n = indices[i];
		for ( j = i; j < result[i].size(); j++)
		{
			m = indices[j];
			result[i][j] = DiracH(basisfuns[m],dualfuns[n],derivativefuns[m], x, NGridPoints, kappa, Z);
		}
	}

	return result;

}

/*
double *CalculateDiracMatrixElements(double **basisfuns, double **derivativefuns, double **dualfuns, const unsigned int *indices, const unsigned int Nindices, const double *x,
		unsigned int NGridPoints, const int kappa, const int Z)
{
	double *result = (double *) malloc(Nindices*Nindices * sizeof(double));
	double h0;

	unsigned int i, j, n, m;
	for ( i = 0; i < Nindices; i++ )
	{
		n = indices[i];
		for ( j = i; j < Nindices; j++)
		{
			m = indices[j];
			h0 = DiracH(basisfuns[m],dualfuns[n],derivativefuns[m], x, NGridPoints, kappa, Z);
			result[Nindices*i+ j] = h0;
			result[Nindices*j+ i] = h0;

		}
	}


	return result;

}
*/


/*
 * for ( i = 0; i < Nindices; i++ )
	{
		n = indices[i];
		for( j = i; j < Nindices; j++ )
		{

			m = indices[j];
			result[(j-i) + i*Nindices - i*(i-1)/2] = DiracH(basisfuns[m],dualfuns[n],derivativefuns[m], x, NGridPoints, kappa, Z);
		}

	}
 */
