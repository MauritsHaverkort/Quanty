/*
 * Old_SlaterIntegrals.h
 *
 *  Created on: Sep 19, 2018
 *      Author: brass
 */

#ifndef MARTIN_OLD_SLATERINTEGRALS_H_
#define MARTIN_OLD_SLATERINTEGRALS_H_

SSlaterIntegral *GetSlaterIntegrals(double **basisfuns, unsigned int NBasisFuns, const double *x,
                                    unsigned int NGridPoints,const int **tupel, unsigned int NIntegrals);

SSlaterIntegral *GetSlaterIntegrals(double **basisfuns, unsigned int NBasisFuns, const double *x,
                                    unsigned int NGridPoints,vector<int> indToK, unsigned int &NIntegrals);

SSlaterIntegral *InitializeIntegrals(const int **tupel, unsigned int NIntegrals,
                                     unsigned int *kMax);

void CalculateSlaterIntegrals(const int NIntegrals, SSlaterIntegral *integrals, const int dim,
                              double **kernels, double **basisfuns, double **dualfuns);

double **CreateDualFunctions(const unsigned int NBasisFuns, const unsigned int NGridPoints,
                             double **basisfuns, const double *x);

void Dual(const int n, const double *f, const double *weights, double *fdual);



#endif /* MARTIN_OLD_SLATERINTEGRALS_H_ */
