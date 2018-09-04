/*
 * CreateOperators.h
 *
 *  Created on: Aug 31, 2018
 *      Author: brass
 */

#ifndef MARTIN_CREATEOPERATORS_H_
#define MARTIN_CREATEOPERATORS_H_

void CreateOperatorCoulombRelativisticOneIntegral(uint32_t const NFermion,int const kappa1,int const kappa2, int const kappa3,int const kappa4,uint16_t const *Jindex1,
                                                  uint16_t const *Jindex2,uint16_t const *Jindex3, uint16_t const *Jindex4,double const *Fk, const unsigned int n1,
												  const unsigned int n2, const unsigned int n3, const unsigned int n4, OperatorType *const op);

uint32_t CreateOperatorFastCoulombRelativistic(uint32_t const NFermion,double **basisfuns, uint32_t NbasisFuns, const double *grid, uint32_t NGridpoints, uint32_t **tuple,
                                               uint32_t NIntegrals,uint16_t **indices, OperatorType *op);

uint32_t CreateOperatorFastCoulombRelativistic(uint32_t const NFermion,double **basisfuns, uint32_t NbasisFuns, const double *grid, uint32_t NGridpoints,vector<int> indToK
                                               ,uint16_t **indices, OperatorType *op);

uint32_t CreateOperatorFullCoulombRelativistic(uint32_t const NFermion,double **basisfuns, uint32_t NbasisFuns, const double *grid, uint32_t NGridpoints, uint32_t **tuple,
                                               uint32_t NIntegrals,uint16_t **indices, OperatorType *op);

uint32_t CreateOperatorDirac(uint32_t const NFermion,const unsigned int NBasisFuns, double **basisfuns, double **derivfuns, const double *x, uint32_t NGridpoints,
		SKappaGroup *groups,uint32_t Ngroups,uint32_t **indices, int Z, OperatorType *op);

double* CreateBaseForContinuum( int NBasisFuns,double **basisfuns,double **derivfuns,double *x,uint32_t NGridpoints,SKappaGroup group,int Ncont,
                                                                double **contfuns, double w);

double *CreateBaseForContinuum( int NBasisFuns, vector<vector<double> > basisfuns,vector<vector<double> > derivfuns,vector<double> x,
		uint32_t NGridpoints,SKappaGroup group,int Ncont, vector<vector<double> > &contfuns, double w);

uint32_t CreateOperatorSingleCoulombRelativistic(uint32_t const NFermion,double **basisfuns, uint32_t NbasisFuns, const double *grid, uint32_t NGridpoints, uint32_t **tuple,
                                               uint32_t NIntegrals,uint16_t **indices, OperatorType *op);

#endif /* MARTIN_CREATEOPERATORS_H_ */
