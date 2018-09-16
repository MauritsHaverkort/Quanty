/*
 * AtomicU.h
 *
 *  Created on: Sep 7, 2018
 *      Author: brass
 */

#ifndef ATOMICOPERATORS_ATOMICU_H_
#define ATOMICOPERATORS_ATOMICU_H_

static int LuaOperatorNewAtomicU (lua_State *L);

uint32_t CreateOperatorAtomicURelativistic(uint32_t const NFermion,vector<vector<double> >basisfuns, vector<double> grid,vector<int> indToK
                                               ,uint16_t **indices, OperatorType *op);

void CreateOperatorAtomicURelativisticOneIntegral(uint32_t const NFermion,int const kappa1,int const kappa2, int const kappa3,int const kappa4,uint16_t const *Jindex1,
                                                  uint16_t const *Jindex2,uint16_t const *Jindex3, uint16_t const *Jindex4,double const *Fk, const unsigned int n1,
												  const unsigned int n2, const unsigned int n3, const unsigned int n4, OperatorType *const op);

uint32_t CreateOperatorAtomicUNonRel(uint32_t const NFermion,vector<Interpolation::InterpolatingFunction *> ipolfuns,vector<int> indToK
                                               ,uint16_t **indices, OperatorType *op);
void CreateOperatorAtomicUNonRelOneIntegral(uint32_t const NFermion,int const l1,int const l2, int const l3,int const l4,
		uint16_t const *Lindex1Up, uint16_t const *Lindex1Dn, uint16_t const *Lindex2Up, uint16_t const *Lindex2Dn,
		uint16_t const *Lindex3Up, uint16_t const *Lindex3Dn, uint16_t const *Lindex4Up, uint16_t const *Lindex4Dn,
		double const *Fk, const unsigned int n1, const unsigned int n2, const unsigned int n3, const unsigned int n4, OperatorType *const op);



#endif /* ATOMICOPERATORS_ATOMICU_H_ */
