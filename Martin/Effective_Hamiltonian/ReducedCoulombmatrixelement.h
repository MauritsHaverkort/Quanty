/*
 * ReducedCoulombmatrixelement.h
 *
 *  Created on: Dec 7, 2021
 *      Author: martin
 */

#ifndef MARTIN_EFFECTIVE_HAMILTONIAN_REDUCEDCOULOMBMATRIXELEMENT_H_
#define MARTIN_EFFECTIVE_HAMILTONIAN_REDUCEDCOULOMBMATRIXELEMENT_H_

class Reduced_Coulomb_matrixelement {
public:
	Reduced_Coulomb_matrixelement();
	Reduced_Coulomb_matrixelement(int L_Start, int L_End);
	virtual ~Reduced_Coulomb_matrixelement();
	double *Fk;
	int L_start, L_end;
};

#endif /* MARTIN_EFFECTIVE_HAMILTONIAN_REDUCEDCOULOMBMATRIXELEMENT_H_ */
