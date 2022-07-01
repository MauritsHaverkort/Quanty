/*
 * ReducedCoulombmatrixelement.cpp
 *
 *  Created on: Dec 7, 2021
 *      Author: martin
 */

#include "ReducedCoulombmatrixelement.h"

Reduced_Coulomb_matrixelement::Reduced_Coulomb_matrixelement() {
	Fk = 0;
	L_start = -1;
	L_end = -2;

}

Reduced_Coulomb_matrixelement::Reduced_Coulomb_matrixelement(int L_Start, int L_End) {
	L_start = L_Start;
	L_end   = L_End;
	if (L_Start <= L_End)
		Fk = (double*) malloc(((L_end-L_start)/2+1)*sizeof(double));
	else
		Fk = 0;
}

Reduced_Coulomb_matrixelement::~Reduced_Coulomb_matrixelement() {
	// TODO Auto-generated destructor stub

}

