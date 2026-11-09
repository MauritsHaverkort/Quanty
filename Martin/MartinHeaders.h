/*
 * MartinHeaders.h
 *
 *  Created on: Aug 31, 2018
 *      Author: brass
 */

#ifndef MARTIN_MARTINHEADERS_H_
#define MARTIN_MARTINHEADERS_H_

/* Speed of light in Hartree atomic units for the continuum code, which is reachable
 * only as Martin.<Name>. It is the compiled CODATA value; the Dirac operator itself
 * (NewOperator("Dirac"), Martin.NewDirac) takes c from PhysicalConstants.cAU.value. */
#define speedoflight SpeedOfLightAU

#include <vector>
#include <map>

using std::vector;
using std::map;

#include "Typedefs_Martin.h"
#include "Quadrature.h"
#include "Eigendifferential.h"
#include "Old_SlaterIntegrals.h"
#include "ContinuumStates.h"
#include "CreateOperators.h"
#include "DiracHamiltonian.h"

#include "functions_martin.h"
#include "Math/Quad.h"



#endif /* MARTIN_MARTINHEADERS_H_ */
