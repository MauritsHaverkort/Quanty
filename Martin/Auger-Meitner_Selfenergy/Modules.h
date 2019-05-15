/*
 * Modules.h
 *
 *  Created on: May 14, 2019
 *      Author: brass
 */

#ifndef MARTIN_AUGER_MEITNER_SELFENERGY_MODULES_H_
#define MARTIN_AUGER_MEITNER_SELFENERGY_MODULES_H_

void SlaterIntegralsToList(const SSlaterIntegral *integrals, const uint32_t N, vector<double> &Fk);

void CreateContinuumSlaterIntegrals(vector<Interpolation::InterpolatingFunction *> &WavefunsBoundG, vector<Interpolation::InterpolatingFunction *> &WavefunsBoundF,
							vector<Interpolation::InterpolatingFunction *> &WavefunsUnboundG, vector<Interpolation::InterpolatingFunction *> &WavefunsUnboundF,
							vector<int> &indToK, const int kappa, vector<vector<double> > &Fks);

static int LuaCreateContinuumSlaterIntegrals (lua_State *L);



#endif /* MARTIN_AUGER_MEITNER_SELFENERGY_MODULES_H_ */
