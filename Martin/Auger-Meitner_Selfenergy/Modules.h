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
void CreateContinuumSlaterIntegrals(vector<Interpolation::InterpolatingFunction *> &WavefunsBoundG, vector<Interpolation::InterpolatingFunction *> &WavefunsBoundF,
							vector<Interpolation::InterpolatingFunction *> &WavefunsUnboundG, vector<Interpolation::InterpolatingFunction *> &WavefunsUnboundF,
							vector<int> &indToK, const int kappa, vector<vector<vector<double> > > &Fks);

void CreateAugerMeitnerU(vector<int> &indToK, uint32_t NF, uint16_t **indices, vector<OperatorType> &UList, bool QRelativistic=true);

void CalcCNSE(const vector<WaveFunctionType> &psiCont, const vector<WaveFunctionType> &phiBound, const vector<OperatorType> &Uk, double wMin, double wMax, double dw,
              double gammaSE, const vector<vector<double> > &Fk, const vector<double> &Econt, const vector<double> &EAuger, vector<vector<Complex16> > &SelfEnergy);

static int LuaCreateContinuumSlaterIntegrals (lua_State *L);
static int LuaCreateAugerMeitnerU(lua_State *L);



#endif /* MARTIN_AUGER_MEITNER_SELFENERGY_MODULES_H_ */
