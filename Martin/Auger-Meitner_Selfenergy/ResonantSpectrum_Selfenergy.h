/*
 * ResonantSpectrum_Selfenergy.h
 *
 *  Created on: May 23, 2019
 *      Author: brass
 */

#ifndef MARTIN_AUGER_MEITNER_SELFENERGY_RESONANTSPECTRUM_SELFENERGY_H_
#define MARTIN_AUGER_MEITNER_SELFENERGY_RESONANTSPECTRUM_SELFENERGY_H_

static int LuaCreateResonantSpectrum(lua_State *L);

uint32_t ResonantSpectrumWithSelfEnergy(const OperatorType OppH1, const OperatorType OppH2, const OperatorType T1, const OperatorType *T2,
										const vector<vector<OperatorType> >&U, const vector<vector<vector<double> > > &Fk, uint32_t NT2,
										const WaveFunctionType psi, TriDiagonalMatrixType **SpecTri2, SpectraType **Spec, double energy1,
										double energy2, ConvergenceSettingsType *settings, const vector<vector<double> > &ECont, double gammaSE);



#endif /* MARTIN_AUGER_MEITNER_SELFENERGY_RESONANTSPECTRUM_SELFENERGY_H_ */
