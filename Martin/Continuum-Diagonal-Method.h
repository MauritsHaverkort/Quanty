/*
 * Continuum-Diagonal-Method.h
 *
 *  Created on: Apr 1, 2019
 *      Author: brass
 */

#ifndef MARTIN_CONTINUUM_DIAGONAL_METHOD_H_
#define MARTIN_CONTINUUM_DIAGONAL_METHOD_H_

/*
 * Continuum-Diagonal-Method.cpp
 *
 *  Created on: Mar 12, 2019
 *      Author: brass
 */

typedef struct mypoint
{
	double w;
	Complex16 y;
}mypoint;

void PsiOppListOppListPsi(WaveFunctionType *psi, vector<OperatorType*> U, CompactMatrixType *M);

void PsiOppListOperatorOppListPsi(WaveFunctionType *psi, vector<OperatorType*> U, OperatorType *H, CompactMatrixType *M);

uint32_t GetNumberOfPolesBelowCutoff(const std::vector<double> &poles, double cutoff);

void CalculateResidues(const std::vector<double> &v, CompactMatrixType *M, uint32_t Nresidues, std::vector<double> &residues);

void SelectDominantResiduesAndPoles(const std::vector<double> &poles, const std::vector<double> &residues, double eps, std::vector<double> &selectedPoles, std::vector<double> &selectedResidues);

Complex16 Resolvent(Complex16 z, const std::vector<double> &poles, const std::vector<double> &residues);

Complex16 Resolvent(Complex16 z, const std::vector<double> &poles, const std::vector<double> &residues,const vector<Complex16> &selfenergy);

void PlotSpectrumFromPoles(double eps,double gamma, double shift, const std::vector<double> &poles, const std::vector<double> &v, CompactMatrixType *M, std::vector<mypoint> &spectrum);

void PlotSpectrumFromPoles(double eps,double gamma, double shift, const std::vector<double> &poles, const std::vector<double> &v, CompactMatrixType *M, std::vector<mypoint> &spectrum,const vector<vector<Complex16> > &selfenergy);

void lua_push(lua_State * L, const mypoint p);

void AugerMeitnerSelfEnergy(vector<Interpolation::InterpolatingFunction *> WavefunsBoundG,vector<Interpolation::InterpolatingFunction *> WavefunsBoundF,
							vector<Interpolation::InterpolatingFunction *> WavefunsUnboundG,vector<Interpolation::InterpolatingFunction *> WavefunsUnboundF,
							vector<int> indToK, uint16_t **indices, uint32_t NF, OperatorType *H, WaveFunctionType * psi, double eps, double gamma,
							vector<double> &EnergiesUnbound, vector<mypoint> &spectrum, const uint32_t NPointsToEvaluateInterpolation = 0);

void CreateAugerMeitnerOperators(SSlaterIntegral *integrals, uint32_t NIntegrals,uint32_t NF, uint16_t **indices, vector<OperatorType> &UList);

void CalculateEnergyIndependentSelfEnergy(const OperatorType &H, const WaveFunctionType  &psi, const vector<OperatorType> &UList, CompactMatrixType &T, vector<double> &poles);

void FlattenSlaterIntegrals(const SSlaterIntegral *integrals, const uint32_t N, vector<double> &Fk);

void Transpose(const vector<vector<double> > &A, vector<vector<double> > &AT);

void InterpolateSlaterIntegrals(const int NIntegrals, SSlaterIntegral *integrals, const int n,double **kernels,
								vector<vector<double> > &WavefunsBound, vector<vector<double> > &dualfunsBound,
								const vector<vector<double> > &WavefunsUnbound,const vector<vector<double> > &dualfunsUnbound,
								const vector<double> &EnergiesUnbound, vector<Interpolation::CubicSpline> &FkInterpol);

#endif /* MARTIN_CONTINUUM_DIAGONAL_METHOD_H_ */
