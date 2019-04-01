/*
 * Continuum-Diagonal-Method.cpp
 *
 *  Created on: Mar 12, 2019
 *      Author: brass
 */

#include "Continuum-Diagonal-Method.h"

void PsiOppListOppListPsi(WaveFunctionType *psi, vector<OperatorType*> U, CompactMatrixType *M)
{
	OperatorType UU;
	double Re, Im;
	for (uint32_t i = 0; i < M->Ni; i++)
	{
		for (uint32_t j = i; j < M->Nj; j++)
		{
			InitOperatorForMult(&UU,2,U[i],U[j]);
			OperatorMultiply2NormalOrder(&UU, *U[j], *U[i], 'D', 'N', &QDetNotOrderedFcFaBcBa);
			PsiOperatorPsi(UU,*psi,&Re,&Im);
			M->VRe[i*M->Nj + j] = Re;
			M->VRe[j*M->Nj + i] = Re;
			FreeOperator(&UU);
		}
	}
}

void PsiOppListOperatorOppListPsi(WaveFunctionType *psi, vector<OperatorType*> U, OperatorType *H, CompactMatrixType *M)
{
	OperatorType UHU, HU;
	double Re, Im;
	if(M->QComplex)
	{
		for (uint32_t i = 0; i < M->Ni; i++)
		{
			InitOperatorForMult(&HU,2,U[i],H);
			OperatorMultiply2NormalOrder(&HU,*H,*U[i],'N','N',&QDetNotOrderedFcFaBcBa);
			for (uint32_t j = i; j < M->Nj; j++)
			{
				//InitOperatorForMult(&UHU,3,U[i],U[j],H);
				//OperatorMultiply3NormalOrder(&UHU, *U[j], *H, *U[i], 'D', 'N', 'N', &QDetNotOrderedFcFaBcBa);
				InitOperatorForMult(&UHU,2,U[j],HU);
				OperatorMultiply2NormalOrder(&UHU,*U[j],HU,'D','N',&QDetNotOrderedFcFaBcBa);
				PsiOperatorPsi(UHU,*psi,&Re,&Im);
				M->V[i*M->Nj + j].Re = Re;
				M->V[i*M->Nj + j].Im = Im;
				M->V[j*M->Nj + i] = M->V[i*M->Nj + j];
				FreeOperator(&UHU);
			}
			FreeOperator(&HU);
		}
	}
	else
	{
#pragma omp parallel for
		for (uint32_t i = 0; i < M->Ni; i++)
		{
			OperatorType hu;
			InitOperatorForMult(&hu,2,H,U[i]);
			OperatorMultiply2NormalOrder(&hu,*H,*U[i],'N','N',&QDetNotOrderedFcFaBcBa);
			for (uint32_t j = i; j < M->Nj; j++)
			{
				//InitOperatorForMult(&UHU,3,U[i],U[j],H);
				//OperatorMultiply3NormalOrder(&UHU, *U[j], *H, *U[i], 'D', 'N', 'N', &QDetNotOrderedFcFaBcBa);
				OperatorType uhu;
				InitOperatorForMult(&uhu,2,U[j],&hu);
				OperatorMultiply2NormalOrder(&uhu,*U[j],hu,'D','N',&QDetNotOrderedFcFaBcBa);
				double re, im;
				PsiOperatorPsi(uhu,*psi,&re,&im);
				M->VRe[i*M->Nj + j] = re;
				M->VRe[j*M->Nj + i] = re;
				FreeOperator(&uhu);
			}
			FreeOperator(&hu);
		}
	}

}

uint32_t GetNumberOfPolesBelowCutoff(const std::vector<double> &poles, double cutoff)
{
	uint32_t i = 0;
	while (i < (uint32_t) poles.size() && poles[i] < cutoff)
	{
		i++;
	}
	return i;
}

void CalculateResidues(const std::vector<double> &v, CompactMatrixType *M, uint32_t Nresidues, std::vector<double> &residues)
{
	if(M->Ni < Nresidues || M->Nj < (uint32_t) v.size())
	{
		printf("error in CalculateResidues: CompactMatrix has wrong dimensions: Ni = %d Nj = %d\n",M->Ni,M->Nj);
		fflush(stdout);
		return;
	}
	if(M->QComplex)
	{
		printf("error in CalculateResidues: CompactMatrix is complex. I should implement that\n");
		fflush(stdout);
		return;
	}
	residues = std::vector<double>(Nresidues);
	for (uint32_t i = 0; i < Nresidues; i++)
	{
		residues[i] = 0;
		for (uint32_t j = 0; j < (uint32_t) v.size(); j++)
		{
			residues[i] += M->VRe[i*M->Nj + j] * v[j];
		}
		residues[i] *= residues[i];
	}
}

void SelectDominantResiduesAndPoles(const std::vector<double> &poles, const std::vector<double> &residues, double eps, std::vector<double> &selectedPoles, std::vector<double> &selectedResidues)
{
	if (poles.size() < residues.size())
	{
		printf("error in SelectDominantResiduesAndPoles: poles.size() = %d < %d = residues.size()\n",poles.size(),residues.size());
		fflush(stdout);
		return;
	}
	for (uint32_t i = 0; i < (uint32_t) residues.size(); i++)
	{
		if (residues[i] > eps)
		{
			selectedPoles.push_back(poles[i]);
			selectedResidues.push_back(residues[i]);
		}
	}
}

Complex16 Resolvent(Complex16 z, const std::vector<double> &poles, const std::vector<double> &residues)
{
	Complex16 result;
	result.Re = 0;
	result.Im = 0;
	if (poles.size() < residues.size())
	{
		printf("error in Resolvent: poles.size() = %d < %d = residues.size()",poles.size(),residues.size());
		fflush(stdout);
		return result;
	}
	for (uint32_t i = 0; i < (uint32_t) residues.size(); i++)
	{
		result += residues[i]/(z - poles[i]);
	}
	return result;
}

Complex16 Resolvent(Complex16 z, const std::vector<double> &poles, const std::vector<double> &residues,const vector<Complex16> &selfenergy)
{
	Complex16 result;
	result.Re = 0;
	result.Im = 0;
	if (poles.size() < residues.size() || poles.size() != selfenergy.size())
	{
		printf("error in Resolvent: poles.size() = %d < %d = residues.size()",poles.size(),residues.size());
		fflush(stdout);
		return result;
	}
	for (uint32_t i = 0; i < (uint32_t) residues.size(); i++)
	{
		result += residues[i]/(z - poles[i] - selfenergy[i]);
	}
	return result;
}

void PlotSpectrumFromPoles(double eps,double gamma, double shift, const std::vector<double> &poles, const std::vector<double> &v, CompactMatrixType *M, std::vector<mypoint> &spectrum)
{

	double wMax = spectrum[spectrum.size()-1].w;
	uint32_t Npoles = GetNumberOfPolesBelowCutoff(poles, wMax);
	std::vector<double> residues, selectedPoles, selectedResidues;
	CalculateResidues(v,M,Npoles,residues);
	SelectDominantResiduesAndPoles(poles,residues,eps,selectedPoles,selectedResidues);
#pragma omp parallel for
	for (uint32_t i = 0; i < (uint32_t) spectrum.size(); i++)
	{
		Complex16 z;
		z.Re = spectrum[i].w - shift;
		z.Im = gamma/2;
		spectrum[i].y = Resolvent(z,selectedPoles,selectedResidues);
	}

}

void CalcSpectrumFromPolesAdd(double eps,double gamma, double shift, const std::vector<double> &poles, const std::vector<double> &v, CompactMatrixType *M, std::vector<mypoint> &spectrum)
{

	double wMax = spectrum[spectrum.size()-1].w;
	uint32_t Npoles = GetNumberOfPolesBelowCutoff(poles, wMax);
	std::vector<double> residues, selectedPoles, selectedResidues;
	CalculateResidues(v,M,Npoles,residues);
	SelectDominantResiduesAndPoles(poles,residues,eps,selectedPoles,selectedResidues);
#pragma omp parallel for
	for (uint32_t i = 0; i < (uint32_t) spectrum.size(); i++)
	{
		Complex16 z;
		z.Re = spectrum[i].w - shift;
		z.Im = gamma/2;
		spectrum[i].y += Resolvent(z,selectedPoles,selectedResidues);
	}

}

void PlotSpectrumFromPoles(double eps,double gamma, double shift, const std::vector<double> &poles, const std::vector<double> &v, CompactMatrixType *M, std::vector<mypoint> &spectrum,const vector<vector<Complex16> > &selfenergy)
{

	double wMax = spectrum[spectrum.size()-1].w;
	uint32_t Npoles = GetNumberOfPolesBelowCutoff(poles, wMax);
	std::vector<double> residues, selectedPoles, selectedResidues;
	CalculateResidues(v,M,Npoles,residues);
	SelectDominantResiduesAndPoles(poles,residues,eps,selectedPoles,selectedResidues);

#pragma omp parallel for
	for (uint32_t i = 0; i < (uint32_t) spectrum.size(); i++)
	{
		Complex16 z;
		z.Re = spectrum[i].w - shift;
		z.Im = gamma/2;
		spectrum[i].y = Resolvent(z,selectedPoles,selectedResidues,selfenergy[i]);
	}

}

void lua_push(lua_State * L, const mypoint p)
{
	lua_createtable(L,0,0);
	lua_push(L,p.w);
	lua_rawseti(L,-2,1);
	Complex16 *z = (Complex16 *) lua_newuserdata(L, sizeof(Complex16));
	luaL_setmetatable(L, "Complex_Type");
	*z = p.y;
	lua_rawseti(L,-2,2);
}

void AugerMeitnerSelfEnergy(vector<Interpolation::InterpolatingFunction *> WavefunsBoundG,vector<Interpolation::InterpolatingFunction *> WavefunsBoundF,
							vector<Interpolation::InterpolatingFunction *> WavefunsUnboundG,vector<Interpolation::InterpolatingFunction *> WavefunsUnboundF,
							vector<int> indToK, uint16_t **indices, uint32_t NF, OperatorType *H, WaveFunctionType * psi, double eps, double gamma,
							vector<double> &EnergiesUnbound, vector<mypoint> &spectrum)
{
	SSlaterIntegral *si;
	uint32_t kMax, n, NIntegrals, NBound, NUnbound;
	double **kernels;
	vector<vector<double> > dualfunsBound, dualfunsUnbound;
	vector<OperatorType> UList;
	vector<double> poles, Fk;
	CompactMatrixType T;


	vector<double> x = WavefunsBoundG[0]->GetX();
	vector<vector<double> > WavefunsBound = PutRelativisticSplineOnGrid(WavefunsBoundG,WavefunsBoundF,x,0);
	vector<vector<double> > WavefunsUnbound = PutRelativisticSplineOnGrid(WavefunsUnboundG,WavefunsUnboundF,x,0);
	NUnbound = (uint32_t) WavefunsUnbound.size();
	NBound = (uint32_t) WavefunsBound.size();

	unsigned int NGridPoints = (unsigned int)x.size();
	n = 2* NGridPoints;
	CreateDualFunctions(WavefunsBound,x,dualfunsBound);
	CreateDualFunctions(WavefunsUnbound,x,dualfunsUnbound);
	si = InitializeIntegrals(indToK, NIntegrals,&kMax,true);
	kernels = CreateAllKernels(kMax, NGridPoints, &x[0]);
	CreateAugerMeitnerOperators(si,NIntegrals,NF,indices,UList);
	Fk = vector<double>(UList.size());
	CalculateEnergyIndependentSelfEnergy(*H,*psi,UList,T,poles);
	WavefunsBound.push_back(WavefunsUnbound[0]);
	dualfunsBound.push_back(dualfunsUnbound[0]);

	//PrintArray(&poles[0],poles.size());
	//printf("NUnbound = %d\n",NUnbound);
	for (uint32_t i = 0; i < NUnbound; i++)
	{
		WavefunsBound[NBound] = WavefunsUnbound[i];
		dualfunsBound[NBound] = dualfunsUnbound[i];
		CalculateSlaterIntegrals(NIntegrals, si,n,kernels,WavefunsBound,dualfunsBound);
		FlattenSlaterIntegrals(si,NIntegrals,Fk);
		CalcSpectrumFromPolesAdd(eps,gamma,EnergiesUnbound[i],poles,Fk,&T,spectrum);
	}

	FreeCompactMatrix(&T);
}

void FlattenSlaterIntegrals(const SSlaterIntegral *integrals, const uint32_t N, vector<double> &Fk)
{
	uint32_t count = 0;
	for (uint32_t i = 0; i < N; i++)
	{
		for (uint32_t j = 0; j < integrals[i].Nval; j++)
		{
			Fk[count] = integrals[i].val[j];
			count++;
		}
	}
}

void CalculateOverlapMatrix(const vector<WaveFunctionType> &psi, CompactMatrixType &S)
{
	S.QComplex = 0;
	S.Ni = (uint32_t) psi.size();
	S.Nj = S.Ni;
	InitCompactMatrix(&S);
	for (uint32_t i = 0; i < S.Ni; i++)
	{
		for (uint32_t j = i; j < S.Nj; j++)
		{
			double muelleimer;
			WaveFunctionDotProduct(psi[i],psi[j],&S.VRe[i*S.Nj + j],&muelleimer);
			S.VRe[j*S.Nj + i] = S.VRe[i*S.Nj + j];
		}
	}
}

void CalculateHermitianMatrixElements(const vector<WaveFunctionType> &basis,const OperatorType &H, CompactMatrixType &M)
{
	M.QComplex = H.QComplex;
	M.Ni = (uint32_t) basis.size();
	M.Nj = M.Ni;
	uint32_t i = 0;
	while (!M.QComplex && i < M.Ni)
	{
		M.QComplex = basis[i].QComplex;
		i++;
	}
	InitCompactMatrix(&M);
	for (uint32_t i = 0; i < M.Ni; i++)
	{
		WaveFunctionType Hpsi;
		OperatorPsi(H,basis[i],&Hpsi);
		for (uint32_t j = i; j < M.Nj; j++)
		{
			if (M.QComplex)
			{

			}
			else
			{
				double muelleimer;
				WaveFunctionDotProduct(basis[j],Hpsi,&M.VRe[i*M.Nj+j],&muelleimer);
				M.VRe[j*M.Nj+i] = M.VRe[i*M.Nj+j];
			}
		}

		WaveFunctionFree(&Hpsi);
	}
}

void CompactMatrixMultiplyDiagonal(const vector<double> &D, const CompactMatrixType &A, CompactMatrixType &B)
{
	B.Ni = A.Ni;
	B.Nj = A.Nj;
	B.QComplex = A.QComplex;
	InitCompactMatrix(&B);
	if (B.QComplex)
	{
		for (uint32_t i = 0; i < A.Ni; i++)
		{
			for (uint32_t j = 0; j < A.Nj; j++)
			{
				B.V[i*B.Nj + j] = D[i] * A.V[i*B.Nj + j];
			}
		}
	}
	else
	{
		for (uint32_t i = 0; i < A.Ni; i++)
		{
			for (uint32_t j = 0; j < A.Nj; j++)
			{
				B.VRe[i*B.Nj + j] = D[i] * A.VRe[i*B.Nj + j];
			}
		}
	}

}

template <typename T>
void Map(T (*f)(T), vector<T> &v)
{
	for (uint32_t i = 0; i < v.size(); i++)
	{
		v[i] = f(v[i]);
	}
}

double Inv(double x)
{
	return 1. / x;
}

void TransformMatrixToONB(const CompactMatrixType &S,const CompactMatrixType &H, CompactMatrixType &S12, CompactMatrixType &HTrafo)
{
	vector<double> eigenvals = vector<double>(S.Ni);
	CompactMatrixType T,temp, S12Inv;
	CompactMatrixEigensystem(&S,&T,&eigenvals[0]);
	Map(sqrt,eigenvals);
	CompactMatrixMultiplyDiagonal(eigenvals,T,temp);
	CompactMatrixMultiply('T','N',&T,&temp,&S12); // TODO check if eigenvals row- or coll-major
	FreeCompactMatrix(&temp);
	Map(Inv,eigenvals);
	CompactMatrixMultiplyDiagonal(eigenvals,T,temp);
	CompactMatrixMultiply('T','N',&T,&temp,&S12Inv);
	CompactMatrixMultiply3('N','N','N',&S12Inv,&H,&S12Inv,&HTrafo);

	FreeCompactMatrix(&T);
	FreeCompactMatrix(&temp);
	FreeCompactMatrix(&S12Inv);
}

void CalculateEigensystemOnBasis(const vector<WaveFunctionType> &basis,const OperatorType &H, CompactMatrixType &T, vector<double> &eigenvals)
{
	CompactMatrixType M, S, S12, Hbasis, eigenfuns;
	CalculateHermitianMatrixElements(basis,H,M);
	CalculateOverlapMatrix(basis,S);
	TransformMatrixToONB(S,M,S12,Hbasis);
	eigenvals = vector<double>(M.Ni);
	CompactMatrixEigensystem(&Hbasis,&eigenfuns,&eigenvals[0]); // TODO check if eigenvals row- or coll-major
	CompactMatrixMultiply('N','N',&eigenfuns,&S12,&T);

	FreeCompactMatrix(&M);
	FreeCompactMatrix(&S);
	FreeCompactMatrix(&S12);
	FreeCompactMatrix(&Hbasis);
	FreeCompactMatrix(&eigenfuns);
}

void CalculateEnergyIndependentSelfEnergy(const OperatorType &H, const WaveFunctionType  &psi, const vector<OperatorType> &UList, CompactMatrixType &T, vector<double> &poles)
{
	vector<WaveFunctionType> Upsi = vector<WaveFunctionType>(UList.size());
	for (uint32_t i = 0; i < UList.size(); i++)
	{
		OperatorPsi(UList[i],psi,&Upsi[i]);
	}
	CalculateEigensystemOnBasis(Upsi,H,T,poles);


}

void CreateAugerMeitnerOperators(SSlaterIntegral *integrals, uint32_t NIntegrals,uint32_t NF, uint16_t **indices, vector<OperatorType> &UList)
{
	for(uint32_t i = 0; i < NIntegrals; i++)
	{
		int *k = integrals[i].tau.kappa;
		uint32_t *n = integrals[i].tau.n;
		double *Fk = integrals[i].val;
		for (uint32_t j = 0; j < integrals[i].Nval; j++)
		{
			Fk[j] = 0;
		}
		for (uint32_t j = 0; j < integrals[i].Nval; j++)
		{
			OperatorType optemp;

			Fk[j] = 1;
			CreateOperatorAtomicURelativisticOneIntegral(NF,k[0],k[1],k[2],k[3],indices[n[0]],indices[n[1]],indices[n[2]],indices[n[3]],Fk,n[0],n[1],n[2],n[3],&optemp);
			UList.push_back(optemp);
			Fk[j] = 0;
		}

	}
}

void InitializeSpectrum(lua_State *L, int arg, vector<mypoint> &spectrum)
{
	double wMin, wMax, dw;
	lua_check(L,arg-1+1,wMin);
	lua_check(L,arg-1+2,wMax);
	lua_check(L,arg-1+3,dw);

	uint32_t Npoints = floor((wMax-wMin)/dw) + 1;
	spectrum = std::vector<mypoint>(Npoints);
	for (uint32_t i = 0; i < Npoints; i++)
	{
		spectrum[i].w = wMin + i * dw;
		spectrum[i].y.Re = 0.;
		spectrum[i].y.Im = 0.;
	}
}




